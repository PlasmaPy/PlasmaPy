"""
Tests for helper function `merge_voltage_clusters` contained in
`plasmapy.analysis.swept_langmuir.helpers`.
"""

from contextlib import nullcontext as does_not_raise
from unittest import mock

import astropy.units as u
import numpy as np
import pytest

from plasmapy.analysis.swept_langmuir.helpers import check_sweep, merge_voltage_clusters
from plasmapy.utils.exceptions import PlasmaPyWarning


@pytest.mark.parametrize(
    ("voltage", "current", "kwargs", "with_context", "expected"),
    [
        # raises
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(-10.0, 30, 100),
            {"voltage_step_size": "wrong type"},
            pytest.raises(TypeError),
            None,
        ),
        # warnings (and values)
        (  # voltage spacing is regular and not step size given
            np.linspace(-40.0, 40, 100),
            np.linspace(-10.0, 30, 100),
            {"voltage_step_size": None},
            pytest.warns(PlasmaPyWarning),
            None,  # same as inputs
        ),
        (  # non-zero voltage_step_size < any step size
            np.array([1.1, 1.5,  2, 3.8, 6.0, 6.2], dtype=float),
            np.array([-5, -5.2, -2, 0, 4.2, 5], dtype=float),
            {"voltage_step_size": 0.1},
            pytest.warns(PlasmaPyWarning),
            None,
        ),
        (
            np.array([1, 1.5, 2.1, 4, 5.5, 6], dtype=float),
            np.array([-5, -2, -2, 0, 5, 5], dtype=float),
            {"voltage_step_size": 0.2},
            pytest.warns(PlasmaPyWarning),
            None,
        ),
        # values
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(-10.0, 30, 100),
            {"voltage_step_size": 0.0},
            does_not_raise(),
            None,  # same as inputs
        ),
        (
            np.array([1, 2, 2, 4, 6, 6], dtype=float),
            np.array([-5, -2, -2.1, 0, 5, 4.9], dtype=float),
            {"voltage_step_size": 0.0},
            does_not_raise(),
            (
                np.array([1, 2, 4, 6], dtype=float),
                np.array([-5, -2.05, 0, 4.95], dtype=float),
            ),
        ),
        (
            np.array([1, 2, 2, 4, 6, 6], dtype=float),
            np.array([-5, -2, -2.1, 0, 5, 4.9], dtype=float),
            {"voltage_step_size": 0},  # integer value
            does_not_raise(),
            (
                np.array([1, 2, 4, 6], dtype=float),
                np.array([-5, -2.05, 0, 4.95], dtype=float),
            ),
        ),
        (  # voltage_step_size is negative
            np.array([1, 1.5, 1.56, 2.1, 4, 5.5, 5.52], dtype=float),
            np.array([-5, -2, -2.1, 0, 5, 7, 7.1], dtype=float),
            {"voltage_step_size": -0.1},
            does_not_raise(),
            (
                np.array([1, 1.53, 2.1, 4, 5.51], dtype=float),
                np.array([-5, -2.05, 0, 5, 7.05], dtype=float),
            ),
        ),
        (  # voltage array starts with a cluster
            np.array([1, 1.1, 1.2, 1.5, 2, 4, 5.5], dtype=float),
            np.array([-5, -5.1, -4.9, -4, 0, 5, 7], dtype=float),
            {"voltage_step_size": 0.1},
            does_not_raise(),
            (
                np.array([1., 1.15, 1.5, 2, 4, 5.5], dtype=float),
                np.array([-5., -5, -4, 0, 5, 7], dtype=float),
            ),
        ),
        (  # voltage array starts with a cluster, spans < 2 * voltage_step_size
            np.array([1, 1.1, 1.2, 1.5, 2, 4, 5.5], dtype=float),
            np.array([-5, -5.1, -4.9, -4, 0, 5, 7], dtype=float),
            {"voltage_step_size": 0.12},
            does_not_raise(),
            (
                np.array([1.1, 1.5, 2, 4, 5.5], dtype=float),
                np.array([-5., -4, 0, 5, 7], dtype=float),
            ),
        ),
    ],
)
def test_merge_voltage_clusters(
    voltage, current, kwargs, with_context, expected
) -> None:
    with (
        with_context,
        mock.patch(
            "plasmapy.analysis.swept_langmuir.helpers.check_sweep",
            wraps=check_sweep,
        ) as mock_sweep,
    ):
        rtn_voltage, rtn_current = merge_voltage_clusters(
            voltage=voltage,
            current=current,
            **kwargs,
        )

        if expected is not None:
            assert np.allclose(rtn_voltage, expected[0])

            rtn_nan_mask = np.isnan(rtn_current)
            expected_nan_mask = np.isnan(expected[1])
            assert np.array_equal(rtn_nan_mask, expected_nan_mask)
            assert np.allclose(
                rtn_current[np.logical_not(rtn_nan_mask)],
                expected[1][np.logical_not(expected_nan_mask)],
            )
        else:
            assert np.allclose(rtn_voltage, voltage)
            assert np.allclose(rtn_current, current)

        mock_sweep.assert_called_once()
        mock_sweep.reset_mock()
