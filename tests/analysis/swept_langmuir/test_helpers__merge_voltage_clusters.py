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
            {"force_regular_spacing": None},  # not a bool
            pytest.raises(TypeError),
            None,
        ),
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
        (  # voltage_step_size = 0 and spacing is irregular
            np.array([1.1, 1.1,  2, 3.8, 6.1, 6.1], dtype=float),
            np.array([-5, -5.2, -2, 0, 4.2, 5], dtype=float),
            {"voltage_step_size": 0, "force_regular_spacing": True},
            pytest.warns(PlasmaPyWarning),
            (
                np.array([1.1, 2, 3.8, 6.1], dtype=float),
                np.array([-5.1, -2, 0, 4.6], dtype=float),
            ),
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
            np.array([1, 2, 4, 6], dtype=float),
            np.array([-5, -2, 0, 5], dtype=float),
            {"voltage_step_size": 0.0, "force_regular_spacing": True},
            does_not_raise(),
            (
                np.array([1, 2, 3, 4, 5, 6], dtype=float),
                np.array([-5, -2, np.nan, 0, np.nan, 5], dtype=float),
            ),
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
            {"voltage_step_size": 0, "force_regular_spacing": True},
            does_not_raise(),
            (
                np.array([1, 2, 3, 4, 5, 6], dtype=float),
                np.array([-5, -2.05, np.nan, 0, np.nan, 4.95], dtype=float),
            ),
        ),
        # voltage_step_size < all actual diffs
        (
            np.array([1, 1.5, 2.1, 4, 5.5, 6], dtype=float),
            np.array([-5, -2, -2, 0, 5, 5], dtype=float),
            {"voltage_step_size": 0.2, "force_regular_spacing": False},
            does_not_raise(),
            None,
        ),
        (
            np.array([1, 1.5, 2.1, 4, 5.5, 6], dtype=float),
            np.array([-5, -2, -2, 0, 5, 5], dtype=float),
            {"voltage_step_size": 0.4, "force_regular_spacing": True},
            does_not_raise(),
            (
                0.4 * np.arange(13, dtype=float) + 1.0,
                np.array(
                    [
                        -5,
                        -2.6,
                        -2,
                        -1.894736,
                        -1.473684,
                        -1.052631,
                        -0.631578,
                        -0.210526,
                        0.666667,
                        2.0,
                        3.333333,
                        4.666667,
                        5.0,
                    ],
                    dtype=float,
                ),  # generated using np.interp
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
