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
        # values
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(-10.0, 30, 100),
            {"voltage_step_size": None},
            pytest.warns(PlasmaPyWarning),
            None,  # same as inputs
        ),
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
            ),  # same as inputs
        ),
    ],
)
def test_merge_voltage_clusters(
    voltage, current, kwargs, with_context, expected
) -> None:
    with (
        with_context,
        mock.patch(
            "plasmapy.analysis.swept_langmuir.helpers.merge_voltage_clusters",
            wraps=merge_voltage_clusters,
        ) as mock_merge,
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

        # mock_merge.assert_called_once()
        # mock_merge.reset_mock()

        mock_sweep.assert_called_once()
        mock_sweep.reset_mock()
