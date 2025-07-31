"""
Tests for helper function `merge_voltage_clusters` contained in
`plasmapy.analysis.swept_langmuir.helpers`.
"""

from contextlib import nullcontext as does_not_raise
from unittest import mock

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
        (  # voltage array has np.nan values (hence not monotonic)
            np.array([1.0, 1.2, np.nan, 2.5]),
            np.array([-1, 0, np.nan, 4.3]),
            {"filter_nan": False},
            pytest.raises(ValueError),
            None,
        ),
        (  # current array has np.nan values (hence not monotonic)
            np.array([1.0, 1.2, 1.8, 2.5]),
            np.array([-1, 0, np.nan, 4.3]),
            {"filter_nan": False},
            pytest.raises(ValueError),
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
            np.array([1.1, 1.5, 2, 3.8, 6.0, 6.2], dtype=float),
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
        (  # no merging needed
            np.linspace(-40.0, 40, 100),
            np.linspace(-10.0, 30, 100),
            {"voltage_step_size": 0.0},
            does_not_raise(),
            None,  # same as inputs
        ),
        (  # merge identical voltages
            np.array([1, 2, 2, 4, 6, 6], dtype=float),
            np.array([-5, -2, -2.1, 0, 5, 4.9], dtype=float),
            {"voltage_step_size": 0.0},
            does_not_raise(),
            (
                np.array([1, 2, 4, 6], dtype=float),
                np.array([-5, -2.05, 0, 4.95], dtype=float),
            ),
        ),
        (  # step size given as an integer
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
                np.array([1.0, 1.15, 1.5, 2, 4, 5.5], dtype=float),
                np.array([-5.0, -5, -4, 0, 5, 7], dtype=float),
            ),
        ),
        (  # voltage array starts with a cluster, spans < 2 * voltage_step_size
            np.array([1, 1.1, 1.2, 1.5, 2, 4, 5.5], dtype=float),
            np.array([-5, -5.1, -4.9, -4, 0, 5, 7], dtype=float),
            {"voltage_step_size": 0.12},
            does_not_raise(),
            (
                np.array([1.1, 1.5, 2, 4, 5.5], dtype=float),
                np.array([-5.0, -4, 0, 5, 7], dtype=float),
            ),
        ),
        (  # multiple voltage clusters, including at beginning and end
            np.array(
                [1, 1.01, 1.05, 1.5, 2, 3.88, 3.9, 3.92, 4.8, 5.5, 5.55],
                dtype=float,
            ),
            np.array(
                [-5, -5.1, -4.9, -4, -3, 0, 0.1, -0.1, 5, 7, 6.9],
                dtype=float,
            ),
            {"voltage_step_size": 0.1},
            does_not_raise(),
            (
                np.array([1.02, 1.5, 2, 3.9, 4.8, 5.525], dtype=float),
                np.array([-5.0, -4, -3, 0, 5, 6.95], dtype=float),
            ),
        ),
        (  # self determine voltage_step_size
            np.array(
                [1, 1.01, 1.05, 1.5, 2, 3.88, 3.9, 3.92, 4.8, 5.5, 5.55],
                dtype=float,
            ),
            np.array(
                [-5, -5.1, -4.9, -4, -3, 0, 0.1, -0.1, 5, 7, 6.9],
                dtype=float,
            ),
            {"voltage_step_size": None},
            does_not_raise(),
            (
                np.array([1.02, 1.5, 2, 3.9, 4.8, 5.525], dtype=float),
                np.array([-5.0, -4, -3, 0, 5, 6.95], dtype=float),
            ),
        ),
        (  # test case for auto filtering NaN values
            np.array(
                [1, 1.01, np.nan, 1.5, 2, 3.88, 3.9, np.nan, 4.8, 5.5, 5.55],
                dtype=float,
            ),
            np.array(
                [-5, -5.1, np.nan, -4, -3, 0, 0.1, np.nan, np.nan, 7, 6.9],
                dtype=float,
            ),
            {"filter_nan": True, "voltage_step_size": 0.12},
            does_not_raise(),
            (
                np.array([1.005, 1.5, 2, 3.89, 5.525], dtype=float),
                np.array([-5.05, -4, -3, 0.05, 6.95], dtype=float),
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
            assert np.allclose(rtn_current, expected[1])
        else:
            assert np.allclose(rtn_voltage, voltage)
            assert np.allclose(rtn_current, current)

        # ensure check_sweep() is called...then we do not have to
        # explicitly add test cases covered by check_sweep()
        assert 0 < mock_sweep.call_count <= 2
        mock_sweep.reset_mock()
