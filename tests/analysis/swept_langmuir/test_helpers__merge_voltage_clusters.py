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
        # (
        #     np.linspace(-40.0, 40, 100),
        #     np.linspace(-10.0, 30, 100),
        #     {"voltage_order": "ascending"},
        #     does_not_raise(),
        #     None,
        # ),
        # (
        #     np.linspace(-40.0, 40, 100),
        #     np.linspace(-10.0, 30, 100),
        #     {"voltage_order": "descending"},
        #     does_not_raise(),
        #     (np.linspace(40.0, -40, 100), np.linspace(30.0, -10, 100)),
        # ),
        # (
        #     np.linspace(40.0, -40, 100),
        #     np.linspace(-10.0, 30, 100),
        #     {"voltage_order": "descending"},
        #     does_not_raise(),
        #     None,
        # ),
        # (
        #     np.linspace(40.0, -40, 100),
        #     np.linspace(-10.0, 30, 100),
        #     {"voltage_order": "ascending"},
        #     does_not_raise(),
        #     (np.linspace(-40.0, 40, 100), np.linspace(30.0, -10, 100)),
        # ),
        # (
        #     np.array([-40.0, 20.0, -22.0, 40.0]),
        #     np.array([-10.0, 10.0, -5.0, 30.0]),
        #     {"voltage_order": "ascending"},
        #     does_not_raise(),
        #     (np.array([-40.0, -22.0, 20.0, 40.0]), np.array([-10.0, -5.0, 10.0, 30.0])),
        # ),
        # (
        #     np.array([-40.0, 20.0, -22.0, 40.0]),
        #     np.array([-10.0, 10.0, -5.0, 30.0]),
        #     {"voltage_order": "descending"},
        #     does_not_raise(),
        #     (np.array([40.0, 20.0, -22.0, -40.0]), np.array([30.0, 10.0, -5.0, -10.0])),
        # ),
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
            assert np.allclose(rtn_current, expected[1])
        else:
            assert np.allclose(rtn_voltage, voltage)
            assert np.allclose(rtn_current, current)

        # mock_merge.assert_called_once()
        # mock_merge.reset_mock()

        mock_sweep.assert_called_once()
        mock_sweep.reset_mock()
