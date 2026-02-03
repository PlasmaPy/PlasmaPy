"""
Tests for
`plasmapy.analysis.swept_langmuir.helpers._condition_voltage_window`.
"""

import numpy as np
import pytest

from plasmapy.analysis.swept_langmuir.helpers import _condition_voltage_window


class TestConditionVoltageWindow:
    @pytest.mark.parametrize(
        ("_raises", "voltage", "voltage_window"),
        [
            # voltage_window is not a list-like or None
            (pytest.raises(TypeError), np.arange(10.0), {"one": 1, "two": 2}),
            (pytest.raises(TypeError), np.arange(10.0), "invalid window"),
            (pytest.raises(TypeError), np.arange(10.0), 5.0),
            # voltage_window is wrong size
            (pytest.raises(ValueError), np.arange(10.0), []),
            (pytest.raises(ValueError), np.arange(10.0), [5.0]),
            (pytest.raises(ValueError), np.arange(10.0), (-1.0, 2.0, 5.0)),
            # voltage_window does not have Real numbers or None
            (pytest.raises(TypeError), np.arange(10.0), ["one", 2]),
            (pytest.raises(TypeError), np.arange(10.0), ["one", "two"]),
            (pytest.raises(TypeError), np.arange(10.0), [(1.0, 2.0), 2.0]),
            # voltage_window is out of range
            (pytest.raises(ValueError), np.arange(10.0), [-5, -2]),
            (pytest.raises(ValueError), np.arange(10.0), [11, 18]),
        ],
    )
    def test_raises(self, _raises, voltage, voltage_window):
        with _raises:
            _condition_voltage_window(voltage, voltage_window)

    @pytest.mark.parametrize(
        ("voltage", "voltage_window", "expected"),
        [
            # voltage_window as list and numpy array
            (np.arange(10), [3.3, 4.5], slice(4, 5, 1)),
            (np.arange(10), np.array([3.3, 4.5]), slice(4, 5, 1)),
            (np.arange(10), [3.3, 7.2], slice(4, 8, 1)),
            (np.arange(10), np.array([3.3, 7.2]), slice(4, 8, 1)),
            # voltage_window has None values
            (np.arange(10), None, slice(None, None, 1)),
            (np.arange(10), [None, None], slice(None, None, 1)),
            (np.arange(10), [3.3, None], slice(4, None, 1)),
            (np.arange(10), [None, 7.2], slice(None, 8, 1)),
            # voltage_window is unsorted
            (np.arange(10), [7.2, 3.3], slice(4, 8, 1)),
            # voltage_window has one index out of bounds
            (np.arange(10), [-5, 7.2], slice(None, 8, 1)),
            (np.arange(10), [3.3, 20], slice(4, None, 1)),
        ],
    )
    def test_expected(self, voltage, voltage_window, expected):
        result = _condition_voltage_window(voltage, voltage_window)
        assert result == expected
