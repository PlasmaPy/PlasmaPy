"""
Tests for
`plasmapy.analysis.swept_langmuir.plasma_potential._condition_smooth_fractions`.
"""

import numpy as np
import pytest

from plasmapy.analysis.swept_langmuir.plasma_potential import (
    _condition_smooth_fractions,
)


class TestConditionSmoothFractions:
    @pytest.mark.parametrize(
        ("_raises", "smooth_fractions", "data_size"),
        [
            # smooth_fractions is not a list-like or None
            (pytest.raises(TypeError), {"one": 1, "two": 2}, 10),
            (pytest.raises(TypeError), "invalid fraction", 10),
            (pytest.raises(TypeError), 5.0, 10),
            # smooth_fractions is not 1D
            (pytest.raises(ValueError), [[1, 2], [1, 2]], 10),
            (pytest.raises(ValueError), np.zeros((2, 5)), 10),
            # smooth_fractions elements are NOT floats
            (pytest.raises(ValueError), [1, 2, 3], 10),
            (pytest.raises(ValueError), ["one", 2.0, "three"], 10),
            # smooth_fractions is out of range.. [0, 1]
            (pytest.raises(ValueError), [-5.0], 10),
            (pytest.raises(ValueError), [5.0], 10),
            (pytest.raises(ValueError), [-0.5, 0.1, 0.9, 1.0], 10),
            (pytest.raises(ValueError), [-0.5, 0, 0.9, 1.0], 10),
            (pytest.raises(ValueError), [0.1, 0.9, 1.1], 10),
            (pytest.raises(ValueError), [0.1, 0.9, 20.0], 10),
            (pytest.raises(ValueError), [0, 0.1, 0.9, 1.1], 10),
            (pytest.raises(ValueError), [-20.0, 0, 0.1, 0.9, 1.1], 10),
            # data_size not integer
            (pytest.raises(TypeError), [0.2, 0.5, 0.8], "hello"),
            (pytest.raises(TypeError), [0.2, 0.5, 0.8], 90.5),
            # data_size not positive, non-zero integer
            (pytest.raises(ValueError), [0.2, 0.5, 0.8], 0),
            (pytest.raises(ValueError), [0.2, 0.5, 0.8], -20),
            # # computed savgol_windows is null
            (pytest.raises(ValueError), [0.2, 0.5, 0.8], 1),
            (pytest.raises(ValueError), [0.02], 20),
        ],
    )
    def test_raises(self, _raises, smooth_fractions, data_size):
        with _raises:
            _condition_smooth_fractions(smooth_fractions, data_size)

    @pytest.mark.parametrize(
        ("smooth_fractions", "data_size", "expected"),
        [
            ([0.2, 0.5, 0.8], 20, np.array([4, 10, 16], dtype=int)),
            ([0.02, 0.2, 0.5, 0.8], 20, np.array([4, 10, 16], dtype=int)),
            ([0.0, 0.02, 0.2, 0.5, 0.8], 20, np.array([0, 4, 10, 16], dtype=int)),
            (None, 50, np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=int)),
        ],
    )
    def test_expected(self, smooth_fractions, data_size, expected):
        result = _condition_smooth_fractions(smooth_fractions, data_size)
        assert np.array_equal(result, expected)
