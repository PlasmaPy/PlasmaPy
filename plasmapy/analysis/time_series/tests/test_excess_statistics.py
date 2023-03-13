"""Tests for excess_statistics.py"""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.analysis.time_series.excess_statistics import ExcessStatistics


@pytest.mark.parametrize(
    "signal, thresholds, time_step, pdf, bins, expected",
    [
        (
            [0, 2, 0],
            1,
            1,
            False,
            32,
            ([1], [1], [1], [0]),
        ),
        (
            [0, 2, 0],
            1.5,
            1,
            False,
            32,
            ([1], [1], [1], [0]),
        ),
        (
            [0, 2, 2, 0, 4, 0],
            [1, 3, 5],
            1,
            False,
            32,
            ([3, 1, 0], [2, 1, 0], [1.5, 1, 0], [0.5, 0, 0]),
        ),
        (
            np.array([2, 2, 0, 4, 0]),
            [1, 3, 5],
            1,
            False,
            32,
            ([3, 1, 0], [1, 1, 0], [1.5, 1, 0], [0.5, 0, 0]),
        ),
        (
            [0, 2, 0, 4, 0] * u.eV,
            [1, 3, 5],
            1,
            False,
            32,
            ([2, 1, 0], [2, 1, 0], [1, 1, 0], [0, 0, 0]),
        ),
        (
            [0, 2, 0, 4, 0],
            [1, 3, 5],
            1,
            True,
            2,
            (
                [2, 1, 0],
                [2, 1, 0],
                [1, 1, 0],
                [0, 0, 0],
                [[0, 2], [0, 2], [0, 0]],
                [[0.75, 1.25], [0.75, 1.25], [0, 0]],
            ),
        ),
    ],
)
def test_ExcessStatistics(signal, thresholds, time_step, pdf, bins, expected):
    """test excess_stat function"""
    tmp = ExcessStatistics(signal, thresholds, time_step)
    assert tmp.total_time_above_threshold == expected[0]
    assert tmp.number_of_crossings == expected[1]
    assert tmp.average_times == expected[2]
    assert tmp.rms_times == expected[3]
    if pdf:
        hist, bin_centers = tmp.hist(bins)
        assert np.allclose(hist, expected[4])
        assert np.allclose(bin_centers, expected[5])


@pytest.mark.parametrize(
    "signal, thresholds, time_step, bins",
    [([1, 2], 1, -1, 32), ([1, 2], 1, 1, 1.5)],
)
def test_ExcessStatistics_exception(signal, thresholds, time_step, bins):
    """test whether exception is risen"""
    with pytest.raises((TypeError, ValueError)):
        tmp = ExcessStatistics(signal, thresholds, time_step)
        tmp.hist(bins)
