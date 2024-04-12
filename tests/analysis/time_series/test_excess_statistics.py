"""Tests for excess_statistics.py"""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.analysis.time_series.excess_statistics import ExcessStatistics


@pytest.mark.parametrize(
    ("signal", "thresholds", "time_step", "pdf", "bins", "expected"),
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
            1,
            1 * u.s,
            False,
            32,
            ([1] * u.s, [1], [1] * u.s, [0] * u.s),
        ),
        (
            [0, 2, 0],
            1,
            0.5,
            False,
            32,
            ([0.5], [1], [0.5], [0]),
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
def test_ExcessStatistics(signal, thresholds, time_step, pdf, bins, expected) -> None:
    """Test ExcessStatistics class"""
    excess_stats = ExcessStatistics(signal, thresholds, time_step)
    assert excess_stats.total_time_above_threshold == expected[0]
    assert excess_stats.number_of_crossings == expected[1]
    assert excess_stats.average_times == expected[2]
    assert excess_stats.rms_times == expected[3]
    if pdf:
        hist, bin_centers = excess_stats.hist(bins)
        assert np.allclose(hist, expected[4])
        assert np.allclose(bin_centers, expected[5])


@pytest.mark.parametrize(
    ("signal", "thresholds", "time_step", "bins", "exception"),
    [([1, 2], 1, -1, 32, ValueError), ([1, 2], 1, 1, 1.5, TypeError)],
)
def test_ExcessStatistics_exception(
    signal, thresholds, time_step, bins, exception
) -> None:
    """Test whether exception is risen"""
    with pytest.raises(exception):
        tmp = ExcessStatistics(signal, thresholds, time_step)
        tmp.hist(bins)
