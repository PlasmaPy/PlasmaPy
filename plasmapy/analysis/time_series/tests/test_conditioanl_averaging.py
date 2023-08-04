"""Tests for excess_averaging.py"""

from math import exp
import astropy.units as u
import numpy as np
from numpy.linalg import cond
import pytest

from plasmapy.analysis.time_series.conditional_averaging import ConditionalEvents


@pytest.mark.parametrize(
    (
        "signal",
        "time",
        "lower_threshold",
        "upper_threshold",
        "reference_signal",
        "length_of_return",
        "distance",
        "exception",
    ),
    [
        ([1, 2], [1, 2], 1.5, None, None, None, -1, ValueError),
        ([1, 2], [1, 2, 3], 1.5, None, None, None, 0, ValueError),
        ([1, 2], [1, 2], 1.5, None, [1, 2, 3], None, 0, ValueError),
        ([1, 2], [1, 2], 1.5, None, None, 5, 0, ValueError),
        ([1, 2], [1, 2], 1.5, None, None, -5, 0, ValueError),
        ([1, 2], [1, 2], 1.5, 1, None, None, 0, ValueError),
    ],
)
def test_ConditionalEvents_ValueErrors(
    signal,
    time,
    lower_threshold,
    upper_threshold,
    reference_signal,
    length_of_return,
    distance,
    exception,
):
    """Test whether exception is risen"""
    with pytest.raises(exception):
        tmp = ConditionalEvents(
            signal,
            time,
            lower_threshold,
            upper_threshold,
            reference_signal,
            length_of_return,
            distance,
        )


@pytest.mark.parametrize(
    (
        "signal",
        "time",
        "lower_threshold",
        "upper_threshold",
        "reference_signal",
        "length_of_return",
        "distance",
        "expected",
    ),
    [
        (
            np.array([1, 2, 1, 1, 2, 1]),
            np.array([1, 2, 3, 4, 5, 6]),
            1.5,
            None,
            None,
            None,
            0,
            (
                [-1.0, 0.0, 1.0],
                [1.0, 2.0, 1.0],
                [1.0, 1.0, 1.0],
                [2.0, 2.0],
                [3.0],
                [2.0, 5.0],
                2,
            ),
        ),
        (
            [1, 2, 1, 1, 2, 1],
            [1, 2, 3, 4, 5, 6],
            1.5,
            None,
            None,
            None,
            0,
            [
                [-1.0, 0.0, 1.0],
                [1.0, 2.0, 1.0],
                [1.0, 1.0, 1.0],
                [2.0, 2.0],
                [3.0],
                [2.0, 5.0],
                2,
            ],
        ),
        (
            [1, 2, 1, 1, 2, 1, 1, 5, 1],
            [1, 2, 3, 4, 5, 6, 7, 8, 9],
            1.5,
            4,
            None,
            3,
            0,
            [
                [-1.0, 0.0, 1.0],
                [1.0, 2.0, 1.0],
                [1.0, 1.0, 1.0],
                [2.0, 2.0],
                [3.0],
                [2.0, 5.0],
                2,
            ],
        ),
        (
            [1, 3, 2, 3, 1, 1, 1, 2, 1, 1, 4, 1],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            1.5,
            None,
            None,
            3,
            0,
            [
                [-1.0, 0.0, 1.0],
                [1.0, 3.0, 1.333333],
                [1.0, 0.93103448, 0.88888888],
                [3.0, 2.0, 4.0],
                [6.0, 3.0],
                [2.0, 8.0, 11.0],
                3,
            ],
        ),
    ],
)
def test_ConditionalEvents_exception(
    signal,
    time,
    lower_threshold,
    upper_threshold,
    reference_signal,
    length_of_return,
    distance,
    expected,
):
    """Test ConditionalEvents class"""
    cond_events = ConditionalEvents(
        signal,
        time,
        lower_threshold,
        upper_threshold,
        reference_signal,
        length_of_return,
        distance,
    )
    assert np.allclose(cond_events.time, expected[0])
    assert np.allclose(cond_events.average, expected[1])
    assert np.allclose(cond_events.variance, expected[2])
    assert np.allclose(cond_events.peaks, expected[3])
    assert np.allclose(cond_events.waiting_times, expected[4])
    assert np.allclose(cond_events.arrival_times, expected[5])
    assert np.allclose(cond_events.number_of_events, expected[6])
