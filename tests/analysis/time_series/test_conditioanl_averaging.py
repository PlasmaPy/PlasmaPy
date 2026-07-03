"""Tests for conditional_averaging.py"""

import astropy.units as u
import numpy as np
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
        ([1, 2] * u.eV, [1, 2], 1.5, None, None, None, 0, u.UnitsError),
        ([1, 2], [1, 2], 1.5 * u.eV, None, None, None, 0, u.UnitsError),
        ([1, 2] * u.eV, [1, 2], 1.5 * u.m, None, None, None, 0, u.UnitsError),
        ([1, 2] * u.eV, [1, 2], 1.5 * u.eV, 4.0, None, None, 0, u.UnitsError),
        ([1, 2], [1, 2], 1.5, None, [1, 2] * u.eV, None, 0, u.UnitsError),
        ([1, 2], [1, 2], 1.5 * u.eV, 4.0, [1, 2] * u.eV, None, 0, u.UnitsError),
        ([1, 2], [1, 2] * u.s, 1.5, 4.0, [1, 2], None, 0, u.UnitsError),
        ([1, 2], [1, 2], 1.5, 4.0, [1, 2], None, 0 * u.s, u.UnitsError),
        ([1, 2], [1, 2] * u.s, 1.5, 4.0, [1, 2], 1, 0 * u.s, u.UnitsError),
        ([1, 2], [1, 2, 3], 1.5, None, None, None, 0, ValueError),
        ([1, 2], [1, 2], 1.5, None, [1, 2, 3], None, 0, ValueError),
        ([1, 2], [1, 2], 1.5, None, None, 5, 0, ValueError),
        ([1, 2], [1, 2], 1.5, None, None, -5, 0, ValueError),
        ([1, 2], [1, 2], 1.5, 1, None, None, 0, ValueError),
    ],
)
def test_ConditionalEvents_Errors(
    signal,
    time,
    lower_threshold,
    upper_threshold,
    reference_signal,
    length_of_return,
    distance,
    exception,
) -> None:
    """Test whether exception is risen"""
    with pytest.raises(exception):
        ConditionalEvents(
            signal,
            time,
            lower_threshold,
            upper_threshold=upper_threshold,
            reference_signal=reference_signal,
            length_of_return=length_of_return,
            distance=distance,
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
            [1, 2, 3, 4, 5, 6, 7, 8, 9] * u.s,
            1.5,
            4,
            None,
            3 * u.s,
            0 * u.s,
            [
                [-1.0, 0.0, 1.0] * u.s,
                [1.0, 2.0, 1.0],
                [1.0, 1.0, 1.0],
                [2.0, 2.0],
                [3.0] * u.s,
                [2.0, 5.0] * u.s,
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
        (
            [1, 3, 2, 3, 1, 1, 1, 2, 1, 1, 4, 1],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            1.5,
            None,
            [1, 1, 1, 5, 1, 1, 9, 1, 1, 1, 1, 1],
            3,
            0,
            [
                [-1.0, 0.0, 1.0],
                [1.5, 2.0, 1.5],
                [0.9, 0.8, 0.9],
                [3.0, 1.0],
                [3.0],
                [4.0, 7.0],
                2,
            ],
        ),
        (
            [1, 3, 2, 3, 1, 1, 1, 2, 1, 1, 4, 1],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            1.5 * u.eV,
            None,
            [1, 1, 1, 5, 1, 1, 9, 1, 1, 1, 1, 1] * u.eV,
            3,
            0,
            [
                [-1.0, 0.0, 1.0],
                [1.5, 2.0, 1.5],
                [0.9, 0.8, 0.9],
                [3.0, 1.0],
                [3.0],
                [4.0, 7.0],
                2,
            ],
        ),
        (
            [1, 3, 2, 3, 1, 1, 1, 2, 1, 1, 4, 1] * u.K,
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            1.5 * u.eV,
            None,
            [1, 1, 1, 5, 1, 1, 9, 1, 1, 1, 1, 1] * u.eV,
            3,
            0,
            [
                [-1.0, 0.0, 1.0],
                [1.5, 2.0, 1.5] * u.K,
                [0.9, 0.8, 0.9],
                [3.0, 1.0] * u.K,
                [3.0],
                [4.0, 7.0],
                2,
            ],
        ),
    ],
)
def test_ConditionalEvents_class(
    signal,
    time,
    lower_threshold,
    upper_threshold,
    reference_signal,
    length_of_return,
    distance,
    expected,
) -> None:
    """Tests for ConditionalEvents class"""
    cond_events = ConditionalEvents(
        signal,
        time,
        lower_threshold,
        upper_threshold=upper_threshold,
        reference_signal=reference_signal,
        length_of_return=length_of_return,
        distance=distance,
    )
    assert np.allclose(cond_events.time, expected[0])
    assert np.allclose(cond_events.average, expected[1])
    assert np.allclose(cond_events.variance, expected[2])
    assert np.allclose(cond_events.peaks, expected[3])
    assert np.allclose(cond_events.waiting_times, expected[4])
    assert np.allclose(cond_events.arrival_times, expected[5])
    assert np.allclose(cond_events.number_of_events, expected[6])


@pytest.mark.parametrize(
    (
        "signal",
        "time",
        "lower_threshold",
        "upper_threshold",
        "reference_signal",
        "length_of_return",
        "distance",
        "remove_non_max_peaks",
        "expected",
    ),
    [
        (
            [1, 2, 5, 3, 1, 2, 1, 1, 1, 1, 1, 1],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            1.5,
            None,
            None,
            5,
            0,
            True,
            [
                [-2.0, -1.0, 0.0, 1.0, 2.0],
                [1.0, 2.0, 5.0, 3.0, 1.0],
                [1, 1, 1, 1, 1],
                [5],
                [3.0],
                [3.0],
                1,
            ],
        ),
        (
            [1, 3, 2, 3, 1, 3, 1, 2, 1, 1, 1, 1],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            1.5,
            None,
            [1, 2, 5, 3, 1, 2, 1, 1, 1, 1, 1, 1],
            5,
            0,
            True,
            [
                [-2.0, -1.0, 0.0, 1.0, 2.0],
                [1.0, 3.0, 2.0, 3.0, 1.0],
                [1, 1, 1, 1, 1],
                [2],
                [3.0],
                [3.0],
                1,
            ],
        ),
    ],
)
def test_peak_not_max_value(
    signal,
    time,
    lower_threshold,
    upper_threshold,
    reference_signal,
    length_of_return,
    distance,
    remove_non_max_peaks,
    expected,
) -> None:
    """Tests for ConditionalEvents class"""
    cond_events = ConditionalEvents(
        signal,
        time,
        lower_threshold,
        upper_threshold=upper_threshold,
        reference_signal=reference_signal,
        length_of_return=length_of_return,
        distance=distance,
        remove_non_max_peaks=remove_non_max_peaks,
    )
    assert np.allclose(cond_events.time, expected[0])
    assert np.allclose(cond_events.average, expected[1])
    assert np.allclose(cond_events.variance, expected[2])
    assert np.allclose(cond_events.peaks, expected[3])
    assert np.allclose(cond_events.waiting_times, expected[4])
    assert np.allclose(cond_events.arrival_times, expected[5])
    assert np.allclose(cond_events.number_of_events, expected[6])
