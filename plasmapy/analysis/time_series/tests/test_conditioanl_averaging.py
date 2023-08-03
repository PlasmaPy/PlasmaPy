"""Tests for excess_averaging.py"""

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
        "weight",
        "exception",
    ),
    [
        ([1, 2], [1, 2], 1.5, None, None, None, 0, "false_input", ValueError),
        ([1, 2], [1, 2], 1.5, None, None, None, -1, "amplitude", ValueError),
        ([1, 2], [1, 2, 3], 1.5, None, None, None, 0, "amplitude", ValueError),
        ([1, 2], [1, 2], 1.5, None, [1, 2, 3], None, 0, "amplitude", ValueError),
        ([1, 2], [1, 2], 1.5, None, None, 5, 0, "amplitude", ValueError),
        ([1, 2], [1, 2], 1.5, None, None, -5, 0, "amplitude", ValueError),
        ([1, 2], [1, 2], 1.5, 1, None, None, 0, "amplitude", ValueError),
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
    weight,
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
            weight,
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
        "weight",
        "expected",
    ),
    [
        (
            [1, 2],
            [1, 2],
            1.5,
            None,
            None,
            None,
            0,
            "amplitude",
            [0, 1, 0],
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
    weight,
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
        weight,
    )
