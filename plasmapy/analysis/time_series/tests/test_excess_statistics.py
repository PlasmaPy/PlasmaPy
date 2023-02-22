"""Tests for excess_statistics.py"""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.analysis.time_series.excess_statistics import excess_stat


@pytest.mark.parametrize(
    "signal, thresholds, time_step, pdf, bins, expected",
    [
        (
            np.array([0, 2, 0, 4, 0]),
            [1, 3, 5],
            1,
            False,
            32,
            ([2, 1, 0], [2, 1, 0], [1, 1, 0], [0, 0, 0]),
        ),
        (
            np.array([0, 2, 0, 4, 0]),
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
def test_excess_stat(signal, thresholds, time_step, pdf, bins, expected):
    """test excess_stat function"""
    results = excess_stat(
        signal=signal, thresholds=thresholds, time_step=time_step, pdf=pdf, bins=bins
    )
    assert np.allclose(results, expected)
