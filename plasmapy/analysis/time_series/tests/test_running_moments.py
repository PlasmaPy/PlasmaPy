"""Tests for running_moments.py"""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.analysis.time_series.running_moments import running_mean


@pytest.mark.parametrize(
    "signal, radius, expected",
    [
        ([1, 2, 3], 1, [2]),
        ([1, 2, 3, 4], 1, [2, 3]),
        ([1, 2, 3, 4] * u.eV, 1, [2, 3] * u.eV),
        (np.array([-0.5, -0.0, 0.5, 1.0, 1.5]), 2, [0.5]),
    ],
)
def test_runnine_mean(signal, radius, expected):
    """test running_mean function"""
    result = running_mean(signal=signal, radius=radius)
    assert np.allclose(result, expected)
