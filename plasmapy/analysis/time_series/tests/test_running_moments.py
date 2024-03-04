"""Tests for running_moments.py"""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.analysis.time_series.running_moments import running_mean, running_moment


@pytest.mark.parametrize(
    ("signal", "radius", "expected"),
    [
        ([1, 2, 3], 1, [2]),
        ([1, 2, 3, 4], 1, [2, 3]),
        ([1, 2, 3, 4] * u.eV, 1, [2, 3] * u.eV),
        (np.array([-0.5, -0.0, 0.5, 1.0, 1.5]), 2, [0.5]),
    ],
)
def test_running_mean(signal, radius: int, expected) -> None:
    """Test running_mean function"""
    result = running_mean(signal=signal, radius=radius)
    assert np.allclose(result, expected)


@pytest.mark.parametrize(("signal", "radius"), [([1, 2], 1), ([1, 2, 3, 4], 1.2)])
def test_running_mean_exception(signal, radius: int) -> None:
    """Test whether exception is risen"""
    with pytest.raises((ValueError, TypeError)):
        running_mean(signal, radius)


@pytest.mark.parametrize(
    ("signal", "radius", "moment", "time", "expected"),
    [
        ([1, 2, 3], 1, 1, [1, 2, 3], ([2], [2])),
        ([1, 2, 3] * u.eV, 1, 1, [1, 2, 3], ([2] * u.eV, [2])),
        ([1, 2, 3, 2, 1], 1, 2, [1, 2, 3, 4, 5], ([2 / 27**0.5], [3])),
        ([1, 2, 3, 2, 1] * u.eV, 1, 2, [1, 2, 3, 4, 5], ([2 / 27**0.5] * u.eV, [3])),
        (
            [1, 2, 3, 2, 1] * u.eV,
            1,
            3,
            [1, 2, 3, 4, 5],
            ([1.7320508075688772] * u.eV, [3]),
        ),
        ([1, 2, 3, 2, 1] * u.eV, 1, 4, [1, 2, 3, 4, 5], ([3] * u.eV, [3])),
    ],
)
def test_running_moment(signal, radius: int, moment, time, expected) -> None:
    """Test running_moment_function"""
    result = running_moment(signal, radius, moment, time)
    assert np.allclose(result, expected)


@pytest.mark.parametrize(
    ("signal", "radius", "moment", "time"),
    [
        ([1, 2, 3, 4, 5], 1, 0, None),
        ([1, 2, 3, 4, 5], 1, 6, None),
        ([1, 2, 3, 4], 1, 2, None),
        ([1, 2, 3, 4, 5], 1, 2, [1, 2, 3, 4]),
    ],
)
def test_running_moment_exception(signal, radius: int, moment, time) -> None:
    """Test whether exception is risen"""
    with pytest.raises(ValueError):
        running_moment(signal, radius, moment, time)
