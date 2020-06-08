import astropy.units as u
import numpy as np

from plasmapy.utils.decorators.converter import angular_freq_to_hz


def test_conversion():
    @angular_freq_to_hz
    def test_function():
        return 2 * np.pi * u.rad / u.s

    assert test_function().unit == (
        u.rad / u.s
    ), f"Unit expected is {(u.rad / u.s)} instead of {test_function().unit}"
    assert (
        test_function().value == 2 * np.pi
    ), f"Value expected is {2 * np.pi} instead of {test_function().value}"
    assert (
        test_function(to_hz=True).unit == u.Hz
    ), f"Unit expected is {u.Hz} instead of {test_function(to_hz=True).unit}"
    assert (
        test_function(to_hz=True).value == 1
    ), f"Value expected is {1} instead of {test_function(to_hz=True).value}"
