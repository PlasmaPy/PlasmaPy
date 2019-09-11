import astropy.units as u
import numpy as np

from plasmapy.utils.decorators.converter import angular_freq_to_hz

def test_conversion():
    @angular_freq_to_hz
    def test_function():
        return 2 * np.pi * u.rad / u.s

    assert test_function().unit == (u.rad / u.s)
    assert test_function().value == 2 * np.pi
    assert test_function(to_hz=True).unit == u.Hz
    assert test_function(to_hz=True).value == 1

