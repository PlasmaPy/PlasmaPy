"""Tests for methods relating to quantities."""

import numpy as np
from astropy import units as u
import pytest

from ...constants import c
from ..checks import (
    _check_quantity, _check_relativistic, check_relativistic
)


def test__check_quantity():

    # exceptions associated with the units keyword

    with pytest.raises(TypeError):
        _check_quantity(5*u.T, 'arg', 'funcname', 5*u.T)

    with pytest.raises(TypeError):
        _check_quantity(5*u.T, 'arg', 'funcname', 5)

    with pytest.raises(TypeError):
        _check_quantity(5*u.T, 'arg', 'funcname', [u.T, 1])

    with pytest.raises(TypeError):
        _check_quantity(5*u.T, 'arg', 'funcname', [1, u.m])

    with pytest.raises(TypeError):
        _check_quantity(u.T, 'arg', 'funcname', u.J)

    with pytest.raises(UserWarning):
        _check_quantity(5.0, 'arg', 'funcname', u.m)

    with pytest.raises(u.UnitConversionError):
        _check_quantity(3*u.m/u.s, 'V', 'f', u.m)

    # check basic functionality

    _check_quantity(5*u.T, 'B', 'f', u.T)
    _check_quantity(3*u.m/u.s, 'V', 'f', u.m/u.s)
    _check_quantity(3*u.m/u.s, 'V', 'f', [u.m/u.s])
    _check_quantity(3*u.m/u.s**2, 'a', 'f', [u.m/u.s, u.m/(u.s**2)])
    _check_quantity(3*u.km/u.yr, 'V', 'f', u.m/u.s)

    # check temperature in units of energy per particle (e.g., eV)

    _check_quantity(5*u.eV, 'T', 'f', u.K)
    _check_quantity(5*u.K, 'T', 'f', u.eV)
    _check_quantity(5*u.keV, 'T', 'f', [u.m, u.K])

    # check keywords relating to numerical values

    _check_quantity(5j*u.m, 'L', 'f', u.m, can_be_complex=True)
    _check_quantity(np.inf*u.T, 'B', 'f', u.T)

    with pytest.raises(ValueError):
        _check_quantity(-5*u.K, 'T', 'f', u.K, can_be_negative=False)

    with pytest.raises(ValueError):
        _check_quantity(5j*u.K, 'T', 'f', u.K)

    with pytest.raises(ValueError):
        _check_quantity(np.inf*u.K, 'T', 'f', u.K, can_be_inf=False)


non_relativistic_speeds = [
    (0*u.m/u.s, 0.1),
    (0.099999*c, 0.1),
    (-0.09*c, 0.1),
    (5*u.AA/u.Gyr, 0.1)
]
relativisitc_error_inputs = [
    (0.11*c, 0.1, UserWarning),
    (1.0*c, 0.1, UserWarning),
    (1.1*c, 0.1, UserWarning),
    (np.inf*u.cm/u.s, 0.1, UserWarning),
    (-0.11*c, 0.1, UserWarning),
    (-1.0*c, 0.1, UserWarning),
    (-1.1*c, 0.1, UserWarning),
    (-np.inf*u.cm/u.s, 0.1, UserWarning),
    (2997924581*u.cm/u.s, 0.1, UserWarning),
    (0.02*c, 0.01, UserWarning),
    (u.m/u.s, 0.1, TypeError),
    (51513.35, 0.1, TypeError),
    (5*u.m, 0.1, u.UnitConversionError),
    (np.nan*u.m/u.s, 0.1, ValueError)
]


@pytest.mark.parametrize("speed, betafrac", non_relativistic_speeds)
def test__check_relativisitc_valid(speed, betafrac):
    _check_relativistic(speed, 'f', betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac, error", relativisitc_error_inputs)
def test__check_relativistic_errors(speed, betafrac, error):
    with pytest.raises(error):
        _check_relativistic(speed, 'f', betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac", non_relativistic_speeds)
def test_check_relativistic_decorator(speed, betafrac):

    @check_relativistic(betafrac=betafrac)
    def speed_func():
        return speed

    speed_func()


@pytest.mark.parametrize(
    "speed",
    [item[0] for item in non_relativistic_speeds])
def test_check_relativistic_decorator_no_args(speed):

    @check_relativistic
    def speed_func():
        return speed

    speed_func()

@pytest.mark.parametrize("speed, betafrac, error", relativisitc_error_inputs)
def test_check_relativistic_decorator_errors(
    speed, betafrac, error):

    @check_relativistic(betafrac=betafrac)
    def speed_func():
        return speed

    with pytest.raises(error):
        speed_func()
