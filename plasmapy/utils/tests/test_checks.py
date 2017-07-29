"""Tests for methods relating to quantities."""

import numpy as np
from astropy import units as u
import pytest

from ...constants import c
from ..checks import (
    _check_quantity, _check_relativistic, check_relativistic
)


# (value, units, can_be_negative, can_be_complex, can_be_inf, error)
quantity_error_examples = [
    # exceptions associated with the units keyword
    (5*u.T, 5*u.T, True, False, True, TypeError),
    (5*u.T, 5, True, False, True, TypeError),
    (5*u.T, [u.T, 1], True, False, True, TypeError),
    (5*u.T, [1, u.m], True, False, True, TypeError),
    (u.T, u.J, True, False, True, TypeError),
    (5.0, u.m, True, False, True, UserWarning),
    (3*u.m/u.s, u.m, True, False, True, u.UnitConversionError),
    (-5*u.K, u.K, False, False, True, ValueError),
    (5j*u.K, u.K, True, False, True, ValueError),
    (np.inf*u.K, u.K, True, False, False, ValueError)
]

# (value, units, can_be_negative, can_be_complex, can_be_inf)
quantity_valid_examples = [
    # check basic functionality
    (5*u.T, u.T, True, False, True),
    (3*u.m/u.s, u.m/u.s, True, False, True),
    (3*u.m/u.s, [u.m/u.s], True, False, True),
    (3*u.m/u.s**2, [u.m/u.s, u.m/(u.s**2)], True, False, True),
    (3*u.km/u.yr, u.m/u.s, True, False, True),
    # check temperature in units of energy per particle (e.g., eV)
    (5*u.eV, u.K, True, False, True),
    (5*u.K, u.eV, True, False, True),
    (5*u.keV, [u.m, u.K], True, False, True),
    # check keywords relating to numerical values
    (5j*u.m, u.m, True, True, True),
    (np.inf*u.T, u.T, True, False, True)
]


@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf, error",
    quantity_error_examples)
def test__check_quantity_errors(
        value, units, can_be_negative, can_be_complex, can_be_inf, error):
    with pytest.raises(error):
        _check_quantity(value, 'arg', 'funcname', units,
                        can_be_negative=can_be_negative,
                        can_be_complex=can_be_complex,
                        can_be_inf=can_be_inf)


@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf",
    quantity_valid_examples)
def test__check_quantity(value, units, can_be_negative, can_be_complex, can_be_inf):

    _check_quantity(value, 'arg', 'funcname', units,
                    can_be_negative=can_be_negative,
                    can_be_complex=can_be_complex,
                    can_be_inf=can_be_inf)

# (speed, betafrac)
non_relativistic_speed_examples = [
    (0*u.m/u.s, 0.1),
    (0.099999*c, 0.1),
    (-0.09*c, 0.1),
    (5*u.AA/u.Gyr, 0.1)
]

# (speed, betafrac, error)
relativisitc_error_examples = [
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


@pytest.mark.parametrize("speed, betafrac", non_relativistic_speed_examples)
def test__check_relativisitc_valid(speed, betafrac):
    _check_relativistic(speed, 'f', betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac, error", relativisitc_error_examples)
def test__check_relativistic_errors(speed, betafrac, error):
    with pytest.raises(error):
        _check_relativistic(speed, 'f', betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac", non_relativistic_speed_examples)
def test_check_relativistic_decorator(speed, betafrac):

    @check_relativistic(betafrac=betafrac)
    def speed_func():
        return speed

    speed_func()


@pytest.mark.parametrize(
    "speed",
    [item[0] for item in non_relativistic_speed_examples])
def test_check_relativistic_decorator_no_args(speed):

    @check_relativistic
    def speed_func():
        return speed

    speed_func()


@pytest.mark.parametrize("speed, betafrac, error", relativisitc_error_examples)
def test_check_relativistic_decorator_errors(speed, betafrac, error):

    @check_relativistic(betafrac=betafrac)
    def speed_func():
        return speed

    with pytest.raises(error):
        speed_func()
