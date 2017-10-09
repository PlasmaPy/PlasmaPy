"""Tests for methods relating to quantities."""

import numpy as np
from astropy import units as u
import pytest

from ...constants import c
from ..checks import (_check_relativistic, check_relativistic, )


# (value, units, error)
quantity_error_examples_default = [
    # exceptions associated with the units keyword
    (5*u.T, 5*u.T, TypeError),
    (5*u.T, 5, TypeError),
    (5*u.T, [u.T, 1], TypeError),
    (5*u.T, [1, u.m], TypeError),
    (u.T, u.J, TypeError),
    (5.0, u.m, UserWarning),
    (3*u.m/u.s, u.m, u.UnitConversionError),
    (5j*u.K, u.K, ValueError),
]

# (value, units, can_be_negative, can_be_complex, can_be_inf, error)
quantity_error_examples_non_default = [
    (-5*u.K, u.K, False, False, True, ValueError),
    (np.inf*u.K, u.K, True, False, False, ValueError)
]

# (value, units)
quantity_valid_examples_default = [
    # check basic functionality
    (5*u.T, u.T),
    (3*u.m/u.s, u.m/u.s),
    (3*u.m/u.s, [u.m/u.s]),
    (3*u.m/u.s**2, [u.m/u.s, u.m/(u.s**2)]),
    (3*u.km/u.yr, u.m/u.s),
    # check temperature in units of energy per particle (e.g., eV)
    (5*u.eV, u.K),
    (5*u.K, u.eV),
    (5*u.keV, [u.m, u.K]),
    # check keywords relating to numerical values
    (np.inf*u.T, u.T)
]

# (value, units, can_be_negative, can_be_complex, can_be_inf)
quantity_valid_examples_non_default = [
    (5j*u.m, u.m, True, True, True)
]

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


# Tests for _check_relativistic
@pytest.mark.parametrize("speed, betafrac", non_relativistic_speed_examples)
def test__check_relativisitc_valid(speed, betafrac):
    _check_relativistic(speed, 'f', betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac, error", relativisitc_error_examples)
def test__check_relativistic_errors(speed, betafrac, error):
    with pytest.raises(error):
        _check_relativistic(speed, 'f', betafrac=betafrac)


# Tests for check_relativistic decorator
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


@pytest.mark.parametrize("speed",
    [item[0] for item in non_relativistic_speed_examples])
def test_check_relativistic_decorator_no_args_parentheses(speed):

    @check_relativistic()
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
