"""Tests for methods relating to quantities."""

import numpy as np
from astropy import units as u
import pytest

from ...utils.exceptions import RelativityWarning, RelativityError
from ...constants import c
from ..checks import (
    _check_quantity, _check_relativistic, check_relativistic,
    check_quantity
)


# (value, units, error)
quantity_error_examples_default = [
    # exceptions associated with the units keyword
    (5 * u.T, 5 * u.T, TypeError),
    (5 * u.T, 5, TypeError),
    (5 * u.T, [u.T, 1], TypeError),
    (5 * u.T, [1, u.m], TypeError),
    (u.T, u.J, TypeError),
    (3 * u.m / u.s, u.m, u.UnitConversionError),
    (5j * u.K, u.K, ValueError),
]

# (value, units, can_be_negative, can_be_complex, can_be_inf, error)
quantity_error_examples_non_default = [
    (-5 * u.K, u.K, False, False, True, ValueError),
    (np.inf * u.K, u.K, True, False, False, ValueError)
]

# (value, units)
quantity_valid_examples_default = [
    # check basic functionality
    (5 * u.T, u.T),
    (3 * u.m / u.s, u.m / u.s),
    (3 * u.m / u.s, [u.m / u.s]),
    (3 * u.m / u.s ** 2, [u.m / u.s, u.m / (u.s ** 2)]),
    (3 * u.km / u.yr, u.m / u.s),
    # check temperature in units of energy per particle (e.g., eV)
    (5 * u.eV, u.K),
    (5 * u.K, u.eV),
    (5 * u.keV, [u.m, u.K]),
    # check keywords relating to numerical values
    (np.inf * u.T, u.T)
]

# (value, units, can_be_negative, can_be_complex, can_be_inf)
quantity_valid_examples_non_default = [
    (5j * u.m, u.m, True, True, True)
]


# Tests for _check_quantity
@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf, error",
    quantity_error_examples_non_default)
def test__check_quantity_errors_non_default(
        value, units, can_be_negative, can_be_complex, can_be_inf, error):
    with pytest.raises(error):
        _check_quantity(value, 'arg', 'funcname', units,
                        can_be_negative=can_be_negative,
                        can_be_complex=can_be_complex,
                        can_be_inf=can_be_inf)


def test__check_quantity_warns_on_casting():
    with pytest.warns(u.UnitsWarning):
        _check_quantity(5, 'arg', 'funcname', u.m,)

@pytest.mark.parametrize(
    "value, units, error", quantity_error_examples_default)
def test__check_quantity_errors_default(value, units, error):
    with pytest.raises(error):
        _check_quantity(value, 'arg', 'funcname', units)


@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf",
    quantity_valid_examples_non_default)
def test__check_quantity_non_default(
        value, units, can_be_negative, can_be_complex, can_be_inf):
    _check_quantity(value, 'arg', 'funcname', units,
                    can_be_negative=can_be_negative,
                    can_be_complex=can_be_complex,
                    can_be_inf=can_be_inf)


@pytest.mark.parametrize("value, units", quantity_valid_examples_default)
def test__check_quantity_default(value, units):
    _check_quantity(value, 'arg', 'funcname', units)


# Tests for check_quantity decorator
@pytest.mark.parametrize(
    "value, units, error", quantity_error_examples_default)
def test_check_quantity_decorator_errors_default(value, units, error):

    @check_quantity({
        "x": {"units": units}
    })
    def func(x):
        return x

    with pytest.raises(error):
        func(value)


@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf, error",
    quantity_error_examples_non_default)
def test_check_quantity_decorator_errors_non_default(
        value, units, can_be_negative, can_be_complex, can_be_inf, error):

    @check_quantity({
        "x": {"units": units, "can_be_negative": can_be_negative,
              "can_be_complex": can_be_complex, "can_be_inf": can_be_inf}
    })
    def func(x):
        return x

    with pytest.raises(error):
        func(value)


@pytest.mark.parametrize("value, units", quantity_valid_examples_default)
def test_check_quantity_decorator_default(value, units):

    @check_quantity({
        "x": {"units": units}
    })
    def func(x):
        return x

    func(value)


@pytest.mark.parametrize(
    "value, units, can_be_negative, can_be_complex, can_be_inf",
    quantity_valid_examples_non_default)
def test_check_quantity_decorator_non_default(
        value, units, can_be_negative, can_be_complex, can_be_inf):

    @check_quantity({
        "x": {"units": units, "can_be_negative": can_be_negative,
              "can_be_complex": can_be_complex, "can_be_inf": can_be_inf}
    })
    def func(x):
        return x

    func(value)


def test_check_quantity_decorator_missing_validated_params():

    @check_quantity({
        "x": {"units": u.m},
        "y": {"units": u.s}
    })
    def func(x):
        return x

    with pytest.raises(TypeError) as e:
        func(1 * u.m)

    assert "Call to func is missing validated params y" == str(e.value)


def test_check_quantity_decorator_two_args_default():

    @check_quantity({
        "x": {"units": u.m},
        "y": {"units": u.s}
    })
    def func(x, y):
        return x / y

    func(1 * u.m, 1 * u.s)


def test_check_quantity_decorator_two_args_not_default():

    @check_quantity({
        "x": {"units": u.m, "can_be_negative": False},
        "y": {"units": u.s}
    })
    def func(x, y):
        return x / y

    with pytest.raises(ValueError):
        func(-1 * u.m, 2 * u.s)


def test_check_quantity_decorator_two_args_one_kwargs_default():

    @check_quantity({
        "x": {"units": u.m},
        "y": {"units": u.s},
        "z": {"units": u.eV}
    })
    def func(x, y, another, z=10 * u.eV):
        return x * y * z

    func(1 * u.m, 1 * u.s, 10 * u.T)

def test_check_quantity_decorator_two_args_one_kwargs_not_default():

    @check_quantity({
        "x": {"units": u.m},
        "y": {"units": u.s, "can_be_negative": False},
        "z": {"units": u.eV, "can_be_inf": False}
    })
    def func(x, y, z=10 * u.eV):
        return x * y * z

    with pytest.raises(ValueError):
        func(1 * u.m, 1 * u.s, z=np.inf * u.eV)


class Test_check_quantity_none_shall_pass:
    @check_quantity({
        "x": {"units": u.m, "none_shall_pass": True},
    })
    def func(self, x = None):
        if x is None:
            return 0 * u.m
        return x

    def test_incorrect_units(self):
        with pytest.raises(u.UnitConversionError):
            self.func(1*u.s)

    def test_none_to_zero(self):
        assert self.func(None) == 0*u.m

# (speed, betafrac)
non_relativistic_speed_examples = [
    (0 * u.m / u.s, 0.1),
    (0.0099999 * c, 0.1),
    (-0.009 * c, 0.1),
    (5 * u.AA / u.Gyr, 0.1)
]

# (speed, betafrac, error)
relativistic_error_examples = [
    (u.m / u.s, 0.1, TypeError),
    (51513.35, 0.1, TypeError),
    (5 * u.m, 0.1, u.UnitConversionError),
    (np.nan * u.m / u.s, 0.1, ValueError),
    (1.0 * c, 0.1, RelativityError),
    (1.1 * c, 0.1, RelativityError),
    (np.inf * u.cm / u.s, 0.1, RelativityError),
    (-1.0 * c, 0.1, RelativityError),
    (-1.1 * c, 0.1, RelativityError),
    (-np.inf * u.cm / u.s, 0.1, RelativityError),
]

# (speed, betafrac, warning)
relativistic_warning_examples = [
    (0.11 * c, 0.1),
    (-0.11 * c, 0.1),
    (2997924581 * u.cm / u.s, 0.1),
    (0.02 * c, 0.01),
]


# Tests for _check_relativistic
@pytest.mark.parametrize("speed, betafrac", non_relativistic_speed_examples)
def test__check_relativisitc_valid(speed, betafrac):
    _check_relativistic(speed, 'f', betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac, error", relativistic_error_examples)
def test__check_relativistic_errors(speed, betafrac, error):
    with pytest.raises(error):
        _check_relativistic(speed, 'f', betafrac=betafrac)


@pytest.mark.parametrize("speed, betafrac",
                         relativistic_warning_examples)
def test__check_relativistic_warnings(speed, betafrac):
    with pytest.warns(RelativityWarning):
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


@pytest.mark.parametrize(
    "speed",
    [item[0] for item in non_relativistic_speed_examples])
def test_check_relativistic_decorator_no_args_parentheses(speed):

    @check_relativistic()
    def speed_func():
        return speed

    speed_func()


@pytest.mark.parametrize("speed, betafrac, error", relativistic_error_examples)
def test_check_relativistic_decorator_errors(speed, betafrac, error):

    @check_relativistic(betafrac=betafrac)
    def speed_func():
        return speed

    with pytest.raises(error):
        speed_func()
