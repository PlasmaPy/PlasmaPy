"""Tests for methods relating to quantities."""

import numpy as np
from astropy import units as u
import pytest

from ...constants import c
from ..checks import (check_quantity, check_relativistic)


def test__check_quantity():
    def run_test(arg_value, units, can_be_negative=True, can_be_complex=False, can_be_inf=True):
        """ Tests the check_quantity decorator
        Tests that the given value is an astropy Quantity with correct units
        and valid numerical values
        Builds a method, decorates it, and immediately calls it
        """
        @check_quantity({
            'arg': {'units': units, 'can_be_negative': can_be_negative, 'can_be_complex': can_be_complex, 'can_be_inf': can_be_inf}
        })
        def test_method(arg): pass
        test_method(arg_value)

    # exceptions associated with the units keyword
    with pytest.raises(TypeError):
        run_test(5*u.T, 5*u.T)

    with pytest.raises(TypeError):
        run_test(5*u.T, 5)

    with pytest.raises(TypeError):
        run_test(5*u.T, [u.T, 1])

    with pytest.raises(TypeError):
        run_test(5*u.T, [1, u.m])

    with pytest.raises(TypeError):
        run_test(u.T, u.J)

    with pytest.raises(UserWarning):
        run_test(5.0, u.m)

    with pytest.raises(u.UnitConversionError):
        run_test(3*u.m/u.s, u.m)

    # check basic functionality

    run_test(5*u.T, u.T)
    run_test(3*u.m/u.s, u.m/u.s)
    run_test(3*u.m/u.s, [u.m/u.s])
    run_test(3*u.m/u.s**2, [u.m/u.s, u.m/(u.s**2)])
    run_test(3*u.km/u.yr, u.m/u.s)

    # check temperature in units of energy per particle (e.g., eV)

    run_test(5*u.eV, u.K)
    run_test(5*u.K, u.eV)
    run_test(5*u.keV, [u.m, u.K])

    # check keywords relating to numerical values

    run_test(5j*u.m, u.m, can_be_complex=True)
    run_test(np.inf*u.T, u.T)

    with pytest.raises(ValueError):
        run_test(-5*u.K, u.K, can_be_negative=False)

    with pytest.raises(ValueError):
        run_test(5j*u.K, u.K)

    with pytest.raises(ValueError):
        run_test(np.inf*u.K, u.K, can_be_inf=False)


def test__check_relativistic():

    def run_test(value, betafrac=0.1):
        @check_relativistic(betafrac=betafrac)
        def test_method():
            return value
        test_method()

    run_test(0*u.m/u.s)
    run_test(0.099999*c)
    run_test(-0.09*c)
    run_test(5*u.AA/u.Gyr)

    with pytest.raises(UserWarning):
        run_test(0.11*c)

    with pytest.raises(UserWarning):
        run_test(1.0*c)

    with pytest.raises(UserWarning):
        run_test(1.1*c)

    with pytest.raises(UserWarning):
        run_test(np.inf*u.cm/u.s)

    with pytest.raises(UserWarning):
        run_test(-0.11*c)

    with pytest.raises(UserWarning):
        run_test(-1.0*c)

    with pytest.raises(UserWarning):
        run_test(-1.1*c)

    with pytest.raises(UserWarning):
        run_test(-np.inf*u.cm/u.s)

    with pytest.raises(UserWarning):
        run_test(2997924581*u.cm/u.s)

    with pytest.raises(UserWarning):
        run_test(0.02*c, betafrac=0.01)

    with pytest.raises(TypeError):
        run_test(u.m/u.s)

    with pytest.raises(TypeError):
        run_test(51513.35)

    with pytest.raises(u.UnitConversionError):
        run_test(5*u.m)

    with pytest.raises(ValueError):
        run_test(np.nan*u.m/u.s)
