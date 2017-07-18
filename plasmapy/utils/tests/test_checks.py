"""Tests for methods relating to quantities."""

import numpy as np
from astropy import units as u
import pytest

from ...constants import c
from ..checks import (_check_quantity, _check_relativistic)


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


def test__check_relativistic():

    _check_relativistic(0*u.m/u.s, 'f')
    _check_relativistic(0.099999*c, 'f')
    _check_relativistic(-0.09*c, 'f')
    _check_relativistic(5*u.AA/u.Gyr, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(0.11*c, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(1.0*c, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(1.1*c, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(np.inf*u.cm/u.s, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(-0.11*c, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(-1.0*c, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(-1.1*c, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(-np.inf*u.cm/u.s, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(2997924581*u.cm/u.s, 'f')

    with pytest.raises(UserWarning):
        _check_relativistic(0.02*c, 'f', betafrac=0.01)

    with pytest.raises(TypeError):
        _check_relativistic(u.m/u.s, 'f')

    with pytest.raises(TypeError):
        _check_relativistic(51513.35, 'f')

    with pytest.raises(u.UnitConversionError):
        _check_relativistic(5*u.m, 'f')

    with pytest.raises(ValueError):
        _check_relativistic(np.nan*u.m/u.s, 'f')
