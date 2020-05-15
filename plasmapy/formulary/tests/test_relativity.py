"""Tests for functions in relativity.py."""

import pytest
import numpy as np
from astropy import units as u
from astropy.constants import c
from plasmapy.formulary.relativity import Lorentz_factor, relativistic_energy
from plasmapy.utils.exceptions import RelativityError


def test_Lorentz_factor():
    r"""Test Lorentz_factor in relativity.py"""

    V = 123456789 * u.m / u.s
    assert np.isclose(Lorentz_factor(V), (1 / np.sqrt(1 - V ** 2 / c ** 2)).value)
    assert Lorentz_factor(-V) == Lorentz_factor(V)

    assert np.isclose(Lorentz_factor(0 * u.m / u.s), 1.0)
    assert Lorentz_factor(c) == np.inf

    V_arr = np.array([987532.0, 299792458]) * u.m / u.s
    gamma_arr = Lorentz_factor(V_arr)
    assert np.isclose(gamma_arr[0], (1 / np.sqrt(1 - V_arr[0] ** 2 / c ** 2)).value)
    assert gamma_arr[1] == np.inf

    assert (
        Lorentz_factor(3 * u.m / u.s) * u.dimensionless_unscaled
    ).unit == u.dimensionless_unscaled

    with pytest.raises(RelativityError):
        Lorentz_factor(1.0000000001 * c)

    with pytest.raises(ValueError), pytest.warns(u.UnitsWarning):
        Lorentz_factor(299792459)

    with pytest.warns(u.UnitsWarning):
        Lorentz_factor(2.2)

    with pytest.raises(u.UnitTypeError):
        Lorentz_factor(4 * u.kg)


def test_relativistic_energy():
    r"""Test relativistic_energy in relativity.py"""

    V = 123456789 * u.m / u.s
    m = 1 * u.kg
    assert np.isclose(relativistic_energy(V, m).value, ((1 / np.sqrt(1 - V ** 2 / c ** 2)) * m * c ** 2).value)
    assert relativistic_energy(-V, m) == relativistic_energy(V, m)

    assert np.isclose(relativistic_energy(0 * u.m / u.s, m).value, (m * c ** 2).value)
    assert relativistic_energy(c, m) == np.inf

    V_arr = np.array([987532.0, 299792458]) * u.m / u.s
    gamma_arr = Lorentz_factor(V_arr)
    assert np.isclose(gamma_arr[0], (1 / np.sqrt(1 - V_arr[0] ** 2 / c ** 2)).value)
    assert gamma_arr[1] == np.inf

    assert (
                   relativistic_energy(3 * u.m / u.s, 2 * u.kg)
           ).unit == u.J

    with pytest.raises(RelativityError):
        relativistic_energy(1.0000000001 * c, m)

    with pytest.raises(ValueError), pytest.warns(u.UnitsWarning):
        relativistic_energy(299792459, 1)

    with pytest.warns(u.UnitsWarning):
        relativistic_energy(2.2, m)

    with pytest.raises(u.UnitTypeError):
        relativistic_energy(4 * u.kg, m)
