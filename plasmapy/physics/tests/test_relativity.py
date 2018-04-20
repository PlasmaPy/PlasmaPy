"""Tests for functions in relativity.py."""

import pytest
import numpy as np
from astropy import units as u
from ...constants import c
from ..relativity import Lorentz_factor
from ...utils.exceptions import RelativityError


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

    assert (Lorentz_factor(3 * u.m / u.s) * u.dimensionless_unscaled).unit == \
        u.dimensionless_unscaled

    with pytest.raises(RelativityError):
        Lorentz_factor(1.0000000001 * c)

    with pytest.raises(ValueError), pytest.warns(u.UnitsWarning):
            Lorentz_factor(299792459)

    with pytest.warns(u.UnitsWarning):
        Lorentz_factor(2.2)

    with pytest.raises(u.UnitConversionError):
        Lorentz_factor(4 * u.kg)
