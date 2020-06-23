import astropy.units as u
import numpy as np
import pytest
from astropy.constants import c

from plasmapy.formulary.ionization import Z_bal, Saha


def test_Z_bal():
    n = 1e19 * u.m ** -3
    T_e = 5000 * u.K

    assert (Z_bal(n, T_e) * u.dimensionless_unscaled).unit == u.dimensionless_unscaled

    with pytest.warns(u.UnitsWarning):
        Z_bal(1e19, T_e)

    with pytest.raises(u.UnitTypeError):
        Z_bal(2e19 * u.kg, T_e)


def test_Saha():
    g_j = 2
    g_k = 1
    n_e = 1e19 * u.m ** -3
    E_jk = 10 * u.Ry
    T_e = 5000 * u.K

    assert (
        Saha(g_j, g_k, n_e, E_jk, T_e) * u.dimensionless_unscaled).unit == u.dimensionless_unscaled

    with pytest.warns(u.UnitsWarning):
        Saha(g_j, g_k, 1e19, E_jk, T_e)

    with pytest.raises(u.UnitTypeError):
        Saha(g_j, g_k, 1e19 * u.kg, E_jk, T_e)
