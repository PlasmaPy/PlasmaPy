import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.ionization import ionization_balance, Saha, Z_bal_


def test_ionization_balance():
    n = 1e19 * u.m**-3
    T_e = 5000 * u.K

    val = ionization_balance(n, T_e)
    assert val.unit == u.dimensionless_unscaled
    assert np.isclose(val.value, 0.27478417234035346)

    # aliases must be the same function
    assert Z_bal_ is ionization_balance

    with pytest.warns(u.UnitsWarning):
        ionization_balance(1e19, T_e)

    with pytest.raises(u.UnitTypeError):
        ionization_balance(2e19 * u.kg, T_e)


def test_Saha():
    g_j = 2
    g_k = 1
    n_e = 1e19 * u.m**-3
    E_jk = 10 * u.Ry
    T_e = 5000 * u.K

    val = Saha(g_j, g_k, n_e, E_jk, T_e)
    assert val.unit == u.dimensionless_unscaled
    assert np.isclose(val.value, 2.4775608756800435e-129)

    with pytest.warns(u.UnitsWarning):
        Saha(g_j, g_k, 1e19, E_jk, T_e)

    with pytest.raises(u.UnitTypeError):
        Saha(g_j, g_k, 1e19 * u.kg, E_jk, T_e)
