"""Tests for functionality contained in `plasmapy.formulary.dimensionless`."""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.dimensionless import (
    beta,
    betaH_,
    Debye_number,
    Hall_parameter,
    Mag_Reynolds,
    nD_,
    quantum_theta,
    Re_,
    Reynolds_number,
    Rm_,
)
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

Z = 1

B = 1.0 * u.T

n = 5e19 * u.m**-3
n_e = Z * 5e19 * u.m**-3

T = 1e6 * u.K
T_e = 1e6 * u.K


@pytest.mark.parametrize(
    "alias, parent",
    [
        (nD_, Debye_number),
        (betaH_, Hall_parameter),
        (Re_, Reynolds_number),
        (Rm_, Mag_Reynolds),
    ],
)
def test_aliases(alias, parent):
    """Test all aliases defined in dimensionless.py"""
    assert alias is parent


def test_beta_dimensionless():
    # Check that beta is dimensionless
    float(beta(T, n, B))


def test_quantum_theta_dimensionless():
    # Check that quantum theta is dimensionless
    float(quantum_theta(T, n))


def test_beta_nan():
    # Check that nans are passed through properly
    B = np.array([1, np.nan]) * u.T
    n = np.array([1, 1]) * u.cm**-3
    T = np.array([1, 1]) * u.K
    out = beta(T, n, B)
    assert np.isnan(out[1])
    assert out[1].unit == u.dimensionless_unscaled


def test_Reynolds_number():
    r"""Test Reynolds_number in dimensionless.py"""
    rho = 1490 * u.kg / u.m**3
    U = 0.1 * u.m / u.s
    L = 0.05 * u.m
    mu = 10 * u.kg / (u.m * u.s)

    assert (
        Reynolds_number(rho, U, L, mu) * u.dimensionless_unscaled
    ).unit == u.dimensionless_unscaled

    with pytest.warns(u.UnitsWarning):
        Reynolds_number(rho, 2.2, L, mu)

    with pytest.raises(u.UnitTypeError):
        Reynolds_number(rho, 4 * u.kg, L, mu)


def test_Mag_Reynolds():
    r"""Test Mag_Reynolds in dimensionless.py"""

    sigma = 1e8 * u.S / u.m
    U = 0.1 * u.m / u.s
    L = 0.05 * u.m

    assert (
        Mag_Reynolds(U, L, sigma) * u.dimensionless_unscaled
    ).unit == u.dimensionless_unscaled

    with pytest.warns(u.UnitsWarning):
        Mag_Reynolds(2.2, L, sigma)

    with pytest.raises(u.UnitTypeError):
        Mag_Reynolds(2.2 * u.kg, L, sigma)


def test_Debye_number():
    r"""Test the Debye_number function in dimensionless.py."""

    assert Debye_number(T_e, n_e).unit.is_equivalent(u.dimensionless_unscaled)

    T_e_eV = T_e.to(u.eV, equivalencies=u.temperature_energy())
    assert np.isclose(Debye_number(T_e, n_e).value, Debye_number(T_e_eV, n_e).value)

    assert np.isclose(Debye_number(1 * u.eV, 1 * u.cm**-3).value, 1720862385.43342)

    with pytest.warns(u.UnitsWarning):
        Debye_number(T_e, 4)

    with pytest.raises(ValueError):
        Debye_number(None, n_e)

    with pytest.raises(u.UnitTypeError):
        Debye_number(5 * u.m, 5 * u.m**-3)

    with pytest.raises(u.UnitTypeError):
        Debye_number(5 * u.K, 5 * u.m**3)

    with pytest.raises(ValueError):
        Debye_number(5j * u.K, 5 * u.cm**-3)

    Tarr2 = np.array([1, 2]) * u.K
    narr3 = np.array([1, 2, 3]) * u.m**-3
    with pytest.raises(ValueError):
        Debye_number(Tarr2, narr3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_number(1.1, 1.1) == Debye_number(1.1 * u.K, 1.1 * u.m**-3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_number(1.1 * u.K, 1.1) == Debye_number(1.1, 1.1 * u.m**-3)

    assert_can_handle_nparray(Debye_number)
