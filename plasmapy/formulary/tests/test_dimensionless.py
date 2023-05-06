"""Tests for functionality contained in `plasmapy.formulary.dimensionless`."""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.dimensionless import (
    beta,
    betaH_,
    Debye_number,
    Hall_parameter,
    Lundquist_number,
    Mag_Reynolds,
    nD_,
    Re_,
    Reynolds_number,
    Rm_,
)
from plasmapy.particles import Particle
from plasmapy.utils import RelativityWarning
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


def test_Hall_parameter():
    r"""Test Hall_parameter in dimensionless.py"""

    ion = Particle("He-4 +1")
    particle = Particle("e-")

    assert Hall_parameter(n, T, B, ion, particle).unit.is_equivalent(
        u.dimensionless_unscaled
    )

    assert np.isclose(Hall_parameter(n, T, B, ion, particle).value, 70461.38821149625)

    with pytest.warns(u.UnitsWarning):
        Hall_parameter(n, T, 1.0, ion, particle)

    with pytest.raises(u.UnitTypeError):
        Hall_parameter(n, T, 1.0 * u.kg, ion, particle)

    with pytest.raises(TypeError):
        Hall_parameter(n, T, B, None, particle)

    with pytest.raises(ValueError):
        Hall_parameter(n, T, B, ion, particle, coulomb_log_method="test")

    with pytest.warns(u.UnitsWarning):
        Hall_parameter(n, T, B, ion, particle, V=100)

    with pytest.raises(TypeError):
        Hall_parameter(n, T, B, ion, particle, coulomb_log="test")

    with pytest.warns(RelativityWarning):
        Hall_parameter(1e10 * u.m**-3, 5.8e3 * u.eV, 2.3 * u.T, ion, particle)


def test_Lundquist_number():
    r"""Test the Lundquist_number function in dimensionless.py."""
    L = 0.05 * u.m
    rho = 1490 * u.kg / u.m**3
    sigma = 1e8 * u.S / u.m

    Lundquist_number(L, B, rho, sigma).unit.is_equivalent(u.dimensionless_unscaled)

    with pytest.warns(u.UnitsWarning):
        Lundquist_number(3.3, B, rho, sigma)

    with pytest.raises(u.UnitTypeError):
        Lundquist_number(3.3 * u.kg, B, rho, sigma)
