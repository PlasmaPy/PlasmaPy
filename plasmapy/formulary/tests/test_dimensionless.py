import astropy.units as u
import numpy as np
import pytest
from astropy.constants import c

from plasmapy.formulary.dimensionless import (
    Mag_Reynolds,
    Reynolds_number,
    beta,
    quantum_theta,
    Re_,
    Rm_,
)
from plasmapy.utils import RelativityError

B = 1.0 * u.T
n = 5e19 * u.m ** -3
T = 1e6 * u.K


def test_beta_dimensionless():
    # Check that beta is dimensionless
    float(beta(T, n, B))


def test_quantum_theta_dimensionless():
    # Check that quantum theta is dimensionless
    float(quantum_theta(T, n))


def test_beta_nan():
    # Check that nans are passed through properly
    B = np.array([1, np.nan]) * u.T
    n = np.array([1, 1]) * u.cm ** -3
    T = np.array([1, 1]) * u.K
    out = beta(T, n, B)
    assert np.isnan(out[1])
    assert out[1].unit == u.dimensionless_unscaled


def test_Reynolds_number():
    r"""Test Reynolds_number in dimensionless.py"""
    rho = 1490 * u.kg / u.m ** 3
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


def test_dimensionless_aliases():
    r"""Test all aliases defined in dimensionless.py"""

    assert Re_ is Reynolds_number
    assert Rm_ is Mag_Reynolds
