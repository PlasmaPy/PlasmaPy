import pytest

from plasmapy.formulary.dimensionless import beta, quantum_theta, Reynolds_number

import astropy.units as u
import numpy as np

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
    v = 0.1 * u.m / u.s
    L = 0.05 * u.m
    mu = 10 * u.kg / (u.m * u.s)

    assert (
        Reynolds_number(rho, v, L, mu) * u.dimensionless_unscaled
           ).unit == u.dimensionless_unscaled

    with pytest.raises(RelativityError):
        Reynolds_number(rho, 1.0000000001 * c, L, mu)

    with pytest.raises(RelativityError), pytest.warns(u.UnitsWarning):
        Reynolds_number(1, 299792459, 0.05, 10)

    with pytest.warns(u.UnitsWarning):
        Reynolds_number(rho, 2.2, L, mu)

    with pytest.raises(u.UnitTypeError):
        Reynolds_number(rho, 4 * u.kg, L, mu)
