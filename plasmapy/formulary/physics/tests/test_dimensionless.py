from plasmapy.formulary.physics.dimensionless import (beta, quantum_theta)

import astropy.units as u
import numpy as np

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
    n = np.array([1, 1]) * u.cm**-3
    T = np.array([1, 1]) * u.K
    out = beta(T, n, B)
    assert np.isnan(out[1])
    assert out[1].unit == u.dimensionless_unscaled
