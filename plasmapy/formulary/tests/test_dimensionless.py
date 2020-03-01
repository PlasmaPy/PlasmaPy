from plasmapy.formulary.dimensionless import (
     beta,
     quantum_theta,
     prandtl_number,
     magnetic_prandtl_number
)

import astropy.units as u
import numpy as np

B = 1.0 * u.T
n = 5e19 * u.m ** -3
T = 1e6 * u.K
v = 1e4 * (u.m ** 2 / u.s)

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


def test_prandtl_number():
    assert1 = prandtl_number(v, 1e4 * (u.m ** 2 / u.s))
    error1 = f"""Prandtl number should be 1 when the momentum diffusivity(v) and
              thermal diffusivity are equal(alpha) instead of {assert1}"""
    assert2 = prandtl_number(v, 1e6 * (u.m ** 2 / u.s))
    error2 = f"""Prandtl number should be 2 when the momentum diffusivity(v)
             and thermal diffusivity(alpha) are 10**4 and 10**6 instead of
             {assert2}"""
    assert assert1 ==  1.0 , error1
    assert assert2 == 0.01 , error2


def test_magnetic_prandtl_number():
    assert1 = magnetic_prandtl_number(v, 1e4 * (u.m ** 2 / u.s))
    error1 = f"""Magnetic prandtl number should be 1 when the momentum
             diffusivity(v) and magnetic diffusivity(n) are equal instead of
             {assert1}"""
    assert2 = magnetic_prandtl_number(v, 1e6 * (u.m ** 2 / u.s))
    error2 = f"""Magnetic prandtl number should be 2 when the momentum
             diffusivity(v) and magnetic diffusivity(n) are 10**4 and 10**6
             instead of {assert2}"""
    assert assert1 ==  1.0 , error1
    assert assert2 == 0.01 , error2
