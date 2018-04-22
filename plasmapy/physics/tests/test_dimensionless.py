from plasmapy.physics.dimensionless import (beta)

import astropy.units as u

B = 1.0 * u.T
n_i = 5e19 * u.m ** -3
T_e = 1e6 * u.K


def test_beta():
    # Check that beta is dimensionless
    float(beta(T_e, n_i, B))
