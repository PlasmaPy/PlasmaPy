"""Tests for functions that calculate plasma dielectric parameters in
dielectry.py"""

import numpy as np
from astropy import units as u

from ..dielectric import (cold_plasma_permittivity_LRP,
                          cold_plasma_permittivity_SDP)

B = 1.0 * u.T
n = [1e18/u.m**3]
omega = 55e6*u.rad/u.s

single_species = ['e']
two_species = ['e', 'D+']
three_species = ['e', 'D+', 'H+']


def test_cold_plasma_permittivity_SDP():
    r"""Test the cold plasma dielectric tensor (S, D, P) function"""

    # test with one species
    S, D, P = cold_plasma_permittivity_SDP(B, single_species, n, omega)
    assert np.isclose(S, 1.10288)
    assert np.isclose(D, 329)
    assert np.isclose(P, -1052100)

    # test with two species and same density
    n_2 = np.array([1, 1])*1e19/u.m**3
    S, D, P = cold_plasma_permittivity_SDP(B, two_species, n_2, omega)
    assert np.isclose(S, -11894.2)
    assert np.isclose(D, 13654.4)
    assert np.isclose(P, -10523881)

    # test with three species and a 5% H minority fraction in a D plasma
    n_3 = np.array([1, 1, 5/100])*1e19/u.m**3
    S, D, P = cold_plasma_permittivity_SDP(B, three_species, n_3, omega)
    assert np.isclose(S, -11753.3)
    assert np.isclose(D, 13408.99181054283)
    assert np.isclose(P, -10524167.9)


def test_cold_plasma_permittivity_LRP():
    r"""Test the cold plasma dielectric tensor (L, R, P) function"""
    # test with one species
    L, R, P = cold_plasma_permittivity_LRP(B, single_species, n, omega)
    assert np.isclose(L, -327.9)
    assert np.isclose(R, 330.105)
    assert np.isclose(P, -1052100.6)


def test_SDP_LRP_relationships():
    r"""
    Test the relationships between (S, D, P) notation in Stix basis and
    (L, R, P) notation in the rotating basis, ie :
    S = (R+L)/2 and D = (R-L)/2
    and
    Checks for R=S+D and L=S-D
    """
    # test with one species
    S, D, P = cold_plasma_permittivity_SDP(B, single_species, n, omega)
    L, R, P = cold_plasma_permittivity_LRP(B, single_species, n, omega)

    assert np.isclose(S, (R+L)/2)
    assert np.isclose(D, (R-L)/2)
    assert np.isclose(R, S+D)
    assert np.isclose(L, S-D)
