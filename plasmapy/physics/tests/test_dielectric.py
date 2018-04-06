"""Tests for functions that calculate plasma dielectric parameters in
dielectry.py"""

import numpy as np
from astropy import units as u

from ..dielectric import (cold_plasma_permittivity_LRP,
                          cold_plasma_permittivity_SDP)

from ..parameters import (plasma_frequency, gyrofrequency)

B = 1.0 * u.T
n = [1e18 / u.m ** 3]
omega = 55e6 * u.rad / u.s

single_species = ['e']
two_species = ['e', 'D+']
three_species = ['e', 'D+', 'H+']


class Test_ColdPlasmaPermittivity(object):
    def test_proton_electron_plasma(self):
        """
        Test proton-electron plasma against the (approximate)
        analytical formulas
        """
        B = 1 * u.T
        n = [1, 1] * 1 / u.m ** 3
        omega = 1 * u.rad / u.s
        omega_ce = gyrofrequency(B, particle='e', signed=True)
        omega_pe = plasma_frequency(n[0], particle='e')
        omega_cp = abs(omega_ce) / 1860
        omega_pp = omega_pe / 43

        S_analytical = 1 \
            - omega_pe ** 2 / (omega ** 2 - omega_ce ** 2) \
            - omega_pp ** 2 / (omega ** 2 - omega_cp ** 2)

        D_analytical = \
            + omega_ce / omega * omega_pe ** 2 / (omega ** 2 - omega_ce ** 2) \
            + omega_cp / omega * omega_pp ** 2 / (omega ** 2 - omega_cp ** 2)

        P_analytical = 1 - (omega_pe ** 2 + omega_pp ** 2) / omega ** 2

        species = ['e', 'p']
        S, D, P = cold_plasma_permittivity_SDP(B, species, n, omega)

        assert np.isclose(S, S_analytical)
        assert np.isclose(D, D_analytical)
        assert np.isclose(P, P_analytical)

    def test_three_species(self):
        """
        Test with three species (2 ions): D plasma with 5%H minority fraction
        """
        n_3 = np.array([1, 1, 5 / 100]) * 1e19 / u.m ** 3
        S, D, P = cold_plasma_permittivity_SDP(B, three_species, n_3, omega)
        assert np.isclose(S, -11753.3)
        assert np.isclose(D, 13408.99181054283)
        assert np.isclose(P, -10524167.9)

    def test_SD_to_LR_relationships(self):
        """
        Test the relationships between (S, D, P) notation in Stix basis and
        (L, R, P) notation in the rotating basis, ie test:
         S = (R+L)/2 and D = (R-L)/2
        and
         R = S+D and L = S-D
        """
        # test with a single species
        S, D, _ = cold_plasma_permittivity_SDP(B, single_species, n, omega)
        L, R, _ = cold_plasma_permittivity_LRP(B, single_species, n, omega)

        assert np.isclose(R, S + D)
        assert np.isclose(L, S - D)
        assert np.isclose(S, (R + L) / 2)
        assert np.isclose(D, (R - L) / 2)
