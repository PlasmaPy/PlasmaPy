"""Tests for functions that calculate plasma dielectric parameters in
dielectry.py"""

import numpy as np
from astropy import units as u

from ..dielectric import (
    RotatingTensorElements,
    StixTensorElements,
    cold_plasma_permittivity_LRP,
    cold_plasma_permittivity_SDP,
    permittivity_1D_Maxwellian,
)
from ..parameters import gyrofrequency, plasma_frequency, thermal_speed

B = 1.0 * u.T
n = [1e18 / u.m ** 3]
omega = 55e6 * u.rad / u.s

single_species = ["e"]
two_species = ["e", "D+"]
three_species = ["e", "D+", "H+"]


class Test_ColdPlasmaPermittivity(object):
    def test_proton_electron_plasma(self):
        """
        Test proton-electron plasma against the (approximate)
        analytical formulas
        """
        B = 1 * u.T
        n = [1, 1] * 1 / u.m ** 3
        omega = 1 * u.rad / u.s
        omega_ce = gyrofrequency(B, particle="e", signed=True)
        omega_pe = plasma_frequency(n[0], particle="e")
        omega_cp = abs(omega_ce) / 1860
        omega_pp = omega_pe / 43

        S_analytical = (
            1
            - omega_pe ** 2 / (omega ** 2 - omega_ce ** 2)
            - omega_pp ** 2 / (omega ** 2 - omega_cp ** 2)
        )

        D_analytical = +omega_ce / omega * omega_pe ** 2 / (
            omega ** 2 - omega_ce ** 2
        ) + omega_cp / omega * omega_pp ** 2 / (omega ** 2 - omega_cp ** 2)

        P_analytical = 1 - (omega_pe ** 2 + omega_pp ** 2) / omega ** 2

        species = ["e", "p"]
        S, D, P = tuple_result = cold_plasma_permittivity_SDP(B, species, n, omega)

        assert tuple_result.sum is S
        assert tuple_result.difference is D
        assert tuple_result.plasma is P
        assert isinstance(tuple_result, StixTensorElements)

        assert np.isclose(S, S_analytical)
        assert np.isclose(D, D_analytical)
        assert np.isclose(P, P_analytical)

        L, R, P = rotating_tuple_result = cold_plasma_permittivity_LRP(
            B, species, n, omega
        )
        assert rotating_tuple_result.left is L
        assert rotating_tuple_result.right is R
        assert rotating_tuple_result.plasma is P
        assert isinstance(rotating_tuple_result, RotatingTensorElements)

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

    def test_numpy_array_workflow(self):
        """
        As per @jhillairet at:
        https://github.com/PlasmaPy/PlasmaPy/issues/539#issuecomment-425337810
        """
        ns = np.logspace(17, 19, 50) / u.m ** 3
        B0 = 4 * u.T
        omega_RF = 2 * np.pi * 50e6 * (u.rad / u.s)

        S, D, P = cold_plasma_permittivity_SDP(
            B=B0, species=["e", "D+"], n=[ns, ns], omega=omega_RF
        )
        assert S.shape == D.shape == P.shape == (50,)


class Test_permittivity_1D_Maxwellian:
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T = 30 * 11600 * u.K
        self.n = 1e18 * u.cm ** -3
        self.particle = "Ne"
        self.z_mean = 8 * u.dimensionless_unscaled
        self.vTh = thermal_speed(self.T, self.particle, method="most_probable")
        self.omega = 5.635e14 * 2 * np.pi * u.rad / u.s
        self.kWave = self.omega / self.vTh
        self.True1 = (
            -6.728092569241431e-08 + 5.760379561405176e-07j
        ) * u.dimensionless_unscaled

    def test_known1(self):
        """
        Tests Fermi_integral for expected value.
        """
        methodVal = permittivity_1D_Maxwellian(
            self.omega, self.kWave, self.T, self.n, self.particle, self.z_mean
        )
        testTrue = np.isclose(methodVal, self.True1, rtol=1e-6, atol=0.0)
        errStr = f"Permittivity value should be {self.True1} and not " f"{methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 + 1e-15
        methodVal = permittivity_1D_Maxwellian(
            self.omega, self.kWave, self.T, self.n, self.particle, self.z_mean
        )
        testTrue = not np.isclose(methodVal, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"Permittivity value test gives {methodVal} "
            f"and should not be equal to {fail1}."
        )
        assert testTrue, errStr
