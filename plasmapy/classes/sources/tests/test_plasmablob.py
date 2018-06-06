import pytest
import numpy as np
import astropy.units as u

from plasmapy.classes.sources import plasmablob
from plasmapy.utils.exceptions import InvalidParticleError


class Test_PlasmaBlobRegimes:
    def test_intermediate_coupling(self):
        r"""
        Method to test for coupling parameter for a plasma.

        Tests against expected value for coupling parameter for a
        plasma in the intermediate coupling regime.

        The input values in this case have no special significance
        and are just to get the desired output.
        """

        T_e = 25 * 15e3 * u.K
        n_e = 1e26 * u.cm ** -3
        Z = 2.0
        particle = 'p'
        blob = plasmablob.PlasmaBlob(T_e=T_e,
                                     n_e=n_e,
                                     Z=Z,
                                     particle=particle)

        expect_regime = 'Intermediate coupling regime: Gamma = 10.585076050938532.'
        regime, _ = blob.regimes()
        testTrue = regime == expect_regime

        errStr = f"Regime should be {expect_regime}, but got {regime} instead."
        assert testTrue, errStr

    def test_strongly_coupled(self):
        r"""
        Method to test for coupling parameter for a plasma.

        Tests against expected value for coupling parameter for a
        plasma in the strongly coupled regime.

        The input values in this case have no special significance
        and are just to get the desired output.
        """

        T_e = 5 * 15e3 * u.K
        n_e = 1e26 * u.cm ** -3
        Z = 3.0
        particle = 'p'
        blob = plasmablob.PlasmaBlob(T_e=T_e,
                                     n_e=n_e,
                                     Z=Z,
                                     particle=particle)

        expect_regime = 'Strongly coupled regime: Gamma = 104.02780112828943.'

        regime, _ = blob.regimes()
        testTrue = regime == expect_regime

        errStr = f"Regime should be {expect_regime}, but got {regime} instead."
        assert testTrue, errStr

    def test_weakly_coupled(self):
        r"""
        Method to test for coupling parameter for a plasma.

        Tests against expected value for coupling parameter for a
        plasma in the weakly coupled regime.

        The input values in this case have no special significance
        and are just to get the desired output.
        """

        T_e = 15 * 11e3 * u.K
        n_e = 1e15 * u.cm ** -3
        Z = 2.5
        particle = 'p'
        blob = plasmablob.PlasmaBlob(T_e=T_e,
                                     n_e=n_e,
                                     Z=Z,
                                     particle=particle)

        expect_regime = 'Weakly coupled regime: Gamma = 0.0075178096952688445.'

        regime, _ = blob.regimes()
        testTrue = regime == expect_regime

        errStr = f"Regime should be {expect_regime}, but got {regime} instead."
        assert testTrue, errStr

    def test_thermal_kinetic_energy_dominant(self):
        r"""
        Method to test for degeneracy parameter for a plasma.

        Tests against expected value for degeneracy parameter for a
        plasma in the thermal degenerate regime.

        The input values in this case have no special significance
        and are just to get the desired output.
        """

        T_e = 10 * 11e3 * u.K
        n_e = 1e20 * u.cm ** -3
        Z = 2.5
        particle = 'p'
        blob = plasmablob.PlasmaBlob(T_e=T_e,
                                     n_e=n_e,
                                     Z=Z,
                                     particle=particle)

        expect_regime = 'Thermal kinetic energy dominant: Theta = 120.65958493847927'

        _, regime = blob.regimes()
        testTrue = regime == expect_regime

        errStr = f"Regime should be {expect_regime}, but got {regime} instead."
        assert testTrue, errStr

    def test_fermi_quantum_energy_dominant(self):
        r"""
        Method to test for degeneracy parameter for a plasma.

        Tests against expected value for degeneracy parameter for a
        plasma in the Fermi degenerate regime.

        The input values in this case have no special significance
        and are just to get the desired output.
        """

        T_e = 6 * 15e3 * u.K
        n_e = 1e26 * u.cm ** -3
        Z = 3.0
        particle = 'p'
        blob = plasmablob.PlasmaBlob(T_e=T_e,
                                     n_e=n_e,
                                     Z=Z,
                                     particle=particle)

        expect_regime = 'Fermi quantum energy dominant: Theta = 0.009872147858602853'

        _, regime = blob.regimes()
        testTrue = regime == expect_regime

        errStr = f"Regime should be {expect_regime}, but got {regime} instead."
        assert testTrue, errStr

    def test_both_fermi_and_thermal_energy_important(self):
        r"""
        Method to test for degeneracy parameter for a plasma.

        Tests against expected value for degeneracy parameter for a
        plasma whose both Fermi and thermal energy are important.

        The input values in this case have no special significance
        and are just to get the desired output.
        """

        T_e = 5 * 15e3 * u.K
        n_e = 1e25 * u.cm ** -3
        Z = 2.0
        particle = 'p'
        blob = plasmablob.PlasmaBlob(T_e=T_e,
                                     n_e=n_e,
                                     Z=Z,
                                     particle=particle)

        expect_regime = 'Both Fermi and thermal energy important: Theta = 0.03818537605355442'

        _, regime = blob.regimes()
        testTrue = regime == expect_regime

        errStr = f"Regime should be {expect_regime}, but got {regime} instead."
        assert testTrue, errStr


class Test_PlasmaBlob:
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T_e = 5 * 11e3 * u.K
        self.n_e = 1e23 * u.cm ** -3
        self.Z = 2.5
        self.particle = 'p'
        self.blob = plasmablob.PlasmaBlob(T_e=self.T_e,
                                          n_e=self.n_e,
                                          Z=self.Z,
                                          particle=self.particle)
        self.couplingVal = 10.468374460435724
        self.thetaVal = 0.6032979246923964

    def test_invalid_particle(self):
        """
        Checks if function raises error for invalid particle.
        """
        with pytest.raises(InvalidParticleError):
            plasmablob.PlasmaBlob(T_e=self.T_e,
                                  n_e=self.n_e,
                                  Z=self.Z,
                                  particle="cupcakes")

    def test_electron_temperature(self):
        """Testing if we get the same electron temperature we put in """
        testTrue = self.T_e == self.blob.electron_temperature
        errStr = (f"Input electron temperature {self.T_e} should be equal to "
                  f"electron temperature of class "
                  f"{self.blob.electron_temperature}.")
        assert testTrue, errStr

    def test_electron_density(self):
        """Testing if we get the same electron density we put in """
        testTrue = self.n_e == self.blob.electron_density
        errStr = (f"Input electron density {self.n_e} should be equal to "
                  f"electron density of class "
                  f"{self.blob.electron_density}.")
        assert testTrue, errStr

    def test_ionization(self):
        """Testing if we get the same ionization we put in """
        testTrue = self.Z == self.blob.ionization
        errStr = (f"Input ionization {self.Z} should be equal to "
                  f"ionization of class "
                  f"{self.blob.ionization}.")
        assert testTrue, errStr

    def test_composition(self):
        """Testing if we get the same composition (particle) we put in """
        testTrue = self.particle == self.blob.composition
        errStr = (f"Input particle {self.particle} should be equal to "
                  f"composition of class "
                  f"{self.blob.composition}.")
        assert testTrue, errStr

    def test_coupling(self):
        """
        Tests if coupling  method value meets expected value.
        """
        methodVal = self.blob.coupling()
        errStr = (f"Coupling parameter should be {self.couplingVal} "
                  f"and not {methodVal.si.value}.")
        testTrue = np.isclose(methodVal.value,
                              self.couplingVal,
                              rtol=1e-8,
                              atol=0.0)
        assert testTrue, errStr

    def test_quantum_theta(self):
        """
        Tests if degeneracy parameter method value meets expected value.
        """
        methodVal = self.blob.quantum_theta()
        errStr = (f"Degeneracy parameter should be {self.thetaVal} "
                  f"and not {methodVal.si.value}.")
        testTrue = np.isclose(methodVal.value,
                              self.thetaVal,
                              rtol=1e-8,
                              atol=0.0)
        assert testTrue, errStr
