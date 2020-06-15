"""
Defines the core Plasma class used by PlasmaPy to represent plasma properties.
"""
import warnings

import astropy.units as u
from plasmapy.formulary.collisions import coupling_parameter
from plasmapy.formulary.dimensionless import quantum_theta
from plasmapy.formulary.parameters import _grab_charge
from plasmapy.particles import particle_mass
from plasmapy.plasma.plasma_base import GenericPlasma
from plasmapy.utils import CouplingWarning, call_string

__all__ = ["PlasmaBlob"]


class PlasmaBlob(GenericPlasma):
    """
    Class for describing and calculating plasma parameters without
    spatial/temporal description.
    """

    @u.quantity_input(T_e=u.K, n_e=u.m ** -3)
    def __init__(self, T_e, n_e, Z=None, particle="p"):
        """
        Initialize plasma paramters.
        The most basic description is composition (ion), temperature,
        density, and ionization.
        """
        self.T_e = T_e
        self.n_e = n_e
        self.particle = particle
        self.Z = _grab_charge(particle, Z)
        # extract mass from particle
        self.ionMass = particle_mass(self.particle)

    def __str__(self):
        """
        Fetches regimes for easy printing

        Examples
        --------
        >>> print(PlasmaBlob(1e4*u.K, 1e20/u.m**3, particle='p'))
        PlasmaBlob(T_e=10000.0*u.K, n_e=1e+20*u.m**-3, particle='p', Z=1)
        Intermediate coupling regime: Gamma = 0.01250283...
        Thermal kinetic energy dominant: Theta = 109690.5...

        """
        return self.__repr__() + "\n" + "\n".join(self.regimes())

    def __repr__(self):
        """

        Returns
        -------
        str

        Examples
        --------
        >>> from astropy import units as u
        >>> PlasmaBlob(1e4*u.K, 1e20/u.m**3, particle='p')
        PlasmaBlob(T_e=10000.0*u.K, n_e=1e+20*u.m**-3, particle='p', Z=1)
        """
        argument_dict = {
            "T_e": self.T_e,
            "n_e": self.n_e,
            "particle": self.particle,
            "Z": self.Z,
        }

        return call_string(PlasmaBlob, (), argument_dict)

    @property
    def electron_temperature(self):
        return self.T_e

    @property
    def electron_density(self):
        return self.n_e

    @property
    def ionization(self):
        return self.Z

    @property
    def composition(self):
        return self.particle

    def regimes(self):
        """
        Generate a comprehensive description of the plasma regimes
        based on plasma properties and consequent plasma parameters.
        """
        # getting dimensionless parameters
        coupling = self.coupling()
        quantum_theta = self.quantum_theta()

        # determining regimes based off dimensionless parameters
        # coupling
        if coupling <= 0.01:
            # weakly coupled
            coupling_str = f"Weakly coupled regime: Gamma = {coupling}."
        elif coupling >= 100:
            # strongly coupled
            coupling_str = f"Strongly coupled regime: Gamma = {coupling}."
        else:
            # intermediate regime
            coupling_str = f"Intermediate coupling regime: Gamma = {coupling}."
        # quantum_theta
        if quantum_theta <= 0.01:
            # Fermi energy dominant
            quantum_theta_str = (
                f"Fermi quantum energy dominant: Theta = " f"{quantum_theta}"
            )
        elif quantum_theta >= 100:
            # thermal kinetic energy dominant
            quantum_theta_str = (
                f"Thermal kinetic energy dominant: Theta = " f"{quantum_theta}"
            )
        else:
            # intermediate regime
            quantum_theta_str = (
                f"Both Fermi and thermal energy important: " f"Theta = {quantum_theta}"
            )

        # summarizing and printing/returning regimes
        aggregateStrs = [coupling_str, quantum_theta_str]
        return aggregateStrs

    def coupling(self):
        """
        Ion-ion coupling parameter to determine if quantum/coupling effects
        are important. This compares Coulomb potential energy to thermal
        kinetic energy.
        """
        couple = coupling_parameter(
            self.T_e, self.n_e, (self.particle, self.particle), self.Z
        )
        if couple < 0.01:
            warnings.warn(
                f"Coupling parameter is {couple}, you might have strong coupling effects",
                CouplingWarning,
            )

        return couple

    def quantum_theta(self):
        """
        Quantum theta parameter, which compares Fermi kinetic energy to
        thermal kinetic energy to check if quantum effects are important.
        """
        theta = quantum_theta(self.T_e, self.n_e)
        return theta

    @classmethod
    def is_datasource_for(cls, **kwargs):
        match = "T_e" in kwargs.keys() and "n_e" in kwargs.keys()
        return match
