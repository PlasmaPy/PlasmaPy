"""
Defines the core Plasma class used by PlasmaPy to represent plasma properties.
"""
import warnings

import numpy as np
import astropy.units as u

from plasmapy.constants import (m_p,
                                m_e,
                                c,
                                mu0,
                                k_B,
                                e,
                                eps0,
                                pi,
                                )
from plasmapy.physics.parameters import _grab_charge
from plasmapy.physics.dimensionless import (quantum_theta,
                                            )
from plasmapy.physics.transport import (coupling_parameter,
                                        )
from plasmapy.atomic import particle_mass

from plasmapy.utils import call_string, CouplingWarning

__all__ = [
    "Plasma3D",
    "PlasmaBlob"
]


class Plasma3D:
    """
    Core class for describing and calculating plasma parameters with
    spatial dimensions.

    Attributes
    ----------
    x : `astropy.units.Quantity`
        x-coordinates within the plasma domain. Equal to the
        `domain_x` input parameter.
    y : `astropy.units.Quantity`
        y-coordinates within the plasma domain. Equal to the
        `domain_y` input parameter.
    z : `astropy.units.Quantity`
        z-coordinates within the plasma domain. Equal to the
        `domain_z` input parameter.
    grid : `astropy.units.Quantity`
        (3, x, y, z) array containing the values of each coordinate at
        every point in the domain.
    domain_shape : tuple
        Shape of the plasma domain.
    density : `astropy.units.Quantity`
        (x, y, z) array of mass density at every point in the domain.
    momentum : `astropy.units.Quantity`
        (3, x, y, z) array of the momentum vector at every point in
        the domain.
    pressure : `astropy.units.Quantity`
        (x, y, z) array of pressure at every point in the domain.
    magnetic_field : `astropy.units.Quantity`
        (3, x, y, z) array of the magnetic field vector at every point
        in the domain.

    Parameters
    ----------
    domain_x : `astropy.units.Quantity`
        1D array of x-coordinates for the plasma domain. Must have
        units convertable to length.
    domain_y : `astropy.units.Quantity`
        1D array of y-coordinates for the plasma domain. Must have
        units convertable to length.
    domain_z : `astropy.units.Quantity`
        1D array of z-coordinates for the plasma domain. Must have
        units convertable to length.

    """
    @u.quantity_input(domain_x=u.m, domain_y=u.m, domain_z=u.m)
    def __init__(self, domain_x, domain_y, domain_z):
        # Define domain sizes
        self.x = domain_x
        self.y = domain_y
        self.z = domain_z

        self.grid = np.array(np.meshgrid(self.x, self.y, self.z,
                                         indexing='ij'))
        self.domain_shape = (len(self.x), len(self.y), len(self.z))

        # Initiate core plasma variables
        self.density = np.zeros(self.domain_shape) * u.kg / u.m**3
        self.momentum = np.zeros((3, *self.domain_shape)) * u.kg / (u.m ** 2 * u.s)
        self.pressure = np.zeros(self.domain_shape) * u.Pa
        self.magnetic_field = np.zeros((3, *self.domain_shape)) * u.T
        self.electric_field = np.zeros((3, *self.domain_shape)) * u.V / u.m

    @property
    def velocity(self):
        return self.momentum / self.density

    @property
    def magnetic_field_strength(self):
        B = self.magnetic_field
        return np.sqrt(np.sum(B * B, axis=0))

    @property
    def electric_field_strength(self):
        E = self.electric_field
        return np.sqrt(np.sum(E * E, axis=0))

    @property
    def alfven_speed(self):
        B = self.magnetic_field
        rho = self.density
        return np.sqrt(np.sum(B * B, axis=0) / (mu0 * rho))


class PlasmaBlob:
    """
    Class for describing and calculating plasma parameters without
    spatial/temporal description.
    """
    @u.quantity_input(T_e=u.K, n_e=u.m**-3)
    def __init__(self, T_e, n_e, Z=None, particle='p'):
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
        Intermediate coupling regime: Gamma = 0.012502837623108332.
        Thermal kinetic energy dominant: Theta = 109690.53176225389

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
        argument_dict = {'T_e': self.T_e,
                         'n_e': self.n_e,
                         'particle': self.particle,
                         'Z': self.Z}

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
            quantum_theta_str = (f"Fermi quantum energy dominant: Theta = "
                                 f"{quantum_theta}")
        elif quantum_theta >= 100:
            # thermal kinetic energy dominant
            quantum_theta_str = (f"Thermal kinetic energy dominant: Theta = "
                                 f"{quantum_theta}")
        else:
            # intermediate regime
            quantum_theta_str = (f"Both Fermi and thermal energy important: "
                                 f"Theta = {quantum_theta}")

        # summarizing and printing/returning regimes
        aggregateStrs = [coupling_str,
                         quantum_theta_str,
                         ]
        return aggregateStrs

    def coupling(self):
        """
        Ion-ion coupling parameter to determine if quantum/coupling effects
        are important. This compares Coulomb potential energy to thermal
        kinetic energy.
        """
        couple = coupling_parameter(self.T_e,
                                    self.n_e,
                                    (self.particle, self.particle),
                                    self.Z)
        if couple < 0.01:
            warnings.warn(f"Coupling parameter is {couple}, you might have strong coupling effects",
                          CouplingWarning)

        return couple
    def quantum_theta(self):
        """
        Quantum theta parameter, which compares Fermi kinetic energy to
        thermal kinetic energy to check if quantum effects are important.
        """
        theta = quantum_theta(self.T_e, self.n_e)
        return theta
