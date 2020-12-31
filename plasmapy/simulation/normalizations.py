"""Classes that describes normalizations for different systems of equations."""

import abc

from astropy import units as u

from plasmapy.particles import Particle


class AbstractNormalizations(abc.ABC):  # coverage: ignore
    """
    An abstract base class to represent the normalizations of systems of
    equations describing plasmas.

    Warnings
    --------
    This interface is unstable and is subject to change.
    """

    @abc.abstractmethod
    def current_density(self) -> u.A / u.m ** 2:
        """The current density normalization."""
        raise NotImplementedError

    @abc.abstractmethod
    def diffusivity(self) -> u.m ** 2 / u.s:
        """The normalization for diffusivity."""
        raise NotImplementedError

    @property
    def dynamic_viscosity(self) -> u.Pa * u.s:
        """The normalization for dynamic viscosity."""
        raise NotImplementedError

    @abc.abstractmethod
    def electric_field(self) -> u.V / u.m:
        """The electric field normalization."""
        raise NotImplementedError

    @abc.abstractmethod
    def heat_flux(self) -> u.J * u.m ** -2 * u.s ** -1:
        """The normalization for heat flux."""
        raise NotImplementedError

    @abc.abstractmethod
    def ion(self) -> Particle:
        """
        The ion in the plasma, which could be a positron in the case
        of a pair plasma.
        """
        raise NotImplementedError

    @abc.abstractmethod
    def length(self) -> u.m:
        """The normalization for length."""
        raise NotImplementedError

    @abc.abstractmethod
    def magnetic_field(self) -> u.T:
        """The magnetic field normalization."""
        raise NotImplementedError

    @abc.abstractmethod
    def magnetic_flux(self) -> u.T * u.m:
        """The normalization for the magnetic flux or vector potential."""
        raise NotImplementedError

    @abc.abstractmethod
    def mass(self) -> u.kg:
        """The normalization for mass."""
        raise NotImplementedError

    @abc.abstractmethod
    def mass_density(self) -> u.kg / u.m ** 3:
        """The normalization for mass density."""
        raise NotImplementedError

    @abc.abstractmethod
    def number_density(self) -> u.m ** -3:
        """The normalization for number density."""
        raise NotImplementedError

    @abc.abstractmethod
    def pressure(self) -> u.Pa:
        """The normalization for pressure."""
        raise NotImplementedError

    @abc.abstractmethod
    def temperature(self) -> u.K:
        """The normalization for temperature."""
        raise NotImplementedError

    @property
    def thermal_conductivity(self) -> u.W / (u.K * u.m):
        """The normalization for thermal conductivity."""
        raise NotImplementedError

    @abc.abstractmethod
    def time(self) -> u.s:
        """The normalization for time."""
        raise NotImplementedError

    @abc.abstractmethod
    def wavenumber(self) -> u.m ** -1:
        """The normalization for inverse length."""
        raise NotImplementedError

    @abc.abstractmethod
    def velocity(self) -> u.m / u.s:
        """The normalization for velocity."""
        raise NotImplementedError

    @abc.abstractmethod
    def volumetric_heating_rate(self) -> u.J * u.m ** -3 * u.s ** -1:
        """The normalization for the volumetric heating rate."""
        raise NotImplementedError

    @abc.abstractmethod
    def volumetric_rate(self) -> u.m ** -3 * u.s ** -1:
        """
        The normalization for a volumetric rate.

        This normalization is applicable to, for example, the number
        of collisions per unit volume per unit time.
        """
        raise NotImplementedError
