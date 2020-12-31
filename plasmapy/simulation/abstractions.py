"""Abstract classes for numerical simulations."""

__all__ = ["AbstractSimulation", "AbstractTimeDependentSimulation"]

from abc import ABC, abstractmethod
from typing import NoReturn


class AbstractSimulation(ABC):
    """
    A prototype abstract interface for numerical simulations.

    Notes
    -----
    This interface is incomplete and unstable, and is thus subject to
    change at any time.
    """

    @abstractmethod
    def summarize(self):
        """
        Print out a summary of the simulation parameters and status.
        """

        pass

    @abstractmethod
    def initialize(self) -> NoReturn:
        """Prepare the simulation to be run."""

        pass

    @abstractmethod
    def simulate(self):
        """Perform the actual simulation."""

        pass

    @abstractmethod
    def finalize(self) -> NoReturn:
        """Perform the steps to close the simulation and output data."""

        pass


class AbstractTimeDependentSimulation(AbstractSimulation):
    """
    A prototype abstract interface for time-dependent numerical simulations.

    Notes
    -----
    This interface is incomplete and unstable, and is thus subject to
    change at any time.
    """

    pass


class AbstractNormalizations(ABC):  # coverage: ignore
    """
    An abstract base class to represent the normalizations of systems of
    equations describing plasmas.

    Warnings
    --------
    This interface is unstable and is subject to change.
    """

    @abstractmethod
    def current_density(self) -> u.A / u.m ** 2:
        """The current density normalization."""
        pass

    @abstractmethod
    def diffusivity(self) -> u.m ** 2 / u.s:
        """The normalization for diffusivity."""
        pass

    @abstractmethod
    def dynamic_viscosity(self) -> u.Pa * u.s:
        """The normalization for dynamic viscosity."""
        pass

    @abstractmethod
    def electric_field(self) -> u.V / u.m:
        """The electric field normalization."""
        pass

    @abstractmethod
    def heat_flux(self) -> u.J * u.m ** -2 * u.s ** -1:
        """The normalization for heat flux."""
        pass

    @abstractmethod
    def ion(self) -> Particle:
        """
        The ion in the plasma, which could be a positron in the case
        of a pair plasma.
        """
        pass

    @abstractmethod
    def length(self) -> u.m:
        """The normalization for length."""
        pass

    @abstractmethod
    def magnetic_field(self) -> u.T:
        """The magnetic field normalization."""
        pass

    @abstractmethod
    def magnetic_flux(self) -> u.T * u.m:
        """The normalization for the magnetic flux or vector potential."""
        pass

    @abstractmethod
    def mass(self) -> u.kg:
        """The normalization for mass."""
        pass

    @abstractmethod
    def mass_density(self) -> u.kg / u.m ** 3:
        """The normalization for mass density."""
        pass

    @abstractmethod
    def number_density(self) -> u.m ** -3:
        """The normalization for number density."""
        pass

    @abstractmethod
    def pressure(self) -> u.Pa:
        """The normalization for pressure."""
        pass

    @abstractmethod
    def temperature(self) -> u.K:
        """The normalization for temperature."""
        pass

    @abstractmethod
    def thermal_conductivity(self) -> u.W / (u.K * u.m):
        """The normalization for thermal conductivity."""
        pass

    @abstractmethod
    def time(self) -> u.s:
        """The normalization for time."""
        pass

    @abstractmethod
    def wavenumber(self) -> u.m ** -1:
        """The normalization for inverse length."""
        pass

    @abstractmethod
    def velocity(self) -> u.m / u.s:
        """The normalization for velocity."""
        pass

    @abstractmethod
    def volumetric_heating_rate(self) -> u.J * u.m ** -3 * u.s ** -1:
        """The normalization for the volumetric heating rate."""
        pass

    @abstractmethod
    def volumetric_rate(self) -> u.m ** -3 * u.s ** -1:
        """
        The normalization for a volumetric rate.

        This normalization is applicable to, for example, the number
        of collisions per unit volume per unit time.
        """
        pass
