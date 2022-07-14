"""Abstract classes for numerical simulations."""

__all__ = [
    "AbstractNormalizations",
    "AbstractSimulation",
    "AbstractTimeDependentSimulation",
]

import astropy.units as u

from abc import ABC, abstractmethod
from typing import NoReturn


class AbstractSimulation(ABC):
    """
    A prototype abstract interface for numerical simulations.

    .. warning::
        This interface is unstable and subject to change.
    """

    @abstractmethod
    def summarize(self):
        """
        Print a summary of the simulation parameters and status.
        """
        ...

    @abstractmethod
    def initialize(self) -> NoReturn:
        """Prepare the simulation to be run."""
        ...

    @abstractmethod
    def simulate(self):
        """Perform the simulation."""
        ...

    @abstractmethod
    def finalize(self) -> NoReturn:
        """Perform the steps to close the simulation."""
        ...


class AbstractTimeDependentSimulation(AbstractSimulation):
    """
    A prototype abstract interface for time-dependent numerical simulations.

    .. warning::
        This interface is unstable and is subject to change.
    """

    ...


class AbstractNormalizations(ABC):
    """
    An abstract base class to represent the normalization constants for
    systems of equations describing plasmas.

    .. warning::
        This interface is unstable and is subject to change.
    """

    @property
    @abstractmethod
    def current_density(self) -> u.A / u.m**2:
        """The current density normalization."""
        ...

    @property
    @abstractmethod
    def diffusivity(self) -> u.m**2 / u.s:
        """The normalization for diffusivity."""
        ...

    @property
    @abstractmethod
    def dynamic_viscosity(self) -> u.Pa * u.s:
        """The normalization for dynamic viscosity."""
        ...

    @property
    @abstractmethod
    def electric_field(self) -> u.V / u.m:
        """The electric field normalization."""
        ...

    @property
    @abstractmethod
    def heat_flux(self) -> u.J * u.m**-2 * u.s**-1:
        """The normalization for heat flux."""
        ...

    @property
    @abstractmethod
    def length(self) -> u.m:
        """The normalization for length."""
        ...

    @property
    @abstractmethod
    def magnetic_field(self) -> u.T:
        """The magnetic field normalization."""
        ...

    @property
    @abstractmethod
    def magnetic_flux(self) -> u.T * u.m:
        """The normalization for the magnetic flux or vector potential."""
        ...

    @property
    @abstractmethod
    def mass(self) -> u.kg:
        """The normalization for mass."""
        ...

    @property
    @abstractmethod
    def mass_density(self) -> u.kg / u.m**3:
        """The normalization for mass density."""
        ...

    @property
    @abstractmethod
    def number_density(self) -> u.m**-3:
        """The normalization for number density."""
        ...

    @property
    @abstractmethod
    def pressure(self) -> u.Pa:
        """The normalization for pressure."""
        ...

    @property
    @abstractmethod
    def temperature(self) -> u.K:
        """The normalization for temperature."""
        ...

    @property
    @abstractmethod
    def thermal_conductivity(self) -> u.W / (u.K * u.m):
        """The normalization for thermal conductivity."""
        ...

    @property
    @abstractmethod
    def time(self) -> u.s:
        """The normalization for time."""
        ...

    @property
    @abstractmethod
    def wavenumber(self) -> u.m**-1:
        """The normalization for inverse length."""
        ...

    @property
    @abstractmethod
    def velocity(self) -> u.m / u.s:
        """The normalization for velocity."""
        ...

    @property
    @abstractmethod
    def volumetric_heating_rate(self) -> u.J * u.m**-3 * u.s**-1:
        """The normalization for the volumetric heating rate."""
        ...

    @property
    @abstractmethod
    def volumetric_rate(self) -> u.m**-3 * u.s**-1:
        """
        The normalization for a volumetric rate.

        This normalization is applicable to, for example, the number
        of collisions per unit volume per unit time.
        """
        ...
