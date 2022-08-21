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
    An abstract base class to represent the :term:`normalization`
    constants for systems of equations describing plasmas.

    .. warning::
        This interface is unstable and is subject to change.
    """

    @property
    @abstractmethod
    def acceleration(self) -> u.m / u.s**2:
        r"""The acceleration :term:`normalization`\ , :math:`a_⭑`\ ."""

    @property
    @abstractmethod
    def current_density(self) -> u.A / u.m**2:
        r"""The current density :term:`normalization`\ ."""
        ...

    @property
    @abstractmethod
    def diffusivity(self) -> u.m**2 / u.s:
        """The :term:`normalization` for diffusivity."""
        ...

    @property
    @abstractmethod
    def dynamic_viscosity(self) -> u.Pa * u.s:
        """The :term:`normalization` for dynamic viscosity."""
        ...

    @property
    @abstractmethod
    def electric_field(self) -> u.V / u.m:
        r"""The electric field :term:`normalization`\ ."""
        ...

    @property
    @abstractmethod
    def frequency(self) -> u.s**-1:
        r"""The frequency :term:`normalization`\ ."""

    @property
    @abstractmethod
    def heat_flux(self) -> u.J * u.m**-2 * u.s**-1:
        """The :term:`normalization` for heat flux."""
        ...

    @property
    @abstractmethod
    def length(self) -> u.m:
        """The :term:`normalization` for length."""
        ...

    @property
    @abstractmethod
    def magnetic_field(self) -> u.T:
        r"""The magnetic field :term:`normalization`\ ."""
        ...

    @property
    @abstractmethod
    def magnetic_flux(self) -> u.T * u.m:
        """
        The :term:`normalization` for the magnetic flux or vector
        potential.
        """
        ...

    @property
    @abstractmethod
    def mass(self) -> u.kg:
        """The :term:`normalization` for mass."""
        ...

    @property
    @abstractmethod
    def mass_density(self) -> u.kg / u.m**3:
        """The :term:`normalization` for mass density."""
        ...

    @property
    @abstractmethod
    def number_density(self) -> u.m**-3:
        """The :term:`normalization` for number density."""
        ...

    @property
    @abstractmethod
    def pressure(self) -> u.Pa:
        """The :term:`normalization` for pressure."""
        ...

    @property
    @abstractmethod
    def temperature(self) -> u.K:
        """The :term:`normalization` for temperature."""
        ...

    @property
    @abstractmethod
    def thermal_conductivity(self) -> u.W / (u.K * u.m):
        """The :term:`normalization` for thermal conductivity."""
        ...

    @property
    @abstractmethod
    def time(self) -> u.s:
        """The :term:`normalization` for time."""
        ...

    @property
    @abstractmethod
    def velocity(self) -> u.m / u.s:
        """The :term:`normalization` for velocity."""
        ...

    @property
    @abstractmethod
    def volumetric_heating_rate(self) -> u.J * u.m**-3 * u.s**-1:
        """The :term:`normalization` for volumetric heating rate."""
        ...

    @property
    @abstractmethod
    def volumetric_rate(self) -> u.m**-3 * u.s**-1:
        """The :term:`normalization` for a volumetric rate."""
        ...

    @property
    @abstractmethod
    def wavenumber(self) -> u.m**-1:
        """The :term:`normalization` for inverse length."""
        ...
