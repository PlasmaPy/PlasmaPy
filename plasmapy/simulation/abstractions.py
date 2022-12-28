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
    An abstract base class to represent the |normalization constants|
    for systems of equations describing plasmas.

    For a list of physical types, see `astropy.units.physical`.

    .. warning::

        This interface is unstable and is subject to change.
    """

    @property
    @abstractmethod
    def acceleration(self) -> u.m / u.s**2:
        r"""The |normalization constant| for acceleration."""

    def area(self) -> u.m**2:
        r"""The |normalization constant| for area."""

    @property
    @abstractmethod
    def current_density(self) -> u.A / u.m**2:
        """
        The |normalization constant| for :wikipedia:`current
        density`.
        """

    @property
    @abstractmethod
    def diffusivity(self) -> u.m**2 / u.s:
        """
        The |normalization constant| for :wikipedia`diffusivity`.

        .. note::

           Use this normalization for :wikipedia:`kinematic viscosity`,
           :wikipedia:`magnetic diffusivity`, :wikipedia:`mass
           diffusivity`, and :wikipedia:`thermal diffusivity`.
        """

    @property
    @abstractmethod
    def dynamic_viscosity(self) -> u.Pa * u.s:
        """
        The |normalization constant| for :wikipedia:`dynamic
        viscosity`.

        .. note::

           For `kinematic viscosity`_, use
           ~`plasmapy.simulation.abstractions.AbstractNormalizations.diffusivity`.

        .. _kinematic viscosity: https://en.wikipedia.org/wiki/Viscosity#Kinematic_viscosity
        """

    @property
    @abstractmethod
    def electric_field(self) -> u.V / u.m:
        """
        The |normalization constant| for :wikipedia`electric
        field`.
        """

    @property
    @abstractmethod
    def energy(self) -> u.J:
        """The |normalization constant| for energy."""

    @property
    @abstractmethod
    def frequency(self) -> u.s**-1:
        """The |normalization constant| for frequency."""

    @property
    @abstractmethod
    def heat_flux(self) -> u.J * u.m**-2 * u.s**-1:
        """
        The |normalization constant| for :wikipedia:`heat flux`.

        .. note::

           For :wikipedia:`thermal diffusivity`, use the normalization
           constant for
           ~`plasmapy.simulation.abstractions.AbstractNormalizations.diffusivity`.
        """

    @property
    @abstractmethod
    def length(self) -> u.m:
        """The |normalization constant| for length."""

    @property
    @abstractmethod
    def magnetic_field(self) -> u.T:
        """The |normalization constant| for :wikipedia:`magnetic field`."""

    @property
    @abstractmethod
    def magnetic_flux(self) -> u.T * u.m:
        """
        The |normalization constant| for :wikipedia:`magnetic
        flux` or :wikipedia:`magnetic vector potential`.
        """

    @property
    @abstractmethod
    def mass(self) -> u.kg:
        """The |normalization constant| for mass."""

    @property
    @abstractmethod
    def mass_density(self) -> u.kg / u.m**3:
        """The |normalization constant| for mass density."""

    @property
    @abstractmethod
    def number_density(self) -> u.m**-3:
        """The |normalization constant| for :wikipedia:`number density`."""

    @property
    @abstractmethod
    def pressure(self) -> u.Pa:
        """The |normalization constant| for :wikipedia:`pressure`."""

    @property
    @abstractmethod
    def temperature(self) -> u.K:
        """
        The |normalization constant| for :term:`temperature` in
        units of kelvin.
        """

    @property
    @abstractmethod
    def thermal_conductivity(self) -> u.W / (u.K * u.m):
        """
        The |normalization constant| for :wikipedia:`thermal
        conductivity`.
        """

    @property
    @abstractmethod
    def time(self) -> u.s:
        """The |normalization constant| for time."""

    @property
    @abstractmethod
    def velocity(self) -> u.m / u.s:
        """The |normalization constant| for velocity."""

    @property
    @abstractmethod
    def volume(self) -> u.m**3:
        """The |normalization constant| for volume."""

    @property
    @abstractmethod
    def wavenumber(self) -> u.m**-1:
        """The |normalization constant| for :wikipedia:`wavenumber`."""
