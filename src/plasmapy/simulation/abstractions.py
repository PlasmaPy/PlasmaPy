"""
Abstract classes for numerical simulations.

.. attention::

   |expect-api-changes|
"""

__all__ = [
    "AbstractNormalizations",
    "AbstractSimulation",
    "AbstractTimeDependentSimulation",
]

from abc import ABC, abstractmethod

import astropy.units as u


class AbstractSimulation(ABC):
    """
    A prototype abstract interface for numerical simulations.

    .. attention::

       |expect-api-changes|
    """

    @abstractmethod
    def summarize(self):
        """
        Print a summary of the simulation parameters and status.
        """
        ...

    @abstractmethod
    def initialize(self) -> None:
        """Prepare the simulation to be run."""
        ...

    @abstractmethod
    def simulate(self):
        """Perform the simulation."""
        ...

    @abstractmethod
    def finalize(self) -> None:
        """Perform the steps to close the simulation."""
        ...


class AbstractTimeDependentSimulation(AbstractSimulation):
    """
    A prototype abstract interface for time-dependent numerical simulations.

    .. warning::
        This interface is unstable and is subject to change.
    """


class AbstractNormalizations(ABC):
    """
    An abstract base class to represent the normalization constants for
    systems of equations describing plasmas.

    .. warning::
        This interface is unstable and is subject to change.
    """

    @property
    @abstractmethod
    def current_density(self) -> u.Quantity[u.A / u.m**2]:  # ty: ignore[not-subscriptable]
        """The current density normalization."""
        ...

    @property
    @abstractmethod
    def diffusivity(self) -> u.Quantity[u.m**2 / u.s]:  # ty: ignore[not-subscriptable]
        """The normalization for diffusivity."""
        ...

    @property
    @abstractmethod
    def dynamic_viscosity(self) -> u.Quantity[u.Pa * u.s]:  # ty: ignore[not-subscriptable]
        """The normalization for dynamic viscosity."""
        ...

    @property
    @abstractmethod
    def electric_field(self) -> u.Quantity[u.V / u.m]:  # ty: ignore[not-subscriptable]
        """The electric field normalization."""
        ...

    @property
    @abstractmethod
    def heat_flux(self) -> u.Quantity[u.J * u.m**-2 * u.s**-1]:  # ty: ignore[not-subscriptable]
        """The normalization for heat flux."""
        ...

    @property
    @abstractmethod
    def length(self) -> u.Quantity[u.m]:  # ty: ignore[not-subscriptable]
        """The normalization for length."""
        ...

    @property
    @abstractmethod
    def magnetic_field(self) -> u.Quantity[u.T]:  # ty: ignore[not-subscriptable]
        """The magnetic field normalization."""
        ...

    @property
    @abstractmethod
    def magnetic_flux(self) -> u.Quantity[u.T * u.m]:  # ty: ignore[not-subscriptable]
        """The normalization for the magnetic flux or vector potential."""
        ...

    @property
    @abstractmethod
    def mass(self) -> u.Quantity[u.kg]:  # ty: ignore[not-subscriptable]
        """The normalization for mass."""
        ...

    @property
    @abstractmethod
    def mass_density(self) -> u.Quantity[u.kg / u.m**3]:  # ty: ignore[not-subscriptable]
        """The normalization for mass density."""
        ...

    @property
    @abstractmethod
    def number_density(self) -> u.Quantity[u.m**-3]:  # ty: ignore[not-subscriptable]
        """The normalization for number density."""
        ...

    @property
    @abstractmethod
    def pressure(self) -> u.Quantity[u.Pa]:  # ty: ignore[not-subscriptable]
        """The normalization for pressure."""
        ...

    @property
    @abstractmethod
    def temperature(self) -> u.Quantity[u.K]:  # ty: ignore[not-subscriptable]
        """The normalization for temperature."""
        ...

    @property
    @abstractmethod
    def thermal_conductivity(self) -> u.Quantity[u.W / u.m / u.K]:  # ty: ignore[not-subscriptable]
        """The normalization for thermal conductivity."""
        ...

    @property
    @abstractmethod
    def time(self) -> u.Quantity[u.s]:  # ty: ignore[not-subscriptable]
        """The normalization for time."""
        ...

    @property
    @abstractmethod
    def wavenumber(self) -> u.Quantity[u.m**-1]:  # ty: ignore[not-subscriptable]
        """The normalization for inverse length."""
        ...

    @property
    @abstractmethod
    def velocity(self) -> u.Quantity[u.m / u.s]:  # ty: ignore[not-subscriptable]
        """The normalization for velocity."""
        ...

    @property
    @abstractmethod
    def volumetric_heating_rate(self) -> u.Quantity[u.J * u.m**-3 * u.s**-1]:  # ty: ignore[not-subscriptable]
        """The normalization for the volumetric heating rate."""
        ...

    @property
    @abstractmethod
    def volumetric_rate(self) -> u.Quantity[u.m**-3 * u.s**-1]:  # ty: ignore[not-subscriptable]
        """
        The normalization for a volumetric rate.

        This normalization is applicable to, for example, the number
        of collisions per unit volume per unit time.
        """
        ...
