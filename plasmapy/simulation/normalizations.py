"""Classes that describes normalizations for different systems of equations."""

import abc

from astropy import constants as const
from astropy import units as u
from typing import Optional

from plasmapy import formulary
from plasmapy.particles import Particle, particle_input
from plasmapy.utils.decorators import validate_quantities


class AbstractNormalizations(abc.ABC):
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
    def thermal_conductivity(self) -> :
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
        of collisions per cubic meter per second.
        """
        raise NotImplementedError

# Below is a first draft, which needs a lot of revising.

class MHDNormalizations(AbstractNormalizations):
    """
    Class to represent normalizations of the equations of ideal
    magnetohydrodynamics (MHD).

    Parameters
    ----------
    magnetic_field
        The normalization for the magnetic field.

    length
        The normalization for length.

    number_density
        The normalization for number density.

    ion : particle-like representation of an ion
        The ion species in the plasma.

    Raises
    ------

    Examples
    --------

    >>> import astropy.units as u
    >>> normalizations = IdealMHDNormalizations(
    ...    magnetic_field = 0.0250 * u.T,
    ...    length = 0.2 * u.m,
    ...    number_density = 1e19 * u.m ** -3,
    ...    ion = 'He-4 1+',
    ... )
    >>> normalizations.magnetic_field
    <Quantity 0.025 T>
    >>> normalizations.length
    <Quantity 0.2 m>
    >>> normalizations.number_density
    <Quantity 1.e+19 1 / m3>
    >>> normalizations.velocity
    <Quantity 86510.5392384 m / s>

    Notes
    -----
    <include a description of the math behind the normalizations here>

    """

    @particle_input
    @validate_quantities
    def __init__(
        self, magnetic_field: u.T, length: u.m, number_density: u.m ** -3, ion: Particle
    ):

        self._data = {}
        self.magnetic_field = magnetic_field
        self.length = length
        self.number_density = number_density
        self.ion = ion

    @property
    def ion(self) -> Optional[Particle]:
        """The ion species of the plasma."""
        return self._data["ion"]

    @ion.setter
    @particle_input
    def ion(self, ion: Particle = None):
        self._data["ion"] = ion

    @property
    def magnetic_field(self) -> u.T:
        """The normalization for the magnetic field."""
        return self._data["B"]

    @magnetic_field.setter
    def magnetic_field(self, B: u.T):
        self._data["B"] = B.to("T")

    @property
    def length(self) -> u.m:
        """The normalization for length."""
        return self._data["L"]

    @length.setter
    def length(self, L: u.m):
        self._data["L"] = L.to("m")

    @property
    def number_density(self) -> u.m ** -3:
        """The normalization for number density."""
        return self._data["n"]

    @number_density.setter
    def number_density(self, n):
        self._data["n"] = n.to(u.m ** -3)

    @property
    def mass_density(self) -> u.kg * u.m ** -3:
        """The normalization for mass density."""
        return self.ion.mass * self.number_density

    @property
    def time(self) -> u.s:
        """The normalization for time."""
        raise self.length / self.velocity

    @property
    def velocity(self) -> u.m / u.s:
        """
        The normalization for velocity.

        This normalization is equal to the Alfvén velocity calculated
        using the normalizations for the magnetic field and mass density.
        """
        return formulary.Alfven_speed(B=self.magnetic_field, density=self.mass_density)

    @property
    def current_density(self) -> u.A / u.m ** 2:
        """The normalization for current density."""
        raise self.magnetic_field / (const.mu0 * self.length)

    @property
    def pressure(self) -> u.Pa:
        """
        The normalization for pressure.

        This normalization is equal to two times the magnetic pressure
        calculated using the normalization for the magnetic field.
        """
        return 2 * formulary.magnetic_pressure(self.magnetic_field)

    @property
    def electric_field(self) -> u.V / u.m:
        """The normalization for electric field."""
        return (self.velocity * self.magnetic_field).to(u.V / u.m)

    @property
    def temperature(self) -> u.K:
        """The normalization for temperature."""
        return self.magnetic_field ** 2 / (const.k_B * const.mu0 * self.number_density)

    @property
    def magnetic_flux(self) -> u.T * u.m:
        """The normalization for magnetic flux."""
        return self.magnetic_field * self.length

    @property
    def resistivity(self) -> u.ohm * u.m:
        """The normalization for the resistivity."""
        return const.mu0 * self.length * self.velocity

    @property
    def diffusivity(self) -> u.m ** 2 / u.s:  #
        """The normalization for a diffusivity."""
        return self.length ** 2 / self.time

    @property
    def kinematic_viscosity(self) -> u.m ** 2 / u.s:
        """
        The normalization for kinematic viscosity.

        This normalization is equal to the dynamic viscosity divided by
        the mass density normalization.
        """
        # This is the same as the diffusivity normalization, so is it
        # worth having both in there?  If only the diffusivity normalization
        # is kept, then the dynamic viscosity docstring should say that
        # the kinematic viscosity normalization equals the diffusivity
        # normalization.  I put this in here explicitly since that will
        # make the intention of future code more clear.
        return self.diffusivity

    @property
    def dynamic_viscosity(self) -> u.Pa * u.s:
        """
        The normalization for dynamic viscosity.

        This normalization is equal to the product of the normalizations
        for time and pressure.
        """
        return self.pressure * self.time

    @property
    def thermal_conductivity(self):
        raise NotImplementedError


def normalization_from_unit(normalizations: AbstractNormalizations, unit: u.UnitBase) -> u.Quantity:
    """
    First draft!  Takes a (probably composite) unit, decomposes it
    to the base units used for the normalization, gets the normalizations
    associated with the relevant base units, takes it to the appropriate
    power, and then multiplies them together.

    """

    unit_decomposition = unit.si.decompose(bases=self._base_units)
    normalization = 1.0 * u.dimensionless_unscaled
    for base_unit, power in unit_decomposition.bases, unit_decomposition.powers:
        normalization_for_base = self._base_units_to_physical_type[base_unit]
        normalization *= normalization_for_base ** power
    return normalization.to(unit)




# Behavior to implement in the future: if we have a classes called, say,
# UniformPlasma and DimensionlessPlasma, then we should be able to do
# operations like:
#
#   uniform_plasma / normalizations -> dimensionless_plasma
#
#   dimensionless_plasma * normalizations = uniform_plasma
#
# More generally:
#
#   dimensional_plasma / normalizations -> dimensionless_plasma
#
#   dimensionless_plasma * normalizations -> dimensional_plasma

_unit_for_physical_type = {
    u.T: "magnetic_field",
    u.m: "length",
    u.s: "time",
    u.K: "temperature",
    u.J: "energy",
    u.m / u.s: "velocity",
    u.A / u.m ** 2: "current_density",
    u.Pa: "pressure",
    u.V / u.m: "electric_field",
    u.m ** -3: "number_density",
}

_irreducible_units_to_physical_type = {
    u.A: "current",
    u.K: "temperature",
    u.kg: "mass",
    u.m: "length",
    u.s: "time",
}  # probably not needed

_base_units_to_physical_type = {
    u.kg: "mass",
    u.m: "length",
    u.s: "time",
    u.T: "magnetic_field",
    u.K: "temperature",
}  # to use with .decompose(bases=...)

# Maybe make the base units the ones associated with the normalizations
# passed to the Normalization class upon instantiation.

_base_units = [u.kg, u.K, u.m, u.s, u.T]
