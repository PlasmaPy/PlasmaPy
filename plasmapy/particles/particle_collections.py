"""Collections of `~plasmapy.particles.particle_class.Particle` objects."""

__all__ = ["ParticleList"]

import astropy.units as u
import collections
import numpy as np

from typing import *

from plasmapy.particles.exceptions import *
from plasmapy.particles.particle_class import (
    CustomParticle,
    DimensionlessParticle,
    Particle,
    ParticleLike,
)


class ParticleList(collections.UserList):
    """
    A `list` like collection of `Particle` and/or `CustomParticle` objects.

    Parameters
    ----------
    particles : iterable
        An iterable that provides a sequence of `ParticleLike` objects.
        Objects that are not a `Particle` or `CustomParticle` instance
        will be cast into a `Particle` instance.

    Raises
    ------
    `~plasmapy.particles.exceptions.InvalidParticleError`
        If an object supplied to `ParticleList` is not `ParticleLike`.

    TypeError
        If a `DimensionlessParticle` is provided.

    Examples
    --------
    A `ParticleList` can be created by calling it with a `list`, `tuple`,
    or other iterable that provides `~plasmapy.particles.ParticleLike`
    objects.

    >>> from plasmapy.particles import ParticleList
    >>> particle_list = ParticleList(["e-", "e+"])
    >>> particle_list[0]
    Particle("e-")

    Attributes such as `~plasmapy.particles.ParticleList.mass`
    and `~plasmapy.particles.ParticleList.charge` will return a
    `~astropy.units.Quantity` array containing the values of the
    corresponding attribute for each particle in the `ParticleList`.

    >>> particle_list.mass
    <Quantity [9.1093...e-31, 9.1093...e-31] kg>
    >>> particle_list.charge
    <Quantity [-1.60217663e-19,  1.60217663e-19] C>
    >>> particle_list.symbols
    ['e-', 'e+']

    `ParticleList` instances can also be created through addition and
    multiplication with `Particle`, `CustomParticle`, and `ParticleList`
    instances.

    >>> from plasmapy.particles import Particle, CustomParticle
    >>> import astropy.units as u
    >>> proton = Particle("p+")
    >>> custom_particle = CustomParticle(mass=1e-26*u.kg, charge=6e-19*u.C)
    >>> 2 * proton + custom_particle
    ParticleList(['p+', 'p+', 'CustomParticle(mass=1e-26 kg, charge=6e-19 C)'])

    These operations may also be performed using `ParticleLike` objects.

    >>> particle_list + "deuteron"
    ParticleList(['e-', 'e+', 'D 1+'])

    Normal `list` methods may also be used on `ParticleList` objects.
    When a `ParticleLike` object is appended to a `ParticleList`, that
    object will be cast into a `Particle`.

    >>> noble_gases = ParticleList(["He", "Ar", "Kr", "Xe", "Rn"])
    >>> noble_gases.append("Og")
    >>> noble_gases[-1]
    Particle("Og")

    The ``>`` operator may be used with `Particle` and `ParticleList`
    instances to access the nuclear reaction energy.

    >>> reactants = ParticleList(["deuterium", "tritium"])
    >>> products = ParticleList(["alpha", "neutron"])
    >>> energy = reactants > products
    >>> energy.to("MeV")
    <Quantity 17.58925... MeV>
    """

    @staticmethod
    def _list_of_particles_and_custom_particles(
        particles: Iterable[ParticleLike],
    ) -> List[Union[Particle, CustomParticle]]:  # TODO #687
        """
        Convert an iterable that provides `ParticleLike` objects into a
        `list` containing `Particle` and `CustomParticle` instances.
        """
        new_particles = []
        for obj in particles:
            if isinstance(obj, (Particle, CustomParticle)):
                new_particles.append(obj)
            elif isinstance(obj, DimensionlessParticle):
                raise TypeError(
                    "ParticleList instances cannot include dimensionless particles."
                )
            else:
                try:
                    new_particles.append(Particle(obj))
                except (TypeError, InvalidParticleError) as exc:
                    raise InvalidParticleError(
                        f"The object {obj} supplied to ParticleList is not a "
                        f"particle-like object."
                    ) from exc

        return new_particles

    def __init__(self, particles: Iterable):
        self._data = self._list_of_particles_and_custom_particles(particles)

    @staticmethod
    def _cast_other_as_particle_list(other):
        if isinstance(other, ParticleList):
            return other

        try:
            return ParticleList(other)
        except (InvalidParticleError, TypeError):
            pass

        try:
            return ParticleList([other])
        except (InvalidParticleError, TypeError):
            raise InvalidParticleError(f"Cannot cast {other} into a ParticleList")

    def __add__(self, other):
        try:
            other_as_particle_list = self._cast_other_as_particle_list(other)
        except (TypeError, InvalidParticleError) as exc:
            raise InvalidParticleError(
                f"Cannot add {repr(other)} to a ParticleList."
            ) from exc
        return ParticleList(self.data + other_as_particle_list.data)

    def __radd__(self, other):
        other_as_particle_list = self._cast_other_as_particle_list(other)
        return other_as_particle_list.__add__(self)

    def __repr__(self):
        return f"ParticleList({repr(self.symbols)})"

    def __gt__(self, other):
        from plasmapy.particles.nuclear import nuclear_reaction_energy

        other_as_particle_list = self._cast_other_as_particle_list(other)
        return nuclear_reaction_energy(
            reactants=self.symbols, products=other_as_particle_list.symbols
        )

    def __str__(self):
        return self.__repr__()

    def _get_particle_attribute(self, attr, unit=None, default=None):
        """
        Get the values of a particular attribute from all of the particles.

        If a ``unit`` is provided, then this function will return a
        `~astropy.units.Quantity` array with that unit.
        """
        values = [getattr(particle, attr, default) for particle in self.data]
        if unit:
            values = u.Quantity(values)
        return values

    def append(self, particle: ParticleLike):
        """Append a particle to the end of the `ParticleList`."""
        # TODO: use particle_input when it works with CustomParticle and ParticleLike
        if not isinstance(particle, (Particle, CustomParticle)):
            particle = Particle(particle)
        self.data.append(particle)

    @property
    def charge(self) -> u.C:
        """
        A `~astropy.units.Quantity` array of the electric charges
        of the particles.
        """
        return self._get_particle_attribute("charge", unit=u.C, default=np.nan * u.C)

    @property
    def data(self) -> List[Union[Particle, CustomParticle]]:
        """
        A `list` containing the particles contained in the
        `ParticleList` instance.

        The `~plasmapy.particles.ParticleList.data` attribute should not
        be modified directly.
        """
        return self._data

    def extend(self, iterable: Iterable[ParticleLike]):
        """
        Extend the sequence by appending `ParticleLike` elements
        from ``iterable``.
        """
        if isinstance(iterable, ParticleList):
            self.data.extend(iterable)
        else:
            for obj in iterable:
                self.append(obj)

    @property
    def half_life(self) -> u.s:
        """
        A `~astropy.units.Quantity` array of the half-lives of the
        particles.
        """
        return self._get_particle_attribute("half_life", unit=u.s, default=np.nan * u.s)

    def insert(self, index, particle: ParticleLike):
        """Insert a particle before an index."""
        # TODO: use particle_input when it works with CustomParticle and ParticleLike
        if not isinstance(particle, (Particle, CustomParticle)):
            particle = Particle(particle)
        self.data.insert(index, particle)

    @property
    def integer_charge(self) -> np.array:
        """
        An array of the quantized charges of the particles, as
        multiples of the elementary charge.
        """
        return np.array(self._get_particle_attribute("integer_charge", default=np.nan))

    @property
    def mass(self) -> u.kg:
        """A `~astropy.units.Quantity` array of the masses of the particles."""
        return self._get_particle_attribute("mass", unit=u.kg, default=np.nan * u.J)

    @property
    def mass_energy(self) -> u.J:
        """
        A `~astropy.units.Quantity` array of the mass energies of the
        particles.

        If the particle is an isotope or nuclide, return the mass energy
        of the nucleus only.
        """
        return self._get_particle_attribute(
            "mass_energy", unit=u.J, default=np.nan * u.J,
        )

    def sort(self, key: Callable = None, reverse: bool = False):
        """
        Sort the `ParticleList` in-place.

        For more information, refer to the documentation for `list.sort`.
        """
        if key is None:
            raise TypeError("Unable to sort a ParticleList without a key.")
        else:
            self._data.sort(key=key, reverse=reverse)

    @property
    def symbols(self) -> List[str]:
        """A `list` of the symbols of the particles."""
        return self._get_particle_attribute("symbol")


# Override the docstrings for the parent class

ParticleList.clear.__doc__ = """Remove all items from the `ParticleList`."""

ParticleList.copy.__doc__ = """Return a shallow copy of the `ParticleList`."""

ParticleList.count.__doc__ = """
Return the number of occurrences of ``item``.  Here, ``item`` may be a
`Particle`, `CustomParticle`, or `ParticleLike` representation of a particle.
"""

ParticleList.extend.__doc__ = """
Extend `ParticleList` by casting `ParticleLike` items from ``iterable``
into `Particle` or `CustomParticle` instances.
"""

ParticleList.index.__doc__ = """
Return first index of a `ParticleLike` value. Raise `ValueError` if
the value is not present.
"""

ParticleList.pop.__doc__ = """
Remove and return item at index (default last).  Raise `IndexError` if
the `ParticleList` is empty or the index is out of range.
"""

ParticleList.remove.__doc__ = """
Remove the first occurrence of a `ParticleLike` item.  Raise
`ValueError` if the value is not present.
"""

ParticleList.reverse.__doc__ = """Reverse the `ParticleList` in place."""
