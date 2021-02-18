"""Collections of `~plasmapy.particles.particle_class.Particle` objects."""

__all__ = ["ParticleList"]

import astropy.units as u
import collections
import numpy as np

from typing import *

from plasmapy.particles.decorators import particle_input
from plasmapy.particles.exceptions import *
from plasmapy.particles.particle_class import (
    AbstractParticle,
    CustomParticle,
    DimensionlessParticle,
    Particle,
    ParticleLike,
)


class ParticleList(collections.UserList):
    """
    A list-like collection of `Particle` and/or `CustomParticle` objects.

    Parameters
    ----------
    particles : iterable
        An iterable that provides a sequence of `ParticleLike` objects.

    Examples
    --------
    >>> from plasmapy.particles import ParticleList
    >>> particle_list = ParticleList(["e-", "e+"])
    >>> particle_list[0]
    Particle("e-")
    >>> particle_list.mass
    <Quantity [9.1093...e-31, 9.1093...e-31] kg>
    >>> particle_list.charge
    <Quantity [-1.60217663e-19,  1.60217663e-19] C>
    >>> particle_list.symbols
    ['e-', 'e+']
    >>> particle_list + [Particle("p+")]
    ParticleList(['e-', 'e+', 'p+'])
    """

    @staticmethod
    def _list_of_particles_and_custom_particles(
        particles: Iterable[ParticleLike],
    ) -> List[AbstractParticle]:
        """"""

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

    def __repr__(self):
        return f"ParticleList({repr(self.symbols)})"

    def __str__(self):
        return str(self.data)

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

    @particle_input
    def append(self, particle: Particle):
        """Append a particle to the end of the `ParticleList`."""
        self.data.append(particle)

    @property
    def charge(self) -> u.C:
        """An array of the electric charges of the particles."""
        return self._get_particle_attribute("charge", unit=u.C, default=np.nan * u.C)

    @property
    def data(self) -> List:
        """
        A regular `list` containing the particles contained in the
        `ParticleList` instance.

        The ``data`` attribute should not be modified directly.
        """
        return self._data

    def extend(self, other):
        if isinstance(other, ParticleList):
            self.data.extend(other)
        else:
            for obj in other:
                self.append(obj)

    @property
    def half_life(self) -> u.s:
        """An array of the half-lives of the particles."""
        return self._get_particle_attribute("half_life", unit=u.s, default=np.nan * u.s)

    @particle_input
    def insert(self, index, particle: Particle):
        """Insert a particle before an index."""
        self.data.insert(index, particle)

    @property
    def integer_charge(self) -> np.array:
        """
        An array of the quantized charges of the particles, as multiples
        of the elementary charge.
        """
        return np.array(self._get_particle_attribute("integer_charge", default=np.nan))

    @property
    def mass(self) -> u.kg:
        """An array of the masses of the particles."""
        return self._get_particle_attribute("mass", unit=u.kg, default=np.nan * u.J)

    @property
    def mass_energy(self) -> u.J:
        """
        An array of the mass energies of the particles in joules.

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
            self._data.sort(key=key, reverse=False)

    @property
    def symbols(self) -> List[str]:
        """A list of the symbols of the particles."""
        return self._get_particle_attribute("symbol")
