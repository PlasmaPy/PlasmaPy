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
)


class ParticleList(collections.UserList):
    """
    A list-like collection of `Particle` and/or `CustomParticle` objects.

    Parameters
    ----------
    *particles : `particle-like`
        A series of particle-like objects.

    Examples
    --------
    >>> from plasmapy.particles import Particle, ParticleList
    >>> ParticleList("e-", "e+")
    ParticleList([Particle("e-"), Particle("e+")])

    """

    def _reset_attribute_values(self) -> NoReturn:
        """
        Delete the cache of attribute values that were stored.


        """
        self._attribute_values = dict()

    @staticmethod
    def _turn_into_particles(particles) -> List[AbstractParticle]:
        """"""
        if len(particles) == 1 and isinstance(particles[0], (tuple, list)):
            particles = particles[0]

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

    def __init__(self, *particles):
        self.data = self._turn_into_particles(particles)
        self._reset_attribute_values()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f"ParticleList({repr(self.data)})"

    def _get_particle_attribute(self, attr, unit=None, default=None):
        """
        Get the values of a particular attribute from all of the particles.

        If a ``unit`` is provided, then this function will return a
        `~astropy.units.Quantity` in that unit.
        """

        if attr in self._attribute_values:
            return self._attribute_values[attr]

        values = []
        for particle in self.data:
            try:
                value = getattr(particle, attr)
            except Exception:
                value = default
            finally:
                values.append(value)

        if unit:
            values = u.Quantity(values)

        self._attribute_values[attr] = values

        return values

    @particle_input
    def append(self, particle: Particle):
        """Append a particle to the end of the list."""
        self.data.append(particle)
        self._reset_attribute_values()

    @property
    def charge(self) -> u.C:
        """An array of the electric charges of the particles."""
        return self._get_particle_attribute("charge", unit=u.C)

    @property
    def half_life(self) -> u.s:
        """An array of the half-lives of the particles."""
        return self._get_particle_attribute("half_life", unit=u.s)

    @particle_input
    def insert(self, index, particle: Particle):
        """Insert a particle before an index."""
        self.data.insert(index, particle)
        self._reset_attribute_values()

    @property
    def integer_charge(self):
        """
        An array of the quantized charges of the particles, as multiples
        of the elementary charge.
        """
        return self._get_particle_attribute("integer_charge", expected_type=float)

    @property
    def mass(self) -> u.kg:
        """An array of the masses of the particles."""
        return self._get_particle_attribute("mass", unit=u.kg)

    @property
    def mass_energy(self) -> u.J:
        """
        An array of the mass energies of the particles in joules.

        If the particle is an isotope or nuclide, return the mass energy
        of the nucleus only.
        """
        return self._get_particle_attribute("mass_energy", unit=u.J)

    @property
    def sort(self):
        raise RuntimeError("Unable to sort a ParticleList.")

    @property
    def symbols(self) -> List[str]:
        """A list of the symbols of the particles."""
        return self._get_particle_attribute("symbol")
