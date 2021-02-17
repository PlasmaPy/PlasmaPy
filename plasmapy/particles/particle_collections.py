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
    *particles : `particle_like`
        A series of particle-like objects.

    Examples
    --------
    >>> from plasmapy.particles import Particle, ParticleList
    >>> ParticleList("e-", "e+")
    ParticleList([Particle("e-"), Particle("e+")])

    """

    @staticmethod
    def _get_list_of_particle_like_instances(particles) -> List[AbstractParticle]:
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
        self._data = self._get_list_of_particle_like_instances(particles)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f"ParticleList({repr(self.data)})"

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
        """Append a particle to the end of the list."""
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
            "mass_energy", unit=u.J, default=np.nan * u.J
        )

    def sort(self):
        # TODO: Enable sorting if a key is provided.
        raise RuntimeError("Unable to sort a ParticleList.")

    @property
    def symbols(self) -> List[str]:
        """A list of the symbols of the particles."""
        return self._get_particle_attribute("symbol")
