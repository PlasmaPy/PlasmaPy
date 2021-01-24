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
    particle_like,
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
        Delete the cache of attribute values that were stored to prevent
        them from needing to be recalculated.
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
                        f"The object {obj} supplied to ParticleList is not particle-like."
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

    def _get_particle_attribute(self, attr, unit=None, expected_type=None):

        if attr in self._attribute_values:
            return self._attribute_values[attr]

        if unit:
            default_value = np.nan * unit
        elif expected_type is str:
            default_value = ""
        elif expected_type in (float, int):
            default_value = np.nan
        else:
            default_value = None

        values = []
        for particle in self.data:
            try:
                values.append(getattr(particle, attr))
            except Exception:
                values.append(default_value)

        if unit:
            values = u.Quantity(values)
        elif expected_type in (int, float):
            values = np.array(values)

        self._attribute_values[attr] = values

        return values

    @property
    def antiparticle(self):
        """A `ParticleList` with the antiparticles of the particles."""
        return ParticleList(self._get_particle_attribute("antiparticle"))

    @particle_input
    def append(self, particle: Particle):
        """Append a particle to the end of the list."""
        self.data.append(particle)
        self._reset_attribute_values()

    @property
    def baryon_number(self) -> List[int]:
        """An array of the baryon number of the particles."""
        return self._get_particle_attribute("baryon_number", expected_type=int)

    @property
    def binding_energy(self) -> u.J:
        """An array of the binding energies of the particles."""
        return self._get_particle_attribute("binding_energy", unit=u.J)

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

    def is_category(self) -> list:
        raise NotImplementedError

    def isotopic_abundance(self) -> float:
        """An array of the isotopic abundances of the particles."""
        return self._get_particle_attribute("isotopic_abundance", expected_type=float)

    @property
    def lepton_number(self) -> list:
        """An array of the lepton numbers of the particles."""
        return self._get_particle_attribute("lepton_number")

    @property
    def mass(self) -> u.kg:
        """An array of the masses of the particles."""
        return self._get_particle_attribute("mass", unit=u.kg)

    @property
    def mass_number(self) -> list:
        """An array of the mass numbers of the particles."""
        return self._get_particle_attribute("mass_number")

    @property
    def nuclide_mass(self) -> u.kg:
        """An array of the nuclide masses of the particles."""
        return self._get_particle_attribute("nuclide_mass", u.kg)

    @property
    def roman_symbol(self) -> List[str]:
        """A list containing the particle symbol in Roman notation."""
        return self._get_particle_attribute("roman_symbol", expected_type=str)

    @property
    def sort(self):
        raise NotImplementedError

    @property
    def spin(self) -> np.ndarray:
        """An array containing the spins of the particles."""
        return self._get_particle_attribute("spin")

    @property
    def standard_atomic_weight(self):
        """
        An array containing the standard atomic weights of the
        particles.
        """
        return self._get_particle_attribute("standard_atomic_weight")
