"""Collections of `~plasmapy.particles.particle_class.Particle` objects."""

__all__ = ["ParticleList"]

import astropy.units as u
import numpy as np
from typing import *
from plasmapy.particles.decorators import particle_input
from plasmapy.particles.particle_class import Particle, AbstractParticle, particle_like, CustomParticle, DimensionlessParticle
from plasmapy.particles.exceptions import *
import collections


class ParticleList(collections.UserList):

    def __init__(self, *particles):
        if len(particles) == 1 and isinstance(particles[0], (tuple, list)):
            particles = particles[0]
        self.data = list()
        for particle in particles:
            if isinstance(particle, DimensionlessParticle):
                raise (
                    f"ParticleList instances cannot include dimensionless particles."
                )
            if not isinstance(particle, (Particle, CustomParticle)):
                particle = Particle(particle)
            self.data.append(particle)

    @particle_input
    def append(self, particle: Particle):
        self.data.append(particle)

    @particle_input
    def insert(self, index, particle: Particle):
        self.data.insert(index, particle)




