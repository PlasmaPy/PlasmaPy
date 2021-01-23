"""Collections of `~plasmapy.particles.particle_class.Particle` objects."""

__all__ = ["ParticleList"]

import astropy.units as u
import numpy as np
from typing import *
from plasmapy.particles.particle_class import Particle, particle_like
from plasmapy.particles.exceptions import InvalidParticleError
import collections

class ParticleList(collections.UserList):

    pass









