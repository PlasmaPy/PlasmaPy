"""Tests for particle collections."""

import pytest

import astropy.units as u

from plasmapy.particles.particle_class import Particle, DimensionlessParticle, CustomParticle
from plasmapy.particles.particle_list import ParticleList
from plasmapy.particles.exceptions import *
from plasmapy.particles import electron, proton, alpha


custom_particle = CustomParticle(mass=1e-25 * u.kg, charge=1e-18 * u.C)
dimensionless_particle = DimensionlessParticle(mass=1.25, charge=1.58)

particle_list_arguments = [
    (),
    (electron,),
    (electron, proton, alpha),
    ("e-", "e+"),
    (electron, "e-"),
    (custom_particle,),
    (custom_particle, electron, "e-"),
]

@pytest.mark.parametrize("args", particle_list_arguments)
def test_particle_list_creation_membership(args):
    particle_list = ParticleList(*args)
    for arg, particle in zip(args, particle_list):
        assert particle == arg
        assert isinstance(particle, (Particle, CustomParticle))


