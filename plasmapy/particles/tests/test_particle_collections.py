"""Tests for particle collections."""

import astropy.units as u
import pytest

from plasmapy.particles import alpha, electron, neutron, proton
from plasmapy.particles.exceptions import *
from plasmapy.particles.particle_class import (
    CustomParticle,
    DimensionlessParticle,
    Particle,
)
from plasmapy.particles.particle_collections import ParticleList

custom_particle = CustomParticle(mass=1e-25 * u.kg, charge=1e-18 * u.C)
dimensionless_particle = DimensionlessParticle(mass=1.25, charge=1.58)

particle_list_arguments = [
    (),
    (electron,),
    (electron, proton),
    (electron, proton, alpha),
    ("e-", "e+"),
    (electron, "e-"),
    (custom_particle,),
    (custom_particle, electron, "e-"),
]

numerical_attributes = [
    "charge",
    "mass",
    "mass_energy",
]


@pytest.mark.parametrize("args", particle_list_arguments)
def test_particle_list_creation_membership(args):
    particle_list = ParticleList(*args)
    for arg, particle in zip(args, particle_list):
        assert particle == arg
        assert isinstance(particle, (Particle, CustomParticle))


@pytest.mark.parametrize("attribute", numerical_attributes)
def test_particle_list_quantity_attributes(attribute):
    """
    Test that the attributes of ParticleList correspond to the
    attributes of the listed particles.

    This class does not test ParticleList instances that include
    CustomParticle instances inside of them.
    """
    particle_list_arguments = (electron, "e+", proton, neutron, alpha)
    particle_list = ParticleList(particle_list_arguments)
    expected_particles = [Particle(arg) for arg in particle_list_arguments]
    actual = getattr(particle_list, attribute)
    expected = [getattr(particle, attribute) for particle in expected_particles]
    assert u.allclose(actual, expected, equal_nan=True)
