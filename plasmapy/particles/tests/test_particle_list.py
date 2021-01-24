"""Tests for particle collections."""

import pytest

import astropy.units as u

from plasmapy.particles.particle_class import Particle, DimensionlessParticle, CustomParticle
from plasmapy.particles.particle_list import ParticleList
from plasmapy.particles.exceptions import *
from plasmapy.particles import electron, proton, alpha, neutron


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

attributes = [
    "antiparticle",
    "baryon_number",
    "binding_energy",
    "charge",
    "half_life",
    "isotopic_abundance",
    "lepton_number",
    "mass",
    "mass_number",
    "roman_symbol",
    "spin",
    "standard_atomic_weight",
]

@pytest.mark.parametrize("args", particle_list_arguments)
def test_particle_list_creation_membership(args):
    particle_list = ParticleList(*args)
    for arg, particle in zip(args, particle_list):
        assert particle == arg
        assert isinstance(particle, (Particle, CustomParticle))


@pytest.mark.parametrize("attribute", attributes)
def test_particle_list_attributes(attribute):
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






