"""Tests for particle collections."""

import astropy.units as u
import pytest

from plasmapy.particles import alpha, electron, neutron, proton
from plasmapy.particles.atomic import atomic_number
from plasmapy.particles.exceptions import InvalidParticleError
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
    ([electron]),
    ([electron, proton]),
    ([electron, proton, alpha]),
    (["e-", "e+"]),
    ([electron, "e-"]),
    ([custom_particle]),
    ([custom_particle, electron, "e-"]),
]

attributes = [
    "charge",
    "mass",
    "mass_energy",
    "half_life",
]


def _everything_is_particle_or_custom_particle(iterable):
    return all([isinstance(p, (Particle, CustomParticle)) for p in iterable])


@pytest.mark.parametrize("args", particle_list_arguments)
def test_particle_list_creation_membership(args):
    particle_list = ParticleList(args)
    for arg, particle in zip(args, particle_list):
        assert particle == arg
    assert _everything_is_particle_or_custom_particle(particle_list)
    assert _everything_is_particle_or_custom_particle(particle_list.data)


@pytest.mark.parametrize("attribute", attributes)
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
    assert actual == expected or u.allclose(actual, expected, equal_nan=True)


@pytest.mark.parametrize("attr", ["mass", "charge", "integer_charge", "symbols"])
def test_that_particle_list_attrs_cannot_be_redefined(attr):
    """
    Test that attributes of ParticleList cannot be manually redefined."""
    particle_list = ParticleList(["D+", "p+", "n"])
    with pytest.raises(expected_exception=AttributeError):
        setattr(particle_list, attr, 42)


def test_particle_list_len():
    original_list = ["n", "p", "e-"]
    particle_list = ParticleList(original_list)
    assert len(particle_list) == len(original_list)


@pytest.fixture
def various_particles():
    return ParticleList(
        [
            "H",
            "He",
            "e-",
            "alpha",
            "tau neutrino",
            CustomParticle(mass=3 * u.kg, charge=5 * u.C),
            CustomParticle(),
            CustomParticle(mass=7 * u.kg),
            CustomParticle(charge=11 * u.C),
        ]
    )


valid_particles = (1, CustomParticle(), Particle("Fe"))
invalid_particles = (0, "not a particle", DimensionlessParticle())


def test_append_particle_like(various_particles):
    original_length = len(various_particles)
    various_particles.append("Li")
    appended_item = various_particles[-1]
    assert len(various_particles) == original_length + 1
    assert isinstance(appended_item, Particle)
    assert appended_item == Particle("Li")
    assert various_particles.data[-1] is appended_item


def test_pop(various_particles):
    expected = various_particles[:-1]
    various_particles.pop()
    assert various_particles == expected


def test_extend(various_particles):
    new_particles = ["Fe", Particle("e-")]
    various_particles.extend(new_particles)
    assert all([isinstance(p, (Particle, CustomParticle)) for p in various_particles])


def test_instantiating_with_invalid_particle():
    list_with_a_nonparticle = ["H", -1, "e+"]
    with pytest.raises(InvalidParticleError):
        ParticleList(list_with_a_nonparticle)


@pytest.mark.parametrize("invalid_particle", invalid_particles)
def test_appending_invalid_particle(various_particles, invalid_particle):
    with pytest.raises((InvalidParticleError, TypeError)):
        various_particles.append(invalid_particle)


def test_particle_list_extend_with_invalid_particles(various_particles):
    with pytest.raises(InvalidParticleError):
        various_particles.extend(invalid_particles)


@pytest.mark.parametrize("invalid_particle", invalid_particles)
def test_particle_list_insert_invalid_particle(various_particles, invalid_particle):
    with pytest.raises((InvalidParticleError, TypeError)):
        various_particles.insert(1, invalid_particle)


def test_particle_list_sort_with_key_and_reverse():
    elements = ["He", "H", "Fe", "U"]
    particle_list = ParticleList(elements)
    particle_list.sort(key=atomic_number, reverse=True)
    assert particle_list.symbols == ["U", "Fe", "He", "H"]


def test_particle_list_sort_without_key(various_particles):
    with pytest.raises(TypeError):
        various_particles.sort()
