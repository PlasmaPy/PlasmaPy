"""Tests for particle collections."""
import astropy.units as u
import numpy as np
import pytest

from typing import Dict

from plasmapy.particles import alpha, electron, neutron, proton
from plasmapy.particles.atomic import atomic_number
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.particles.nuclear import nuclear_reaction_energy
from plasmapy.particles.particle_class import (
    CustomParticle,
    DimensionlessParticle,
    Particle,
    ParticleLike,
)
from plasmapy.particles.particle_collections import ParticleList

custom_particle = CustomParticle(mass=1e-25 * u.kg, charge=1e-18 * u.C)
dimensionless_particle = DimensionlessParticle(mass=1.25, charge=1.58)


attributes = [
    "charge",
    "half_life",
    "charge_number",
    "mass",
    "mass_energy",
]


@pytest.fixture
def various_particles():
    """A sample `ParticleList` with several different valid particles."""
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


def _everything_is_particle_or_custom_particle(iterable):
    """
    Test that every object in an iterable is either a `Particle` instance
    or a `CustomParticle` instance.
    """
    return all(isinstance(p, (Particle, CustomParticle)) for p in iterable)


@pytest.mark.parametrize(
    "args",
    [
        (),
        ([electron]),
        ([electron, proton]),
        ([electron, proton, alpha]),
        (["e-", "e+"]),
        ([electron, "e-"]),
        ([custom_particle]),
        ([custom_particle, electron, "e-"]),
    ],
)
def test_particle_list_membership(args):
    """
    Test that the particles in the `ParticleList` match the particles
    (or particle-like objects) that are passed to it.
    """
    particle_list = ParticleList(args)
    for arg, particle in zip(args, particle_list):
        assert particle == arg
    assert _everything_is_particle_or_custom_particle(particle_list)
    assert _everything_is_particle_or_custom_particle(particle_list.data)


@pytest.mark.parametrize("attribute", attributes)
def test_particle_list_attributes(attribute, various_particles):
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


@pytest.mark.parametrize("attribute", attributes)
def test_particle_list_no_redefining_attributes(various_particles, attribute):
    """
    Test that attributes of `ParticleList` cannot be manually redefined.

    This test may fail if `@cached_property` is used instead of `@property`
    because `@cached_property` allows reassignment while `@property` does
    not.
    """
    with pytest.raises(AttributeError):
        various_particles.__setattr__(attribute, 42)


def test_particle_list_len():
    """Test that using `len` on a `ParticleList` returns the expected number."""
    original_list = ["n", "p", "e-"]
    particle_list = ParticleList(original_list)
    assert len(particle_list) == len(original_list)


def test_particle_list_append(various_particles):
    """Test that a particle-like object can get appended to a `ParticleList`."""
    original_length = len(various_particles)
    various_particles.append("Li")
    appended_item = various_particles[-1]
    assert len(various_particles) == original_length + 1
    assert isinstance(appended_item, Particle)
    assert appended_item == Particle("Li")
    assert various_particles.data[-1] is appended_item


def test_particle_list_pop(various_particles):
    """Test that the last item in the `ParticleList` is removed when"""
    expected = various_particles[:-1]
    particle_that_should_be_removed = various_particles[-1]
    removed_particle = various_particles.pop()
    assert various_particles == expected
    assert removed_particle == particle_that_should_be_removed


def test_particle_list_extend(various_particles):
    """
    Test that a `ParticleList` can be extended when provided with an
    iterable that yields particle-like objects.
    """
    new_particles = ["Fe", Particle("e-"), CustomParticle()]
    various_particles.extend(new_particles)
    assert _everything_is_particle_or_custom_particle(various_particles)
    assert various_particles[-3:] == new_particles


def test_particle_list_extended_with_particle_list(various_particles):
    """Test that a `ParticleList` can be extended with another `ParticleList`."""
    particle_list = ParticleList(["D", "T", CustomParticle()])
    various_particles.extend(particle_list)
    assert various_particles[-3:] == particle_list


def test_particle_list_insert(various_particles):
    """Test insertion of particle-like objects into"""
    various_particles.insert(0, "tau neutrino")
    assert various_particles[0] == "tau neutrino"
    assert _everything_is_particle_or_custom_particle(various_particles)


invalid_particles = (0, "not a particle", DimensionlessParticle())


def test_particle_list_instantiate_with_invalid_particles():
    """
    Test that a `ParticleList` instance cannot be created when it is
    provided with invalid particles.
    """
    with pytest.raises(InvalidParticleError):
        ParticleList(invalid_particles)


@pytest.mark.parametrize("invalid_particle", invalid_particles)
def test_particle_list_append_invalid_particle(various_particles, invalid_particle):
    """
    Test that objects that are not particle-like cannot be appended to
    a `ParticleList` instance.
    """
    with pytest.raises((InvalidParticleError, TypeError)):
        various_particles.append(invalid_particle)


def test_particle_list_extend_with_invalid_particles(various_particles):
    """
    Test that a `ParticleList` instance cannot be extended with any
    objects that are not particle-like.
    """
    with pytest.raises(InvalidParticleError):
        various_particles.extend(invalid_particles)


@pytest.mark.parametrize("invalid_particle", invalid_particles)
def test_particle_list_insert_invalid_particle(various_particles, invalid_particle):
    """
    Test that objects that are not particle-like cannot be inserted into
    a `ParticleList` instance.
    """
    with pytest.raises((InvalidParticleError, TypeError)):
        various_particles.insert(1, invalid_particle)


def test_particle_list_sort_with_key_and_reverse():
    """
    Test that a `ParticleList` instance can be sorted if a key is
    provided, and that the ``reverse`` keyword argument works too.
    """
    elements = ["He", "H", "Fe", "U"]
    particle_list = ParticleList(elements)
    particle_list.sort(key=atomic_number, reverse=True)
    assert particle_list.symbols == ["U", "Fe", "He", "H"]


def test_particle_list_sort_without_key(various_particles):
    """Test that a `ParticleList` cannot be sorted if a key is not provided."""
    with pytest.raises(TypeError):
        various_particles.sort()


def test_particle_list_dimensionless_particles():
    """
    Test that a `ParticleList` cannot be instantiated with a
    `DimensionlessParticle`.
    """
    with pytest.raises(TypeError):
        ParticleList([DimensionlessParticle()])


def test_particle_list_adding_particle_list(various_particles):
    """Test that a `ParticleList` can be added to another `ParticleList`."""
    extra_particles = ParticleList(["H", "D", "T"])
    new_particles_list = various_particles + extra_particles
    assert new_particles_list[-3:] == extra_particles
    assert isinstance(new_particles_list, ParticleList)


def test_add_particle_list_and_particle(various_particles):
    """
    Test that a `ParticleList` can be added to a `Particle` on the right
    and then return a `ParticleList`.
    """
    new_particle_list = various_particles + electron
    assert new_particle_list[-1] == electron
    assert new_particle_list[:-1] == various_particles
    assert isinstance(new_particle_list, ParticleList)


def test_add_particle_and_particle_list(various_particles):
    """
    Test that a `Particle` can be added to a `ParticleList` on the right
    and then return a `ParticleList`.
    """
    new_particle_list = electron + various_particles
    assert new_particle_list[0] == electron
    assert new_particle_list[1:] == various_particles
    assert isinstance(new_particle_list, ParticleList)


def test_add_particle_and_particle_like():
    """
    Test that a `Particle` can be added to a particle-like object on the
    right and then return a `ParticleList`.
    """
    heavy_isotopes_of_hydrogen = Particle("D") + "T"
    assert isinstance(heavy_isotopes_of_hydrogen, ParticleList)
    assert heavy_isotopes_of_hydrogen[0] == "D"
    assert heavy_isotopes_of_hydrogen[1] == "T"


def test_add_particle_like_and_particle():
    """
    Test that a particle-like object on the left can be added to a
    `Particle` instance on the right and then return a `ParticleList`.
    """
    heavy_isotopes_of_hydrogen = "D" + Particle("T")
    assert isinstance(heavy_isotopes_of_hydrogen, ParticleList)
    assert heavy_isotopes_of_hydrogen[0] == "D"
    assert heavy_isotopes_of_hydrogen[1] == "T"


def test_particle_list_gt_as_nuclear_reaction_energy():
    """
    Test that `ParticleList.__gt__` can be used to get the same result
    as `nuclear_reaction_energy`.
    """
    reactants = ParticleList(["D+", "T+"])
    products = ParticleList(["alpha", "n"])
    expected_energy = nuclear_reaction_energy("D + T --> alpha + n")
    actual_energy = reactants > products
    assert u.allclose(expected_energy, actual_energy)


def test_particle_gt_as_radioactive_decay():
    """
    Test a nuclear reaction where a `Particle` instance is the sole
    reactant on the left side.
    """
    tritium = Particle("T")
    expected_energy = nuclear_reaction_energy("T -> He-3 + e")
    actual_energy = tritium > Particle("He-3") + "e"
    assert u.allclose(expected_energy, actual_energy)


@pytest.mark.parametrize("method", ["append", "__add__", "__radd__"])
@pytest.mark.parametrize("invalid_particle", invalid_particles)
def test_particle_list_invalid_ops(various_particles, invalid_particle, method):
    """
    Test that operations with invalid particles raise the appropriate
    exceptions.
    """
    with pytest.raises((InvalidParticleError, TypeError)):
        getattr(various_particles, method)(invalid_particle)


@pytest.mark.parametrize("particle", [electron, CustomParticle()])
@pytest.mark.parametrize("method", ["__mul__", "__rmul__"])
def test_particle_multiplication(method, particle):
    """
    Test that multiplying a `Particle` or `CustomParticle` with an
    integer returns a `ParticleList` that contains the correct values.
    """
    particle_list = getattr(particle, method)(3)
    assert particle_list == [particle, particle, particle]


@pytest.mark.parametrize(
    "particles, args, kwargs, expected",
    [
        [
            ["electron", "proton", "neutron"],
            ["lepton"],
            {},
            [True, False, False],
        ],
        [
            ["electron", "proton", "neutron"],
            [],
            {"require": "lepton"},
            [True, False, False],
        ],
        [
            ["electron", "proton", "neutron"],
            [],
            {"exclude": "lepton"},
            [False, True, True],
        ],
        [
            ["electron", "proton", "neutron"],
            [],
            {"any_of": {"lepton", "charged"}},
            [True, True, False],
        ],
    ],
)
def test_particle_list_is_category(particles, args, kwargs, expected):
    """
    Test that ``ParticleList.is_category()`` behaves as expected.
    """
    sample_list = ParticleList(particles)
    assert sample_list.is_category(*args, **kwargs) == expected


def test_mean_particle():
    """
    Test that ``ParticleList.average_particle()`` returns a particle with
    the mean mass and mean charge of a |ParticleList|.
    """
    massless_uncharged_particle = CustomParticle(mass=0 * u.kg, charge=0 * u.C)
    particle_list = ParticleList([proton, electron, alpha, massless_uncharged_particle])
    expected_mass = (proton.mass + electron.mass + alpha.mass) / 4
    expected_charge = (proton.charge + electron.charge + alpha.charge) / 4
    average_particle = particle_list.average_particle()
    assert u.isclose(average_particle.mass, expected_mass, rtol=1e-14)
    assert u.isclose(average_particle.charge, expected_charge, rtol=1e-14)


def test_weighted_mean_particle():
    """
    Test that ``ParticleList.average_particle()`` returns a particle with
    the weighted mean.
    """
    custom_proton = CustomParticle(mass=proton.mass, charge=proton.charge)
    particle_list = ParticleList([proton, electron, alpha, custom_proton])
    abundances = [1, 2, 0, 1]
    expected_mass = (proton.mass + electron.mass) / 2
    expected_charge = 0 * u.C
    average_particle = particle_list.average_particle(abundances=abundances)
    assert u.isclose(average_particle.mass, expected_mass, rtol=1e-14)
    assert u.isclose(average_particle.charge, expected_charge, rtol=1e-14)


boolean_pairs = [(False, False), (True, False), (False, True), (True, True)]


@pytest.mark.parametrize("use_rms_charge, use_rms_mass", boolean_pairs)
def test_root_mean_square_particle(use_rms_charge, use_rms_mass):
    """
    Test that ``ParticleList.average_particle`` returns the mean or root
    mean square of the charge and mass, as appropriate.
    """

    particle_list = ParticleList(["p+", "e-"])
    average_particle = particle_list.average_particle(
        use_rms_charge=use_rms_charge, use_rms_mass=use_rms_mass
    )

    expected_average_charge = (1 if use_rms_charge else 0) * proton.charge
    assert u.isclose(average_particle.charge, expected_average_charge, rtol=1e-14)

    if use_rms_mass:
        expected_average_mass = np.sqrt((proton.mass**2 + electron.mass**2) / 2)
    else:
        expected_average_mass = (proton.mass + electron.mass) / 2

    assert u.isclose(average_particle.mass, expected_average_mass, atol=1e-35 * u.kg)


particle_multiplicities = [
    {"e-": 1},
    {"p+": 5, "e-": 11, "Fe-56 5+": 2},
    {"p+": 4},
    {CustomParticle(mass=1 * u.kg, charge=1 * u.C): 1},
    {"p+": 5, CustomParticle(mass=1 * u.kg, charge=1 * u.C): 1},
    {CustomParticle(): 1},
    {"p+": 2, "p+": 1},
]


@pytest.mark.parametrize("particle_multiplicities", particle_multiplicities)
@pytest.mark.parametrize("use_rms_charge, use_rms_mass", boolean_pairs)
def test_weighted_averages_of_particles(
    particle_multiplicities: Dict[ParticleLike, int],
    use_rms_charge,
    use_rms_mass,
):
    """
    Compare the mass and charge of the average particle for two |ParticleList|
    instances.

    The first |ParticleList| contains repeated particles.

    The second |ParticleList| contains only one of each kind of particle
    present in the first list, with the number of each particle recorded
    in a separate array.

    The unweighted averages of the first |ParticleList| should equal the
    weighted averages of the second |ParticleList|, with the number of
    each particle provided as the abundances.
    """
    all_particles = ParticleList([])
    for particle, multiplicity in particle_multiplicities.items():
        all_particles.extend(ParticleList(multiplicity * [particle]))

    unique_particles = ParticleList(particle_multiplicities.keys())
    number_of_each_particle = list(particle_multiplicities.values())

    unweighted_mean_of_all_particles = all_particles.average_particle(
        use_rms_charge=use_rms_charge,
        use_rms_mass=use_rms_mass,
    )

    weighted_mean_of_unique_particles = unique_particles.average_particle(
        use_rms_charge=use_rms_charge,
        use_rms_mass=use_rms_mass,
        abundances=number_of_each_particle,
    )

    assert u.isclose(
        unweighted_mean_of_all_particles.mass,
        weighted_mean_of_unique_particles.mass,
        rtol=1e-14,
        equal_nan=True,
    )

    assert u.isclose(
        unweighted_mean_of_all_particles.charge,
        weighted_mean_of_unique_particles.charge,
        rtol=1e-14,
        equal_nan=True,
    )

    if len(unique_particles) == 1 and isinstance(unique_particles[0], Particle):
        assert isinstance(unweighted_mean_of_all_particles, Particle)
        assert isinstance(weighted_mean_of_unique_particles, Particle)


def test_particle_list_with_no_arguments():
    """Test that `ParticleList()` returns an empty `ParticleList`."""
    empty_particle_list = ParticleList()
    assert isinstance(empty_particle_list, ParticleList)
    assert len(empty_particle_list) == 0
