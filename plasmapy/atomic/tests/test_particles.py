import pytest

from ..particles import (
    _Particles,
    _special_particles,
    _leptons,
    _antileptons,
    _baryons,
    _antibaryons,
    _particles,
    _antiparticles,
    _fermions,
    _bosons,
    _neutrinos,
    _antineutrinos,
    _everything,
)

particle_antiparticle_pairs = [
    ('e-', 'e+'),
    ('mu-', 'mu+'),
    ('tau-', 'tau+'),
    ('p+', 'p-'),
    ('n', 'antineutron'),
    ('nu_e', 'anti_nu_e'),
    ('nu_mu', 'anti_nu_mu'),
    ('nu_tau', 'anti_nu_tau'),
]


@pytest.mark.parametrize("particle,antiparticle", particle_antiparticle_pairs)
def test_particle_antiparticle_pairs(particle, antiparticle):
    """Test that particles and antiparticles have the same or exact
    opposite properties in the _Particles dictionary."""

    assert not _Particles[particle]['antimatter'], \
        f"{particle} is incorrectly marked as antimatter."

    assert _Particles[antiparticle]['antimatter'], \
        f"{antiparticle} is incorrectly marked as matter."

    identical_keys = ['half-life', 'spin']

    if 'nu' not in particle:
        identical_keys.append('mass')

    if particle in ['e-', 'mu-', 'tau-'] or 'nu' in particle:
        identical_keys.append('generation')

    opposite_keys = ['charge', 'lepton number', 'baryon number']

    for key in identical_keys:
        assert _Particles[particle][key] == _Particles[antiparticle][key], \
            f"{particle} and {antiparticle} do not have identical {key}."

    for key in opposite_keys:
        assert _Particles[particle][key] == -_Particles[antiparticle][key], \
            f"{particle} and {antiparticle} do not have exact opposite {key}."

    if particle not in ['e-', 'n']:
        assert _Particles[particle]['name'] == \
            _Particles[antiparticle]['name'].replace('anti', ''), \
            (f"{particle} and {antiparticle} do not have same name except "
             "for 'anti'.")


required_keys = [
    'name',
    'spin',
    'class',
    'lepton number',
    'baryon number',
    'charge',
    'half-life',
    'mass',
    'antimatter',
]


@pytest.mark.parametrize("particle", _everything)
def test__Particles_required_keys(particle):
    r"""Test that required keys are present for all particles."""

    missing_keys = []

    for key in required_keys:
        try:
            _Particles[particle][key]
        except KeyError:
            missing_keys.append(key)

    if missing_keys:
        raise KeyError(
            "The following keys are missing from "
            f"_Particles['{particle}']:\n{missing_keys}")
