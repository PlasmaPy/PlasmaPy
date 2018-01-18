import pytest
from ..particles import _Particles

particle_antiparticle_pairs = [
    ('e-', 'e+'),
    ('mu-', 'mu+'),
    ('tau-', 'tau+'),
    ('p', 'p-'),
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
