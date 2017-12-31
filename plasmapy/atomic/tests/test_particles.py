import pytest

from ..particles import (
    _get_standard_symbol,
    _Particles,
    _case_sensitive_aliases,
    _case_insensitive_aliases)


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


aliases_and_symbols = [
    ('electron', 'e-'),
    ('beta-', 'e-'),
    ('beta+', 'e+'),
    ('positron', 'e+'),
    ('proton', 'p'),
    ('', ''),
    (5, 5),
    ('deuterium+', 'D 1+'),
    ('deuterium 1+', 'D 1+'),
    ('tritium +1', 'T 1+'),
    ('alpha', 'He-4 2+'),
    ('D+', 'D 1+'),
    ('Deuterium', 'D'),
    ('deuteron', 'D 1+'),
    ('triton', 'T 1+'),
    ('muon', 'mu-'),
    ('antimuon', 'mu+'),
    ('tau particle', 'tau-'),
    ('antitau', 'tau+'),
    ('p+', 'p'),
]


@pytest.mark.parametrize("alias,symbol", aliases_and_symbols)
def test_get_standard_symbol(alias, symbol):
    """Test that _get_standard_symbol correctly takes in aliases and
    returns the corresponding symbols, and returns the original argument
    if the argument does not correspond to an alias."""
    result = _get_standard_symbol(alias)
    assert result == symbol, \
        (f"_get_standard_symbol({alias}) returns {result}, which differs "
         f"from the expected symbol of {symbol}.\n\n"
         f"_case_insensitive_aliases:\n{_case_insensitive_aliases}\n\n"
         f"_case_sensitive_aliases:\n{_case_sensitive_aliases}")


alias_dictionaries = [_case_sensitive_aliases, _case_insensitive_aliases]


@pytest.mark.parametrize("alias_dict", alias_dictionaries)
def test_alias_dict_properties(alias_dict):
    """Test properties of the alias dictionaries."""
    for key in alias_dict.keys():
        assert isinstance(key, str), f"{key}\n{alias_dict}"
    for value in alias_dict.values():
        assert isinstance(value, str), f"{value}\n{alias_dict}"
