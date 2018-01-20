import pytest

from ...constants import m_p

from ..particles import (
    _leptons,
    _antileptons,
    _baryons,
    _antibaryons,
    _everything,
    _particles,
    _antiparticles,
    _fermions,
    _bosons,
    _neutrinos,
    _antineutrinos,
    _special_particles,
)

from ..classes import Particle


def test_Particle():
    r"""Original, temporary tests of Particle class."""
    p = Particle('H', mass_numb=1, Z=1)
    assert p._generation is None


# (arg, kwargs, results_dict
test_Particle_table = [
    ('H+', {},
     {'m': m_p}
     )
]

@pytest.mark.parametrize("symbol", _neutrinos + _antineutrinos)
def test_Particle_neutrinos(symbol):
    r"""Test the properties of neutrinos in the Particle class."""

    nu = Particle(symbol)

    assert 'nu' in nu._particle_symbol


    # lepton_number = 1 for neutrinos and -1 for antineutrinos


@pytest.mark.parametrize("symbol", _special_particles)
def test_Particle_special(symbol):
    r"""Test the properties of special particles."""

    particle = Particle(symbol)

    assert particle._atomic_symbol is None
    assert particle._atomic_number is None
    assert particle._isotope_symbol is None
    assert particle._element_name is None
    assert particle._mass_number is None
