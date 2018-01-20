import pytest

from ...constants import m_p

from ...utils import (
    MissingAtomicDataError,
    InvalidParticleError,
    InvalidElementError,
    AtomicError,
)

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


# (arg, kwargs, results_dict
test_Particle_table = [
    ('H+', {},
     {'m': m_p}
     )
]


@pytest.mark.parametrize("symbol", _everything)
def test_Particle_everything(symbol):
    r"""Test required properties of items in _Particles dictionary."""

    particle = Particle(symbol)

    if not particle.element:
        pass


@pytest.mark.parametrize("symbol", _special_particles)
def test_Particle_special(symbol):
    r"""Test the properties of special particles that do not
    correspond to elements."""

    particle = Particle(symbol)

    assert particle._atomic_symbol is None, \
        f"Particle('symbol')._atomic_symbol is not None"
    assert particle._atomic_number is None, \
        f"Particle('symbol')._atomic_number is not None"
    assert particle._isotope_symbol is None, \
        f"Particle('symbol')._isotope_symbol is not None"
    assert particle._element_name is None, \
        f"Particle('symbol')._element_name is not None"
    assert particle._mass_number is None, \
        f"Particle('symbol')._mass_number is not None"


@pytest.mark.parametrize("symbol", _neutrinos + _antineutrinos)
def test_Particle_neutrinos(symbol):
    r"""Test the properties of neutrinos in the Particle class."""

    nu = Particle(symbol)
    assert 'nu' in nu._particle_symbol

    assert nu._mass is None, \
        f"Particle('{symbol}')._mass should be None."

    with pytest.raises(MissingAtomicDataError, message=(
            f"Particle('{symbol}').m is not raising an exception")):
        nu.m


@pytest.mark.parametrize("symbol", _leptons + _antileptons)
def test_Particle_leptons(symbol):

    lepton = Particle(symbol)

    assert lepton.spin == 1/2
