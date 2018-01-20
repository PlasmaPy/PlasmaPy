import pytest
import numpy as np
from astropy import units as u
import inspect

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
    ('p+', {},
     {'particle': 'p+',
      'element': 'H',
      'isotope': 'H-1',
      'ion': 'p+',
      'm': m_p,
      'Z': 1,
      'spin': 1/2,
      'half-life': np.inf * u.s,
      'atomic_number': 1,
      'mass_number': 1,
      'lepton_number': 0,
      'baryon_number': 1,
      }
     )
]


@pytest.mark.parametrize("arg, kwargs, expected_dict", test_Particle_table)
def test_Particle_class(arg, kwargs, expected_dict):
    r"""Test required properties of items in _Particles dictionary."""

    particle = Particle(arg)

    errmsg = ""

    for key in expected_dict.keys():
        expected = expected_dict[key]

        if inspect.isclass(expected) and issubclass(expected, Exception):
            try:
                with pytest.raises(expected):
                    exec(f"particle.{key}")
            except pytest.fail.Exception as exc_failed_fail:
                errmsg += f"{key} {expected}\n"
            except Exception as exc_bad:
                errmsg += f"{key} {expected} {exc_bad}\n"
        else:
            try:
                result = eval(f"particle.{key}")
                assert result == expected
            except AssertionError as exc_assert:
                errmsg += f"{key} {result} != {expected}\n"
            except Exception as exc_general:
                errmsg += f"{key} {exc_general}\n"

    if len(errmsg) > 0:
        raise Exception(errmsg)


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
