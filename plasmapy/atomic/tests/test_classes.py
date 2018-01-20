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
      'half_life': np.inf * u.s,
      'atomic_number': 1,
      'mass_number': 1,
      'lepton_number': 0,
      'baryon_number': 1,
      }),

    ('e-', {},
        {'particle': 'e-',
         'element': InvalidElementError,
        })
]


@pytest.mark.parametrize("arg, kwargs, expected_dict", test_Particle_table)
def test_Particle_class(arg, kwargs, expected_dict):
    r"""Test required properties of items in _Particles dictionary."""

    # To allow tests of a Particle class to continue after the first
    # exception, error messages will be appended as new lines to errmsg.
    errmsg = ""

    cls = f"Particle('{arg}')"

    try:
        particle = Particle(arg)
    except Exception as exc:
        raise Exception(f"Unable to create {cls}.") from exc

    for key in expected_dict.keys():
        expected = expected_dict[key]

        if inspect.isclass(expected) and issubclass(expected, Exception):
            # Exceptions are expected to be raised when accessing certain
            # attributes for some particles.  For example, accessing a
            # neutrino's mass should raise a MissingAtomicDataError.
            # If expected_dict[key] is an exception, then check to make
            # sure that this exception is raised.
            try:
                with pytest.raises(expected):
                    exec(f"particle.{key}")
            except pytest.fail.Exception as exc_failed_fail:
                errmsg += f"\n{cls}.{key} does not raise {expected}."
            except Exception as exc_bad:
                errmsg += (f"\n{cls}.{key} does not raise {expected} and " \
                           f"instead raises a different exception.")
        else:
            try:
                result = eval(f"particle.{key}")
                assert result == expected
            except AssertionError as exc_assert:
                errmsg += f"\n{cls}.{key} does not equal {expected}."
            except Exception as exc_general:
                errmsg += f"\n{cls}.{key} raises an unexpected exception."

    if len(errmsg) > 0:
        raise Exception("The following problems were found for "
                        f"Particle('{key}'):" + errmsg)


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
