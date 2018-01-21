import pytest
import numpy as np
from astropy import units as u
import inspect

from ...constants import m_p, m_e

from ...utils import (
    MissingAtomicDataError,
    InvalidParticleError,
    InvalidElementError,
    InvalidIsotopeError,
    InvalidIonError,
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
      'mass': m_p,
      'integer_charge': 1,
      'spin': 1/2,
      'half_life': np.inf * u.s,
      'atomic_number': 1,
      'mass_number': 1,
      'lepton_number': 0,
      'baryon_number': 1,
      }),

    ('p-', {},
     {'particle': 'p-',
      'element': InvalidElementError,
      'isotope': InvalidIsotopeError,
      'ion': InvalidIonError,
      'mass': m_p,
      'integer_charge': -1,
      'spin': 1/2,
      'half_life': np.inf * u.s,
      'atomic_number': InvalidElementError,
      'mass_number': InvalidIsotopeError,
      'lepton_number': 0,
      'baryon_number': -1,
      }),

    ('e-', {},
        {'particle': 'e-',
         'element': InvalidElementError,
         'isotope': InvalidIsotopeError,
         'ion': InvalidIonError,
         'mass': m_e,
         'integer_charge': -1,
         'spin': 1/2,
         'half_life': np.inf * u.s,
         'atomic_number': InvalidElementError,
         'lepton_number': 1,
         'baryon_number': 0,
         }),

    ('e+', {},
        {'particle': 'e+',
         'element': InvalidElementError,
         'isotope': InvalidIsotopeError,
         'ion': InvalidIonError,
         'mass': m_e,
         'integer_charge': 1,
         'spin': 1/2,
         'half_life': np.inf * u.s,
         'atomic_number': InvalidElementError,
         'lepton_number': -1,
         'baryon_number': 0,
         }),
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
                errmsg += (f"\n{cls}.{key} does not raise {expected} and "
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


@pytest.mark.parametrize("symbol", _neutrinos + _antineutrinos)
def test_Particle_neutrinos(symbol):
    r"""Test the properties of neutrinos in the Particle class."""

    nu = Particle(symbol)
    assert 'nu' in nu._particle_symbol

    assert nu._mass is None, \
        f"Particle('{symbol}')._mass should be None."

    with pytest.raises(MissingAtomicDataError, message=(
            f"Particle('{symbol}').m is not raising an exception")):
        nu.mass


@pytest.mark.parametrize("symbol", _fermions)
def test_Particle_fermions(symbol):
    r"""Test that fermions have spin = 1/2."""
    fermion = Particle(symbol)
    assert np.isclose(fermion.spin, 1/2, atol=1e-15)


equivalent_particles_table = [
    ['H', 'hydrogen', 'hYdRoGeN'],
    ['p+', 'proton', 'H-1+', 'H-1 1+', 'H-1 +1'],
    ['D', 'H-2', 'Hydrogen-2', 'deuterium'],
    ['T', 'H-3', 'Hydrogen-3', 'tritium'],
    ['alpha', 'He-4++', 'He-4 2+', 'He-4 +2'],
    ['e-', 'electron', 'e'],
    ['e+', 'positron'],
    ['p-', 'antiproton'],
    ['n', 'n-1', 'neutron', 'NEUTRON'],
    ['muon', 'mu-', 'muon-'],
    ['tau', 'tau-'],
]


@pytest.mark.parametrize("equivalent_particles", equivalent_particles_table)
def test_Particle_equivalent_cases(equivalent_particles):
    r"""Test that all instances of a list of particles are equivalent,
    except for the _original_* private variables which will differ."""

    equivalent_Particle_classes = []

    for particle in equivalent_particles:
        equivalent_Particle_classes.append(Particle(particle))

    for Q in equivalent_Particle_classes:
        del(Q._original_argument)
        del(Q._original_mass_number)
        del(Q._original_integer_charge)

    for Q in equivalent_Particle_classes:
        assert Q == equivalent_Particle_classes[0], \
            f"{equivalent_particles}"
