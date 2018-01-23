import pytest
import numpy as np
from astropy import units as u
import inspect

from ...constants import m_p, m_e

from ...utils import (
    AtomicWarning,
    AtomicError,
    MissingAtomicDataError,
    InvalidParticleError,
    InvalidElementError,
    InvalidIsotopeError,
    InvalidIonError,
    ChargeError,
)

from typing import Union, Dict

from ..particles import (
    _fermions,
    _neutrinos,
    _antineutrinos,
)

from ..classes import Particle


def _call_string(arg: Union[str, int], kwargs: Dict) -> str:
    r"""Return a string that recreates the call to create a particular
    particle from """
    if kwargs != {}:
        keyword_string = ", " \
            + str(kwargs).strip(r"}{'").replace("'", "").replace(":", " =")
    else:
        keyword_string = ""
    return f"Particle({repr(arg)}{keyword_string})"


# (arg, kwargs, results_dict
test_Particle_table = [

    ('p+', {},
     {'particle': 'p+',
      'element': 'H',
      'isotope': 'H-1',
      'ion': 'p+',
      'mass': m_p,
      'integer_charge': 1,
      'spin': 1 / 2,
      'half_life': np.inf * u.s,
      'atomic_number': 1,
      'mass_number': 1,
      'lepton_number': 0,
      'baryon_number': 1,
      'reduced_mass(Particle("p"))': m_p / 2,
      'reduced_mass(m_p)': m_p / 2,
      '__str__()': 'p+',
      '__repr__()': 'Particle("p+")'}),

    ('p-', {},
     {'particle': 'p-',
      'element': InvalidElementError,
      'isotope': InvalidIsotopeError,
      'ion': InvalidIonError,
      'mass': m_p,
      'integer_charge': -1,
      'spin': 1 / 2,
      'half_life': np.inf * u.s,
      'atomic_number': InvalidElementError,
      'mass_number': InvalidIsotopeError,
      'lepton_number': 0,
      'baryon_number': -1,
      '__str__()': 'p-',
      '__repr__()': 'Particle("p-")'}),

    ('e-', {},
     {'particle': 'e-',
      'element': InvalidElementError,
      'isotope': InvalidIsotopeError,
      'ion': InvalidIonError,
      'mass': m_e,
      'integer_charge': -1,
      'spin': 1 / 2,
      'half_life': np.inf * u.s,
      'atomic_number': InvalidElementError,
      'lepton_number': 1,
      'baryon_number': 0,
      'reduced_mass(Particle("e+"))': m_e / 2,
      'reduced_mass("e-")': m_e / 2,
      '__str__()': 'e-',
      '__repr__()': 'Particle("e-")'}),

    ('e+', {},
     {'particle': 'e+',
      'element': InvalidElementError,
      'isotope': InvalidIsotopeError,
      'ion': InvalidIonError,
      'mass': m_e,
      'integer_charge': 1,
      'spin': 1 / 2,
      'half_life': np.inf * u.s,
      'atomic_number': InvalidElementError,
      'lepton_number': -1,
      'baryon_number': 0,
      'is_stable': True,
      'is_antimatter': True,
      '__str__()': 'e+',
      '__repr__()': 'Particle("e+")'}),

    ('Fe', {'Z': 17, 'mass_numb': 56},
     {'particle': 'Fe-56 17+',
      'element': 'Fe',
      'isotope': 'Fe-56',
      'ion': 'Fe-56 17+',
      'integer_charge': 17,
      'atomic_number': 26,
      'mass_number': 56,
      'baryon_number': 56,
      '__str__()': 'Fe-56 17+',
      '__repr__()': 'Particle("Fe-56 17+")'}),

    ('alpha', {},
     {'particle': 'He-4 2+',
      'element': 'He',
      'isotope': 'He-4',
      'ion': 'He-4 2+',
      'integer_charge': 2,
      'atomic_number': 2,
      'mass_number': 4,
      'baryon_number': 4,
      'lepton_number': 0,
      'half_life': np.inf * u.s,
      'is_stable': True}),

    ('D+', {},
     {'particle': 'D 1+',
      'element': 'H',
      'isotope': 'D',
      'ion': 'D 1+',
      'integer_charge': 1,
      'atomic_number': 1,
      'mass_number': 2,
      'baryon_number': 2,
      'lepton_number': 0}),

    ('tritium', {'Z': 1},
     {'particle': 'T 1+',
      'element': 'H',
      'isotope': 'T',
      'ion': 'T 1+',
      'integer_charge': 1,
      'atomic_number': 1,
      'mass_number': 3,
      'baryon_number': 3,
      'lepton_number': 0}),

    ('muon', {},
     {'particle': 'mu-',
      'element': InvalidElementError,
      'isotope': InvalidIsotopeError,
      'ion': InvalidIonError,
      'integer_charge': -1,
      'atomic_number': InvalidElementError,
      'mass_number': InvalidIsotopeError,
      'baryon_number': 0,
      'lepton_number': 1,
      'is_antimatter': False,
      'is_stable': False}),

    ('neutron', {},
     {'particle': 'n',
      'element': InvalidElementError,
      'isotope': InvalidIsotopeError,
      'ion': InvalidIonError,
      'integer_charge': 0,
      'atomic_number': InvalidElementError,
      'mass_number': InvalidIsotopeError,
      'baryon_number': 1,
      'lepton_number': 0,
      'is_stable': False,
      'is_antimatter': False}),

    ('H', {},
     {'particle': 'H',
      'element': 'H',
      'isotope': InvalidIsotopeError,
      'ion': InvalidIonError,
      'charge': ChargeError,
      'integer_charge': ChargeError,
      'mass_number': InvalidIsotopeError,
      'baryon_number': AtomicError,
      'lepton_number': 0,
      'is_stable': InvalidIsotopeError,  # or  a different error?
      'half_life': InvalidIsotopeError,
      'is_antimatter': False,
      'standard_atomic_weight': (1.008 * u.u).to(u.kg),
      'mass': (1.008 * u.u).to(u.kg),
      })
]


@pytest.mark.parametrize("arg, kwargs, expected_dict", test_Particle_table)
def test_Particle_class(arg, kwargs, expected_dict):
    r"""Test that Particle objects for different subatomic particles,
    elements, isotopes, and ions return the expected properties.  Provide
    a detailed error message that lists all of the inconsistencies with
    the expected results."""

    call = _call_string(arg, kwargs)
    errmsg = ""

    try:
        particle = Particle(arg, **kwargs)
    except Exception as exc:
        raise AtomicError(f"Problem creating {call}") from exc

    for key in expected_dict.keys():
        expected = expected_dict[key]

        if inspect.isclass(expected) and issubclass(expected, Exception):

            # Exceptions are expected to be raised when accessing certain
            # attributes for some particles.  For example, accessing a
            # neutrino's mass should raise a MissingAtomicDataError since
            # only upper limits of neutrino masses are presently available.
            # If expected_dict[key] is an exception, then check to make
            # sure that this exception is raised.

            try:
                with pytest.raises(expected):
                    exec(f"particle.{key}")
            except pytest.fail.Exception as exc_failed_fail:
                errmsg += f"\n{call}[{key}] does not raise {expected}."
            except Exception as exc_bad:
                errmsg += (f"\n{call}[{key}] does not raise {expected} but "
                           f"raises a different exception.")

        else:

            try:
                result = eval(f"particle.{key}")
                assert result == expected
            except AssertionError as exc_assert:
                errmsg += f"\n{call}.{key} does not equal {expected}."
            except Exception as exc_general:
                errmsg += f"\n{call}.{key} raises an unexpected exception."

    if len(errmsg) > 0:
        raise Exception(f"Problems with {call}:{errmsg}")


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

    for Q in equivalent_Particle_classes[1:]:
        assert Q == equivalent_Particle_classes[0], \
            f"{equivalent_particles}"


# arg, kwargs, attribute, exception
test_Particle_error_table = [
    ('a', {}, "", InvalidParticleError),
    ('d+', {'mass_numb': 9}, "", InvalidParticleError),
    ('H', {'mass_numb': 99}, "", InvalidParticleError),
    ('e-', {'Z': -1}, "", InvalidParticleError),
    ('nu_e', {}, '.mass', MissingAtomicDataError),
    ('e-', {}, '.element', InvalidElementError),
    ('H', {'Z': 0}, '.isotope', InvalidIsotopeError),
    ('He', {'mass_numb': 3}, '.ion', InvalidIonError),
    ('He', {'mass_numb': 4}, '.charge', ChargeError),
    ('He', {'mass_numb': 4}, '.integer_charge', ChargeError),
    ('tau+', {}, '.element', InvalidElementError),
    ('neutron', {}, '.atomic_number', InvalidElementError),
    ('neutron', {}, '.mass_number', InvalidIsotopeError),
    ('Fe', {}, '.spin', MissingAtomicDataError),
    ('Og', {}, '.standard_atomic_weight', MissingAtomicDataError),
    ('alpha', {}, '.standard_atomic_weight', InvalidElementError),
    ('Fe-56', {}, '.standard_atomic_weight', InvalidElementError),
]


@pytest.mark.parametrize(
    "arg, kwargs, attribute, exception", test_Particle_error_table)
def test_Particle_errors(arg, kwargs, attribute, exception):
    r"""Test that the appropriate exceptions are raised during the creation
    and use of a Particle object."""
    call = _call_string(arg, kwargs)
    with pytest.raises(exception, message=(
            f"The following command: "
            f"\n\n >>> {_call_string(arg, kwargs)}{attribute}\n\n"
            f"did not raise a {exception.__name__} as expected")):
        exec(f'Particle(arg, **kwargs){attribute}')


# arg, kwargs, attribute, exception
test_Particle_warning_table = [
    ('H----', {}, "", AtomicWarning),
    ('alpha', {'mass_numb': 4}, "", AtomicWarning),

]


@pytest.mark.parametrize(
    "arg, kwargs, attribute, warning", test_Particle_warning_table)
def test_Particle_warnings(arg, kwargs, attribute, warning):
    r"""Test that the appropriate warnings are issued during the creation
    and use of a Particle object."""
    with pytest.warns(warning, message=(
            f"The following command: "
            f"\n\n >>> {_call_string(arg, kwargs)}{attribute}\n\n"
            f"did not issue a {warning.__name__} as expected")):
        exec(f'Particle(arg, **kwargs){attribute}')


def test_Particle_cmp():
    proton1 = Particle('p+')
    proton2 = Particle('proton')
    electron = Particle('e-')

    assert proton1 == proton2, "Particle('p+') == Particle('proton') is False."
    assert proton1 != electron, "Particle('p+') == Particle('e-') is True."
    assert not proton1 == 1, "Particle('p+') == 1 is True."
