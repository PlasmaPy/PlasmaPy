import pytest
import numpy as np
from astropy import units as u
import inspect

from ...constants import m_p, m_e, m_n, e

from ...utils import (
    AtomicWarning,
    AtomicError,
    MissingAtomicDataError,
    MissingAtomicDataWarning,
    InvalidParticleError,
    InvalidElementError,
    InvalidIsotopeError,
    ChargeError,
)

from ..atomic import known_isotopes
from ..isotopes import _Isotopes
from ..particle_class import Particle
from ..parsing import _call_string

# (arg, kwargs, results_dict
test_Particle_table = [

    ('neutron', {},
     {'particle': 'n',
      'element': None,
      'isotope': None,
      'ion': None,
      'integer_charge': 0,
      'atomic_number': InvalidElementError,
      'mass_number': InvalidIsotopeError,
      'baryon_number': 1,
      'lepton_number': 0,
      'mass': m_n,
      'nuclide_mass': m_n,
      'binding_energy': 0 * u.J,
      'periodic_table.group': InvalidElementError,
      }),

    ('p+', {},
     {'particle': 'p+',
      'element': 'H',
      'element_name': 'hydrogen',
      'isotope': 'H-1',
      'ion': 'p+',
      'mass': m_p,
      'nuclide_mass': m_p,
      'integer_charge': 1,
      'charge.value': e.si.value,
      'spin': 1 / 2,
      'half_life': np.inf * u.s,
      'atomic_number': 1,
      'mass_number': 1,
      'lepton_number': 0,
      'baryon_number': 1,
      '__str__()': 'p+',
      '__repr__()': 'Particle("p+")',
      'is_category("fermion")': True,
      'is_category(["fermion"])': True,
      'is_category({"fermion"})': True,
      'is_category(any_of=("boson", "fermion"))': True,
      'is_category(require=["boson", "fermion"])': False,
      'is_category(("element", "isotope", "ion"))': True,
      'is_category("charged")': True,
      'periodic_table.group': 1,
      'periodic_table.block': 's',
      'periodic_table.period': 1,
      'periodic_table.category': 'nonmetal',
      'binding_energy': 0 * u.J,
      }),

    ('p-', {},
     {'particle': 'p-',
      'element': None,
      'element_name': InvalidElementError,
      'isotope': None,
      'ion': None,
      'mass': m_p,
      'integer_charge': -1,
      'spin': 1 / 2,
      'half_life': np.inf * u.s,
      'atomic_number': InvalidElementError,
      'mass_number': InvalidIsotopeError,
      'lepton_number': 0,
      'baryon_number': -1,
      '__str__()': 'p-',
      '__repr__()': 'Particle("p-")',
      'periodic_table.group': InvalidElementError,
      }),

    ('e-', {},
     {'particle': 'e-',
      'element': None,
      'element_name': InvalidElementError,
      'isotope': None,
      'ion': None,
      'mass': m_e,
      'integer_charge': -1,
      'spin': 1 / 2,
      'half_life': np.inf * u.s,
      'atomic_number': InvalidElementError,
      'lepton_number': 1,
      'baryon_number': 0,
      '__str__()': 'e-',
      '__repr__()': 'Particle("e-")',
      'binding_energy': InvalidIsotopeError,
      'periodic_table.group': InvalidElementError,
      'periodic_table.block': InvalidElementError,
      'periodic_table.period': InvalidElementError,
      'periodic_table.category': InvalidElementError,
      }),

    ('e+', {},
     {'particle': 'e+',
      'element': None,
      'isotope': None,
      'ion': None,
      'mass': m_e,
      'nuclide_mass': InvalidIsotopeError,
      'integer_charge': 1,
      'spin': 1 / 2,
      'half_life': np.inf * u.s,
      'atomic_number': InvalidElementError,
      'lepton_number': -1,
      'baryon_number': 0,
      'is_category(require="positron")': True,
      'is_category(any_of={"positron"})': True,
      'is_category(exclude="positron")': False,
      '__str__()': 'e+',
      '__repr__()': 'Particle("e+")',
      'periodic_table.group': InvalidElementError,
      'periodic_table.block': InvalidElementError,
      'periodic_table.period': InvalidElementError,
      'periodic_table.category': InvalidElementError,
      }),

    ('H', {},
     {'particle': 'H',
      'element': 'H',
      'isotope': None,
      'ion': None,
      'charge': ChargeError,
      'integer_charge': ChargeError,
      'mass_number': InvalidIsotopeError,
      'baryon_number': AtomicError,
      'lepton_number': 0,
      'half_life': InvalidIsotopeError,
      'standard_atomic_weight': (1.008 * u.u).to(u.kg),
      'mass': (1.008 * u.u).to(u.kg),
      'nuclide_mass': InvalidIsotopeError,
      'is_category("charged")': False,
      'is_category("nonmetal")': True,
      'is_category("proton")': False,
      }),

    ('D+', {},
     {'particle': 'D 1+',
      'element': 'H',
      'element_name': 'hydrogen',
      'isotope': 'D',
      'ion': 'D 1+',
      'integer_charge': 1,
      'atomic_number': 1,
      'mass_number': 2,
      'baryon_number': 2,
      'lepton_number': 0,
      'neutron_number': 1,
      'is_category(require=("ion", "isotope"))': True,
      'periodic_table.group': 1,
      'periodic_table.block': 's',
      'periodic_table.period': 1,
      'periodic_table.category': 'nonmetal',
      }),

    ('tritium', {'Z': 1},
     {'particle': 'T 1+',
      'element': 'H',
      'isotope': 'T',
      'ion': 'T 1+',
      'integer_charge': 1,
      'atomic_number': 1,
      'mass_number': 3,
      'baryon_number': 3,
      'lepton_number': 0,
      'neutron_number': 2,
      'is_category("ion", "isotope")': True,
      'is_category(require="uncharged")': False,
      'periodic_table.group': 1,
      }),

    ('Fe', {'Z': 17, 'mass_numb': 56},
     {'particle': 'Fe-56 17+',
      'element': 'Fe',
      'element_name': 'iron',
      'isotope': 'Fe-56',
      'ion': 'Fe-56 17+',
      'integer_charge': 17,
      'atomic_number': 26,
      'mass_number': 56,
      'baryon_number': 56,
      '__str__()': 'Fe-56 17+',
      '__repr__()': 'Particle("Fe-56 17+")',
      'is_category("element")': True,
      'is_category("ion")': True,
      'is_category("isotope")': True,
      }),

    ('alpha', {},
     {'particle': 'He-4 2+',
      'element': 'He',
      'element_name': 'helium',
      'isotope': 'He-4',
      'ion': 'He-4 2+',
      'integer_charge': 2,
      'atomic_number': 2,
      'mass_number': 4,
      'baryon_number': 4,
      'lepton_number': 0,
      'half_life': np.inf * u.s,
      }),

    ('Li', {'mass_numb': 7},
     {'particle': 'Li-7',
      'element': 'Li',
      'element_name': 'lithium',
      'isotope': 'Li-7',
      'ion': None,
      'integer_charge': ChargeError,
      'atomic_number': 3,
      'mass_number': 7,
      'neutron_number': 4,
      'baryon_number': 7,
      'half_life': np.inf * u.s,
      'nuclide_mass': 1.1647614796180463e-26 * u.kg,
      }),

    ('Cn-276', {"Z": 22},
     {'particle': 'Cn-276 22+',
      'element': 'Cn',
      'isotope': 'Cn-276',
      'ion': 'Cn-276 22+',
      'element_name': 'copernicium',
      'integer_charge': 22,
      'atomic_number': 112,
      'mass_number': 276,
      'neutron_number': 164,
      'baryon_number': 276,
      'lepton_number': 0}),

    ('muon', {},
     {'particle': 'mu-',
      'element': None,
      'isotope': None,
      'ion': None,
      'integer_charge': -1,
      'atomic_number': InvalidElementError,
      'mass_number': InvalidIsotopeError,
      'baryon_number': 0,
      'lepton_number': 1,
      }),

    ('nu_tau', {},
     {'particle': 'nu_tau',
      'element': None,
      'isotope': None,
      'mass': MissingAtomicDataError,
      'integer_charge': 0,
      'mass_number': InvalidIsotopeError,
      'element_name': InvalidElementError,
      'baryon_number': 0,
      'lepton_number': 1,
      'half_life': np.inf * u.s,
      'is_category("fermion")': True,
      'is_category("neutrino")': True,
      'is_category("boson")': False,
      'is_category("matter", exclude={"antimatter"})': True,
      'is_category("matter", exclude=["antimatter"])': True,
      'is_category("matter", exclude="antimatter")': True,
      'is_category(any_of={"matter", "boson"})': True,
      'is_category(any_of=["antimatter", "boson", "charged"])': False,
      'is_category(["fermion", "lepton"], exclude="matter")': False,
      'is_category("lepton", "invalid")': AtomicError,
      'is_category(["boson"], exclude=["lepton", "invalid"])': AtomicError,
      'is_category("boson", exclude="boson")': AtomicError,
      'is_category(any_of="boson", exclude="boson")': AtomicError,
      }),
]


@pytest.mark.parametrize("arg, kwargs, expected_dict", test_Particle_table)
def test_Particle_class(arg, kwargs, expected_dict):
    """Test that `~plasmapy.atomic.Particle` objects for different
    subatomic particles, elements, isotopes, and ions return the
    expected properties.  Provide a detailed error message that lists
    all of the inconsistencies with the expected results.
    """

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
                errmsg += (f"\n{call}.{key} returns {result} instead "
                           f"of the expected value of {expected}.")
            except Exception as exc_general:
                errmsg += f"\n{call}.{key} raises an unexpected exception."

    if len(errmsg) > 0:
        raise Exception(f"Problems with {call}:{errmsg}")


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
    """Test that all instances of a list of particles are equivalent."""

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
    ('Au-818', {}, "", InvalidParticleError),
    ('Au-12', {}, "", InvalidParticleError),
    ('Au', {'mass_numb': 13}, "", InvalidParticleError),
    ('Au', {'mass_numb': 921}, "", InvalidParticleError),
    ('e-', {'Z': -1}, "", InvalidParticleError),
    ('e-', {}, '.atomic_number', InvalidElementError),
    ('alpha', {}, '.standard_atomic_weight', InvalidElementError),
    ('Fe-56', {}, '.standard_atomic_weight', InvalidElementError),
    ('e-', {}, '.standard_atomic_weight', InvalidElementError),
    ('tau-', {}, '.element_name', InvalidElementError),
    ('tau+', {}, '.atomic_number', InvalidElementError),
    ('neutron', {}, '.atomic_number', InvalidElementError),
    ('H', {'Z': 0}, '.mass_number', InvalidIsotopeError),
    ('neutron', {}, '.mass_number', InvalidIsotopeError),
    ('He', {'mass_numb': 4}, '.charge', ChargeError),
    ('He', {'mass_numb': 4}, '.integer_charge', ChargeError),
    ('Fe', {}, '.spin', MissingAtomicDataError),
    ('nu_e', {}, '.mass', MissingAtomicDataError),
    ('Og', {}, '.standard_atomic_weight', MissingAtomicDataError),
    ([], {}, "", TypeError),
]


@pytest.mark.parametrize(
    "arg, kwargs, attribute, exception", test_Particle_error_table)
def test_Particle_errors(arg, kwargs, attribute, exception):
    """Test that the appropriate exceptions are raised during the creation
    and use of a `~plasmapy.atomic.Particle` object.
    """
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
    ('alpha', {'Z': 2}, "", AtomicWarning)
]


@pytest.mark.parametrize(
    "arg, kwargs, attribute, warning", test_Particle_warning_table)
def test_Particle_warnings(arg, kwargs, attribute, warning):
    """Test that the appropriate warnings are issued during the creation
    and use of a `~plasmapy.atomic.Particle` object."""
    with pytest.warns(warning, message=(
            f"The following command: "
            f"\n\n >>> {_call_string(arg, kwargs)}{attribute}\n\n"
            f"did not issue a {warning.__name__} as expected")):
        exec(f'Particle(arg, **kwargs){attribute}')


def test_Particle_cmp():
    """Test __eq__ and __ne__ in the Particle class."""
    proton1 = Particle('p+')
    proton2 = Particle('proton')
    electron = Particle('e-')

    assert proton1 == proton2, "Particle('p+') == Particle('proton') is False."
    assert proton1 != electron, "Particle('p+') == Particle('e-') is True."
    assert not proton1 == 1, "Particle('p+') == 1 is True."


nuclide_mass_and_mass_equiv_table = [
    ('n', 'neutron'),
    ('p+', 'proton'),
    ('H-1', 'p+'),
    ('D', 'D+'),
    ('T', 'T+'),
    ('He-4', 'alpha'),
    ('Fe-56', 'Fe-56 26+'),
]


@pytest.mark.parametrize('isotope, ion', nuclide_mass_and_mass_equiv_table)
def test_particle_class_mass_nuclide_mass(isotope: str, ion: str):
    """Test that the `mass` and `nuclide_mass` attributes return
    equivalent values when appropriate.  The inputs should generally be
    an isotope with no charge information, and a fully ionized ion of
    that isotope, in order to make sure that the nuclide mass of the
    isotope equals the mass of the fully ionized ion.  This method may
    also check neutrons and protons.
    """

    Isotope = Particle(isotope)
    Ion = Particle(ion)

    if Isotope.particle == Ion.particle and Isotope.particle in ('n', 'p+'):
        particle = Isotope.particle

        assert Isotope.nuclide_mass == Ion.mass, (
            f"Particle({repr(particle)}).nuclide_mass does not equal "
            f"Particle({repr(particle)}).mass")

    else:

        inputerrmsg = (f"isotope = {repr(isotope)} and ion = {repr(ion)} are "
                       f"not valid inputs to this test. The inputs should be "
                       f"an isotope with no charge information, and a fully "
                       f"ionized ion of that isotope, in order to make sure "
                       f"that the nuclide mass of the isotope equals the mass "
                       f"of the ion.")

        assert Isotope.isotope and not Isotope.ion, inputerrmsg
        assert Isotope.isotope == Ion.isotope, inputerrmsg
        assert Ion.integer_charge == Ion.atomic_number, inputerrmsg

        assert Isotope.nuclide_mass == Ion.mass, (
            f"The nuclide mass of {isotope} does not equal the mass of {ion} "
            f"which is the fully ionized ion of that isotope. The results of "
            f"the test are:\n\n"
            f"Particle({repr(ion)}).mass = {Ion.mass}\n"
            f"Particle({repr(isotope)}).nuclide_mass = {Isotope.nuclide_mass}"
            "\n")


def test_particle_half_life_string():
    """Find the first isotope where the half-life is stored as a string
    (because the uncertainties are too great), and tests that requesting
    the half-life of that isotope causes a MissingAtomicDataWarning
    whilst returning a string.
    """

    for isotope in known_isotopes():
        half_life = _Isotopes[isotope].get('half-life', None)
        if isinstance(half_life, str):
            break

    with pytest.warns(MissingAtomicDataWarning):
        assert isinstance(Particle(isotope).half_life, str)


@pytest.mark.parametrize("p, is_one", [(Particle('e-'), True), (Particle('p+'), False)])
def test_particle_is_electron(p, is_one):
    assert p.is_electron() == is_one
