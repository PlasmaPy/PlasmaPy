import pytest
import numpy as np
from itertools import product
from astropy import units as u, constants as const

from ..atomic import (atomic_symbol,
                      isotope_symbol,
                      atomic_number,
                      mass_number,
                      element_name,
                      standard_atomic_weight,
                      isotope_mass,
                      ion_mass,
                      is_isotope_stable,
                      half_life,
                      known_isotopes,
                      common_isotopes,
                      stable_isotopes,
                      isotopic_abundance,
                      charge_state,
                      electric_charge,
                      Isotopes,
                      _is_neutron,
                      _is_hydrogen,
                      _is_electron,
                      _is_positron,
                      _is_antiproton,
                      _is_alpha,
                      _extract_charge_state,
                      _is_proton)
from ..nuclear import (nuclear_binding_energy, nuclear_reaction_energy)
from ...utils import (AtomicWarning,
                      ElementError,
                      IsotopeError,
                      IonError,
                      ChargeError)

# (argument, expected)
atomic_symbol_table = [
    (1, 'H'),
    ('H', 'H'),
    ('p', 'H'),
    ('T', 'H'),
    ('deuterium', 'H'),
    ('deuteron', 'H'),
    ('Tritium', 'H'),
    ('triton', 'H'),
    ('H-2', 'H'),
    ('D', 'H'),
    ('T', 'H'),
    ('H-3', 'H'),
    ('Hydrogen-3', 'H'),
    ('helium', 'He'),
    (2, 'He'),
    ('alpha', 'He'),
    ('gold', 'Au'),
    ('Gold', 'Au'),
    (79, 'Au'),
    ('79', 'Au'),
    ('P', 'P'),
    (118, 'Og'),
    ('N-14', 'N'),
    ('N', 'N'),
    ('H +1', 'H'),
    ('H 1+', 'H'),
    ('hydrogen 1+', 'H'),
    ('deuterium 1+', 'H'),
    ('Fe 24+', 'Fe'),
    ('Fe +24', 'Fe'),
    ('Fe 2-', 'Fe'),
    ('Fe -2', 'Fe'),
    ('Fe+', 'Fe'),
    ('Fe++', 'Fe'),
    ('Fe-', 'Fe'),
    ('Fe++++++++++++++', 'Fe')]


@pytest.mark.parametrize(
    'argument, expected', atomic_symbol_table)
def test_atomic_symbol(argument, expected):
    """Test that atomic_symbol returns the expected result."""
    assert atomic_symbol(argument) == expected, \
        (f"atomic_symbol({argument}) is returning {atomic_symbol(argument)} "
         f"which differs from the expected value of {expected}.")


# (argument, expected_error)
atomic_symbol_error_table = [
    ('H-0', IsotopeError),
    (3.14159, TypeError),
    ('Og-294b', IsotopeError),
    ('H-934361079326356530741942970523610389', IsotopeError),
    ('Fe 2+4', ChargeError),
    ('Fe+24', ElementError),
    ('Fe +59', IonError),
    ('C++++++++++++++++', IonError),
    ('C-++++', ChargeError),
    ('neutron', ElementError),
    ('n', ElementError),
    ('n-1', ElementError),
    ('h', ElementError),
    ('d', ElementError),
    ('he', ElementError),
    ('au', ElementError),
    ('p-', ElementError),
    ('antiproton', ElementError)]


@pytest.mark.parametrize(
    'argument, expected_error', atomic_symbol_error_table)
def test_atomic_symbol_error(argument, expected_error):
    """Test that atomic_symbol raises the expected exceptions."""
    with pytest.raises(expected_error, message=(
            f"atomic_symbol({argument}) is not raising {expected_error}.")):
        atomic_symbol(argument)


# (argument, expected)
isotope_symbol_table = [
    (('He', 4), 'He-4'),
    (('helium-4',), 'He-4'),
    (('H-2',), 'D'),
    (('Deuterium',), 'D'),
    (('deuterium',), 'D'),
    (('deuteron',), 'D'),
    (('tritium',), 'T'),
    (('triton',), 'T'),
    (('Hydrogen-3',), 'T'),
    (('hydrogen-3',), 'T'),
    (('H-3',), 'T'),
    ((1, 2), 'D'),
    (('Hydrogen', 3), 'T'),
    (('tritium',), 'T'),
    (('H', 2), 'D'),
    (('Alpha',), 'He-4'),
    (('alpha',), 'He-4'),
    ((79, 197), 'Au-197'),
    (('p',), 'H-1'),
    (('beryllium-8',), 'Be-8'),
    (('N-13',), 'N-13'),
    (('p',), 'H-1'),
    (('proton',), 'H-1'),
    (('protium',), 'H-1'),
    (('N-13 2+',), 'N-13'),
    (('Hydrogen-3 +1',), 'T'),
    (('neutron',), 'n'),
    (('n',), 'n'),
    ((0, 1), 'n'),
    (('neutron',), 'n'),
    (('Neutron',), 'n'),
    (('n-1',), 'n')]


@pytest.mark.parametrize(
    "arguments, expected", isotope_symbol_table)
def test_isotope_symbol(arguments, expected):
    """Test that isotope_symbol returns the expected results."""
    assert isotope_symbol(*arguments) == expected, \
        (f"isotope_symbol is returning {isotope_symbol(*arguments)} "
         f"for arguments of {arguments}, which differs from the "
         f"expected value of {expected}.")


# (argument, kwargs, expected_error)
isotope_symbol_error_table = [
    ('Md-260', {"mass_numb": 261}, IsotopeError),
    ('protium', {"mass_numb": 2}, IsotopeError),
    ('alpha', {"mass_numb": 3}, IsotopeError),
    ('O-18', {"mass_numb": 19}, IsotopeError),
    ('lead-209', {"mass_numb": 511}, IsotopeError),
    ('He-1', {}, IsotopeError),
    (24, {"mass_numb": 23}, IsotopeError),
    ('H', {"mass_numb": 0}, IsotopeError),
    ('H-1', {"mass_numb": 2}, IsotopeError),
    ('P', {}, IsotopeError),
    (1, {}, IsotopeError),
    (4, {}, IsotopeError),
    ('hydrogen-444444', {}, IsotopeError),
    ('Fe', {"mass_numb": 2.1}, TypeError),
    ('He', {"mass_numb": 'c'}, TypeError),
    ('He-3', {"mass_numb": 4}, IsotopeError),
    ('alpha', {"mass_numb": 3}, IsotopeError),
    ('D', {"mass_numb": 3}, IsotopeError),
    ('T', {"mass_numb": 2}, IsotopeError),
    ('Fe', {"mass_numb": None}, IsotopeError),
    ('d', {}, ElementError),
    ('h-3', {}, ElementError),
    ('h', {}, ElementError),
    ('d+', {}, ElementError),
]


@pytest.mark.parametrize(
    "argument, kwargs, expected_error", isotope_symbol_error_table)
def test_isotope_symbol_error(argument, kwargs, expected_error):
    """Test that isotope_symbol raises the expected exceptions."""
    with pytest.raises(expected_error, message=(
            f"isotope_symbol({argument}, **{kwargs}) is not raising a "
            f"{expected_error}.")):
        isotope_symbol(argument, **kwargs)
    pass


# (argument, kwargs, expected_warning)
isotope_symbol_warning_table = [
    ('H-1', {"mass_numb": 1}, AtomicWarning),
    ('H-2', {"mass_numb": 2}, AtomicWarning),
    ('T', {"mass_numb": 3}, AtomicWarning),
    ('Li-6', {"mass_numb": 6}, AtomicWarning),
    ('lithium-6', {"mass_numb": 6}, AtomicWarning),
    ('alpha', {"mass_numb": 4}, AtomicWarning),
    ('p', {"mass_numb": 1}, AtomicWarning)]


@pytest.mark.parametrize(
    "argument, kwargs, expected_warning", isotope_symbol_warning_table)
def test_isotope_symbol_warnings(argument, kwargs, expected_warning):
    """Test that isotope_symbol issues the expected warnings."""
    with pytest.warns(expected_warning, message=(
            f"isotope_symbol({argument}, **{kwargs}) is not issuing a "
            f"{expected_warning}.")):
        isotope_symbol(argument, **kwargs)


# (argument, expected)
atomic_number_table = [
    ('H', 1),
    ('D', 1),
    ('deuterium', 1),
    ('Deuterium', 1),
    ('tritium', 1),
    ('p', 1),
    ('P', 15),
    ('Alpha', 2),
    ('C-12', 6),
    ('Argon', 18),
    ('protium', 1),
    ('H-3', 1),
    ('p+', 1),
    ('Be-8', 4),
    ('N', 7),
    ('N 2+', 7),
    ('N +1', 7),
    ('N+++', 7)]


@pytest.mark.parametrize("argument, expected", atomic_number_table)
def test_atomic_number(argument, expected):
    """Test that atomic_number returns the expected results."""
    assert atomic_number(argument) == expected, \
        (f"atomic_number({argument}) is expecting a result of {expected} but "
         f"is getting a result of {atomic_number(argument)}.")


# (argument, expected_error)
atomic_number_error_table = [
    ('H-3934', IsotopeError),
    ('C-12b', IsotopeError),
    (-1.5, TypeError),
    ('n', ElementError),
    ('n-1', ElementError),
    ('neutron', ElementError),
    ('Neutron', ElementError),
    ('d', ElementError),
    ('t', ElementError),
    ('s-36', ElementError)]


@pytest.mark.parametrize(
    "argument, expected_error", atomic_number_error_table)
def test_atomic_number_error(argument, expected_error):
    """Test that atomic_number raises the expected exceptions."""
    with pytest.raises(expected_error, warning=(
            f"atomic_number({argument}) is not raising a {expected_error}")):
        atomic_number(argument)



# (isotope, expected)
mass_number_table = [
    ('helium-3', 3),
    ('Au-197', 197),
    ('deuterium', 2),
    ('D', 2),
    ('H-2', 2),
    ('tritium', 3),
    ('T', 3),
    ('alpha', 4),
    ('p', 1),
    ('n', 1),
    ('neutron', 1),
    ('n-1', 1),
    ('Be-8', 8),
    ('N-13', 13),
    ('N-13 2+', 13),
    ('N-13 +2', 13),
    ('N-13+++', 13)]


@pytest.mark.parametrize("isotope, expected", mass_number_table)
def test_mass_number(isotope, expected):
    """Test that mass_number returns the expected results."""
    assert mass_number(isotope) == expected, \
        (f"mass_number({isotope}) is returning a value of "
         f"{mass_number(isotope)}, which differs from the expected "
         f"value of {expected}.")


# (argument, expected_error)
mass_number_error_table = [
    ('H-359', IsotopeError),
    ('C-12b', IsotopeError),
    (-1.5, Exception),
    ('N-13+-+-', ChargeError),
    ('h-3', ElementError)]


@pytest.mark.parametrize(
    "argument, expected_error", mass_number_error_table)
def test_mass_number_error(argument, expected_error):
    """Test that mass_number raises the expected exceptions."""
    with pytest.raises(expected_error):
        mass_number(argument)


# (argument, expected)
element_name_table = [
    ('D', 'hydrogen'),
    ('deuterium', 'hydrogen'),
    ('Au', 'gold'),
    ('alpha', 'helium'),
    ('helium-4', 'helium'),
    ('H-2', 'hydrogen'),
    ('Deuterium', 'hydrogen'),
    ('Hydrogen-3', 'hydrogen'),
    ('hydrogen-3', 'hydrogen'),
    ('H-3', 'hydrogen'),
    ('tritium', 'hydrogen'),
    ('Alpha', 'helium'),
    ('alpha', 'helium'),
    (1, 'hydrogen'),
    (26, 'iron'),
    (79, 'gold'),
    ('p', 'hydrogen'),
    ('P', 'phosphorus'),
    ('Be-8', 'beryllium'),
    ('Li-7', 'lithium'),
    ('N', 'nitrogen'),
    ('N+++', 'nitrogen'),
    ('D-', 'hydrogen')]


@pytest.mark.parametrize("argument, expected", element_name_table)
def test_element_name(argument, expected):
    """Test that element_name returns the expected results."""
    assert element_name(argument) == expected, \
        (f"element_name({argument}) is returning a value of "
         f"{element_name(argument)}, which differs from the expected "
         f"value of {expected}.")


# (argument, expected_error)
element_name_error_table = [
    ('vegan cupcakes', ValueError),
    ('C-13-14-15-51698024', ValueError),
    (1.24, TypeError),
    ('n', ValueError),
    ('neutron', ValueError),
    (0, ValueError),
    ('H++', ValueError),
    ('t', ValueError),
    ('pb', ValueError),
    ('d', ValueError),
    ('h-3', ValueError)]


@pytest.mark.parametrize("argument, expected_error", element_name_error_table)
def test_element_name_error(argument, expected_error):
    """Test that element_name raises the expected exceptions."""
#    with pytest.raises(expected_error):
#        element_name(argument)
    pass


def test_standard_atomic_weight_value_between():
    """Test that standard_atomic_weight returns approximately the correct
    value for phosphorus."""
    assert 30.973 < standard_atomic_weight('P').value < 30.974, \
        "Incorrect standard atomic weight for phosphorus."


def test_standard_atomic_weight_unit():
    """Test that standard_atomic_weight returns a Quantity with the
    expected units."""
    assert standard_atomic_weight('Au').unit == u.u, \
        "Incorrect units from standard_atomic_weight for gold."


# (argument, expected)
standard_atomic_weight_table = [
    ('H', 1.008),
    (1, 1.008),
    ('Hydrogen', 1.008)]


@pytest.mark.parametrize("argument, expected", standard_atomic_weight_table)
def test_standard_atomic_weight(argument, expected):
    """Test that standard_atomic_weight returns the expected values for
    hydrogen."""
    assert standard_atomic_weight(argument).value == expected, \
        f"Incorrect standard_atomic_weight for {argument}."


# (argument, expected_error)
standard_atomic_weight_error_table = [
    ('H-1', ValueError),
    ("help i'm trapped in a unit test", ValueError),
    (1.1, TypeError),
    ('n', ValueError),
    ('p', ValueError),
    ('alpha', ValueError),
    ('deuteron', ValueError),
    ('tritium', ValueError),
    ('Au+', ValueError),
    ('Fe -2', ValueError),
    ('Og 2+', ValueError),
    ('h', ValueError),
    ('fe', ValueError)]


@pytest.mark.parametrize("argument, expected_error",
                         standard_atomic_weight_error_table)
def test_standard_atomic_weight_error(argument, expected_error):
    """Test that standard_atomic_weight raises the expected exceptions."""
#    with pytest.raises(expected_error, message=(
#            f"standard_atomic_weight({argument}) is not raising a "
#            "{expected_error}.")):
#        standard_atomic_weight(argument)
    pass


def test_isotope_mass_berkelium_249():
    """Test that isotope_mass returns the correct value for Bk-249."""
    assert np.isclose(isotope_mass('berkelium-249').value, 249.0749877), \
        "Incorrect isotope mass for berkelium."


def test_isotope_mass_n():
    """Test that isotope_mass returns the correct value for neutrons."""
    assert np.isclose(isotope_mass('n') / (1.008664 * u.u), 1, atol=1e-6), \
        "Incorrect isotope mass for neutrons."


def test_isotope_mass_si_30_units():
    """Test that isotope_mass returns a Quantity with the correct unit
    for Si-30."""
    assert isotope_mass('Si-30').unit == u.u, \
        "Incorrect unit for isotope mass for Si-30."


# (arg1, arg2)
isotope_mass_table = [
    (('H-1',), ('protium',)),
    (('H-1',), (1, 1)),
    (('D',), ('H-2',)),
    (('H-2',), ('deuterium',)),
    (('deuterium',), (1, 2)),
    (('T',), ('H-3',)),
    (('H-3',), ('tritium',)),
    (('tritium',), (1, 3))]


@pytest.mark.parametrize("arg1, arg2", isotope_mass_table)
def test_isotope_mass(arg1, arg2):
    """Test that isotope_mass returns equivalent results for equivalent
    arguments."""
    assert isotope_mass(*arg1) == isotope_mass(*arg2), \
        f"isotope_mass(*{arg1}) is not equivalent to isotope_mass(*{arg2})"


# (argument, expected_error)
isotope_mass_error_table = [
    ("H", ValueError),
    (1.1, TypeError),
    ('alpha', ValueError),
    ('He-4 2+', ValueError),
    ('he-4', ValueError),
    ('Fe 2+', ValueError),
    ('Fe -2', ValueError),
    ('deuteron', ValueError),
    ('triton', ValueError),
    ('H-1 +1', ValueError),
    ('H-1+', ValueError)]


@pytest.mark.parametrize("argument, expected_error", isotope_mass_error_table)
def test_isotope_mass_error(argument, expected_error):
    """Test that isotope_mass raises the expected exceptions."""
#    with pytest.raises(expected_error, warning=(
#            f"isotope_mass({argument}) is not raising a {expected_error}")):
#        isotope_mass(argument)
    pass


def test_ion_mass_for_hydrogen_with_no_mass_number():
    """Test that ion_mass does not return the proton mass when no
    mass number is specified for hydrogen.  In this case, the
    standard atomic weight should be used to account for the small
    fraction of deuterium."""
    assert ion_mass('H', Z=1) > const.m_p
    assert ion_mass('hydrogen', Z=1) > const.m_p


def test_ion_mass_unit():
    """Test that ion_mass returns a Quantity with the correct units."""
    assert ion_mass('F-19', Z=3).unit == u.kg


# (arg, kwargs)
inputs_that_should_return_proton_mass = [
    ('proton', {}),
    ('H-1+', {}),
    ('H-1 +1', {}),
    ('H-1 1+', {}),
    ('H-1', {'Z': 1}),
    ('hydrogen-1', {'Z': 1}),
    ('p+', {}),
    ('antiproton', {}),
    ('p-', {}),
]


@pytest.mark.parametrize("arg, kwargs", inputs_that_should_return_proton_mass)
def test_ion_mass_proton_mass(arg, kwargs):
    should_be_proton_mass = ion_mass(arg, **kwargs)
    assert should_be_proton_mass == const.m_p, \
        (f"ion_mass({arg}, **{kwargs}) should be returning the proton mass, "
         f"but is instead returning {should_be_proton_mass}.")


def test_ion_mass_miscellaneous_cases():
    """Test miscellaneous cases for ion_mass."""
    assert np.isclose(ion_mass(9.11e-31 * u.kg).value, 9.10938291e-31,
                      atol=1e-37)
    assert ion_mass(1.67e-27 * u.kg) == 1.67e-27 * u.kg
    assert np.isclose(ion_mass(1 * u.u).value, 1.660538921e-27, atol=1e-35)
    assert ion_mass('alpha') > ion_mass('He-3 2+')


# (arg1, kwargs1, arg2, kwargs2, expected)
equivalent_ion_mass_args = [
    ('e+', {}, 'positron', {}, const.m_e),
    ('alpha', {}, "He-4++", {}, None),
    ('alpha', {}, "helium-4 2+", {}, None),
    ('deuteron', {}, "H", {"Z": 1, "mass_numb": 2}, None),
    ('D+', {}, "H-2+", {}, None),
    ('D+', {}, "D 1+", {}, None),
    ('Deuterium+', {}, "D", {"Z": 1}, None),
    ('triton', {}, "H", {"Z": 1, "mass_numb": 3}, None),
    ('T+', {}, "H-3+", {}, None),
    ('T+', {}, "T 1+", {}, None),
    ('Tritium+', {}, "T", {"Z": 1}, None),
    ('Fe-56 1+', {}, 'Fe', {"mass_numb": 56, "Z": 1},
     ion_mass('Fe-56 1-') - 2 * const.m_e),
    ('Fe-56 +1', {}, 26, {"mass_numb": 56, "Z": 1}, None),
    # The functionality for the following test parameters may be deprecated.
    ('H', {'Z': 1, 'mass_numb': 2}, 1, {'Z': '1', 'mass_numb': '2'}, None)
]


@pytest.mark.parametrize(
    "arg1, kwargs1, arg2, kwargs2, expected", equivalent_ion_mass_args)
def test_ion_mass_equivalent_args(arg1, kwargs1, arg2, kwargs2, expected):
    """Test that """

    result1 = ion_mass(arg1, **kwargs1)
    result2 = ion_mass(arg2, **kwargs2)

    assert result1 == result2, \
        (f"ion_mass({arg1}, **{kwargs1}) = {result1}, whereas "
         f"ion_mass({arg2}, **{kwargs2}) = {result2}.  "
         f"These results are not equivalent as expected.")

    if expected is not None:
        assert result1 == result2 == expected, \
            (f"ion_mass({arg1}, **{kwargs1}) = {result1} and "
             f"ion_mass({arg2}, **{kwargs2}) = {result2}, but  "
             f"these results are not equivalent to {expected} as expected.")


# (argument, kwargs, expected_error)
ion_mass_error_table = [
    ('0g 1+', {}, ValueError),  # since it has no standard atomic weight
    ('Fe-56', {"Z": 1.4}, TypeError),
    ('n', {}, ValueError),
    ('H-1 +1', {"Z": 0}, ValueError),
    (26, {"Z": 1, "mass_numb": 'a'}, TypeError),
    (26, {"Z": 27, "mass_numb": 56}, ValueError),
    ('Og', {"Z": 1}, ValueError),
    ('Og', {"mass_numb": 296, "Z": 1}, ValueError),
    ('n', {}, ValueError),
    ('He 1+', {"mass_numb": 99}, ValueError),
    (1 * u.m, {}, u.UnitConversionError),
    ('Og', {"Z": 1}, ValueError),
    ('fe-56 1+', {}, ValueError)]


@pytest.mark.parametrize("argument, kwargs, expected_error",
                         ion_mass_error_table)
def test_ion_mass_error(argument, kwargs, expected_error):
    """Test errors that should be raised by ion_mass."""
#    with pytest.raises(expected_error, message=(
#            f"ion_mass({argument}, **{kwargs}) is not raising a "
#            f"{expected_error}.")):
#        ion_mass(argument, **kwargs)
    pass


# (argument, kwargs, expected_warning)
ion_mass_warning_table = [
    (1.6e-27 * u.kg, {}, AtomicWarning),
    (8e-25 * u.kg, {}, AtomicWarning)]


@pytest.mark.parametrize("argument, kwargs, expected_warning",
                         ion_mass_warning_table)
def test_ion_mass_warnings(argument, kwargs, expected_warning):
    """Test that ion_mass issues the expected warnings."""
    with pytest.warns(expected_warning, message=(
            f"ion_mass({argument}, **{kwargs}) is not issuing a "
            f"{expected_warning}.")):
        ion_mass(argument, **kwargs)


# (argument)
is_isotope_stable_table = [
    ('H-1',),
    (1, 1),
    ('1', '1'),
    ('N-14',),
    ('N', 14),
    ('P-31',),
    ('P', 31),
    ('p',),
    ('alpha',),
    ('Xe-124',),
    ('Fe', 56),
    ('Fe', '56'),
    ('Fe-56',),
    ('iron-56',),
    ('Iron-56',),
    (26, 56),
]


@pytest.mark.parametrize("argument", is_isotope_stable_table)
def test_is_isotope_stable(argument):
    """Test that is_isotope_stable returns True for stable isotopes."""
    assert is_isotope_stable(*argument), \
        f"is_isotope_stable is not returning True for {argument}"


# (argument)
is_isotope_stable_false_table = [
    ('Be-8',),
    ('n',),
    ('n-1',),
    (0, 1),
    ('U-235',),
    ('uranium-235',),
    ('T',),
    (4, 8),
    ('tritium',),
    ('neutron',),
    ('Pb-209',),
    ('lead-209',),
    ('Lead-209',),
    ('Pb', 209),
    (82, 209),
    ('82', '209')]


@pytest.mark.parametrize("argument", is_isotope_stable_false_table)
def test_is_isotope_stable_false(argument):
    """Test that is_isotope_stable returns False for unstable isotopes."""
    assert not is_isotope_stable(*argument), \
        f"is_isotope_stable is not returning False for {argument}"


# (argument, expected_error)
is_isotope_stable_error_table = [
    (('hydrogen-444444',), ValueError),
    (('hydrogen', 0), ValueError),
    (('',), ValueError),
    (('pb-209',), ValueError),
    (('h',), ValueError)]


@pytest.mark.parametrize("argument, expected_error",
                         is_isotope_stable_error_table)
def test_is_isotope_stable_error(argument, expected_error):
    """Test errors that should be raised by is_isotope_stable."""
#    with pytest.raises(expected_error, message=(
#            f"is_isotope_stable({argument}) is not raising a "
#            f"{expected_error}")):
#        is_isotope_stable(*argument)
    pass


def test_known_common_stable_isotopes():
    """Test that known_isotopes, common_isotopes, and stable_isotopes return
    the correct values for hydrogen."""

    known_should_be = ['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7']
    common_should_be = ['H-1', 'D']
    stable_should_be = ['He-3', 'He-4']

    assert known_isotopes('H') == known_should_be, \
        (f"known_isotopes('H') should return {known_should_be}, but is "
         f"instead returning {known_isotopes('H')}")

    assert common_isotopes('H') == common_should_be, \
        (f"common_isotopes('H') should return {common_should_be}, but is "
         f"instead returning {common_isotopes('H')}")

    assert stable_isotopes('He') == stable_should_be, \
        (f"stable_isotopes('He') should return {stable_should_be}, but is "
         f"instead returning {stable_isotopes('He')}")


def test_half_life():
    """Test that half_life returns the correct values for various isotopes."""
    assert half_life('H-1') == np.inf * u.s, "Incorrect half-life for H-1'."

    assert np.isclose(half_life('tritium').to(u.s).value,
                      (12.32 * u.yr).to(u.s).value, rtol=2e-4), \
        "Incorrect half-life for tritium."

    assert half_life('H-1').unit == 's', "Incorrect unit for H-1."
    assert half_life('tritium').unit == 's', "Incorrect unit for tritium."


def test_half_life_unstable_isotopes():
    """Test that half_life returns None and issues a warning for all isotopes
    that do not yet have half-life data."""
    for isotope in Isotopes.keys():
        if 'half_life' not in Isotopes[isotope].keys() and \
                not Isotopes[isotope].keys():
            with pytest.warns(AtomicWarning, message=(
                    f"No AtomicWarning issued for {isotope}")):
                assert half_life(isotope) is None


def test_half_life_u_220():
    """Test that half_life returns None and issues a warning for an isotope
    without half-life data.

    If half-life data is added for this isotope, then this test should fail
    and a different isotope without half-life data should be chosen instead
    until all isotopes have half-life data."""

    isotope_without_half_life_data = "U-220"

    with pytest.warns(AtomicWarning):

        try:
            half_life_isotope = half_life(isotope_without_half_life_data)
        except Exception:
            raise ValueError(f"half_life is raising an exception instead of "
                             f"issuing a AtomicWarning for an isotope without "
                             f"half-life data")

        assert half_life_isotope is None, \
            (f"half_life should return None for an isotope without half-life "
             f"data, but is returning {half_life_isotope}")


atomic_TypeError_funcs_table = [
    atomic_symbol,
    isotope_symbol,
    atomic_number,
    is_isotope_stable,
    half_life,
    mass_number,
    element_name,
    standard_atomic_weight,
    isotope_mass,
    ion_mass,
    nuclear_binding_energy,
    nuclear_reaction_energy
]

atomic_TypeError_bad_arguments = [1.1, {'cats': 'bats'}, 1 + 1j]


@pytest.mark.parametrize(
    "func, argument",
    product(atomic_TypeError_funcs_table,
            atomic_TypeError_bad_arguments))
def test_atomic_TypeErrors(func, argument):
    """Test that atomic functions raise TypeErrors when arguments of the
    incorrect time are provided."""
#    with pytest.raises(TypeError):
#        func(argument)
    pass


atomic_ValueErrors_funcs_table = [
    atomic_symbol, isotope_symbol, atomic_number,
    is_isotope_stable, half_life, mass_number,
    element_name, standard_atomic_weight]
atomic_ValueError_bad_arguments = [-1, 119, 'grumblemuffins', 'Oj']


@pytest.mark.parametrize(
    "func, argument",
    product(atomic_ValueErrors_funcs_table,
            atomic_ValueError_bad_arguments))
def test_atomic_ValueErrors(func, argument):
    """Test that atomic functions raise ValueErrors when incorrect
    arguments are provided."""
#    with pytest.raises(ValueError):
#        func(argument)
    pass


def test_known_common_stable_isotopes_cases():
    """Test that known_isotopes, common_isotopes, and stable_isotopes
    return certain isotopes that fall into these categories."""
    assert 'H-1' in known_isotopes('H')
    assert 'D' in known_isotopes('H')
    assert 'T' in known_isotopes('H')
    assert 'Be-8' in known_isotopes('Be')
    assert 'Og-294' in known_isotopes(118)
    assert 'H-1' in common_isotopes('H')
    assert 'H-4' not in common_isotopes(1)
    assert 'H-1' in stable_isotopes('H')
    assert 'D' in stable_isotopes('H')
    assert 'T' not in stable_isotopes('H')
    assert 'Fe-56' in common_isotopes('Fe', most_common_only=True)
    assert 'He-4' in common_isotopes('He', most_common_only=True)


def test_known_common_stable_isotopes_len():
    """Test that known_isotopes, common_isotopes, and stable_isotopes return
    lists of the expected lengths.

    The number of common isotopes may change if isotopic composition data
    has any significant changes.

    The number of stable isotopes may decrease slightly if some isotopes are
    discovered to be unstable but with extremely long half-lives.

    The number of known isotopes will increase as new isotopes are discovered,
    so a buffer is included in the test."""

    assert len(common_isotopes()) == 288, \
        ("The length of the list returned by common_isotopes() is "
         f"{len(common_isotopes())}, which is not the expected value.")

    assert len(stable_isotopes()) == 254, \
        ("The length of the list returned by stable_isotopes() is "
         f"{len(stable_isotopes())}, which is not the expected value.")

    assert 3352 <= len(known_isotopes()) <= 3400, \
        ("The length of the list returned by known_isotopes() is "
         f"{len(known_isotopes())}, which is not within the expected range.")


@pytest.mark.parametrize(
    "func", [common_isotopes, stable_isotopes, known_isotopes])
def test_known_common_stable_isotopes_error(func):
    """Test that known_isotopes, common_isotopes, and stable_isotopes
    raise ValueErrors for neutrons."""
#    with pytest.raises(ValueError, message=(
#            f"{func} is not raising a ValueError for neutrons.")):
#        func('n')
    pass


def test_isotopic_abundance():
    """Test that isotopic_abundance returns the appropriate values or
    raises appropriate errors for various isotopes."""
    assert isotopic_abundance('H', 1) == isotopic_abundance('protium')
    assert np.isclose(isotopic_abundance('D'), 0.000115)
    assert isotopic_abundance('Be-8') == 0.0, 'Be-8'
    assert isotopic_abundance('Li-8') == 0.0, 'Li-8'
    assert isotopic_abundance('Og', 294) == 0.0

    with pytest.raises(ValueError, message="No exception raised for neutrons"):
        isotopic_abundance('neutron')

    with pytest.raises(ValueError):
        isotopic_abundance('Og-2')


isotopic_abundance_elements = (
    atomic_number(atomic_numb) for atomic_numb in range(1, 119))

isotopic_abundance_isotopes = (
    common_isotopes(element) for element in isotopic_abundance_elements)

isotopic_abundance_sum_table = (
    (element, isotopes) for element, isotopes in
    zip(isotopic_abundance_elements, isotopic_abundance_isotopes)
    if isotopes)


@pytest.mark.parametrize("element, isotopes", isotopic_abundance_sum_table)
def test_isotopic_abundances_sum(element, isotopes):
    """Test that the sum of isotopic abundances for each element with
    isotopic abundances is one."""
    sum_of_iso_abund = sum(isotopic_abundance(isotope) for isotope in isotopes)
    assert np.isclose(sum_of_iso_abund, 1, atol=1e-6), \
        f"The sum of the isotopic abundances for {element} does not equal 1."


# (argument, expected)
charge_state_table = [
    ('H+', 1),
    ('D +1', 1),
    ('tritium 1+', 1),
    ('H-', -1),
    ('Fe -2', -2),
    ('Fe 2-', -2),
    ('N---', -3),
    ('N++', 2),
    ('alpha', 2),
    ('proton', 1),
    ('deuteron', 1),
    ('triton', 1),
    ('electron', -1),
    ('e-', -1),
    ('e+', 1),
    ('positron', 1),
    ('n', 0),
    ('neutron', 0),
    ('p-', -1),
    ('antiproton', -1),
]


@pytest.mark.parametrize("argument, expected", charge_state_table)
def test_charge_state(argument, expected):
    """Test that charge_state returns the expected results."""
    assert charge_state(argument) == expected, \
        (f"charge_state({argument}) is returning {charge_state(argument)} "
         f"which differs from the expected result of {expected}.")


# (argument, expected_error)
charge_state_error_table = [
    ('fads', ValueError),
    ('H++', ValueError),
    ('h+', ValueError),
    ('fe 1+', ValueError),
    ('d+', ValueError),
    ('Fe 29+', ValueError),
    ('H-1', ValueError),
]


@pytest.mark.parametrize("argument, expected_error", charge_state_error_table)
def test_charge_state_error(argument, expected_error):
    """Test that charge_state raises the expected exceptions."""
#    with pytest.raises(expected_error, message=(
#            f"charge_state({argument} is not raising a {expected_error}.")):
#        charge_state(argument)
    pass


# (argument, expected_warning)
charge_state_warning_table = [
    ('H---', AtomicWarning),
    ('Fe -26', AtomicWarning),
    ('Og 10-', AtomicWarning)]


@pytest.mark.parametrize("argument, expected_warning",
                         charge_state_warning_table)
def test_charge_state_warnings(argument, expected_warning):
    """Test that charge_state issues appropriate warnings."""
    with pytest.warns(expected_warning, message=(
            f"charge_state({argument}) is not issuing a {expected_warning}")):
        charge_state(argument)


def test_electric_charge():
    """Test that the results from electric_charge provide the correct
    values that have the correct characteristics."""
    assert electric_charge('p').value == 1.6021766208e-19
    assert electric_charge('p').unit == 'C'
    assert electric_charge('e').value == -1.6021766208e-19
    assert electric_charge('alpha').value == 3.2043532416e-19
    assert electric_charge('n').value == 0


# (argument, expected_error)
electric_charge_error_table = [
    ('badinput', ValueError),
    (' ', ValueError),
    ('h+', ValueError),
    ('Au 81+', ValueError)]


@pytest.mark.parametrize("argument, expected_error",
                         electric_charge_error_table)
def test_electric_charge_error(argument, expected_error):
    """Test that electric_charge raises the expected exceptions."""
#    with pytest.raises(expected_error, message=(
#            f"electric_charge({argument}) is not raising a "
#            f"{expected_error}.")):
#        electric_charge(argument)
    pass


# (argument, expected_warning)
electric_charge_warning_table = [
    ('Au 81-', AtomicWarning),
    ('H---', AtomicWarning)]


@pytest.mark.parametrize("argument, expected_warning",
                         electric_charge_warning_table)
def test_electric_charge_warning(argument, expected_warning):
    """Test that electric_charge issues the expected warnings."""
    with pytest.warns(expected_warning, message=(
            f"electric_charge({argument}) is not issuing a "
            f"{expected_warning}.")):
        electric_charge(argument)


@pytest.mark.parametrize("test_input,kwargs,expected",
                         [('n', {}, True),
                          ('n-1', {}, True),
                          ('N', {}, False),
                          ('N-1', {}, False),
                          ('N-7', {}, False),
                          ('neutron', {}, True),
                          ('James Chadwick', {}, False),
                          (0, {}, False),
                          (0, {"mass_numb": 1}, True),
                          ('n0', {}, True)])
def test_is_neutron(test_input, kwargs, expected):
    """Test that _is_neutron returns True when the argument corresponds to a
    neutron and False otherwise."""
    assert _is_neutron(test_input, **kwargs) == expected


@pytest.mark.parametrize("test_input,can_be_atom_numb,expected",
                         [('hydrogen', False, True),
                          ('hydrogen+', False, True),
                          ('hydrogen-', False, True),
                          ('hydrogen--', False, True),
                          ('H', False, True),
                          ('H+', False, True),
                          ('H-', False, True),
                          ('proton', False, True),
                          ('protium', False, True),
                          ('deuterium', False, True),
                          ('tritium', False, True),
                          ('triton', False, True),
                          ('deuteron', False, True),
                          ('h', False, False),
                          ('D', False, True),
                          ('D+', False, True),
                          ('H-2', False, True),
                          ('H-2+', False, True),
                          ('H-2 1+', False, True),
                          ('H-2 +1', False, True),
                          ('H-3 -1', False, True),
                          ('He', False, False),
                          ('H-1', False, True),
                          ('H-7', False, True),
                          ('antiproton', False, False),
                          (1, True, True),
                          (1, False, False),
                          ('p-', False, False)])
def test_is_hydrogen(test_input, can_be_atom_numb, expected):
    """Test that _is_hydrogen returns True if the argument corresponds to
    hydrogen and False otherwise."""
    assert _is_hydrogen(test_input,
                        can_be_atomic_number=can_be_atom_numb) == expected


@pytest.mark.parametrize("test_input,kwargs,expected_error",
                         [('H 2+', {}, ValueError),
                          ('D++', {}, ValueError)])
def test_is_hydrogen_errors(test_input, kwargs, expected_error):
    """Test that _is_hydrogen raises the expected exceptions."""
    with pytest.raises(expected_error):
        _is_hydrogen(test_input, **kwargs)


@pytest.mark.parametrize("test_input,expected",
                         [('e-', True),
                          ('e+', False),
                          ('Electron', True),
                          ('electron', True),
                          ('positron', False),
                          ('p', False),
                          ('E', False),
                          ('E-', False),
                          ('beta', False),
                          (-1, False)])
def test_is_electron(test_input, expected):
    """Test that _is_electron returns True if the argument corresponds to
    an electron and False otherwise."""
    assert _is_electron(test_input) == expected


@pytest.mark.parametrize("test_input,expected",
                         [('e-', False),
                          ('e+', True),
                          ('Electron', False),
                          ('electron', False),
                          ('positron', True),
                          ('p', False),
                          ('E', False),
                          ('E-', False),
                          ('beta', False),
                          (1, False)])
def test_is_positron(test_input, expected):
    """Test that _is_positron returns True if the argument corresponds to
    positron and False otherwise."""
    assert _is_positron(test_input) == expected


@pytest.mark.parametrize("test_input,kwargs,expected",
                         [('p', {}, True),
                          ('p+', {}, True),
                          ('hydrogen-1+', {}, True),
                          ('H-1 1+', {}, True),
                          ('H-1', {}, False),
                          ('H', {}, False),
                          ('p-', {}, False),
                          ('antiproton', {}, False),
                          ('Antiproton', {}, False),
                          ('proton', {}, True),
                          ('Proton', {}, True),
                          ('P', {}, False),
                          ('P+', {}, False),
                          (1, {}, False),
                          (1, {"mass_numb": 1, "Z": 1}, True),
                          ('H', {"mass_numb": 1, "Z": 1}, True),
                          ('H-1', {"Z": 1}, True),
                          ('H', {"Z": 1}, False),
                          ('H-1', {"Z": 0}, False)])
def test_is_proton(test_input, kwargs, expected):
    """Test that _is_proton returns True if the argument corresponds to
    a proton and False otherwise."""
    assert _is_proton(test_input, **kwargs) == expected


@pytest.mark.parametrize("test_input,expected",
                         [('e-', False),
                          ('e+', False),
                          ('Electron', False),
                          ('electron', False),
                          ('positron', False),
                          ('p', False),
                          ('E', False),
                          ('E-', False),
                          ('beta', False),
                          ('p-', True),
                          ('Antiproton', True),
                          ('antiproton', True),
                          ('p+', False),
                          ('p--', False),
                          ('P-', False),
                          (57, False)])
def test_is_antiproton(test_input, expected):
    """Test that _is_antiproton returns True if the argument corresponds to
    an antiproton and False otherwise."""
    assert _is_antiproton(test_input) == expected


@pytest.mark.parametrize("test_input,expected",
                         [('e-', False),
                          ('e+', False),
                          ('Electron', False),
                          ('electron', False),
                          ('positron', False),
                          ('p', False),
                          ('E', False),
                          ('E-', False),
                          ('beta', False),
                          ('p-', False),
                          ('Antiproton', False),
                          ('antiproton', False),
                          ('p+', False),
                          ('p--', False),
                          ('P-', False),
                          (57, False),
                          ('alpha', True),
                          ('He-4 2+', True),
                          ('He-4++', True),
                          ('He-3 2+', False),
                          ('He-5 2+', False),
                          ('Helium-4 +2', True),
                          ('Helium-4 -2', False),
                          ('He-4', False),
                          ('helium', False),
                          ('He', False),
                          ('Fe-56', False),
                          ('Fe', False),
                          ('he-4 2+', False)])
def test_is_alpha(test_input, expected):
    """Test that _is_alpha returns True if the argument corresponds to
    an alpha particle and False otherwise."""
    assert _is_alpha(test_input) == expected


@pytest.mark.parametrize("test_input,expected_newarg,expected_Z",
                         [('H', 'H', None),
                          ('H+', 'H', 1),
                          ('D 1+', 'D', 1),
                          ('alpha', 'alpha', 2),
                          ('Fe', 'Fe', None),
                          ('Titanium', 'Titanium', None),
                          ('N-7+++', 'N-7', 3),
                          ('H-1-', 'H-1', -1),
                          ('He-4-', 'He-4', -1)])
def test_extract_charge_state(test_input, expected_newarg, expected_Z):
    """Test that _extract_charge_state returns the expected values."""
    new_symbol, new_Z = _extract_charge_state(test_input)
    assert new_symbol == expected_newarg, \
        (f"_extract_charge_state should return {expected_newarg} as "
         f" its first argument, but is instead returning {new_symbol}")
    assert new_Z == expected_Z, \
        (f"_extract_charge_state should return {expected_Z} as its second"
         f"argument, but is instead returning {new_Z}.")


@pytest.mark.parametrize("test_input,expected_error",
                         [('H-1-+-+', ValueError),
                          ('H ++', ValueError),
                          ('Fe +21+', ValueError)])
def test_extract_charge_state_errors(test_input, expected_error):
    """Test that _extract_charge_state raises the expected exceptions."""
#    with pytest.raises(expected_error):
#        _extract_charge_state(test_input)
    pass


@pytest.mark.parametrize("test_input,expected_warning",
                         [('H-1----', AtomicWarning),
                          ('Fe -4', AtomicWarning),
                          ('lead 4-', AtomicWarning)])
def test_extract_charge_state_warnings(test_input, expected_warning):
    """Test that _extract_charge_state issues the expected warnings."""
    with pytest.warns(expected_warning):
        _extract_charge_state(test_input)
