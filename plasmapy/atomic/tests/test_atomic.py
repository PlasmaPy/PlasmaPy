from itertools import product
from astropy import units as u, constants as const
import numpy as np

from ..atomic import (atomic_symbol, isotope_symbol, atomic_number,
                      mass_number, element_name, standard_atomic_weight,
                      isotope_mass, ion_mass, is_isotope_stable,
                      half_life, known_isotopes, common_isotopes,
                      stable_isotopes, isotopic_abundance, charge_state,
                      electric_charge, Elements, Isotopes,
                      _is_neutron, _is_hydrogen, _is_electron,
                      _is_positron, _is_antiproton, _is_alpha,
                      _extract_charge_state, _is_proton)

from ..nuclear import (nuclear_binding_energy, nuclear_reaction_energy)

import pytest

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
    'argument, expected', atomic_symbol_table
)
def test_atomic_symbol(argument, expected):
    assert atomic_symbol(argument) == expected


# (argument, expected_error)
atomic_symbol_error_table = [
    ('H-0', ValueError),
    (3.141592653589793238462643383279502884, TypeError),
    ('Og-294b', ValueError),
    ('H-934361079326356530741942970523610389', ValueError),
    ('Fe 2+4', ValueError),
    ('Fe+24', ValueError),
    ('Fe +59', ValueError),
    ('C++++++++++++++++', ValueError),
    ('C-++++', ValueError),
    ('neutron', ValueError),
    ('n', ValueError),
    ('n-1', ValueError),
    ('h', ValueError),
    ('d', ValueError),
    ('he', ValueError),
    ('au', ValueError)]


@pytest.mark.parametrize(
    'argument, expected_error', atomic_symbol_error_table)
def test_atomic_symbol_error(argument, expected_error):
    with pytest.raises(expected_error):
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
    "argument, expected", isotope_symbol_table)
def test_isotope_symbol(argument, expected):
    assert isotope_symbol(*argument) == expected


# (argument, kwargs, expected_error)
isotope_symbol_error_table = [
    ('Md-260', {"mass_numb": 261}, ValueError),
    ('protium', {"mass_numb": 2}, ValueError),
    ('alpha', {"mass_numb": 3}, ValueError),
    ('O-18', {"mass_numb": 19}, ValueError),
    ('lead-209', {"mass_numb": 511}, ValueError),
    ('He-1', {}, ValueError),
    (24, {"mass_numb": 23}, ValueError),
    ('H', {"mass_numb": 0}, ValueError),
    ('H-1', {"mass_numb": 2}, ValueError),
    ('P', {}, ValueError),
    (1, {}, ValueError),
    (4, {}, ValueError),
    ('hydrogen-444444', {}, ValueError),
    ('Fe', {"mass_numb": 2.1}, TypeError),
    ('He', {"mass_numb": 'c'}, TypeError),
    ('He-3', {"mass_numb": 4}, ValueError),
    ('alpha', {"mass_numb": 3}, ValueError),
    ('D', {"mass_numb": 3}, ValueError),
    ('T', {"mass_numb": 2}, ValueError),
    ('Fe', {"mass_numb": None}, ValueError),
    ('d', {}, ValueError),
    ('h-3', {}, ValueError),
    ('h', {}, ValueError),
    ('d+', {}, ValueError),
]


@pytest.mark.parametrize(
    "argument, kwargs, expected_error", isotope_symbol_error_table)
def test_isotope_symbol_error(argument, kwargs, expected_error):
    with pytest.raises(expected_error):
        isotope_symbol(argument, **kwargs)


# (argument, kwargs, expected_warning)
isotope_symbol_warning_table = [
    ('H-1', {"mass_numb": 1}, UserWarning),
    ('H-2', {"mass_numb": 2}, UserWarning),
    ('T', {"mass_numb": 3}, UserWarning),
    ('Li-6', {"mass_numb": 6}, UserWarning),
    ('lithium-6', {"mass_numb": 6}, UserWarning),
    ('alpha', {"mass_numb": 4}, UserWarning),
    ('p', {"mass_numb": 1}, UserWarning)]


@pytest.mark.parametrize(
    "argument, kwargs, expected_warning", isotope_symbol_warning_table)
def test_isotope_symbol_error(argument, kwargs, expected_warning):
    with pytest.warns(expected_warning):
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
    assert atomic_number(argument) == expected


# (argument, expected_error)
atomic_number_error_table = [
    ('H-3934', ValueError),
    ('C-12b', ValueError),
    (-1.5, TypeError),
    ('n', ValueError),
    ('n-1', ValueError),
    ('neutron', ValueError),
    ('Neutron', ValueError),
    ('d', ValueError),
    ('t', ValueError),
    ('s-36', ValueError)]


@pytest.mark.parametrize(
    "argument, expected_error", atomic_number_error_table)
def test_atomic_number_error(argument, expected_error):
    with pytest.raises(expected_error):
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
    assert mass_number(isotope) == expected


# (argument, expected_error)
mass_number_error_table = [
    ('H-359', ValueError),
    ('C-12b', ValueError),
    (-1.5, Exception),
    ('N-13+-+-', ValueError),
    ('h-3', ValueError)]


@pytest.mark.parametrize(
    "argument, expected_error", mass_number_error_table)
def test_mass_number_error(argument, expected_error):
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
    assert element_name(argument) == expected


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

    with pytest.raises(expected_error):
        element_name(argument)


def test_standard_atomic_weight_not_none():
    assert standard_atomic_weight('N') is not None


def test_standard_atomic_weight_value_between():
    assert 30.973 < standard_atomic_weight('P').value < 30.974


def test_standard_atomic_weight_unit():
    assert standard_atomic_weight('Au').unit == u.u


# (argument, expected)
standard_atomic_weight_table = [
    ('H', 1.008),
    (1, 1.008),
    ('Hydrogen', 1.008)]


@pytest.mark.parametrize("argument, expected", standard_atomic_weight_table)
def test_standard_atomic_weight(argument, expected):
    assert standard_atomic_weight(argument).value == expected


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
    with pytest.raises(expected_error):
        standard_atomic_weight(argument)


def test_isotope_mass_berkelium_249():
    assert np.isclose(isotope_mass('berkelium-249').value, 249.0749877)


def test_isotope_mass_n():
    assert np.isclose(isotope_mass('n') / (1.008664 * u.u), 1, atol=1e-6)


def test_isotope_mass_si_30_units():
    assert isotope_mass('Si-30').unit == u.u


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
    assert isotope_mass(*arg1) == isotope_mass(*arg2)


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
    with pytest.raises(expected_error):
        isotope_mass(argument)


def test_ion_mass_hydrogen():
    assert ion_mass('H') > const.m_p, "Use standard_atomic_weight of 'H'"
    assert ion_mass('hydrogen') > const.m_p


def test_ion_mass_unit():
    assert ion_mass('F-19', Z=3).unit == u.kg
    assert ion_mass('F-19', Z='3').unit == u.kg


def test_ion_mass():
    assert ion_mass('proton') == const.m_p
    assert ion_mass('e+') == ion_mass('positron') == const.m_e
    assert np.isclose(ion_mass('alpha') / ion_mass('He-4', 2), 1.0)
    assert ion_mass('protium') == const.m_p
    assert ion_mass('Ne-22', 2) == 21.991385114 * u.u - 2 * const.m_e
    assert ion_mass('H-1+') == const.m_p
    assert ion_mass('He+') == ion_mass('He')
    assert ion_mass('He 1+') == ion_mass('He')
    assert ion_mass('He-4 2+') == ion_mass('alpha')
    assert np.isclose(ion_mass('Fe 1-').value,
                      (ion_mass('Fe 1+') + 2 * const.m_e).value, rtol=1e-14)
    assert np.isclose(ion_mass('Fe-56 1-').value,
                      (ion_mass('Fe-56 1+') + 2 * const.m_e).value, rtol=1e-14)
    assert np.isclose((ion_mass('Fe-56 1+')).value,
                      (ion_mass('Fe', Z=1, mass_numb=56)).value)
    assert np.isclose((ion_mass('Fe-56 1+')).value,
                      (ion_mass(26, Z=1, mass_numb=56)).value)
    assert ion_mass(1, Z=1, mass_numb=1) == ion_mass('p')
    assert ion_mass('deuteron') == ion_mass('D +1')
    assert ion_mass('T', Z=1) == ion_mass('T +1')
    assert ion_mass('Fe', mass_numb=56) == ion_mass('Fe', mass_numb='56')
    assert np.isclose(ion_mass(9.11e-31 * u.kg).value, 9.10938291e-31,
                      atol=1e-37)
    assert ion_mass(1.67e-27 * u.kg) == 1.67e-27 * u.kg
    assert np.isclose(ion_mass(1 * u.u).value, 1.660538921e-27, atol=1e-35)
    assert ion_mass('alpha') > ion_mass('He-3 2+')
    assert ion_mass('antiproton') == ion_mass('p-') == ion_mass('p+')


# (argument, kwargs, expected_error)
ion_mass_error_table = [
    ('0g', {}, ValueError),  # since it has no standard atomic weight
    ('Fe-56', {"Z": 1.4}, TypeError),
    ('n', {}, ValueError),
    ('H-1 +1', {"Z": 0}, ValueError),
    (26, {"Z": 1, "mass_numb": 'a'}, TypeError),
    (26, {"Z": 27, "mass_numb": 56}, ValueError),
    ('Og', {"Z": 1}, ValueError),
    ('Og', {"mass_numb": 296, "Z": 1}, ValueError),
    ('n', {}, ValueError),
    ('He', {"mass_numb": 99}, ValueError),
    (1 * u.m, {}, u.UnitConversionError),
    ('Og', {"Z": 1}, ValueError),
    ('fe-56 1+', {}, ValueError)]


@pytest.mark.parametrize("argument, kwargs, expected_error",
                         ion_mass_error_table)
def test_ion_mass_error(argument, kwargs, expected_error):
    with pytest.raises(expected_error):
        ion_mass(argument, **kwargs)


# (argument, kwargs, expected_warning)
ion_mass_warning_table = [
    (1.6e-27 * u.kg, {}, UserWarning),
    (8e-25 * u.kg, {}, UserWarning)]


@pytest.mark.parametrize("argument, kwargs, expected_warning",
                         ion_mass_warning_table)
def test_ion_mass_warnings(argument, kwargs, expected_warning):
    with pytest.warns(expected_warning):
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
    (26, 56)]


@pytest.mark.parametrize("argument", is_isotope_stable_table)
def test_is_isotope_stable(argument):
    assert is_isotope_stable(*argument)


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
    assert not is_isotope_stable(*argument)


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
    with pytest.raises(expected_error):
        is_isotope_stable(*argument)


def test_isotope_calls():
    assert known_isotopes('H') == \
        ['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7']
    assert common_isotopes('H') == ['H-1', 'D']
    assert stable_isotopes('He') == ['He-3', 'He-4']


def test_half_life():
    assert half_life('H-1') == np.inf * u.s
    assert np.isclose(half_life('tritium').to(u.s).value,
                      (12.32 * u.yr).to(u.s).value, rtol=2e-4)
    assert half_life('H-1').unit == 's'
    assert half_life('tritium').unit == 's'


def test_half_life_unstable_isotopes():
    for isotope in Isotopes.keys():  # unstable isotopes w/o half-life data
        if 'half_life' not in Isotopes[isotope].keys() and \
                not Isotopes[isotope].keys():
            with pytest.warns(UserWarning):
                assert half_life(isotope) is None


def test_half_life_u_220():
    with pytest.warns(UserWarning):
        half_life('U-220')
        assert half_life('U-220') is None, \
            ("If half-life data is added for this isotope, then this test "
             "*should* fail and a different isotope without half-life data "
             "should be chosen instead")


atomic_TypeError_funcs_table = [
    atomic_symbol, isotope_symbol, atomic_number,
    is_isotope_stable, half_life, mass_number,
    element_name, standard_atomic_weight, isotope_mass,
    ion_mass, nuclear_binding_energy,
    nuclear_reaction_energy]
atomic_TypeError_bad_arguments = [1.1, {'cats': 'bats'}, 1 + 1j]


@pytest.mark.parametrize(
    "function, argument",
    product(atomic_TypeError_funcs_table,
            atomic_TypeError_bad_arguments))
def test_atomic_TypeErrors(function, argument):
    with pytest.raises(TypeError):
        function(argument)


atomic_ValueErrors_funcs_table = [
    atomic_symbol, isotope_symbol, atomic_number,
    is_isotope_stable, half_life, mass_number,
    element_name, standard_atomic_weight]
atomic_ValueError_bad_arguments = [-1, 119, 'grumblemuffins', 'Oj']


@pytest.mark.parametrize(
    "function, argument",
    product(atomic_ValueErrors_funcs_table,
            atomic_ValueError_bad_arguments))
def test_atomic_ValueErrors(function, argument):
    with pytest.raises(ValueError):
        function(argument)


def test_known_common_stable_isotopes():
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

    assert len(common_isotopes()) == 288, \
        ("The number of common isotopes may change if isotopic composition "
         "data has any significant changes.")

    assert len(stable_isotopes()) == 254, \
        ("The number of stable isotopes may decrese slightly if some are "
         "discovered to be unstable but with extremely long half-lives.")

    assert 3351 <= len(known_isotopes()) <= 3600, \
        ("The number of known isotopes ")

    assert 'Fe-56' in common_isotopes('Fe', most_common_only=True)
    assert 'He-4' in common_isotopes('He', most_common_only=True)


@pytest.mark.parametrize(
    "func", [common_isotopes, stable_isotopes, known_isotopes])
def test_known_common_stable_isotopes_error(func):
    with pytest.raises(ValueError):
        func('n')


def test_isotopic_abundance():

    assert isotopic_abundance('H', 1) == isotopic_abundance('protium')
    assert np.isclose(isotopic_abundance('D'), 0.000115)
    assert isotopic_abundance('Be-8') == 0.0, 'Be-8'
    assert isotopic_abundance('Li-8') == 0.0, 'Li-8'
    assert isotopic_abundance('Og', 294) == 0.0

    with pytest.raises(ValueError):
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
    # Check that the sum of isotopic abundances for each element is one
    sum_of_iso_abund = sum(isotopic_abundance(isotope) for isotope in isotopes)
    assert np.isclose(sum_of_iso_abund, 1, atol=1e-6), \
        "Problem with: " + element


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
    ('neutron', 0)]


@pytest.mark.parametrize("argument, expected", charge_state_table)
def test_charge_state(argument, expected):
    assert charge_state(argument) == expected


# (argument, expected_error)
charge_state_error_table = [
    ('fads', ValueError),
    ('H++', ValueError),
    ('h+', ValueError),
    ('fe 1+', ValueError),
    ('d+', ValueError),
    ('Fe 29+', ValueError)]


@pytest.mark.parametrize("argument, expected_error", charge_state_error_table)
def test_charge_state_error(argument, expected_error):
    with pytest.raises(expected_error):
        charge_state(argument)


# (argument, expected_warning)
charge_state_warning_table = [
    ('H---', UserWarning),
    ('Fe -26', UserWarning),
    ('Og 10-', UserWarning)]


@pytest.mark.parametrize("argument, expected_warning",
                         charge_state_warning_table)
def test_charge_state_warnings(argument, expected_warning):
    with pytest.warns(expected_warning):
        charge_state(argument)


def test_electric_charge():
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
    with pytest.raises(expected_error):
        electric_charge(argument)


# (argument, expected_warning)
electric_charge_warning_table = [
    ('Au 81-', UserWarning),
    ('H---', UserWarning)]


@pytest.mark.parametrize("argument, expected_warning",
                         electric_charge_warning_table)
def test_electric_charge_error(argument, expected_warning):
    with pytest.warns(expected_warning):
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
    assert _is_hydrogen(test_input,
                        can_be_atomic_number=can_be_atom_numb) == expected


@pytest.mark.parametrize("test_input,kwargs,expected_error",
                         [('H 2+', {}, ValueError),
                          ('D++', {}, ValueError)])
def test_is_hydrogen_errors(test_input, kwargs, expected_error):
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
    assert _is_positron(test_input) == expected


@pytest.mark.parametrize("test_input,expected",
                         [('p', True),
                          ('p+', True),
                          ('hydrogen-1+', True),
                          ('H-1 1+', True),
                          ('H-1', False),
                          ('H', False),
                          ('p-', False),
                          ('antiproton', False),
                          ('Antiproton', False),
                          ('proton', True),
                          ('Proton', True),
                          ('P', False),
                          ('P+', False),
                          (1, False),
                          ])
def test_is_proton(test_input, expected):
    assert _is_proton(test_input) == expected


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
    new_symbol, new_Z = _extract_charge_state(test_input)
    assert new_symbol == expected_newarg, str(test_input)
    assert new_Z == expected_Z


@pytest.mark.parametrize("test_input,expected_error",
                         [('H-1-+-+', ValueError),
                          ('H ++', ValueError),
                          ('Fe +21+', ValueError),
                          ('lead-12345', ValueError)])
def test_extract_charge_state_errors(test_input, expected_error):
    with pytest.raises(expected_error):
        _extract_charge_state(test_input)


@pytest.mark.parametrize("test_input,expected_warning",
                         [('H-1----', UserWarning),
                          ('Fe -4', UserWarning),
                          ('lead 4-', UserWarning)])
def test_extract_charge_state_errors(test_input, expected_warning):
    with pytest.warns(expected_warning):
        _extract_charge_state(test_input)
