from itertools import product
from astropy import units as u, constants as const
import numpy as np
from ..atomic import (atomic_symbol, isotope_symbol, atomic_number,
                      mass_number, element_name, standard_atomic_weight,
                      isotope_mass, ion_mass, nuclear_binding_energy,
                      energy_from_nuclear_reaction, is_isotope_stable,
                      half_life, known_isotopes, common_isotopes,
                      stable_isotopes, isotopic_abundance, charge_state,
                      electric_charge, Elements, Isotopes)

import pytest

# (argument, expected)
atomic_symbol_table = [
    (1, 'H'),
    ('H', 'H'),
    ('h', 'H'),
    ('d', 'H'),
    ('T', 'H'),
    ('deuterium', 'H'),
    ('deuteron', 'H'),
    ('Tritium', 'H'),
    ('triton', 'H'),
    ('H-3', 'H'),
    ('Hydrogen-3', 'H'),
    ('helium', 'He'),
    (2, 'He'),
    ('alpha', 'He'),
    ('he', 'He'),
    ('gold', 'Au'),
    ('Gold', 'Au'),
    ('au', 'Au'),
    (79, 'Au'),
    ('79', 'Au'),
    ('p', 'H'),
    ('P', 'P'),
    (118, 'Og'),
    ('neutron', 'n'),
    ('n-1', 'n'),
    ('N-14', 'N'),
    ('n', 'n'),
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
    ('C-++++', ValueError)]


@pytest.mark.parametrize(
    'argument, expected_error', atomic_symbol_error_table
)
def test_atomic_symbol_error(argument, expected_error):
    with pytest.raises(expected_error):
        atomic_symbol(argument)


# (argument, expected)
isotype_symbol_table = [
    (('He', 4), 'He-4'),
    (('helium-4',), 'He-4'),
    (('H-2',), 'D'),
    (('d',), 'D'),
    (('Deuterium',), 'D'),
    (('deuterium',), 'D'),
    (('deuteron',), 'D'),
    (('tritium',), 'T'),
    (('triton',), 'T'),
    (('Hydrogen-3',), 'T'),
    (('hydrogen-3',), 'T'),
    (('h-3',), 'T'),
    (('H-3',), 'T'),
    ((1, 2), 'D'),
    (('h', 2), 'D'),
    (('Hydrogen', 3), 'T'),
    (('tritium',), 'T'),
    (('H', 2), 'D'),
    (('Alpha',), 'He-4'),
    (('alpha',), 'He-4'),
    ((79, 197), 'Au-197'),
    (('p',), 'H-1'),
    (('beryllium-8',), 'Be-8'),
    (('neutron',), 'n'),
    (('n',), 'n'),
    ((0, 1), 'n'),
    (('n-1',), 'n'),
    (('N-13',), 'N-13'),
    (('p',), 'H-1'),
    (('proton',), 'H-1'),
    (('protium',), 'H-1'),
    (('N-13 2+',), 'N-13'),
    (('d+',), 'D'),
    (('Hydrogen-3 +1',), 'T')]


@pytest.mark.parametrize(
    "argument, expected", isotype_symbol_table)
def test_isotope_symbol(argument, expected):
    assert isotope_symbol(*argument) == expected


# (argument, kwargs, expected_error)
isotope_symbol_error_table = [
    ('H-1', {"mass_numb": 1}, UserWarning),
    ('H-2', {"mass_numb": 2}, UserWarning),
    ('T', {"mass_numb": 3}, UserWarning),
    ('Li-6', {"mass_numb": 6}, UserWarning),
    ('lithium-6', {"mass_numb": 6}, UserWarning),
    ('alpha', {"mass_numb": 4}, UserWarning),
    ('p', {"mass_numb": 1}, UserWarning),
    ('n', {"mass_numb": 1}, UserWarning),
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
    ('Fe', {"mass_numb": None}, ValueError)]


@pytest.mark.parametrize(
    "argument, kwargs, expected_error", isotope_symbol_error_table)
def test_isotope_symbol_error(argument, kwargs, expected_error):
    with pytest.raises(expected_error):
        isotope_symbol(argument, **kwargs)


# (argument, expected)
atomic_number_table = [
    ('H', 1),
    ('d', 1),
    ('D', 1),
    ('deuterium', 1),
    ('Deuterium', 1),
    ('t', 1),
    ('tritium', 1),
    ('p', 1),
    ('P', 15),
    ('Alpha', 2),
    ('C-12', 6),
    ('s-36', 16),
    ('Argon', 18),
    ('protium', 1),
    ('H-3', 1),
    ('p+', 1),
    ('Be-8', 4),
    ('n', 0),
    ('n-1', 0),
    ('Neutron', 0),
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
    (-1.5, TypeError)]


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
    ('N-13+-+-', ValueError)]


@pytest.mark.parametrize(
    "argument, expected_error", mass_number_error_table)
def test_mass_number_error(argument, expected_error):
    with pytest.raises(expected_error):
        mass_number(argument)


# (argument, expected)
element_name_table = [
    ('D', 'hydrogen'),
    ('deuterium', 'hydrogen'),
    ('t', 'hydrogen'),
    ('Au', 'gold'),
    ('pb', 'lead'),
    ('alpha', 'helium'),
    ('helium-4', 'helium'),
    ('H-2', 'hydrogen'),
    ('d', 'hydrogen'),
    ('Deuterium', 'hydrogen'),
    ('Hydrogen-3', 'hydrogen'),
    ('hydrogen-3', 'hydrogen'),
    ('h-3', 'hydrogen'),
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
    ('n', 'neutron'),
    (0, 'neutron'),
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
    ('H++', ValueError)]


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
    ('wrong input', ValueError),
    (1.1, TypeError),
    ('n', ValueError),
    ('p', ValueError),
    ('alpha', ValueError),
    ('deuteron', ValueError),
    ('tritium', ValueError),
    ('Au+', ValueError),
    ('Fe -2', ValueError),
    ('Og 2+', ValueError)]


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
    ('p', ValueError),
    ('Fe 2+', ValueError),
    ('Fe -2', ValueError),
    ('deuteron', ValueError),
    ('triton', ValueError),
    ('alpha', ValueError),
    ('p', ValueError)]


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
    (1.6e-27 * u.kg, {}, UserWarning),
    (8e-25 * u.kg, {}, UserWarning),
    (1 * u.m, {}, u.UnitConversionError),
    ('Og', {"Z": 1}, ValueError)]


@pytest.mark.parametrize("argument, kwargs, expected_error",
                         ion_mass_error_table)
def test_ion_mass_error(argument, kwargs, expected_error):
    with pytest.raises(expected_error):
        ion_mass(argument, **kwargs)


def test_nuclear_binding_energy():
    assert nuclear_binding_energy('p') == 0
    assert nuclear_binding_energy('n') == 0
    assert nuclear_binding_energy('He-4') == nuclear_binding_energy('alpha') \
        == nuclear_binding_energy('He', 4)


def test_nuclear_binding_energy_D_T():
    before = nuclear_binding_energy("D") + nuclear_binding_energy("T")
    after = nuclear_binding_energy("alpha")
    E_in_MeV = (after - before).to(u.MeV).value  # D + T --> alpha + n + E
    assert np.isclose(E_in_MeV, 17.58, rtol=0.01)


# (argument, expected_error)
nuclear_binding_energy_table = [
    ("H", ValueError),
    (1.1, TypeError)]


@pytest.mark.parametrize("argument, expected_error",
                         nuclear_binding_energy_table)
def test_nuclear_binding_energy_error(argument, expected_error):
    with pytest.raises(expected_error):
        nuclear_binding_energy(argument)


def test_energy_from_nuclear_reaction():
    reaction1 = 'D + T --> alpha + n'
    reaction2 = 'T + D -> n + alpha'
    released_energy1 = energy_from_nuclear_reaction(reaction1)
    released_energy2 = energy_from_nuclear_reaction(reaction2)
    assert np.isclose(released_energy1.to(u.MeV).value, 17.58, rtol=0.01)
    assert released_energy1 == released_energy2


def test_energy_from_nuclear_reaction_triple_alpha():
    triple_alpha1 = 'alpha + He-4 --> Be-8'
    triple_alpha2 = 'Be-8 + alpha --> carbon-12'
    energy_triplealpha1 = energy_from_nuclear_reaction(triple_alpha1)
    energy_triplealpha2 = energy_from_nuclear_reaction(triple_alpha2)
    assert np.isclose(energy_triplealpha1.to(u.keV).value, -91.8, atol=0.1)
    assert np.isclose(energy_triplealpha2.to(u.MeV).value, 7.367, atol=0.1)


def test_energy_from_nuclear_reaction_alpha_decay():
    alpha_decay_example = 'U-238 --> Th-234 + alpha'
    energy_alpha_decay = energy_from_nuclear_reaction(alpha_decay_example)
    assert np.isclose(energy_alpha_decay.to(u.MeV).value, 4.26975, atol=1e-5)


def test_energy_from_nuclear_reaction_triple_alpha_r():
    triple_alpha1_r = '4He-4 --> 2Be-8'
    energy_triplealpha1_r = energy_from_nuclear_reaction(triple_alpha1_r)
    assert np.isclose(energy_triplealpha1_r.to(u.keV).value,
                      -91.8 * 2, atol=0.1)


# (reaction, expected_error)
energy_from_nuclear_reaction_error_table = [
    ('H + H --> H', ValueError),
    (1, TypeError),
    ('H-1 + H-1 --> H-1', ValueError),
    ("I didn't like unstable eigenfunctions "
     "at first, but then they grew on me", ValueError)
]


@pytest.mark.parametrize(
    "reaction, expected_error", energy_from_nuclear_reaction_error_table)
def test_energy_from_nuclear_reaction_error(reaction, expected_error):
    with pytest.raises(expected_error):
        energy_from_nuclear_reaction(reaction)


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
    ('fe-56',),
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
    ('pb-209',),
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
    (('',), ValueError)]


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
            with pytest.raises(UserWarning):
                assert half_life(isotope) is None


def test_half_lift_u_220():
    with pytest.raises(UserWarning):
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
    energy_from_nuclear_reaction]
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
    ('positron', 1)]


@pytest.mark.parametrize("argument, expected", charge_state_table)
def test_charge_state(argument, expected):
    assert charge_state(argument) == expected


# (argument, expected_error)
charge_state_error_table = [
    ('fads', ValueError),
    ('H++', ValueError),
    ('Fe 29+', ValueError),
    ('H---', UserWarning),
    ('Fe -26', UserWarning),
    ('Og 10-', UserWarning)]


@pytest.mark.parametrize("argument, expected_error", charge_state_error_table)
def test_charge_state_error(argument, expected_error):
    with pytest.raises(expected_error):
        charge_state(argument)


def test_electric_charge():
    assert electric_charge('p').value == 1.6021766208e-19
    assert electric_charge('p').unit == 'C'
    assert electric_charge('e').value == -1.6021766208e-19
    assert electric_charge('alpha').value == 3.2043532416e-19


# (argument, expected_error)
electric_charge_error_table = [
    ('badinput', ValueError),
    (' ', ValueError),
    ('Au 81+', ValueError),
    ('Au 81-', UserWarning),
    ('H---', UserWarning)]


@pytest.mark.parametrize("argument, expected_error",
                         electric_charge_error_table)
def test_electric_charge_error(argument, expected_error):
    with pytest.raises(expected_error):
        electric_charge(argument)
