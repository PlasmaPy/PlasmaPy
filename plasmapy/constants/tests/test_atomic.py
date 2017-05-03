from astropy import units as u, constants as const
import numpy as np

from ..atomic import (element_symbol, isotope_symbol, atomic_number,
                      mass_number, element_name, standard_atomic_weight,
                      isotope_mass, ion_mass, nuclear_binding_energy,
                      energy_from_nuclear_reaction, is_isotope_stable,
                      half_life, known_isotopes, common_isotopes,
                      stable_isotopes, Elements, Isotopes)

import pytest


def test_element_symbol():
    assert element_symbol(1) == 'H'
    assert element_symbol('H') == 'H'
    assert element_symbol('h') == 'H'
    assert element_symbol('d') == 'H'
    assert element_symbol('T') == 'H'
    assert element_symbol('deuterium') == 'H'
    assert element_symbol('Tritium') == 'H'
    assert element_symbol('H-3') == 'H'
    assert element_symbol('Hydrogen-3') == 'H'
    assert element_symbol('helium') == 'He'
    assert element_symbol('Helium') == 'He'
    assert element_symbol(2) == 'He'
    assert element_symbol('alpha') == 'He'
    assert element_symbol('he') == 'He'
    assert element_symbol('gold') == 'Au'
    assert element_symbol('Gold') == 'Au'
    assert element_symbol('au') == 'Au'
    assert element_symbol(79) == 'Au'
    assert element_symbol('79') == 'Au'
    assert element_symbol('p') == 'H'
    assert element_symbol('P') == 'P'
    assert element_symbol(118) == 'Og'
    assert element_symbol('neutron') == 'n'
    assert element_symbol('n-1') == 'n'
    assert element_symbol('N-14') == 'N'
    assert element_symbol('n') == 'n'
    assert element_symbol('N') == 'N'


def test_isotope_symbol():
    assert isotope_symbol('He', 4) == 'He-4'
    assert isotope_symbol('helium-4') == 'He-4'
    assert isotope_symbol('H-2') == 'D'
    assert isotope_symbol('d') == 'D'
    assert isotope_symbol('Deuterium') == 'D'
    assert isotope_symbol('deuterium') == 'D'
    assert isotope_symbol('Hydrogen-3') == 'T'
    assert isotope_symbol('hydrogen-3') == 'T'
    assert isotope_symbol('h-3') == 'T'
    assert isotope_symbol('H-3') == 'T'
    assert isotope_symbol(1, 2) == 'D'
    assert isotope_symbol('h', 2) == 'D'
    assert isotope_symbol('Hydrogen', 3) == 'T'
    assert isotope_symbol('tritium') == 'T'
    assert isotope_symbol('H', 2) == 'D'
    assert isotope_symbol('Alpha') == 'He-4'
    assert isotope_symbol('alpha') == 'He-4'
    assert isotope_symbol(79, 197) == 'Au-197'
    assert isotope_symbol('p') == 'H-1'
    assert isotope_symbol('beryllium-8') == 'Be-8'
    assert isotope_symbol('neutron') == 'n'
    assert isotope_symbol('n') == 'n'
    assert isotope_symbol(0, 1) == 'n'
    assert isotope_symbol('n-1') == 'n'
    assert isotope_symbol('N-13') == 'N-13'
    assert isotope_symbol('p') == 'H-1'
    assert isotope_symbol('proton') == 'H-1'
    assert isotope_symbol('protium') == 'H-1'

    with pytest.raises(UserWarning):
        isotope_symbol('H-1', 1)

    with pytest.raises(UserWarning):
        isotope_symbol('H-2', 2)

    with pytest.raises(UserWarning):
        isotope_symbol('T', 3)

    with pytest.raises(UserWarning):
        isotope_symbol('Li-6', 6)

    with pytest.raises(UserWarning):
        isotope_symbol('lithium-6', 6)

    with pytest.raises(UserWarning):
        isotope_symbol('alpha', 4)

    with pytest.raises(UserWarning):
        isotope_symbol('p', 1)

    with pytest.raises(UserWarning):
        isotope_symbol('n', 1)

    with pytest.raises(ValueError):
        isotope_symbol('Md-260', 261)

    with pytest.raises(ValueError):
        isotope_symbol('protium', 2)

    with pytest.raises(ValueError):
        isotope_symbol('alpha', 3)

    with pytest.raises(ValueError):
        isotope_symbol('O-18', 19)

    with pytest.raises(ValueError):
        isotope_symbol('lead-209', 511)

    with pytest.raises(ValueError):
        isotope_symbol('He-1')

    with pytest.raises(ValueError):
        isotope_symbol(24, 23)

    with pytest.raises(ValueError):
        isotope_symbol('H', 0)

    with pytest.raises(ValueError):
        isotope_symbol('H-1', 2)

    with pytest.raises(ValueError):
        isotope_symbol('P')

    with pytest.raises(ValueError):
        isotope_symbol(1)

    with pytest.raises(ValueError):
        isotope_symbol(4)

    with pytest.raises(UserWarning):
        isotope_symbol('hydrogen-444444')


def test_atomic_number():
    assert atomic_number('H') == 1
    assert atomic_number('d') == 1
    assert atomic_number('D') == 1
    assert atomic_number('deuterium') == 1
    assert atomic_number('Deuterium') == 1
    assert atomic_number('t') == 1
    assert atomic_number('tritium') == 1
    assert atomic_number('p') == 1
    assert atomic_number('P') == 15
    assert atomic_number('Alpha') == 2
    assert atomic_number('C-12') == 6
    assert atomic_number('s-36') == 16
    assert atomic_number('Argon') == 18
    assert atomic_number('protium') == 1
    assert atomic_number('H-3') == 1
    assert atomic_number('p+') == 1
    assert atomic_number('Be-8') == 4
    assert atomic_number('n') == 0
    assert atomic_number('n-1') == 0
    assert atomic_number('Neutron') == 0
    assert atomic_number('N') == 7


def test_mass_number():
    assert mass_number('helium-3') == 3
    assert mass_number('Au-197') == 197
    assert mass_number('deuterium') == 2
    assert mass_number('D') == 2
    assert mass_number('H-2') == 2
    assert mass_number('tritium') == 3
    assert mass_number('T') == 3
    assert mass_number('alpha') == 4
    assert mass_number('p') == 1
    assert mass_number('n') == 1
    assert mass_number('Be-8') == 8
    assert mass_number('N-13') == 13


def test_element_name():
    assert element_name('D') == 'hydrogen'
    assert element_name('deuterium') == 'hydrogen'
    assert element_name('t') == 'hydrogen'
    assert element_name('Au') == 'gold'
    assert element_name('pb') == 'lead'
    assert element_name('alpha') == 'helium'
    assert element_name('helium-4') == 'helium'
    assert element_name('H-2') == 'hydrogen'
    assert element_name('d') == 'hydrogen'
    assert element_name('Deuterium') == 'hydrogen'
    assert element_name('Hydrogen-3') == 'hydrogen'
    assert element_name('hydrogen-3') == 'hydrogen'
    assert element_name('h-3') == 'hydrogen'
    assert element_name('H-3') == 'hydrogen'
    assert element_name('tritium') == 'hydrogen'
    assert element_name('Alpha') == 'helium'
    assert element_name('alpha') == 'helium'
    assert element_name(1) == 'hydrogen'
    assert element_name(26) == 'iron'
    assert element_name(79) == 'gold'
    assert element_name('p') == 'hydrogen'
    assert element_name('P') == 'phosphorus'
    assert element_name('Be-8') == 'beryllium'
    assert element_name('Li-7') == 'lithium'
    assert element_name('n') == 'neutron'
    assert element_name(0) == 'neutron'
    assert element_name('N') == 'nitrogen'


def test_standard_atomic_weight():
    assert standard_atomic_weight('H').value == 1.008
    assert standard_atomic_weight(1).value == 1.008
    assert standard_atomic_weight('Hydrogen').value == 1.008
    assert 30.973 < standard_atomic_weight('P').value < 30.974
    assert standard_atomic_weight('Au').unit == u.u
    assert standard_atomic_weight('N') is not None


def test_isotope_mass():
    assert isotope_mass('H-1') == isotope_mass('protium') == isotope_mass(1, 1)
    assert isotope_mass('D') == isotope_mass('H-2') == \
        isotope_mass('deuterium') == isotope_mass(1, 2)
    assert isotope_mass('T') == isotope_mass('H-3') == \
        isotope_mass('tritium') == isotope_mass(1, 3)
    assert np.isclose(isotope_mass('berkelium-249').value, 249.0749877)
    assert isotope_mass('Si-30').unit == u.u
    assert np.isclose(isotope_mass('n')/(1.008664*u.u), 1, atol=1e-6)


def test_ion_mass():
    assert ion_mass('H') > const.m_p
    assert ion_mass('hydrogen') > const.m_p
    assert ion_mass('proton') == const.m_p
    assert ion_mass('e+') == ion_mass('positron') == const.m_e
    assert np.isclose(ion_mass('alpha')/ion_mass('He-4', 2), 1.0)
    assert ion_mass('protium') == const.m_p
    assert ion_mass('Ne-22', 2) == 21.991385114*u.u - 2*const.m_e
    assert ion_mass('F-19', Z=3).unit == u.kg


def test_nuclear_binding_energy():
    assert nuclear_binding_energy('p') == 0
    assert nuclear_binding_energy('n') == 0
    assert nuclear_binding_energy('He-4') == nuclear_binding_energy('alpha') \
        == nuclear_binding_energy('He', 4)

    before = nuclear_binding_energy("D") + nuclear_binding_energy("T")
    after = nuclear_binding_energy("alpha")
    E_in_MeV = (after - before).to(u.MeV).value  # D + T --> alpha + n + E
    assert np.isclose(E_in_MeV, 17.58, rtol=0.01)


def test_energy_from_nuclear_reaction():
    reaction1 = 'D + T --> alpha + n'
    reaction2 = 'T + D -> n + alpha'
    released_energy1 = energy_from_nuclear_reaction(reaction1)
    released_energy2 = energy_from_nuclear_reaction(reaction2)
    assert np.isclose(released_energy1.to(u.MeV).value, 17.58, rtol=0.01)
    assert released_energy1 == released_energy2

    triple_alpha1 = 'alpha + He-4 --> Be-8'
    triple_alpha2 = 'Be-8 + alpha --> carbon-12'
    energy_triplealpha1 = energy_from_nuclear_reaction(triple_alpha1)
    energy_triplealpha2 = energy_from_nuclear_reaction(triple_alpha2)
    assert np.isclose(energy_triplealpha1.to(u.keV).value, -91.8, atol=0.1)
    assert np.isclose(energy_triplealpha2.to(u.MeV).value, 7.367, atol=0.1)

    alpha_decay_example = 'U-238 --> Th-234 + alpha'
    energy_alpha_decay = energy_from_nuclear_reaction(alpha_decay_example)
    assert np.isclose(energy_alpha_decay.to(u.MeV).value, 4.26975, atol=1e-5)


def test_is_isotope_stable():
    assert is_isotope_stable('H-1') is True
    assert is_isotope_stable(1, 1) is True
    assert is_isotope_stable('1', '1') is True
    assert is_isotope_stable('N-14') is True
    assert is_isotope_stable('N', 14) is True
    assert is_isotope_stable('P-31') is True
    assert is_isotope_stable('P', 31) is True
    assert is_isotope_stable('p') is True
    assert is_isotope_stable('alpha') is True
    assert is_isotope_stable('Xe-124') is True
    assert is_isotope_stable('Fe', 56) is True
    assert is_isotope_stable('Fe', '56') is True
    assert is_isotope_stable('Fe-56') is True
    assert is_isotope_stable('fe-56') is True
    assert is_isotope_stable('iron-56') is True
    assert is_isotope_stable('Iron-56') is True
    assert is_isotope_stable(26, 56) is True

    assert is_isotope_stable('Be-8') is False
    assert is_isotope_stable('n') is False
    assert is_isotope_stable('n-1') is False
    assert is_isotope_stable(0, 1) is False
    assert is_isotope_stable('U-235') is False
    assert is_isotope_stable('uranium-235') is False
    assert is_isotope_stable('T') is False
    assert is_isotope_stable(4, 8) is False
    assert is_isotope_stable('tritium') is False
    assert is_isotope_stable('neutron') is False
    assert is_isotope_stable('Pb-209') is False
    assert is_isotope_stable('pb-209') is False
    assert is_isotope_stable('lead-209') is False
    assert is_isotope_stable('Lead-209') is False
    assert is_isotope_stable('Pb', 209) is False
    assert is_isotope_stable(82, 209) is False
    assert is_isotope_stable('82', '209') is False

    with pytest.raises(ValueError):
        is_isotope_stable('hydrogen-444444')

    with pytest.raises(ValueError):
        is_isotope_stable('hydrogen', 0)

    with pytest.raises(ValueError):
        is_isotope_stable('')


def test_isotope_calls():
    assert known_isotopes('H') == \
        ['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7']
    assert common_isotopes('H') == ['H-1', 'D']
    assert stable_isotopes('He') == ['He-3', 'He-4']


def test_half_life():
    assert half_life('H-1') == np.inf*u.s
    assert np.isclose(half_life('tritium').to(u.s).value,
                      (12.32*u.yr).to(u.s).value, rtol=2e-4)
    assert half_life('H-1').unit == 's'
    assert half_life('tritium').unit == 's'

    for isotope in Isotopes.keys():  # unstable isotopes w/o half-life data
        if 'half_life' not in Isotopes[isotope].keys() and \
                not Isotopes[isotope].keys():
            with pytest.raises(UserWarning):
                assert half_life(isotope) is None


def test_atomic_TypeErrors():

    TypeErrorFunctions = [element_symbol, isotope_symbol, atomic_number,
                          is_isotope_stable, half_life, mass_number,
                          element_name, standard_atomic_weight, isotope_mass,
                          ion_mass, nuclear_binding_energy,
                          energy_from_nuclear_reaction]

    BadArguments = [1.1, True, False, {'cats': 'bats'}, 1+1j]

    for function in TypeErrorFunctions:
        for argument in BadArguments:
            with pytest.raises(TypeError):
                function(argument)

    with pytest.raises(TypeError):
        energy_from_nuclear_reaction(1)


def test_atomic_ValueErrors():

    ValueErrorFunctions = [element_symbol, isotope_symbol, atomic_number,
                           is_isotope_stable, half_life, mass_number,
                           element_name, standard_atomic_weight]

    BadArguments = [-1, 119, 'grumblemuffins', 'Oj']

    for function in ValueErrorFunctions:
        for argument in BadArguments:
            with pytest.raises(ValueError):
                function(argument)

    with pytest.raises(ValueError):
        energy_from_nuclear_reaction("I didn't like unstable eigenfunctions "
                                     "at first, but then they grew on me")
