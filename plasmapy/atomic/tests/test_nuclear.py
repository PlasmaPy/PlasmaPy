from itertools import product
from astropy import units as u, constants as const
import numpy as np
from ..nuclear import (nuclear_binding_energy, nuclear_reaction_energy)

import pytest


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


def test_nuclear_reaction_energy():
    reaction1 = 'D + T --> alpha + n'
    reaction2 = 'T + D -> n + alpha'
    released_energy1 = nuclear_reaction_energy(reaction1)
    released_energy2 = nuclear_reaction_energy(reaction2)
    assert np.isclose((released_energy1.to(u.MeV)).value, 17.58, rtol=0.01)
    assert released_energy1 == released_energy2
    nuclear_reaction_energy('n + p+ --> n + p+ + p- + p+')
    nuclear_reaction_energy('neutron + antineutron --> neutron + antineutron')

    
def test_nuclear_reaction_energy_triple_alpha():
    triple_alpha1 = 'alpha + He-4 --> Be-8'
    triple_alpha2 = 'Be-8 + alpha --> carbon-12'
    energy_triplealpha1 = nuclear_reaction_energy(triple_alpha1)
    energy_triplealpha2 = nuclear_reaction_energy(triple_alpha2)
    assert np.isclose(energy_triplealpha1.to(u.keV).value, -91.8, atol=0.1)
    assert np.isclose(energy_triplealpha2.to(u.MeV).value, 7.367, atol=0.1)

    reactants = ['He-4', 'alpha']
    products = ['Be-8']
    energy = nuclear_reaction_energy(reactants=reactants, products=products)
    assert np.isclose(energy.to(u.keV).value, -91.8, atol=0.1)


def test_nuclear_reaction_energy_alpha_decay():
    alpha_decay_example = 'U-238 --> Th-234 + alpha'
    energy_alpha_decay = nuclear_reaction_energy(alpha_decay_example)
    assert np.isclose(energy_alpha_decay.to(u.MeV).value, 4.26975, atol=1e-5)


def test_nuclear_reaction_energy_triple_alpha_r():
    triple_alpha1_r = '4He-4 --> 2Be-8'
    energy_triplealpha1_r = nuclear_reaction_energy(triple_alpha1_r)
    assert np.isclose(energy_triplealpha1_r.to(u.keV).value,
                      -91.8 * 2, atol=0.1)


def test_nuclear_reaction_energy_beta():
    energy1 = nuclear_reaction_energy(reactants=['n'], products=['p', 'e-'])
    assert np.isclose(energy1.to(u.MeV).value, 0.78, atol=0.01)
    energy2 = nuclear_reaction_energy(
        reactants=['Mg-23'], products=['Na-23', 'e+'])
    assert np.isclose(energy2.to(u.MeV).value, 3.034591, atol=1e-5)


# (reaction, kwargs, expected_error)
nuclear_reaction_energy_error_table = [
    ('H + H --> H', {}, ValueError),
    (1, {}, TypeError),
    ('H-1 + H-1 --> H-1', {}, ValueError),
    ("invalid input", {}, ValueError),
    ('p --> n', {}, ValueError),
    ('p --> p', {'reactants': ['p'], 'products': ['p']}, ValueError),
]


@pytest.mark.parametrize(
    "reaction, kwargs, expected_error", nuclear_reaction_energy_error_table)
def test_nuclear_reaction_energy_error(reaction, kwargs, expected_error):
    with pytest.raises(expected_error):
        nuclear_reaction_energy(reaction, **kwargs)


# (reactants, products, expectedMeV, tol)
nuclear_reaction_energy_kwargs_table = [
    ('H-1', 'p', 0.0, 0.0),
    (['B-10', 'n'], ['Li-7', 'He-4'], 2.8, 0.06),
    (['Li-6', 'D'], ['2alpha'], 22.2, 0.06),
    (['C-12', 'p'], 'N-13', 1.95, 0.006),
    (['N-13'], ['C-13', 'e+'], 1.20, 0.006),
    (['C-13', 'hydrogen-1'], ['Nitrogen-14'], 7.54, 0.006),
    (['N-14', 'H-1'], ['O-15'], 7.35, 0.006),
    (['O-15'], ['N-15', 'e+'], 1.73, 0.006),
    (('N-15', 'H-1'), ('C-12', 'He-4'), 4.96, 0.006),
]


@pytest.mark.parametrize(
    "reactants, products, expectedMeV, tol",
    nuclear_reaction_energy_kwargs_table)
def test_nuclear_reaction_energy_kwargs(reactants, products, expectedMeV, tol):
    energy = nuclear_reaction_energy(reactants=reactants, products=products).si
    expected = (expectedMeV*u.MeV).si
    assert np.isclose(expected.value, energy.value, atol=tol)


# (reactants, products, expected_error)
nuclear_reaction_energy_kwerrors_table = [
    ('n', 3, TypeError),
    ('n', [3], ValueError),
    (['n'], ['p'], ValueError),
    (['n'], ['He-4'], ValueError),
    (['h'], ['H-1'], ValueError),
    (['e-', 'n'], 'p', ValueError),
    (['e+', 'n'], ['p-'], ValueError),
    (['kljsdf'], 'H-3', ValueError),
    (['H'], ['H-1'], ValueError),
    (['p'], ['n', 'n', 'e+'], ValueError),
]


@pytest.mark.parametrize("reactants, products, expected_error",
                         nuclear_reaction_energy_kwerrors_table)
def test_nuclear_reaction_energy_kwerrors(reactants, products, expected_error):
    with pytest.raises(expected_error):
        nuclear_reaction_energy(reactants=reactants, products=products)
