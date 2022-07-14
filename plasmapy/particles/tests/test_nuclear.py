import numpy as np
import pytest

from astropy import constants as const
from astropy import units as u

from plasmapy.particles.nuclear import (
    mass_energy,
    nuclear_binding_energy,
    nuclear_reaction_energy,
)
from plasmapy.utils.pytest_helpers import run_test, run_test_equivalent_calls

test_nuclear_equivalent_calls_table = [
    [nuclear_binding_energy, ["He-4", {}], ["alpha", {}], ["He", {"mass_numb": 4}]]
]


@pytest.mark.parametrize("test_inputs", test_nuclear_equivalent_calls_table)
def test_nuclear_equivalent_calls(test_inputs):
    run_test_equivalent_calls(test_inputs)


def test_nuclear_binding_energy_D_T():
    before = nuclear_binding_energy("D") + nuclear_binding_energy("T")
    after = nuclear_binding_energy("alpha")
    E_in_MeV = (after - before).to(u.MeV).value  # D + T --> alpha + n + E
    assert np.isclose(E_in_MeV, 17.58, rtol=0.01)


def test_nuclear_reaction_energy():
    reaction1 = "D + T --> alpha + n"
    reaction2 = "T + D -> n + alpha"
    released_energy1 = nuclear_reaction_energy(reaction1)
    released_energy2 = nuclear_reaction_energy(reaction2)
    assert np.isclose((released_energy1.to(u.MeV)).value, 17.58, rtol=0.01)
    assert released_energy1 == released_energy2
    assert nuclear_reaction_energy(
        "n + p+ --> n + p+ + p- + p+"
    ) == nuclear_reaction_energy("n + p+ --> n + 2*p+ + p-")
    nuclear_reaction_energy("neutron + antineutron --> neutron + antineutron")


def test_nuclear_reaction_energy_triple_alpha():
    triple_alpha1 = "alpha + He-4 --> Be-8"
    triple_alpha2 = "Be-8 + alpha --> carbon-12"
    energy_triplealpha1 = nuclear_reaction_energy(triple_alpha1)
    energy_triplealpha2 = nuclear_reaction_energy(triple_alpha2)
    assert np.isclose(energy_triplealpha1.to(u.keV).value, -91.8, atol=0.1)
    assert np.isclose(energy_triplealpha2.to(u.MeV).value, 7.367, atol=0.1)
    reactants = ["He-4", "alpha"]
    products = ["Be-8"]
    energy = nuclear_reaction_energy(reactants=reactants, products=products)
    assert np.isclose(energy.to(u.keV).value, -91.8, atol=0.1)


def test_nuclear_reaction_energy_alpha_decay():
    alpha_decay_example = "U-238 --> Th-234 + alpha"
    energy_alpha_decay = nuclear_reaction_energy(alpha_decay_example)
    assert np.isclose(energy_alpha_decay.to(u.MeV).value, 4.26975, atol=1e-5)


def test_nuclear_reaction_energy_triple_alpha_r():
    triple_alpha1_r = "4*He-4 --> 2*Be-8"
    energy_triplealpha1_r = nuclear_reaction_energy(triple_alpha1_r)
    assert np.isclose(energy_triplealpha1_r.to(u.keV).value, -91.8 * 2, atol=0.1)


def test_nuclear_reaction_energy_beta():
    energy1 = nuclear_reaction_energy(reactants=["n"], products=["p", "e-"])
    assert np.isclose(energy1.to(u.MeV).value, 0.78, atol=0.01)
    energy2 = nuclear_reaction_energy(reactants=["Mg-23"], products=["Na-23", "e+"])
    assert np.isclose(energy2.to(u.MeV).value, 3.034591, atol=1e-5)


# (reactants, products, expectedMeV, tol)
nuclear_reaction_energy_kwargs_table = [
    ("H-1", "p", 0.0, 0.0),
    (["B-10", "n"], ["Li-7", "He-4"], 2.8, 0.06),
    (["Li-6", "D"], ["2*alpha"], 22.2, 0.06),
    (["C-12", "p"], "N-13", 1.95, 0.006),
    (["N-13"], ["C-13", "e+"], 1.20, 0.006),
    (["C-13", "hydrogen-1"], ["Nitrogen-14"], 7.54, 0.006),
    (["N-14", "H-1"], ["O-15"], 7.35, 0.006),
    (["O-15"], ["N-15", "e+"], 1.73, 0.006),
    (("N-15", "H-1"), ("C-12", "He-4"), 4.96, 0.006),
]


@pytest.mark.parametrize(
    "reactants, products, expectedMeV, tol", nuclear_reaction_energy_kwargs_table
)
def test_nuclear_reaction_energy_kwargs(reactants, products, expectedMeV, tol):
    energy = nuclear_reaction_energy(reactants=reactants, products=products).si
    expected = (expectedMeV * u.MeV).si
    assert np.isclose(expected.value, energy.value, atol=tol)


def test_nuclear_reaction_energy():
    reactants = ["D", "T"]
    products = ["alpha", "n"]
    expected = 17.6 * u.MeV
    actual = nuclear_reaction_energy(reactants=reactants, products=products)
    assert u.isclose(actual, expected, rtol=1e-3)


table_of_nuclear_tests = [
    [nuclear_binding_energy, ["p"], {}, 0 * u.J],
    [nuclear_binding_energy, ["n"], {}, 0 * u.J],
    [nuclear_binding_energy, ["p"], {}, 0 * u.J],
    [mass_energy, ["e-"], {}, (const.m_e * const.c**2).to(u.J)],
    [mass_energy, ["p+"], {}, (const.m_p * const.c**2).to(u.J)],
    [mass_energy, ["H-1"], {}, (const.m_p * const.c**2).to(u.J)],
    [mass_energy, ["H-1 0+"], {}, (const.m_p * const.c**2).to(u.J)],
    [mass_energy, ["n"], {}, (const.m_n * const.c**2).to(u.J)],
]


@pytest.mark.parametrize(
    ["tested_object", "args", "kwargs", "expected_value"], table_of_nuclear_tests
)
def test_nuclear_table(tested_object, args, kwargs, expected_value):
    run_test(tested_object, args, kwargs, expected_value, rtol=1e-3)
