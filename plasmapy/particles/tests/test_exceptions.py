import numpy as np
from astropy import units as u
from plasmapy.particles.exceptions import InvalidIsotopeError, ParticleError, InvalidParticleError
from plasmapy.particles import IonizationState, IonizationStateCollection, atomic_symbol, nuclear_binding_energy, nuclear_reaction_energy
import pytest 
# IonizationState
from plasmapy.particles.exceptions import (
    ChargeError,
    InvalidElementError,
    InvalidIsotopeError,
    InvalidParticleError,
    MissingParticleDataError,
    ParticleError,
    ParticleWarning,
)
from ..atomic import (
    _is_electron,
    atomic_number,
    common_isotopes,
    electric_charge,
    half_life,
    integer_charge,
    is_stable,
    isotopic_abundance,
    known_isotopes,
    mass_number,
    particle_mass,
    periodic_table_block,
    periodic_table_category,
    periodic_table_group,
    periodic_table_period,
    reduced_mass,
    stable_isotopes,
    standard_atomic_weight,
)
from ..isotopes import _isotopes
from ..nuclear import nuclear_binding_energy, nuclear_reaction_energy
from ..symbols import atomic_symbol, element_name, isotope_symbol

tests_for_exceptions_IonizationState = [
    ("too few nstates", IonizationState, [], {"particle": "H", "ionic_fractions": [1.0]}, ParticleError),
    ( "too many nstates",
        IonizationState,
        [],
        {"particle": "H", "ionic_fractions": [1, 0, 0, 0]},
        ParticleError,
    ),
    (
        "ionic fraction < 0",
        IonizationState,
        [],
        {"particle": "He", "ionic_fractions": [-0.1, 0.1, 1]},
        ParticleError,
    ),
    (
        "ionic fraction > 1",
        IonizationState,
        [],
        {"particle": "He", "ionic_fractions": [1.1, 0.0, 0.0]},
        ParticleError,
    ),
    (
        "invalid ionic fraction",
        IonizationState,
        [],
        {"particle": "He", "ionic_fractions": [1.0, 0.0, "a"]},
        ParticleError,
    ),
    (
        "bad n_elem units",
        IonizationState,
        [],
        {"particle": "H", "ionic_fractions": [0, 1], "n_elem": 3 * u.m ** 3},
        u.UnitTypeError,
    ),
    (
        "bad T_e units",
        IonizationState,
        [],
        {"particle": "H", "ionic_fractions": [0, 1], "T_e": 1 * u.m},
        u.UnitTypeError,
    ),
    (
        "negative n_elem",
        IonizationState,
        [],
        {
            "particle": "He",
            "ionic_fractions": [1.0, 0.0, 0.0],
            "n_elem": -1 * u.m ** -3,
        },
        ParticleError,
    ),
    (
        "negative T_e",
        IonizationState,
        [],
        {"particle": "He", "ionic_fractions": [1.0, 0.0, 0.0], "T_e": -1 * u.K},
        ParticleError,
    ),
    (
        "redundant ndens",
        IonizationState,
        [],
        {
            "particle": "H",
            "ionic_fractions": np.array([3, 4]) * u.m ** -3,
            "n_elem": 4 * u.m ** -3,
        },
        ParticleError,
    ),
]


# IonizationStateCollection

tests_for_exceptions_IonizationStateCollection = [
    ("wrong type",
        IonizationStateCollection,
        [],
     {"inputs": None}, ParticleError),
    ("not normalized", 
        IonizationStateCollection,
        [],
        {"inputs": {"He": [0.4, 0.5, 0.0]}, "tol": 1e-9}, ParticleError
    ),
    ("negative ionfrac",
        IonizationStateCollection,
        [],
     {"inputs": {"H": [-0.1, 1.1]}}, ParticleError),
    ("ion",
        IonizationStateCollection,
        [],
     {"inputs": {"H": [0.1, 0.9], "He+": [0.0, 0.9, 0.1]}}, ParticleError),
    ("repeat elements", 
        IonizationStateCollection,
        [],
        {"inputs": {"H": [0.1, 0.9], "hydrogen": [0.2, 0.8]}}, ParticleError
    ),
    ("isotope of element", 
        IonizationStateCollection,
        [],
        {"inputs": {"H": [0.1, 0.9], "D": [0.2, 0.8]}}, ParticleError
    ),
    ("negative abundance", 
        IonizationStateCollection,
        [],
        {
            "inputs": {"H": [0.1, 0.9], "He": [0.4, 0.5, 0.1]},
            "abundances": {"H": 1, "He": -0.1},
        },
        ParticleError,
    ),
    ("imaginary abundance", 
        IonizationStateCollection,
        [],
        {
            "inputs": {"H": [0.1, 0.9], "He": [0.4, 0.5, 0.1]},
            "abundances": {"H": 1, "He": 0.1j},
        },
        ParticleError,
    ),
    ("wrong density units", 
        IonizationStateCollection,
        [],
        {
            "inputs": {"H": [10, 90] * u.m ** -3, "He": [0.1, 0.9, 0] * u.m ** -2},
            "abundances": {"H": 1, "He": 0.1},
        },
        ParticleError,
    ),
    ("abundance redundance", 
        IonizationStateCollection,
        [],
        {
            "inputs": {"H": [10, 90] * u.m ** -3, "He": [0.1, 0.9, 0] * u.m ** -3},
            "abundances": {"H": 1, "He": 0.1},
        },
        ParticleError,
    ),
    ("abundance contradiction", 
        IonizationStateCollection,
        [],
        {
            "inputs": {"H": [10, 90] * u.m ** -3, "He": [0.1, 0.9, 0] * u.m ** -3},
            "abundances": {"H": 1, "He": 0.11},
        },
        ParticleError,
    ),
    ("kappa too small",
        IonizationStateCollection,
        [],
     {"inputs": ["H"], "kappa": 1.499999}, ParticleError),
    ("negative n",
        IonizationStateCollection,
        [],
     {"inputs": ["H"], "n0": -1 * u.cm ** -3}, ParticleError),
    ("negative T_e",
        IonizationStateCollection,
        [],
     {"inputs": ["H-1"], "T_e": -1 * u.K}, ParticleError),
]


tests_from_nuclear = [
    [nuclear_reaction_energy,
     [],
     {"reactants": ["n"], "products": 3}, TypeError],
    [
        nuclear_reaction_energy,
     [],
        {"reactants": ["n"], "products": ["He-4"]},
        ParticleError,
    ],
    [
        nuclear_reaction_energy,
     [],
        {"reactants": ["h"], "products": ["H-1"]},
        ParticleError,
    ],
    [
        nuclear_reaction_energy,
     [],
        {"reactants": ["e-", "n"], "products": ["p+"]},
        ParticleError,
    ],
    [
        nuclear_reaction_energy,
     [],
        {"reactants": ["e+", "n"], "products": ["p-"]},
        ParticleError,
    ],
    [
        nuclear_reaction_energy,
     [],
        {"reactants": ["ksdf"], "products": ["H-3"]},
        ParticleError,
    ],
    [
        nuclear_reaction_energy,
     [],
        {"reactants": ["H"], "products": ["H-1"]},
        ParticleError,
    ],
    [
        nuclear_reaction_energy,
     [],
        {"reactants": ["p"], "products": ["n", "n", "e-"]},
        ParticleError,
    ],
    [
        nuclear_reaction_energy,
        ["p --> p"],
        {"reactants": "p", "products": "p"},
        ParticleError,
    ],
    [nuclear_binding_energy, ["H"], {}, ParticleError],
    [nuclear_binding_energy, ["He-99"], {}, InvalidParticleError],
    [nuclear_binding_energy, ["He"], {"mass_numb": 99}, InvalidParticleError],
    [nuclear_binding_energy, [3.1415926535j], {}, TypeError],
]
tests_for_exceptions = tests_for_exceptions_IonizationState + tests_for_exceptions_IonizationStateCollection


@pytest.mark.parametrize(["test_name", "tested_object", "args", "kwargs", "expected_exception"], tests_for_exceptions)
def test_named_tests_for_exceptions(test_name, tested_object, args, kwargs, expected_exception):
    """
    Test that appropriate exceptions are raised for inappropriate inputs
    to `IonizationState`.
    """
    print(test_name) # TODO find better ways for this
    with pytest.raises(expected_exception):
        tested_object(*args, **kwargs)

@pytest.mark.parametrize(["tested_object", "args", "kwargs", "expected_exception"], tests_from_nuclear)
def test_unnamed_tests_exceptions(tested_object, args, kwargs, expected_exception):
    """
    Test that appropriate exceptions are raised for inappropriate inputs
    to `IonizationState`.
    """
    with pytest.raises(expected_exception):
        tested_object(*args, **kwargs)


single_arg_only_test = [
    [atomic_symbol, "H-0", InvalidParticleError],
    [atomic_symbol, 3.14159, TypeError],
    [atomic_symbol, "Og-294b", InvalidParticleError],
    [atomic_symbol, "H-934361079326356530741942970523610389", InvalidParticleError],
    [atomic_symbol, "Fe 2+4", InvalidParticleError],
    [atomic_symbol, "Fe+24", InvalidParticleError],
    [atomic_symbol, "Fe +59", InvalidParticleError],
    [atomic_symbol, "C++++++++++++++++", InvalidParticleError],
    [atomic_symbol, "C-++++", InvalidParticleError],
    [atomic_symbol, "neutron", InvalidElementError],
    [atomic_symbol, "n", InvalidElementError],
    [atomic_symbol, "n-1", InvalidElementError],
    [atomic_symbol, "h", InvalidParticleError],
    [atomic_symbol, "d", InvalidParticleError],
    [atomic_symbol, "he", InvalidParticleError],
    [atomic_symbol, "au", InvalidParticleError],
    [atomic_symbol, "p-", InvalidElementError],
    [atomic_symbol, 0, InvalidParticleError],
    [atomic_symbol, 119, InvalidParticleError],
    [atomic_symbol, "antiproton", InvalidElementError],
    [atomic_number, "H-3934", InvalidParticleError],
    [atomic_number, "C-12b", InvalidParticleError],
    [atomic_number, -1.5, TypeError],
    [atomic_number, "n", InvalidElementError],
    [atomic_number, "n-1", InvalidElementError],
    [atomic_number, "neutron", InvalidElementError],
    [atomic_number, "Neutron", InvalidElementError],
    [atomic_number, "d", InvalidParticleError],
    [atomic_number, "t", InvalidParticleError],
    [atomic_number, "s-36", InvalidParticleError],
    [mass_number, "H-359", InvalidParticleError],
    [mass_number, "C-12b", InvalidParticleError],
    [mass_number, -1.5, TypeError],
    [mass_number, "N-13+-+-", InvalidParticleError],
    [mass_number, "h-3", InvalidParticleError],
    [mass_number, "n", InvalidIsotopeError],
    [mass_number, "n-1", InvalidIsotopeError],
    [element_name, "vegancupcakes", InvalidParticleError],
    [element_name, "C-+-", InvalidParticleError],
    [element_name, 1.24, TypeError],
    [element_name, "n", InvalidElementError],
    [element_name, "neutron", InvalidElementError],
    [element_name, 0, InvalidParticleError],
    [element_name, "H++", InvalidParticleError],
    [element_name, "t", InvalidParticleError],
    [element_name, "pb", InvalidParticleError],
    [element_name, "d", InvalidParticleError],
    [element_name, "h-3", InvalidParticleError],
    [element_name, "Pb-9", InvalidParticleError],
    [element_name, "H 2+", InvalidParticleError],
    [standard_atomic_weight, "H-1", ParticleError],
    [standard_atomic_weight, "help i'm trapped in a unit test", InvalidParticleError],
    [standard_atomic_weight, 1.1, TypeError],
    [standard_atomic_weight, "n", InvalidElementError],
    [standard_atomic_weight, "p", ParticleError],
    [standard_atomic_weight, "alpha", ParticleError],
    [standard_atomic_weight, "deuteron", ParticleError],
    [standard_atomic_weight, "tritium", ParticleError],
    [standard_atomic_weight, "Au+", ParticleError],
    [standard_atomic_weight, "Fe -2", ParticleError],
    [standard_atomic_weight, "Og 2+", ParticleError],
    [standard_atomic_weight, "h", InvalidParticleError],
    [standard_atomic_weight, "fe", InvalidParticleError],
    [electric_charge, "badinput", InvalidParticleError],
    [electric_charge, "h+", InvalidParticleError],
    [electric_charge, "Au 81+", InvalidParticleError],
    [electric_charge, "Au 81-", ParticleWarning],
    [electric_charge, "H---", ParticleWarning],
    [integer_charge, "fads", InvalidParticleError],
    [integer_charge, "H++", InvalidParticleError],
    [integer_charge, "h+", InvalidParticleError],
    [integer_charge, "fe 1+", InvalidParticleError],
    [integer_charge, "d+", InvalidParticleError],
    [integer_charge, "Fe 29+", InvalidParticleError],
    [integer_charge, "H-1", ChargeError],
    [integer_charge, "H---", ParticleWarning],
    [integer_charge, "Fe -26", ParticleWarning],
    [integer_charge, "Og 10-", ParticleWarning],
]
@pytest.mark.parametrize(["tested_object", "arg", "expected_exception"], single_arg_only_test)
def test_single_arg_exceptions(tested_object, arg, expected_exception):
    """
    Test that appropriate exceptions are raised for inappropriate inputs
    to `IonizationState`.
    """
    if issubclass(expected_exception, Warning):
        with pytest.warns(expected_exception):
            tested_object(arg)
    elif issubclass(expected_exception, Exception):
        with pytest.raises(expected_exception):
            tested_object(arg)

# from test_atomic.py




# The tables above do not include the function to be tested in order to
# avoid cluttering up the code.  The following block of code prepends
# the correct function to each list containing args, kwargs, and the
# expected outcome prior to being passed through to run_test.


# tables_and_functions = [
#     (atomic_symbol, atomic_symbol_table),
#     (isotope_symbol, isotope_symbol_table),
#     (atomic_number, atomic_number_table),
#     (mass_number, mass_number_table),
#     (element_name, element_name_table),
#     (standard_atomic_weight, standard_atomic_weight_table),
#     (is_stable, is_stable_table),
#     (particle_mass, particle_mass_table),
#     (integer_charge, integer_charge_table),
#     (electric_charge, electric_charge_table),
#     (half_life, half_life_table),
# ]

# all_tests = []

# for func, table in tables_and_functions:
#     for inputs in table:
#         inputs.insert(0, func)
#         if len(inputs) == 3:
#             inputs.insert(2, {})
#     all_tests += table

# # Set up tests for a variety of atomic functions to make sure that bad
# # inputs lead to the expected errors.

# atomic_TypeError_funcs_table = [
#     atomic_symbol,
#     isotope_symbol,
#     atomic_number,
#     is_stable,
#     half_life,
#     mass_number,
#     element_name,
#     standard_atomic_weight,
#     nuclear_binding_energy,
#     nuclear_reaction_energy,
# ]

# atomic_TypeError_badargs = [1.1, {"cats": "bats"}, 1 + 1j]

# atomic_ParticleErrors_funcs_table = [
#     atomic_symbol,
#     isotope_symbol,
#     atomic_number,
#     is_stable,
#     half_life,
#     mass_number,
#     element_name,
#     standard_atomic_weight,
#     particle_mass,
#     known_isotopes,
#     stable_isotopes,
#     common_isotopes,
#     isotopic_abundance,
#     integer_charge,
#     electric_charge,
# ]

# atomic_ParticleError_badargs = [
#     -1,
#     119,
#     "grumblemuffins",
#     "H-0",
#     "Og-294b",
#     "H-9343610",
#     "Fe 2+4",
#     "Fe+24",
#     "Fe +59",
#     "C++++++++++++++++",
#     "C-++++",
#     "h",
#     "d",
#     "he",
#     "au",
#     "alpha 1+",
#     "alpha-4",
# ]

# metatable = [
#     (atomic_TypeError_funcs_table, atomic_TypeError_badargs, TypeError),
#     (
#         atomic_ParticleErrors_funcs_table,
#         atomic_ParticleError_badargs,
#         InvalidParticleError,
#     ),
# ]

# for funcs, badargs, error in metatable:
#     for func in funcs:
#         for badarg in badargs:
#             all_tests += [[func, badarg, error]]


# @pytest.mark.parametrize("inputs", all_tests)
# def test_atomic_functions(inputs):
#     print(inputs)
#     run_test(inputs)




