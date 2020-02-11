from collections import namedtuple
import pytest
import numpy as np
from astropy import units as u, constants as const

from plasmapy.tests.helper import function_test_runner
from plasmapy.particles.symbols import atomic_symbol, isotope_symbol, element_name
from plasmapy.particles.isotopes import _Isotopes

from plasmapy.particles.atomic import (
    atomic_number,
    mass_number,
    standard_atomic_weight,
    particle_mass,
    reduced_mass,
    is_stable,
    half_life,
    known_isotopes,
    common_isotopes,
    stable_isotopes,
    isotopic_abundance,
    integer_charge,
    electric_charge,
    periodic_table_period,
    periodic_table_block,
    periodic_table_category,
    periodic_table_group,
    _is_electron,
)

from plasmapy.particles.exceptions import (
    AtomicError,
    MissingAtomicDataError,
    ChargeError,
    InvalidIsotopeError,
    InvalidElementError,
    InvalidParticleError,
    AtomicWarning,
)

Inputs = namedtuple("Inputs", ["function", "args", "kwargs", "expected"])

nokwargs = {}

test_inputs = [
    Inputs(atomic_symbol, 1, nokwargs, "H"),
    Inputs(atomic_symbol, 1, nokwargs, "H"),
    Inputs(atomic_symbol, "H", nokwargs, "H"),
    Inputs(atomic_symbol, "p", nokwargs, "H"),
    Inputs(atomic_symbol, "T", nokwargs, "H"),
    Inputs(atomic_symbol, "deuterium", nokwargs, "H"),
    Inputs(atomic_symbol, "deuteron", nokwargs, "H"),
    Inputs(atomic_symbol, "Tritium", nokwargs, "H"),
    Inputs(atomic_symbol, "triton", nokwargs, "H"),
    Inputs(atomic_symbol, "H-2", nokwargs, "H"),
    Inputs(atomic_symbol, "D", nokwargs, "H"),
    Inputs(atomic_symbol, "T", nokwargs, "H"),
    Inputs(atomic_symbol, "H-3", nokwargs, "H"),
    Inputs(atomic_symbol, "Hydrogen-3", nokwargs, "H"),
    Inputs(atomic_symbol, "helium", nokwargs, "He"),
    Inputs(atomic_symbol, 2, nokwargs, "He"),
    Inputs(atomic_symbol, "alpha", nokwargs, "He"),
    Inputs(atomic_symbol, "gold", nokwargs, "Au"),
    Inputs(atomic_symbol, "Gold", nokwargs, "Au"),
    Inputs(atomic_symbol, 79, nokwargs, "Au"),
    Inputs(atomic_symbol, "79", nokwargs, "Au"),
    Inputs(atomic_symbol, "P", nokwargs, "P"),
    Inputs(atomic_symbol, 118, nokwargs, "Og"),
    Inputs(atomic_symbol, "N-14", nokwargs, "N"),
    Inputs(atomic_symbol, "N", nokwargs, "N"),
    Inputs(atomic_symbol, "H +1", nokwargs, "H"),
    Inputs(atomic_symbol, "H 1+", nokwargs, "H"),
    Inputs(atomic_symbol, "hydrogen 1+", nokwargs, "H"),
    Inputs(atomic_symbol, "deuterium 1+", nokwargs, "H"),
    Inputs(atomic_symbol, "Fe 24+", nokwargs, "Fe"),
    Inputs(atomic_symbol, "Fe +24", nokwargs, "Fe"),
    Inputs(atomic_symbol, "Fe 2-", nokwargs, "Fe"),
    Inputs(atomic_symbol, "Fe -2", nokwargs, "Fe"),
    Inputs(atomic_symbol, "Fe+", nokwargs, "Fe"),
    Inputs(atomic_symbol, "Fe++", nokwargs, "Fe"),
    Inputs(atomic_symbol, "Fe-", nokwargs, "Fe"),
    Inputs(atomic_symbol, "Fe++++++++++++++", nokwargs, "Fe"),
    Inputs(atomic_symbol, "H-0", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, 3.14159, nokwargs, TypeError),
    Inputs(atomic_symbol, "Og-294b", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "H-93436107932635", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "Fe 2+4", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "Fe+24", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "Fe +59", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "C++++++++++++++++", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "C-++++", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "neutron", nokwargs, InvalidElementError),
    Inputs(atomic_symbol, "n", nokwargs, InvalidElementError),
    Inputs(atomic_symbol, "n-1", nokwargs, InvalidElementError),
    Inputs(atomic_symbol, "h", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "d", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "he", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "au", nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "p-", nokwargs, InvalidElementError),
    Inputs(atomic_symbol, 0, nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, 119, nokwargs, InvalidParticleError),
    Inputs(atomic_symbol, "antiproton", nokwargs, InvalidElementError),
    Inputs(isotope_symbol, ("He", 4), nokwargs, "He-4"),
    Inputs(isotope_symbol, ("helium-4",), nokwargs, "He-4"),
    Inputs(isotope_symbol, ("H-2",), nokwargs, "D"),
    Inputs(isotope_symbol, ("Deuterium",), nokwargs, "D"),
    Inputs(isotope_symbol, ("deuterium",), nokwargs, "D"),
    Inputs(isotope_symbol, ("deuteron",), nokwargs, "D"),
    Inputs(isotope_symbol, ("tritium",), nokwargs, "T"),
    Inputs(isotope_symbol, ("triton",), nokwargs, "T"),
    Inputs(isotope_symbol, ("Hydrogen-3",), nokwargs, "T"),
    Inputs(isotope_symbol, ("hydrogen-3",), nokwargs, "T"),
    Inputs(isotope_symbol, ("H-3",), nokwargs, "T"),
    Inputs(isotope_symbol, (1, 2), nokwargs, "D"),
    Inputs(isotope_symbol, ("Hydrogen", 3), nokwargs, "T"),
    Inputs(isotope_symbol, ("tritium",), nokwargs, "T"),
    Inputs(isotope_symbol, ("H", 2), nokwargs, "D"),
    Inputs(isotope_symbol, ("Alpha",), nokwargs, "He-4"),
    Inputs(isotope_symbol, ("alpha",), nokwargs, "He-4"),
    Inputs(isotope_symbol, (79, 197), nokwargs, "Au-197"),
    Inputs(isotope_symbol, ("p",), nokwargs, "H-1"),
    Inputs(isotope_symbol, ("beryllium-8",), nokwargs, "Be-8"),
    Inputs(isotope_symbol, ("N-13",), nokwargs, "N-13"),
    Inputs(isotope_symbol, ("p",), nokwargs, "H-1"),
    Inputs(isotope_symbol, ("proton",), nokwargs, "H-1"),
    Inputs(isotope_symbol, ("protium",), nokwargs, "H-1"),
    Inputs(isotope_symbol, ("N-13 2+",), nokwargs, "N-13"),
    Inputs(isotope_symbol, ("Hydrogen-3 +1",), nokwargs, "T"),
    Inputs(isotope_symbol, "Md-260", {"mass_numb": 261}, InvalidParticleError),
    Inputs(isotope_symbol, "protium", {"mass_numb": 2}, InvalidParticleError),
    Inputs(isotope_symbol, "alpha", {"mass_numb": 3}, InvalidParticleError),
    Inputs(isotope_symbol, "O-18", {"mass_numb": 19}, InvalidParticleError),
    Inputs(isotope_symbol, "lead-209", {"mass_numb": 511}, InvalidParticleError),
    Inputs(isotope_symbol, "He-1", nokwargs, InvalidParticleError),
    Inputs(isotope_symbol, 24, {"mass_numb": 23}, InvalidParticleError),
    Inputs(isotope_symbol, "H", {"mass_numb": 0}, InvalidParticleError),
    Inputs(isotope_symbol, "H-1", {"mass_numb": 2}, InvalidParticleError),
    Inputs(isotope_symbol, "P", nokwargs, InvalidIsotopeError),
    Inputs(isotope_symbol, 1, nokwargs, InvalidIsotopeError),
    Inputs(isotope_symbol, 4, nokwargs, InvalidIsotopeError),
    Inputs(isotope_symbol, "hydrogen-444444", nokwargs, InvalidParticleError),
    Inputs(isotope_symbol, "Fe", {"mass_numb": 2.1}, TypeError),
    Inputs(isotope_symbol, "He", {"mass_numb": "c"}, TypeError),
    Inputs(isotope_symbol, "He-3", {"mass_numb": 4}, InvalidParticleError),
    Inputs(isotope_symbol, "D", {"mass_numb": 3}, InvalidParticleError),
    Inputs(isotope_symbol, "T", {"mass_numb": 2}, InvalidParticleError),
    Inputs(isotope_symbol, "Fe", {"mass_numb": None}, InvalidIsotopeError),
    Inputs(isotope_symbol, "He", {"mass_numb": 99}, InvalidParticleError),
    Inputs(isotope_symbol, "d", nokwargs, InvalidParticleError),
    Inputs(isotope_symbol, "h-3", nokwargs, InvalidParticleError),
    Inputs(isotope_symbol, "h", nokwargs, InvalidParticleError),
    Inputs(isotope_symbol, "d+", nokwargs, InvalidParticleError),
    Inputs(isotope_symbol, "H-1", {"mass_numb": 1}, AtomicWarning),
    Inputs(isotope_symbol, "H-2", {"mass_numb": 2}, AtomicWarning),
    Inputs(isotope_symbol, "T", {"mass_numb": 3}, AtomicWarning),
    Inputs(isotope_symbol, "Li-6", {"mass_numb": 6}, AtomicWarning),
    Inputs(isotope_symbol, "lithium-6", {"mass_numb": 6}, AtomicWarning),
    Inputs(isotope_symbol, "alpha", {"mass_numb": 4}, AtomicWarning),
    Inputs(isotope_symbol, "p", {"mass_numb": 1}, AtomicWarning),
    Inputs(atomic_number, "H", nokwargs, 1),
    Inputs(atomic_number, "D", nokwargs, 1),
    Inputs(atomic_number, "deuterium", nokwargs, 1),
    Inputs(atomic_number, "Deuterium", nokwargs, 1),
    Inputs(atomic_number, "tritium", nokwargs, 1),
    Inputs(atomic_number, "p", nokwargs, 1),
    Inputs(atomic_number, "P", nokwargs, 15),
    Inputs(atomic_number, "Alpha", nokwargs, 2),
    Inputs(atomic_number, "C-12", nokwargs, 6),
    Inputs(atomic_number, "Argon", nokwargs, 18),
    Inputs(atomic_number, "protium", nokwargs, 1),
    Inputs(atomic_number, "H-3", nokwargs, 1),
    Inputs(atomic_number, "p+", nokwargs, 1),
    Inputs(atomic_number, "Be-8", nokwargs, 4),
    Inputs(atomic_number, "N", nokwargs, 7),
    Inputs(atomic_number, "N 2+", nokwargs, 7),
    Inputs(atomic_number, "N +1", nokwargs, 7),
    Inputs(atomic_number, "N+++", nokwargs, 7),
    Inputs(atomic_number, "H-3934", nokwargs, InvalidParticleError),
    Inputs(atomic_number, "C-12b", nokwargs, InvalidParticleError),
    Inputs(atomic_number, -1.5, nokwargs, TypeError),
    Inputs(atomic_number, "n", nokwargs, InvalidElementError),
    Inputs(atomic_number, "n-1", nokwargs, InvalidElementError),
    Inputs(atomic_number, "neutron", nokwargs, InvalidElementError),
    Inputs(atomic_number, "Neutron", nokwargs, InvalidElementError),
    Inputs(atomic_number, "d", nokwargs, InvalidParticleError),
    Inputs(atomic_number, "t", nokwargs, InvalidParticleError),
    Inputs(atomic_number, "s-36", nokwargs, InvalidParticleError),
    Inputs(mass_number, "helium-3", nokwargs, 3),
    Inputs(mass_number, "Au-197", nokwargs, 197),
    Inputs(mass_number, "deuterium", nokwargs, 2),
    Inputs(mass_number, "D", nokwargs, 2),
    Inputs(mass_number, "H-2", nokwargs, 2),
    Inputs(mass_number, "tritium", nokwargs, 3),
    Inputs(mass_number, "T", nokwargs, 3),
    Inputs(mass_number, "alpha", nokwargs, 4),
    Inputs(mass_number, "p", nokwargs, 1),
    Inputs(mass_number, "Be-8", nokwargs, 8),
    Inputs(mass_number, "N-13", nokwargs, 13),
    Inputs(mass_number, "N-13 2+", nokwargs, 13),
    Inputs(mass_number, "N-13 +2", nokwargs, 13),
    Inputs(mass_number, "N-13+++", nokwargs, 13),
    Inputs(mass_number, "H-359", nokwargs, InvalidParticleError),
    Inputs(mass_number, "C-12b", nokwargs, InvalidParticleError),
    Inputs(mass_number, -1.5, nokwargs, TypeError),
    Inputs(mass_number, "N-13+-+-", nokwargs, InvalidParticleError),
    Inputs(mass_number, "h-3", nokwargs, InvalidParticleError),
    Inputs(mass_number, "n", nokwargs, InvalidIsotopeError),
    Inputs(mass_number, "n-1", nokwargs, InvalidIsotopeError),
    Inputs(element_name, "D", nokwargs, "hydrogen"),
    Inputs(element_name, "deuterium", nokwargs, "hydrogen"),
    Inputs(element_name, "Au", nokwargs, "gold"),
    Inputs(element_name, "alpha", nokwargs, "helium"),
    Inputs(element_name, "helium-4", nokwargs, "helium"),
    Inputs(element_name, "H-2", nokwargs, "hydrogen"),
    Inputs(element_name, "Deuterium", nokwargs, "hydrogen"),
    Inputs(element_name, "Hydrogen-3", nokwargs, "hydrogen"),
    Inputs(element_name, "hydrogen-3", nokwargs, "hydrogen"),
    Inputs(element_name, "H-3", nokwargs, "hydrogen"),
    Inputs(element_name, "tritium", nokwargs, "hydrogen"),
    Inputs(element_name, "Alpha", nokwargs, "helium"),
    Inputs(element_name, "alpha", nokwargs, "helium"),
    Inputs(element_name, 1, nokwargs, "hydrogen"),
    Inputs(element_name, 26, nokwargs, "iron"),
    Inputs(element_name, 79, nokwargs, "gold"),
    Inputs(element_name, "p", nokwargs, "hydrogen"),
    Inputs(element_name, "P", nokwargs, "phosphorus"),
    Inputs(element_name, "Be-8", nokwargs, "beryllium"),
    Inputs(element_name, "Li-7", nokwargs, "lithium"),
    Inputs(element_name, "N", nokwargs, "nitrogen"),
    Inputs(element_name, "N+++", nokwargs, "nitrogen"),
    Inputs(element_name, "D-", nokwargs, "hydrogen"),
    Inputs(element_name, "vegancupcakes", nokwargs, InvalidParticleError),
    Inputs(element_name, "C-+-", nokwargs, InvalidParticleError),
    Inputs(element_name, 1.24, nokwargs, TypeError),
    Inputs(element_name, "n", nokwargs, InvalidElementError),
    Inputs(element_name, "neutron", nokwargs, InvalidElementError),
    Inputs(element_name, 0, nokwargs, InvalidParticleError),
    Inputs(element_name, "H++", nokwargs, InvalidParticleError),
    Inputs(element_name, "t", nokwargs, InvalidParticleError),
    Inputs(element_name, "pb", nokwargs, InvalidParticleError),
    Inputs(element_name, "d", nokwargs, InvalidParticleError),
    Inputs(element_name, "h-3", nokwargs, InvalidParticleError),
    Inputs(element_name, "Pb-9", nokwargs, InvalidParticleError),
    Inputs(element_name, "H 2+", nokwargs, InvalidParticleError),
    Inputs(standard_atomic_weight, "H", nokwargs, (1.008 * u.u).to(u.kg)),
    Inputs(standard_atomic_weight, 1, nokwargs, (1.008 * u.u).to(u.kg)),
    Inputs(standard_atomic_weight, "Hydrogen", nokwargs, (1.008 * u.u).to(u.kg)),
    Inputs(standard_atomic_weight, "Au", nokwargs, u.kg),
    Inputs(standard_atomic_weight, "H-1", nokwargs, AtomicError),
    Inputs(standard_atomic_weight, 1.1, nokwargs, TypeError),
    Inputs(standard_atomic_weight, "n", nokwargs, InvalidElementError),
    Inputs(standard_atomic_weight, "p", nokwargs, AtomicError),
    Inputs(standard_atomic_weight, "alpha", nokwargs, AtomicError),
    Inputs(standard_atomic_weight, "deuteron", nokwargs, AtomicError),
    Inputs(standard_atomic_weight, "tritium", nokwargs, AtomicError),
    Inputs(standard_atomic_weight, "Au+", nokwargs, AtomicError),
    Inputs(standard_atomic_weight, "Fe -2", nokwargs, AtomicError),
    Inputs(standard_atomic_weight, "Og 2+", nokwargs, AtomicError),
    Inputs(standard_atomic_weight, "h", nokwargs, InvalidParticleError),
    Inputs(standard_atomic_weight, "fe", nokwargs, InvalidParticleError),
    Inputs(particle_mass, "proton", nokwargs, const.m_p.to("kg")),
    Inputs(particle_mass, "H-1+", nokwargs, const.m_p.to("kg")),
    Inputs(particle_mass, "H-1 +1", nokwargs, const.m_p.to("kg")),
    Inputs(particle_mass, "H-1 1+", nokwargs, const.m_p.to("kg")),
    Inputs(particle_mass, "H-1", {"Z": 1}, const.m_p.to("kg")),
    Inputs(particle_mass, "hydrogen-1", {"Z": 1}, const.m_p.to("kg")),
    Inputs(particle_mass, "p+", nokwargs, const.m_p.to("kg")),
    Inputs(particle_mass, "F-19", {"Z": 3}, u.kg),
    Inputs(particle_mass, "Og 1+", nokwargs, MissingAtomicDataError),
    Inputs(particle_mass, "Fe-56", {"Z": 1.4}, TypeError),
    Inputs(particle_mass, "H-1 +1", {"Z": 0}, InvalidParticleError),
    Inputs(particle_mass, 26, {"Z": 1, "mass_numb": "a"}, TypeError),
    Inputs(particle_mass, 26, {"Z": 27, "mass_numb": 56}, InvalidParticleError),
    Inputs(particle_mass, "Og", {"Z": 1}, MissingAtomicDataError),
    Inputs(particle_mass, "Og", {"mass_numb": 696, "Z": 1}, InvalidParticleError),
    Inputs(particle_mass, "He 1+", {"mass_numb": 99}, InvalidParticleError),
    Inputs(particle_mass, "fe-56 1+", nokwargs, InvalidParticleError),
    Inputs(particle_mass, "H-1", {"mass_numb": 1, "Z": 1}, AtomicWarning),
    Inputs(particle_mass, "H", nokwargs, standard_atomic_weight("H")),
    Inputs(is_stable, "H-1", nokwargs, True),
    Inputs(is_stable, (1, 1), nokwargs, True),
    Inputs(is_stable, "N-14", nokwargs, True),
    Inputs(is_stable, ("N", 14), nokwargs, True),
    Inputs(is_stable, "P-31", nokwargs, True),
    Inputs(is_stable, ("P", 31), nokwargs, True),
    Inputs(is_stable, "p", nokwargs, True),
    Inputs(is_stable, "alpha", nokwargs, True),
    Inputs(is_stable, "Xe-124", nokwargs, True),
    Inputs(is_stable, "Fe", {"mass_numb": 56}, True),
    Inputs(is_stable, "Fe-56", nokwargs, True),
    Inputs(is_stable, "iron-56", nokwargs, True),
    Inputs(is_stable, "Iron-56", nokwargs, True),
    Inputs(is_stable, (26, 56), nokwargs, True),
    Inputs(is_stable, "Be-8", nokwargs, False),
    Inputs(is_stable, "U-235", nokwargs, False),
    Inputs(is_stable, "uranium-235", nokwargs, False),
    Inputs(is_stable, "T", nokwargs, False),
    Inputs(is_stable, (4, 8), nokwargs, False),
    Inputs(is_stable, "tritium", nokwargs, False),
    Inputs(is_stable, "Pb-209", nokwargs, False),
    Inputs(is_stable, "lead-209", nokwargs, False),
    Inputs(is_stable, "Lead-209", nokwargs, False),
    Inputs(is_stable, "Pb", {"mass_numb": 209}, False),
    Inputs(is_stable, (82, 209), nokwargs, False),
    Inputs(is_stable, ("hydrogen-444444",), nokwargs, InvalidParticleError),
    Inputs(is_stable, ("hydrogen", 0), nokwargs, InvalidParticleError),
    Inputs(is_stable, ("",), nokwargs, InvalidParticleError),
    Inputs(is_stable, ("pb-209",), nokwargs, InvalidParticleError),
    Inputs(is_stable, ("h",), nokwargs, InvalidParticleError),
    Inputs(is_stable, ("He",), nokwargs, InvalidIsotopeError),
    Inputs(is_stable, ("B",), nokwargs, InvalidIsotopeError),
    Inputs(integer_charge, "H+", nokwargs, 1),
    Inputs(integer_charge, "D +1", nokwargs, 1),
    Inputs(integer_charge, "tritium 1+", nokwargs, 1),
    Inputs(integer_charge, "H-", nokwargs, -1),
    Inputs(integer_charge, "Fe -2", nokwargs, -2),
    Inputs(integer_charge, "Fe 2-", nokwargs, -2),
    Inputs(integer_charge, "N--", nokwargs, -2),
    Inputs(integer_charge, "N++", nokwargs, 2),
    Inputs(integer_charge, "alpha", nokwargs, 2),
    Inputs(integer_charge, "proton", nokwargs, 1),
    Inputs(integer_charge, "deuteron", nokwargs, 1),
    Inputs(integer_charge, "triton", nokwargs, 1),
    Inputs(integer_charge, "electron", nokwargs, -1),
    Inputs(integer_charge, "e-", nokwargs, -1),
    Inputs(integer_charge, "e+", nokwargs, 1),
    Inputs(integer_charge, "positron", nokwargs, 1),
    Inputs(integer_charge, "n", nokwargs, 0),
    Inputs(integer_charge, "neutron", nokwargs, 0),
    Inputs(integer_charge, "p-", nokwargs, -1),
    Inputs(integer_charge, "antiproton", nokwargs, -1),
    Inputs(integer_charge, "fads", nokwargs, InvalidParticleError),
    Inputs(integer_charge, "H++", nokwargs, InvalidParticleError),
    Inputs(integer_charge, "h+", nokwargs, InvalidParticleError),
    Inputs(integer_charge, "fe 1+", nokwargs, InvalidParticleError),
    Inputs(integer_charge, "d+", nokwargs, InvalidParticleError),
    Inputs(integer_charge, "Fe 29+", nokwargs, InvalidParticleError),
    Inputs(integer_charge, "H-1", nokwargs, ChargeError),
    Inputs(integer_charge, "H---", nokwargs, AtomicWarning),
    Inputs(integer_charge, "Fe -26", nokwargs, AtomicWarning),
    Inputs(integer_charge, "Og 10-", nokwargs, AtomicWarning),
    Inputs(electric_charge, "p", nokwargs, u.C),
    Inputs(electric_charge, "p", nokwargs, const.e),
    Inputs(electric_charge, "e", nokwargs, -1.6021766208e-19 * u.C),
    Inputs(electric_charge, "alpha", nokwargs, 3.2043532416e-19 * u.C),
    Inputs(electric_charge, "n", nokwargs, 0 * u.C),
    Inputs(electric_charge, "bad input", nokwargs, InvalidParticleError),
    Inputs(electric_charge, "h+", nokwargs, InvalidParticleError),
    Inputs(electric_charge, "Au 81+", nokwargs, InvalidParticleError),
    Inputs(electric_charge, "Au 81-", nokwargs, AtomicWarning),
    Inputs(electric_charge, "H---", nokwargs, AtomicWarning),
    Inputs(half_life, "H-1", nokwargs, u.s),
    Inputs(half_life, "tritium", nokwargs, u.s),
    Inputs(half_life, "H-1", nokwargs, np.inf * u.s),
    Inputs(half_life, "No-248", nokwargs, MissingAtomicDataError),  # isotope with missing data
    Inputs(periodic_table_period, ("Ne", "Na"), nokwargs, TypeError),
    Inputs(periodic_table_block, ("N", "C", "F"), nokwargs, TypeError),
    Inputs(periodic_table_category, ("Rb", "He", "Li"), nokwargs, TypeError),
    Inputs(periodic_table_group, ("B", "Ti", "Ge"), nokwargs, TypeError),
    Inputs(isotopic_abundance, ("H", 1), nokwargs, isotopic_abundance("protium")),
    Inputs(isotopic_abundance, "D", nokwargs, 0.000115),
    Inputs(isotopic_abundance, "Be-8", nokwargs, 0.0),
    Inputs(isotopic_abundance, "Li-8", nokwargs, 0.0),
    Inputs(isotopic_abundance, ("Og", 294), nokwargs, AtomicWarning),
    Inputs(isotopic_abundance, "neutron", nokwargs, InvalidIsotopeError),
    Inputs(isotopic_abundance, "Og-2", nokwargs, InvalidParticleError),
    Inputs(particle_mass, "berkelium-249", nokwargs, (249.0749877 * u.u).to("kg")),
    Inputs(common_isotopes, "n", nokwargs, InvalidElementError),
    Inputs(stable_isotopes, "n", nokwargs, InvalidElementError),
    Inputs(known_isotopes, "n", nokwargs, InvalidElementError),
    Inputs(reduced_mass, ("N", 6e-26 * u.l), nokwargs, u.UnitConversionError),
    Inputs(reduced_mass, ("Og", "H"), nokwargs, MissingAtomicDataError),
    Inputs(known_isotopes, "H", nokwargs, ["H-1", "D", "T", "H-4", "H-5", "H-6", "H-7"]),
    Inputs(common_isotopes, "H", nokwargs, ["H-1", "D"]),
    Inputs(stable_isotopes, "He", nokwargs, ["He-3", "He-4"]),
    Inputs(_is_electron, "e-", nokwargs, True),
    Inputs(_is_electron, "e+", nokwargs, False),
    Inputs(_is_electron, "e", nokwargs, True),
    Inputs(_is_electron, "electron", nokwargs, True),
    Inputs(_is_electron, "ELECTRON", nokwargs, True),
    Inputs(_is_electron, "H", nokwargs, False),
    Inputs(_is_electron, "positron", nokwargs, False),
    Inputs(_is_electron, "Carbon", nokwargs, False),
    Inputs(half_life, "tritium", nokwargs, 3.888e8 * u.s),
]


def append_invalid_particle_tests(list_of_tests):
    """
    Add test cases for bad inputs that should each raise the same type
    of exception for most functions related to atomic data.
    """

    functions_to_test = [
        atomic_symbol,
        isotope_symbol,
        atomic_number,
        is_stable,
        half_life,
        mass_number,
        element_name,
        standard_atomic_weight,
        particle_mass,
        known_isotopes,
        stable_isotopes,
        common_isotopes,
        isotopic_abundance,
        integer_charge,
        electric_charge,
    ]

    invalid_particles_and_exceptions = [
        (-1, InvalidParticleError),
        (119, InvalidParticleError),
        ("...", InvalidParticleError),
        ("H-0", InvalidParticleError),
        ("Og-294b", InvalidParticleError),
        ("H-9343610", InvalidParticleError),
        ("Fe 2+4", InvalidParticleError),
        ("Fe+24", InvalidParticleError),
        ("Fe +59", InvalidParticleError),
        ("C++++++++++++++++", InvalidParticleError),
        ("C-++++", InvalidParticleError),
        ("h", InvalidParticleError),
        ("d", InvalidParticleError),
        ("he", InvalidParticleError),
        ("au", InvalidParticleError),
        ("alpha 1+", InvalidParticleError),
        ("alpha-4", InvalidParticleError),
        (1.1, TypeError),
        ({"...": "..."}, TypeError),
        (1 + 1j, TypeError),
    ]

    for function in functions_to_test:
        for invalid_particle, exception in invalid_particles_and_exceptions:
            list_of_tests.append(Inputs(function, invalid_particle, nokwargs, exception))


append_invalid_particle_tests(test_inputs)


@pytest.mark.parametrize("function, args, kwargs, expected", test_inputs)
def test_particle_functions(function, args, kwargs, expected):
    """
    Test that providing the principal functions in `~plasmapy.particles.atomic`
    result in the expected outcome.
    """

    function_test_runner(function=function, args=args, kwargs=kwargs, expected=expected)



def test_standard_atomic_weight_value_between():
    """Test that `standard_atomic_weight` returns approximately the
    correct value for phosphorus."""
    assert 30.973 < standard_atomic_weight('P').to(u.u).value < 30.974, \
        "Incorrect standard atomic weight for phosphorus."


def test_particle_mass_for_hydrogen_with_no_mass_number():
    """Test that `particle_mass` does not return the proton mass when no
    mass number is specified for hydrogen.  In this case, the
    standard atomic weight should be used to account for the small
    fraction of deuterium."""
    assert particle_mass('H', Z=1) > const.m_p
    assert particle_mass('hydrogen', Z=1) > const.m_p


def test_particle_mass_helium():
    """Test miscellaneous cases for `particle_mass`."""
    assert particle_mass('alpha') > particle_mass('He-3 2+')


# (arg1, kwargs1, arg2, kwargs2, expected)
equivalent_particle_mass_args = [
    ['e+', {}, 'positron', {}, const.m_e],
    ['alpha', {}, 'He-4++', {}, None],
    ['alpha', {}, 'helium-4 2+', {}, None],
    ['deuteron', {}, 'H', {'Z': 1, 'mass_numb': 2}, None],
    ['D+', {}, 'H-2+', {}, None],
    ['D+', {}, 'D 1+', {}, None],
    ['Deuterium+', {}, 'D', {'Z': 1}, None],
    ['triton', {}, 'H', {'Z': 1, 'mass_numb': 3}, None],
    ['T+', {}, 'H-3+', {}, None],
    ['T+', {}, 'T 1+', {}, None],
    ['Tritium+', {}, 'T', {'Z': 1}, None],
    ['Fe-56 1+', {}, 'Fe', {'mass_numb': 56, 'Z': 1},
     particle_mass('Fe-56 1-') - 2 * const.m_e],
    ['Fe-56 +1', {}, 26, {'mass_numb': 56, 'Z': 1}, None],
]


@pytest.mark.parametrize(
    "arg1, kwargs1, arg2, kwargs2, expected", equivalent_particle_mass_args)
def test_particle_mass_equivalent_args(arg1, kwargs1, arg2, kwargs2, expected):
    """Test that `particle_mass` returns equivalent results for
    equivalent positional and keyword arguments."""

    result1 = particle_mass(arg1, **kwargs1)
    result2 = particle_mass(arg2, **kwargs2)

    assert u.isclose(result1, result2), \
        (f"particle_mass({repr(arg1)}, **{kwargs1}) = {repr(result1)}, whereas "
         f"particle_mass({repr(arg2)}, **{kwargs2}) = {repr(result2)}.  "
         f"These results are not equivalent as expected.")

    if expected is not None:
        assert u.isclose(result1, result2) and u.isclose(result2, expected), \
            (f"particle_mass({repr(arg1)}, **{kwargs1}) = {repr(result1)} and "
             f"particle_mass({repr(arg2)}, **{kwargs2}) = {repr(result2)}, but "
             f"these results are not equal to {repr(expected)} as expected.")


def test_half_life_unstable_isotopes():
    """Test that `half_life` returns `None` and raises an exception for
    all isotopes that do not yet have half-life data."""
    for isotope in _Isotopes.keys():
        if 'half_life' not in _Isotopes[isotope].keys() and \
                not _Isotopes[isotope].keys():
            with pytest.raises(MissingAtomicDataError):
                half_life(isotope)


isotope_inputs = [
    (known_isotopes, "H", {}, "H-1", True),
    (known_isotopes, "H", {}, "D", True),
    (known_isotopes, "H", {}, "T", True),
    (known_isotopes, "Be", {}, "Be-8", True),
    (known_isotopes, 118, {}, "Og-294", True),
    (common_isotopes, "H", {}, "H-1", True),
    (stable_isotopes, "H", {}, "H-1", True),
    (stable_isotopes, "H", {}, "D", True),
    (common_isotopes, "Fe", {"most_common_only": True}, "Fe-56", True),
    (common_isotopes, "He", {"most_common_only": True}, "He-4", True),
    (stable_isotopes, "H", {}, "T", False),
    (common_isotopes, 1, {}, "H-4", False),
]


@pytest.mark.parametrize("function, element, kwargs, isotope, should_be_in_list", isotope_inputs)
def test_isotope_contents(function, element, kwargs, isotope, should_be_in_list):
    """
    Test that `known_isotopes`, `common_isotopes`, and `stable_isotopes`
    return certain isotopes that fall into these categories.
    """

    returned_isotopes = function(element, **kwargs)
    is_actually_in_list = isotope in returned_isotopes

    if is_actually_in_list is not should_be_in_list:

        errmsg = (
            f"{isotope} in {function.__name__}({element}, **{kwargs}) "
            f"should return {should_be_in_list} but instead returned "
            f"{is_actually_in_list}."
        )

        pytest.fail(errmsg)


def test_known_common_stable_isotopes_len():
    """Test that `known_isotopes`, `common_isotopes`, and
    `stable_isotopes` each return a `list` of the expected length.

    The number of common isotopes may change if isotopic composition
    data has any significant changes.

    The number of stable isotopes may decrease slightly if some isotopes
    are discovered to be unstable but with extremely long half-lives.

    The number of known isotopes will increase as new isotopes are
    discovered, so a buffer is included in the test.

    """

    assert len(common_isotopes()) == 288, \
        ("The length of the list returned by common_isotopes() is "
         f"{len(common_isotopes())}, which is not the expected value.")

    assert len(stable_isotopes()) == 254, \
        ("The length of the list returned by stable_isotopes() is "
         f"{len(stable_isotopes())}, which is not the expected value.")

    assert 3352 <= len(known_isotopes()) <= 3400, \
        ("The length of the list returned by known_isotopes() is "
         f"{len(known_isotopes())}, which is not within the expected range.")


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
