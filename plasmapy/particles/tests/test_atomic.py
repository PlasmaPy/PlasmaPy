import numpy as np
import pytest
from astropy import constants as const
from astropy import units as u

from plasmapy.particles.atomic import (
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
from plasmapy.particles.exceptions import (
    AtomicError,
    AtomicWarning,
    ChargeError,
    InvalidElementError,
    InvalidIsotopeError,
    InvalidParticleError,
    MissingAtomicDataError,
)
from plasmapy.particles.isotopes import _Isotopes
from plasmapy.particles.symbols import atomic_symbol, element_name, isotope_symbol
from plasmapy.tests.helpers import FunctionTestCase, test_runner

test_cases = [
    FunctionTestCase(function=atomic_symbol, args=1, expected="H"),
    FunctionTestCase(function=atomic_symbol, args=1, expected="H"),
    FunctionTestCase(function=atomic_symbol, args="H", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="p", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="T", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="deuterium", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="deuteron", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="Tritium", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="triton", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="H-2", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="D", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="T", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="H-3", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="Hydrogen-3", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="helium", expected="He"),
    FunctionTestCase(function=atomic_symbol, args=2, expected="He"),
    FunctionTestCase(function=atomic_symbol, args="alpha", expected="He"),
    FunctionTestCase(function=atomic_symbol, args="gold", expected="Au"),
    FunctionTestCase(function=atomic_symbol, args="Gold", expected="Au"),
    FunctionTestCase(function=atomic_symbol, args=79, expected="Au"),
    FunctionTestCase(function=atomic_symbol, args="79", expected="Au"),
    FunctionTestCase(function=atomic_symbol, args="P", expected="P"),
    FunctionTestCase(function=atomic_symbol, args=118, expected="Og"),
    FunctionTestCase(function=atomic_symbol, args="N-14", expected="N"),
    FunctionTestCase(function=atomic_symbol, args="N", expected="N"),
    FunctionTestCase(function=atomic_symbol, args="H +1", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="H 1+", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="hydrogen 1+", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="deuterium 1+", expected="H"),
    FunctionTestCase(function=atomic_symbol, args="Fe 24+", expected="Fe"),
    FunctionTestCase(function=atomic_symbol, args="Fe +24", expected="Fe"),
    FunctionTestCase(function=atomic_symbol, args="Fe 2-", expected="Fe"),
    FunctionTestCase(function=atomic_symbol, args="Fe -2", expected="Fe"),
    FunctionTestCase(function=atomic_symbol, args="Fe+", expected="Fe"),
    FunctionTestCase(function=atomic_symbol, args="Fe++", expected="Fe"),
    FunctionTestCase(function=atomic_symbol, args="Fe-", expected="Fe"),
    FunctionTestCase(function=atomic_symbol, args="Fe++++++++++++++", expected="Fe"),
    FunctionTestCase(function=atomic_symbol, args="H-0", expected=InvalidParticleError),
    FunctionTestCase(function=atomic_symbol, args=3.14159, expected=TypeError),
    FunctionTestCase(
        function=atomic_symbol, args="Og-294b", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=atomic_symbol, args="H-93436107932635", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=atomic_symbol, args="Fe 2+4", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=atomic_symbol, args="Fe+24", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=atomic_symbol, args="Fe +59", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=atomic_symbol, args="C++++++++++++++++", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=atomic_symbol, args="C-++++", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=atomic_symbol, args="neutron", expected=InvalidElementError
    ),
    FunctionTestCase(function=atomic_symbol, args="n", expected=InvalidElementError),
    FunctionTestCase(function=atomic_symbol, args="n-1", expected=InvalidElementError),
    FunctionTestCase(function=atomic_symbol, args="h", expected=InvalidParticleError),
    FunctionTestCase(function=atomic_symbol, args="d", expected=InvalidParticleError),
    FunctionTestCase(function=atomic_symbol, args="he", expected=InvalidParticleError),
    FunctionTestCase(function=atomic_symbol, args="au", expected=InvalidParticleError),
    FunctionTestCase(function=atomic_symbol, args="p-", expected=InvalidElementError),
    FunctionTestCase(function=atomic_symbol, args=0, expected=InvalidParticleError),
    FunctionTestCase(function=atomic_symbol, args=119, expected=InvalidParticleError),
    FunctionTestCase(
        function=atomic_symbol, args="antiproton", expected=InvalidElementError
    ),
    FunctionTestCase(function=isotope_symbol, args=("He", 4), expected="He-4"),
    FunctionTestCase(function=isotope_symbol, args=("helium-4",), expected="He-4"),
    FunctionTestCase(function=isotope_symbol, args=("H-2",), expected="D"),
    FunctionTestCase(function=isotope_symbol, args=("Deuterium",), expected="D"),
    FunctionTestCase(function=isotope_symbol, args=("deuterium",), expected="D"),
    FunctionTestCase(function=isotope_symbol, args=("deuteron",), expected="D"),
    FunctionTestCase(function=isotope_symbol, args=("tritium",), expected="T"),
    FunctionTestCase(function=isotope_symbol, args=("triton",), expected="T"),
    FunctionTestCase(function=isotope_symbol, args=("Hydrogen-3",), expected="T"),
    FunctionTestCase(function=isotope_symbol, args=("hydrogen-3",), expected="T"),
    FunctionTestCase(function=isotope_symbol, args=("H-3",), expected="T"),
    FunctionTestCase(function=isotope_symbol, args=(1, 2), expected="D"),
    FunctionTestCase(function=isotope_symbol, args=("Hydrogen", 3), expected="T"),
    FunctionTestCase(function=isotope_symbol, args=("tritium",), expected="T"),
    FunctionTestCase(function=isotope_symbol, args=("H", 2), expected="D"),
    FunctionTestCase(function=isotope_symbol, args=("Alpha",), expected="He-4"),
    FunctionTestCase(function=isotope_symbol, args=("alpha",), expected="He-4"),
    FunctionTestCase(function=isotope_symbol, args=(79, 197), expected="Au-197"),
    FunctionTestCase(function=isotope_symbol, args=("p",), expected="H-1"),
    FunctionTestCase(function=isotope_symbol, args=("beryllium-8",), expected="Be-8"),
    FunctionTestCase(function=isotope_symbol, args=("N-13",), expected="N-13"),
    FunctionTestCase(function=isotope_symbol, args=("p",), expected="H-1"),
    FunctionTestCase(function=isotope_symbol, args=("proton",), expected="H-1"),
    FunctionTestCase(function=isotope_symbol, args=("protium",), expected="H-1"),
    FunctionTestCase(function=isotope_symbol, args=("N-13 2+",), expected="N-13"),
    FunctionTestCase(function=isotope_symbol, args=("Hydrogen-3 +1",), expected="T"),
    FunctionTestCase(
        function=isotope_symbol,
        args="Md-260",
        kwargs={"mass_numb": 261},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="protium",
        kwargs={"mass_numb": 2},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="alpha",
        kwargs={"mass_numb": 3},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="O-18",
        kwargs={"mass_numb": 19},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="lead-209",
        kwargs={"mass_numb": 511},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol, args="He-1", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args=24,
        kwargs={"mass_numb": 23},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="H",
        kwargs={"mass_numb": 0},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="H-1",
        kwargs={"mass_numb": 2},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(function=isotope_symbol, args="P", expected=InvalidIsotopeError),
    FunctionTestCase(function=isotope_symbol, args=1, expected=InvalidIsotopeError),
    FunctionTestCase(function=isotope_symbol, args=4, expected=InvalidIsotopeError),
    FunctionTestCase(
        function=isotope_symbol, args="hydrogen-444444", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="Fe",
        kwargs={"mass_numb": 2.1},
        expected=TypeError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="He",
        kwargs={"mass_numb": "c"},
        expected=TypeError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="He-3",
        kwargs={"mass_numb": 4},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="D",
        kwargs={"mass_numb": 3},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="T",
        kwargs={"mass_numb": 2},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="Fe",
        kwargs={"mass_numb": None},
        expected=InvalidIsotopeError,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="He",
        kwargs={"mass_numb": 99},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(function=isotope_symbol, args="d", expected=InvalidParticleError),
    FunctionTestCase(
        function=isotope_symbol, args="h-3", expected=InvalidParticleError
    ),
    FunctionTestCase(function=isotope_symbol, args="h", expected=InvalidParticleError),
    FunctionTestCase(function=isotope_symbol, args="d+", expected=InvalidParticleError),
    FunctionTestCase(
        function=isotope_symbol,
        args="H-1",
        kwargs={"mass_numb": 1},
        expected=AtomicWarning,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="H-2",
        kwargs={"mass_numb": 2},
        expected=AtomicWarning,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="T",
        kwargs={"mass_numb": 3},
        expected=AtomicWarning,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="Li-6",
        kwargs={"mass_numb": 6},
        expected=AtomicWarning,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="lithium-6",
        kwargs={"mass_numb": 6},
        expected=AtomicWarning,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="alpha",
        kwargs={"mass_numb": 4},
        expected=AtomicWarning,
    ),
    FunctionTestCase(
        function=isotope_symbol,
        args="p",
        kwargs={"mass_numb": 1},
        expected=AtomicWarning,
    ),
    FunctionTestCase(function=atomic_number, args="H", expected=1),
    FunctionTestCase(function=atomic_number, args="deuterium", expected=1),
    FunctionTestCase(function=atomic_number, args="Deuterium", expected=1),
    FunctionTestCase(function=atomic_number, args="tritium", expected=1),
    FunctionTestCase(function=atomic_number, args="p", expected=1),
    FunctionTestCase(function=atomic_number, args="P", expected=15),
    FunctionTestCase(function=atomic_number, args="Alpha", expected=2),
    FunctionTestCase(function=atomic_number, args="C-12", expected=6),
    FunctionTestCase(function=atomic_number, args="Argon", expected=18),
    FunctionTestCase(function=atomic_number, args="protium", expected=1),
    FunctionTestCase(function=atomic_number, args="H-3", expected=1),
    FunctionTestCase(function=atomic_number, args="p+", expected=1),
    FunctionTestCase(function=atomic_number, args="Be-8", expected=4),
    FunctionTestCase(function=atomic_number, args="N", expected=7),
    FunctionTestCase(function=atomic_number, args="N 2+", expected=7),
    FunctionTestCase(function=atomic_number, args="N +1", expected=7),
    FunctionTestCase(function=atomic_number, args="N+++", expected=7),
    FunctionTestCase(
        function=atomic_number, args="H-3934", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=atomic_number, args="C-12b", expected=InvalidParticleError
    ),
    FunctionTestCase(function=atomic_number, args=-1.5, expected=TypeError),
    FunctionTestCase(function=atomic_number, args="n", expected=InvalidElementError),
    FunctionTestCase(function=atomic_number, args="n-1", expected=InvalidElementError),
    FunctionTestCase(
        function=atomic_number, args="neutron", expected=InvalidElementError
    ),
    FunctionTestCase(
        function=atomic_number, args="Neutron", expected=InvalidElementError
    ),
    FunctionTestCase(function=atomic_number, args="d", expected=InvalidParticleError),
    FunctionTestCase(function=atomic_number, args="t", expected=InvalidParticleError),
    FunctionTestCase(
        function=atomic_number, args="s-36", expected=InvalidParticleError
    ),
    FunctionTestCase(function=mass_number, args="helium-3", expected=3),
    FunctionTestCase(function=mass_number, args="Au-197", expected=197),
    FunctionTestCase(function=mass_number, args="deuterium", expected=2),
    FunctionTestCase(function=mass_number, args="D", expected=2),
    FunctionTestCase(function=mass_number, args="H-2", expected=2),
    FunctionTestCase(function=mass_number, args="tritium", expected=3),
    FunctionTestCase(function=mass_number, args="T", expected=3),
    FunctionTestCase(function=mass_number, args="alpha", expected=4),
    FunctionTestCase(function=mass_number, args="p", expected=1),
    FunctionTestCase(function=mass_number, args="Be-8", expected=8),
    FunctionTestCase(function=mass_number, args="N-13", expected=13),
    FunctionTestCase(function=mass_number, args="N-13 2+", expected=13),
    FunctionTestCase(function=mass_number, args="N-13 +2", expected=13),
    FunctionTestCase(function=mass_number, args="N-13+++", expected=13),
    FunctionTestCase(function=mass_number, args="H-359", expected=InvalidParticleError),
    FunctionTestCase(function=mass_number, args="C-12b", expected=InvalidParticleError),
    FunctionTestCase(function=mass_number, args=-1.5, expected=TypeError),
    FunctionTestCase(
        function=mass_number, args="N-13+-+-", expected=InvalidParticleError
    ),
    FunctionTestCase(function=mass_number, args="h-3", expected=InvalidParticleError),
    FunctionTestCase(function=mass_number, args="n", expected=InvalidIsotopeError),
    FunctionTestCase(function=mass_number, args="n-1", expected=InvalidIsotopeError),
    FunctionTestCase(function=element_name, args="D", expected="hydrogen"),
    FunctionTestCase(function=element_name, args="deuterium", expected="hydrogen"),
    FunctionTestCase(function=element_name, args="Au", expected="gold"),
    FunctionTestCase(function=element_name, args="alpha", expected="helium"),
    FunctionTestCase(function=element_name, args="helium-4", expected="helium"),
    FunctionTestCase(function=element_name, args="H-2", expected="hydrogen"),
    FunctionTestCase(function=element_name, args="Deuterium", expected="hydrogen"),
    FunctionTestCase(function=element_name, args="Hydrogen-3", expected="hydrogen"),
    FunctionTestCase(function=element_name, args="hydrogen-3", expected="hydrogen"),
    FunctionTestCase(function=element_name, args="H-3", expected="hydrogen"),
    FunctionTestCase(function=element_name, args="tritium", expected="hydrogen"),
    FunctionTestCase(function=element_name, args="Alpha", expected="helium"),
    FunctionTestCase(function=element_name, args="alpha", expected="helium"),
    FunctionTestCase(function=element_name, args=1, expected="hydrogen"),
    FunctionTestCase(function=element_name, args=26, expected="iron"),
    FunctionTestCase(function=element_name, args=79, expected="gold"),
    FunctionTestCase(function=element_name, args="p", expected="hydrogen"),
    FunctionTestCase(function=element_name, args="P", expected="phosphorus"),
    FunctionTestCase(function=element_name, args="Be-8", expected="beryllium"),
    FunctionTestCase(function=element_name, args="Li-7", expected="lithium"),
    FunctionTestCase(function=element_name, args="N", expected="nitrogen"),
    FunctionTestCase(function=element_name, args="N+++", expected="nitrogen"),
    FunctionTestCase(function=element_name, args="D-", expected="hydrogen"),
    FunctionTestCase(
        function=element_name, args="vegancupcakes", expected=InvalidParticleError
    ),
    FunctionTestCase(function=element_name, args="C-+-", expected=InvalidParticleError),
    FunctionTestCase(function=element_name, args=1.24, expected=TypeError),
    FunctionTestCase(function=element_name, args="n", expected=InvalidElementError),
    FunctionTestCase(
        function=element_name, args="neutron", expected=InvalidElementError
    ),
    FunctionTestCase(function=element_name, args=0, expected=InvalidParticleError),
    FunctionTestCase(function=element_name, args="H++", expected=InvalidParticleError),
    FunctionTestCase(function=element_name, args="t", expected=InvalidParticleError),
    FunctionTestCase(function=element_name, args="pb", expected=InvalidParticleError),
    FunctionTestCase(function=element_name, args="d", expected=InvalidParticleError),
    FunctionTestCase(function=element_name, args="h-3", expected=InvalidParticleError),
    FunctionTestCase(function=element_name, args="Pb-9", expected=InvalidParticleError),
    FunctionTestCase(function=element_name, args="H 2+", expected=InvalidParticleError),
    FunctionTestCase(
        function=standard_atomic_weight, args="H", expected=(1.008 * u.u).to(u.kg)
    ),
    FunctionTestCase(
        function=standard_atomic_weight, args=1, expected=(1.008 * u.u).to(u.kg)
    ),
    FunctionTestCase(
        function=standard_atomic_weight,
        args="Hydrogen",
        expected=(1.008 * u.u).to(u.kg),
    ),
    FunctionTestCase(function=standard_atomic_weight, args="Au", expected=u.kg),
    FunctionTestCase(function=standard_atomic_weight, args="H-1", expected=AtomicError),
    FunctionTestCase(function=standard_atomic_weight, args=1.1, expected=TypeError),
    FunctionTestCase(
        function=standard_atomic_weight, args="n", expected=InvalidElementError
    ),
    FunctionTestCase(function=standard_atomic_weight, args="p", expected=AtomicError),
    FunctionTestCase(
        function=standard_atomic_weight, args="alpha", expected=AtomicError
    ),
    FunctionTestCase(
        function=standard_atomic_weight, args="deuteron", expected=AtomicError
    ),
    FunctionTestCase(
        function=standard_atomic_weight, args="tritium", expected=AtomicError
    ),
    FunctionTestCase(function=standard_atomic_weight, args="Au+", expected=AtomicError),
    FunctionTestCase(
        function=standard_atomic_weight, args="Fe -2", expected=AtomicError
    ),
    FunctionTestCase(
        function=standard_atomic_weight, args="Og 2+", expected=AtomicError
    ),
    FunctionTestCase(
        function=standard_atomic_weight, args="h", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=standard_atomic_weight, args="fe", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=particle_mass, args="proton", expected=const.m_p.to("kg")
    ),
    FunctionTestCase(function=particle_mass, args="H-1+", expected=const.m_p.to("kg")),
    FunctionTestCase(
        function=particle_mass, args="H-1 +1", expected=const.m_p.to("kg")
    ),
    FunctionTestCase(
        function=particle_mass, args="H-1 1+", expected=const.m_p.to("kg")
    ),
    FunctionTestCase(
        function=particle_mass, args="H-1", kwargs={"Z": 1}, expected=const.m_p.to("kg")
    ),
    FunctionTestCase(
        function=particle_mass,
        args="hydrogen-1",
        kwargs={"Z": 1},
        expected=const.m_p.to("kg"),
    ),
    FunctionTestCase(function=particle_mass, args="p+", expected=const.m_p.to("kg")),
    FunctionTestCase(function=particle_mass, args="F-19", expected=u.kg),
    FunctionTestCase(
        function=particle_mass, args="Og 1+", expected=MissingAtomicDataError
    ),
    FunctionTestCase(
        function=particle_mass, args="Fe-56", kwargs={"Z": 1.4}, expected=TypeError
    ),
    FunctionTestCase(
        function=particle_mass,
        args="H-1 +1",
        kwargs={"Z": 0},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=particle_mass,
        args=26,
        kwargs={"Z": 1, "mass_numb": "a"},
        expected=TypeError,
    ),
    FunctionTestCase(
        function=particle_mass,
        args=26,
        kwargs={"Z": 27, "mass_numb": 56},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=particle_mass,
        args="Og",
        kwargs={"Z": 1},
        expected=MissingAtomicDataError,
    ),
    FunctionTestCase(
        function=particle_mass,
        args="Og",
        kwargs={"mass_numb": 696, "Z": 1},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=particle_mass,
        args="He 1+",
        kwargs={"mass_numb": 99},
        expected=InvalidParticleError,
    ),
    FunctionTestCase(
        function=particle_mass, args="fe-56 1+", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=particle_mass,
        args="H-1",
        kwargs={"mass_numb": 1, "Z": 1},
        expected=AtomicWarning,
    ),
    FunctionTestCase(
        function=particle_mass, args="H", expected=standard_atomic_weight("H")
    ),
    FunctionTestCase(function=is_stable, args="H-1", expected=True),
    FunctionTestCase(function=is_stable, args=(1, 1), expected=True),
    FunctionTestCase(function=is_stable, args="N-14", expected=True),
    FunctionTestCase(function=is_stable, args=("N", 14), expected=True),
    FunctionTestCase(function=is_stable, args="P-31", expected=True),
    FunctionTestCase(function=is_stable, args=("P", 31), expected=True),
    FunctionTestCase(function=is_stable, args="p", expected=True),
    FunctionTestCase(function=is_stable, args="alpha", expected=True),
    FunctionTestCase(function=is_stable, args="Xe-124", expected=True),
    FunctionTestCase(
        function=is_stable, args="Fe", kwargs={"mass_numb": 56}, expected=True
    ),
    FunctionTestCase(function=is_stable, args="Fe-56", expected=True),
    FunctionTestCase(function=is_stable, args="iron-56", expected=True),
    FunctionTestCase(function=is_stable, args="Iron-56", expected=True),
    FunctionTestCase(function=is_stable, args=(26, 56), expected=True),
    FunctionTestCase(function=is_stable, args="Be-8", expected=False),
    FunctionTestCase(function=is_stable, args="U-235", expected=False),
    FunctionTestCase(function=is_stable, args="uranium-235", expected=False),
    FunctionTestCase(function=is_stable, args="T", expected=False),
    FunctionTestCase(function=is_stable, args=(4, 8), expected=False),
    FunctionTestCase(function=is_stable, args="tritium", expected=False),
    FunctionTestCase(function=is_stable, args="Pb-209", expected=False),
    FunctionTestCase(function=is_stable, args="lead-209", expected=False),
    FunctionTestCase(function=is_stable, args="Lead-209", expected=False),
    FunctionTestCase(
        function=is_stable, args="Pb", kwargs={"mass_numb": 209}, expected=False
    ),
    FunctionTestCase(function=is_stable, args=(82, 209), expected=False),
    FunctionTestCase(
        function=is_stable, args=("hydrogen-444444",), expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=is_stable, args=("hydrogen", 0), expected=InvalidParticleError
    ),
    FunctionTestCase(function=is_stable, args=("",), expected=InvalidParticleError),
    FunctionTestCase(
        function=is_stable, args=("pb-209",), expected=InvalidParticleError
    ),
    FunctionTestCase(function=is_stable, args=("h",), expected=InvalidParticleError),
    FunctionTestCase(function=is_stable, args=("He",), expected=InvalidIsotopeError),
    FunctionTestCase(function=is_stable, args=("B",), expected=InvalidIsotopeError),
    FunctionTestCase(function=integer_charge, args="H+", expected=1),
    FunctionTestCase(function=integer_charge, args="D +1", expected=1),
    FunctionTestCase(function=integer_charge, args="tritium 1+", expected=1),
    FunctionTestCase(function=integer_charge, args="H-", expected=-1),
    FunctionTestCase(function=integer_charge, args="Fe -2", expected=-2),
    FunctionTestCase(function=integer_charge, args="Fe 2-", expected=-2),
    FunctionTestCase(function=integer_charge, args="N--", expected=-2),
    FunctionTestCase(function=integer_charge, args="N++", expected=2),
    FunctionTestCase(function=integer_charge, args="alpha", expected=2),
    FunctionTestCase(function=integer_charge, args="proton", expected=1),
    FunctionTestCase(function=integer_charge, args="deuteron", expected=1),
    FunctionTestCase(function=integer_charge, args="triton", expected=1),
    FunctionTestCase(function=integer_charge, args="electron", expected=-1),
    FunctionTestCase(function=integer_charge, args="e-", expected=-1),
    FunctionTestCase(function=integer_charge, args="e+", expected=1),
    FunctionTestCase(function=integer_charge, args="positron", expected=1),
    FunctionTestCase(function=integer_charge, args="n", expected=0),
    FunctionTestCase(function=integer_charge, args="neutron", expected=0),
    FunctionTestCase(function=integer_charge, args="p-", expected=-1),
    FunctionTestCase(function=integer_charge, args="antiproton", expected=-1),
    FunctionTestCase(
        function=integer_charge, args="fads", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=integer_charge, args="H++", expected=InvalidParticleError
    ),
    FunctionTestCase(function=integer_charge, args="h+", expected=InvalidParticleError),
    FunctionTestCase(
        function=integer_charge, args="fe 1+", expected=InvalidParticleError
    ),
    FunctionTestCase(function=integer_charge, args="d+", expected=InvalidParticleError),
    FunctionTestCase(
        function=integer_charge, args="Fe 29+", expected=InvalidParticleError
    ),
    FunctionTestCase(function=integer_charge, args="H-1", expected=ChargeError),
    FunctionTestCase(function=integer_charge, args="H---", expected=AtomicWarning),
    FunctionTestCase(function=integer_charge, args="Fe -26", expected=AtomicWarning),
    FunctionTestCase(function=integer_charge, args="Og 10-", expected=AtomicWarning),
    FunctionTestCase(function=electric_charge, args="p", expected=u.C),
    FunctionTestCase(function=electric_charge, args="p", expected=const.e),
    FunctionTestCase(
        function=electric_charge, args="e", expected=-1.6021766208e-19 * u.C
    ),
    FunctionTestCase(
        function=electric_charge, args="alpha", expected=3.2043532416e-19 * u.C
    ),
    FunctionTestCase(function=electric_charge, args="n", expected=0 * u.C),
    FunctionTestCase(
        function=electric_charge, args="bad input", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=electric_charge, args="h+", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=electric_charge, args="Au 81+", expected=InvalidParticleError
    ),
    FunctionTestCase(function=electric_charge, args="Au 81-", expected=AtomicWarning),
    FunctionTestCase(function=electric_charge, args="H---", expected=AtomicWarning),
    FunctionTestCase(function=half_life, args="H-1", expected=u.s),
    FunctionTestCase(function=half_life, args="tritium", expected=u.s),
    FunctionTestCase(function=half_life, args="H-1", expected=np.inf * u.s),
    FunctionTestCase(
        function=half_life, args="No-248", expected=MissingAtomicDataError
    ),  # isotope with missing data
    FunctionTestCase(
        function=periodic_table_period, args=("Ne", "Na"), expected=TypeError
    ),
    FunctionTestCase(
        function=periodic_table_block, args=("N", "C", "F"), expected=TypeError
    ),
    FunctionTestCase(
        function=periodic_table_category, args=("Rb", "He", "Li"), expected=TypeError
    ),
    FunctionTestCase(
        function=periodic_table_group, args=("B", "Ti", "Ge"), expected=TypeError
    ),
    FunctionTestCase(
        function=isotopic_abundance,
        args=("H", 1),
        expected=isotopic_abundance("protium"),
    ),
    FunctionTestCase(function=isotopic_abundance, args="D", expected=0.000115),
    FunctionTestCase(function=isotopic_abundance, args="Be-8", expected=0.0),
    FunctionTestCase(function=isotopic_abundance, args="Li-8", expected=0.0),
    FunctionTestCase(
        function=isotopic_abundance, args=("Og", 294), expected=AtomicWarning
    ),
    FunctionTestCase(
        function=isotopic_abundance, args="neutron", expected=InvalidIsotopeError
    ),
    FunctionTestCase(
        function=isotopic_abundance, args="Og-2", expected=InvalidParticleError
    ),
    FunctionTestCase(
        function=particle_mass,
        args="berkelium-249",
        expected=(249.0749877 * u.u).to("kg"),
    ),
    FunctionTestCase(function=common_isotopes, args="n", expected=InvalidElementError),
    FunctionTestCase(function=stable_isotopes, args="n", expected=InvalidElementError),
    FunctionTestCase(function=known_isotopes, args="n", expected=InvalidElementError),
    FunctionTestCase(
        function=reduced_mass, args=("N", 6e-26 * u.l), expected=u.UnitConversionError
    ),
    FunctionTestCase(
        function=reduced_mass, args=("Og", "H"), expected=MissingAtomicDataError
    ),
    FunctionTestCase(
        function=known_isotopes,
        args="H",
        expected=["H-1", "D", "T", "H-4", "H-5", "H-6", "H-7"],
    ),
    FunctionTestCase(function=common_isotopes, args="H", expected=["H-1", "D"]),
    FunctionTestCase(function=stable_isotopes, args="He", expected=["He-3", "He-4"]),
    FunctionTestCase(function=_is_electron, args="e-", expected=True),
    FunctionTestCase(function=_is_electron, args="e+", expected=False),
    FunctionTestCase(function=_is_electron, args="e", expected=True),
    FunctionTestCase(function=_is_electron, args="electron", expected=True),
    FunctionTestCase(function=_is_electron, args="ELECTRON", expected=True),
    FunctionTestCase(function=_is_electron, args="H", expected=False),
    FunctionTestCase(function=_is_electron, args="positron", expected=False),
    FunctionTestCase(function=_is_electron, args="Carbon", expected=False),
    FunctionTestCase(function=half_life, args="tritium", expected=3.888e8 * u.s),
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
            list_of_tests.append(
                FunctionTestCase(
                    function=function, args=invalid_particle, expected=exception
                )
            )


append_invalid_particle_tests(test_cases)


@pytest.mark.parametrize("test_case", test_cases)
def test_atomic_functions(test_case):
    """
    Test that providing the principal functions in `~plasmapy.particles.atomic`
    result in the expected outcome.
    """

    test_runner(test_case)


# @pytest.mark.parametrize("function, args, kwargs, expected", test_cases)
# def test_particle_functions(function, args, kwargs, expected):
#    """

#    """

#    function_test_runner(function=function, args=args, kwargs=kwargs, expected=expected)


def test_standard_atomic_weight_value_between():
    """Test that `standard_atomic_weight` returns approximately the
    correct value for phosphorus."""
    assert (
        30.973 < standard_atomic_weight("P").to(u.u).value < 30.974
    ), "Incorrect standard atomic weight for phosphorus."


def test_particle_mass_for_hydrogen_with_no_mass_number():
    """Test that `particle_mass` does not return the proton mass when no
    mass number is specified for hydrogen.  In this case, the
    standard atomic weight should be used to account for the small
    fraction of deuterium."""
    assert particle_mass("H", Z=1) > const.m_p
    assert particle_mass("hydrogen", Z=1) > const.m_p


def test_particle_mass_helium():
    """Test miscellaneous cases for `particle_mass`."""
    assert particle_mass("alpha") > particle_mass("He-3 2+")


# (arg1, kwargs1, arg2, kwargs2, expected)
equivalent_particle_mass_args = [
    ["e+", {}, "positron", {}, const.m_e],
    ["alpha", {}, "He-4++", {}, None],
    ["alpha", {}, "helium-4 2+", {}, None],
    ["deuteron", {}, "H", {"Z": 1, "mass_numb": 2}, None],
    ["D+", {}, "H-2+", {}, None],
    ["D+", {}, "D 1+", {}, None],
    ["Deuterium+", {}, "D", {"Z": 1}, None],
    ["triton", {}, "H", {"Z": 1, "mass_numb": 3}, None],
    ["T+", {}, "H-3+", {}, None],
    ["T+", {}, "T 1+", {}, None],
    ["Tritium+", {}, "T", {"Z": 1}, None],
    [
        "Fe-56 1+",
        {},
        "Fe",
        {"mass_numb": 56, "Z": 1},
        particle_mass("Fe-56 1-") - 2 * const.m_e,
    ],
    ["Fe-56 +1", {}, 26, {"mass_numb": 56, "Z": 1}, None],
]


@pytest.mark.parametrize(
    "arg1, kwargs1, arg2, kwargs2, expected", equivalent_particle_mass_args
)
def test_particle_mass_equivalent_args(arg1, kwargs1, arg2, kwargs2, expected):
    """Test that `particle_mass` returns equivalent results for
    equivalent positional and keyword arguments."""

    result1 = particle_mass(arg1, **kwargs1)
    result2 = particle_mass(arg2, **kwargs2)

    assert u.isclose(result1, result2), (
        f"particle_mass({repr(arg1)}, **{kwargs1}) = {repr(result1)}, whereas "
        f"particle_mass({repr(arg2)}, **{kwargs2}) = {repr(result2)}.  "
        f"These results are not equivalent as expected."
    )

    if expected is not None:
        assert u.isclose(result1, result2) and u.isclose(result2, expected), (
            f"particle_mass({repr(arg1)}, **{kwargs1}) = {repr(result1)} and "
            f"particle_mass({repr(arg2)}, **{kwargs2}) = {repr(result2)}, but "
            f"these results are not equal to {repr(expected)} as expected."
        )


def test_half_life_unstable_isotopes():
    """Test that `half_life` returns `None` and raises an exception for
    all isotopes that do not yet have half-life data."""
    for isotope in _Isotopes.keys():
        if (
            "half_life" not in _Isotopes[isotope].keys()
            and not _Isotopes[isotope].keys()
        ):
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


@pytest.mark.parametrize(
    "function, element, kwargs, isotope, should_be_in_list", isotope_inputs
)
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

    assert len(common_isotopes()) == 288, (
        "The length of the list returned by common_isotopes() is "
        f"{len(common_isotopes())}, which is not the expected value."
    )

    assert len(stable_isotopes()) == 254, (
        "The length of the list returned by stable_isotopes() is "
        f"{len(stable_isotopes())}, which is not the expected value."
    )

    assert 3352 <= len(known_isotopes()) <= 3400, (
        "The length of the list returned by known_isotopes() is "
        f"{len(known_isotopes())}, which is not within the expected range."
    )


isotopic_abundance_elements = (
    atomic_number(atomic_numb) for atomic_numb in range(1, 119)
)

isotopic_abundance_isotopes = (
    common_isotopes(element) for element in isotopic_abundance_elements
)

isotopic_abundance_sum_table = (
    (element, isotopes)
    for element, isotopes in zip(
        isotopic_abundance_elements, isotopic_abundance_isotopes
    )
    if isotopes
)


@pytest.mark.parametrize("element, isotopes", isotopic_abundance_sum_table)
def test_isotopic_abundances_sum(element, isotopes):
    """Test that the sum of isotopic abundances for each element with
    isotopic abundances is one."""
    sum_of_iso_abund = sum(isotopic_abundance(isotope) for isotope in isotopes)
    assert np.isclose(
        sum_of_iso_abund, 1, atol=1e-6
    ), f"The sum of the isotopic abundances for {element} does not equal 1."
