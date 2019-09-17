"""Functions from Atomic to benchmark with aerospeed velocity """

import astropy.constants as const
import astropy.units as u
from plasmapy.atomic import (
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
    periodic_table_group)


class PhysicsSuite:
    """
    Benchmark that times the performance of funcions from Atomic package
    """
    def setup(self):
        pass


    def time_atomic_number(self):
        atomic_number("H")
        atomic_number("tritium")
        atomic_number("alpha")
        atomic_number("oganesson")


    def time_mass_number(self):
        mass_number("H-1")
        mass_number("Pb-208")
        mass_number("tritium")
        mass_number("alpha")


    def time_standard_atomic_weight(self):
        standard_atomic_weight("H")
        standard_atomic_weight("lead")
        standard_atomic_weight("Au")
        standard_atomic_weight("nitrogen")


    def time_particle_mass(self):
        particle_mass("Au")
        particle_mass("xenon")
        particle_mass(2)
        particle_mass("90")


    def time_isotopic_abundance(self):
        isotopic_abundance('Pb-208')
        isotopic_abundance('hydrogen', 1)
        isotopic_abundance('hydrogen', 3)
        isotopic_abundance('hydrogen', 2)


    def time_integer_charge(self):
        integer_charge('Fe-56 2+')
        integer_charge('He -2')
        integer_charge('H+')
        integer_charge('N-14++')


    def time_electric_charge(self):
        electric_charge('p+')
        electric_charge('O-')
        electric_charge('N-')
        electric_charge('He-')


    def time_is_stable(self):
        is_stable("H-1")
        is_stable("tritium")
        is_stable("e-")
        is_stable("tau+")


    def time_half_life(self):
        half_life('T')
        half_life('n')
        half_life('H-1')
        half_life('C-14')


    def time_known_isotopes(self):
        known_isotopes('H')
        known_isotopes('helium 1+')
        known_isotopes()[0:10]
        known_isotopes()[4:15]


    def time_common_isotopes(self):
        common_isotopes('H')
        common_isotopes('Fe')
        common_isotopes('Fe', most_common_only=True)
        common_isotopes()[0:7]


    def time_stable_isotopes(self):
        stable_isotopes('beryllium')
        stable_isotopes('Pb-209')
        stable_isotopes(118)
        stable_isotopes('U', unstable=True)[:5]


    def time_reduced_mass(self):
        reduced_mass('p+', 'e-')
        reduced_mass(5.4e-27 * u.kg, 8.6e-27 * u.kg)
        reduced_mass(6.4e-10 * u.kg, 8.6e-11 * u.kg)
        reduced_mass(1.4e-10 * u.kg, 2.6e-40 * u.kg)


    def time_periodic_table_period(self):
        periodic_table_period(5)
        periodic_table_period("5")
        periodic_table_period("Au")
        periodic_table_period("nitrogen")


    def time_periodic_table_group(self):
        periodic_table_group(18)
        periodic_table_group(24)
        periodic_table_group("Al")
        periodic_table_group("neon")


    def time_periodic_table_block(self):
        periodic_table_block(66)
        periodic_table_block("Tl")
        periodic_table_block("thallium")
        periodic_table_block("francium")


    def time_periodic_table_category(self):
        periodic_table_category(82)
        periodic_table_category("85")
        periodic_table_category("rhodium")
        periodic_table_category(41)
