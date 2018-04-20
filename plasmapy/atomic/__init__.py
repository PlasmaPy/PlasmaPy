"""
The `plasmapy.atomic` subpackage provides access to information about
atoms, isotopes, ions, and other particles.
"""

from .special_particles import ParticleZoo

from .particle_class import Particle
from .particle_input import particle_input

from .symbols import (
    atomic_symbol,
    isotope_symbol,
    ionic_symbol,
    particle_symbol,
    element_name,
)

from .atomic import (
    atomic_number,
    is_stable,
    half_life,
    mass_number,
    standard_atomic_weight,
    particle_mass,
    known_isotopes,
    common_isotopes,
    stable_isotopes,
    isotopic_abundance,
    integer_charge,
    electric_charge,
    reduced_mass,
)

from .nuclear import (
    nuclear_binding_energy,
    nuclear_reaction_energy,
)

proton = Particle("p+")
electron = Particle("e-")
neutron = Particle("n")
positron = Particle("e+")
deuteron = Particle("D 1+")
triton = Particle("T 1+")
alpha = Particle("He-4 2+")
