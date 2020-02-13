"""
The `plasmapy.particles` subpackage provides access to information about
atoms, isotopes, ions, and other particles.
"""

from plasmapy.particles.special_particles import ParticleZoo
from plasmapy.particles.particle_input import particle_input
from plasmapy.particles.particle_class import (
    Particle,
    AbstractParticle,
    DimensionlessParticle,
    CustomParticle,
)

from plasmapy.particles.symbols import (
    atomic_symbol,
    isotope_symbol,
    ionic_symbol,
    particle_symbol,
    element_name,
)

from plasmapy.particles.atomic import (
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

from plasmapy.particles.nuclear import nuclear_binding_energy, nuclear_reaction_energy

from plasmapy.particles.ionization_state import IonizationState, State
from plasmapy.particles.ionization_states import IonizationStates

# Create instances of the most commonly used particles

proton = Particle("p+")
electron = Particle("e-")
neutron = Particle("n")
positron = Particle("e+")
deuteron = Particle("D 1+")
triton = Particle("T 1+")
alpha = Particle("He-4 2+")
