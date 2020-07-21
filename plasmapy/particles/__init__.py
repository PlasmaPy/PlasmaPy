"""
The `plasmapy.particles` subpackage provides access to information about
atoms, isotopes, ions, and other particles.
"""
__all__ = [
    "json_load_particle",
    "json_loads_particle",
    "AbstractParticle",
    "CustomParticle",
    "DimensionlessParticle",
    "Particle",
    "ParticleJSONDecoder",
]

from plasmapy.particles.atomic import (
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
    reduced_mass,
    stable_isotopes,
    standard_atomic_weight,
)
from plasmapy.particles.ionization_state import IonizationState, State
from plasmapy.particles.ionization_states import IonizationStates
from plasmapy.particles.nuclear import nuclear_binding_energy, nuclear_reaction_energy
from plasmapy.particles.particle_class import (
    AbstractParticle,
    CustomParticle,
    DimensionlessParticle,
    Particle,
)
from plasmapy.particles.particle_input import particle_input
from plasmapy.particles.serialization import (
    json_load_particle,
    json_loads_particle,
    ParticleJSONDecoder,
)
from plasmapy.particles.special_particles import ParticleZoo
from plasmapy.particles.symbols import (
    atomic_symbol,
    element_name,
    ionic_symbol,
    isotope_symbol,
    particle_symbol,
)

# Create instances of the most commonly used particles

proton = Particle("p+")
electron = Particle("e-")
neutron = Particle("n")
positron = Particle("e+")
deuteron = Particle("D 1+")
triton = Particle("T 1+")
alpha = Particle("He-4 2+")
