import numpy as np

from astropy import constants
from typing import List

from plasmapy.particles.species import Species


def xab_ratio(a, b):
    return b.thermal_speed() / a.thermal_speed()


def M_matrix(species_a, species_b):
    a, b = species_a, species_b
    xab = xab_ratio(a, b)
    temperature_ratio = a.temperature / b.temperature
    mass_ratio = a.particle.mass / b.particle.mass
    """equations A5a through A5f, Houlberg_1997"""
    M11 = -(1 + mass_ratio) / (1 + xab ** 2) ** (3 / 2)
    M12 = 3 / 2 * (1 + mass_ratio) / (1 + xab ** 2) ** (5 / 2)
    M21 = M12
    M22 = (13 / 4 + 4 * xab ** 2 + 15 / 2 * xab ** 4) / (1 + xab ** 2) ** (5 / 2)
    M13 = -15 / 8 * (1 + mass_ratio) / (1 + xab ** 2) ** (7 / 2)
    M31 = M13
    M23 = (69 / 16 + 6 * xab ** 2 + 63 / 4 * xab ** 4) / (1 + xab ** 2) ** (7 / 2)
    M32 = M23
    M33 = -(433 / 64 + 17 * xab ** 2 + 459 / 8 * xab ** 4 + 175 / 8 * xab ** 6) / (
        1 + xab ** 2
    ) ** (9 / 2)
    M = np.array([[M11, M12, M13], [M21, M22, M23], [M31, M32, M33]])
    return M


def N_matrix(species_a, species_b):
    """equations A6a through A6f, Houlberg_1997"""
    a, b = species_a, species_b
    xab = xab_ratio(a, b)
    temperature_ratio = a.temperature / b.temperature
    mass_ratio = a.particle.mass / b.particle.mass
    N11 = (1 + mass_ratio) / (1 + xab ** 2) ** (3 / 2)
    N21 = -3 / 2 * (1 + mass_ratio) / (1 + xab ** 2) ** (5 / 2)
    N31 = 15 / 8 * (1 + mass_ratio) / (1 + xab ** 2) ** (7 / 2)
    M12 = 3 / 2 * (1 + mass_ratio) / (1 + xab ** 2) ** (5 / 2)
    N12 = -(xab ** 2) * M12
    N22 = 27 / 4 * (temperature_ratio) ** (1 / 2) * xab ** 2 / (1 + xab ** 2) ** (5 / 2)
    M13 = -15 / 8 * (1 + mass_ratio) / (1 + xab ** 2) ** (7 / 2)
    N13 = -(xab ** 4) * M13
    N23 = -225 / 16 * temperature_ratio * xab ** 4 / (1 + xab ** 2) ** (7 / 2)
    N32 = N23 / temperature_ratio
    N33 = (
        2625 / 64 * temperature_ratio ** (1 / 2) * xab ** 4 / (1 + xab ** 2) ** (9 / 2)
    )
    N = np.array([[N11, N12, N13], [N21, N22, N23], [N31, N32, N33]])
    return N


from plasmapy.formulary.collisions import Coulomb_logarithm

CL = lambda ai, bj: Coulomb_logarithm(
    bj.temperature, bj.number_density, (ai.particle, bj.particle)
)


def collision_frequency_ai_bj(ai, bj):
    return (
        4
        / (3 * np.sqrt(np.pi))
        * (
            4
            * np.pi
            * ai.particle.charge ** 2
            * bj.particle.charge ** 2
            * bj.number_density
            * CL(ai, bj)
        )
        / (
            (4 * np.pi * constants.eps0) ** 2
            * ai.particle.mass ** 2
            * ai.thermal_speed() ** 3
        )
    ).si


def effective_momentum_relaxation_rate(
    charge_states_a: List[Species], charge_states_b: List[Species]
):
    def contributions():
        for ai in charge_states_a:
            for bj in charge_states_b:
                # Eq. A4, Houlberg_1997
                # Eq. A3, Houlberg_1997
                yield ai.mass_density() * collision_frequency_ai_bj(ai, bj)

    return sum(contributions())


def charge_weighting_factor(i, a_states):
    """Charge weighting factor; Houlberg_1997, eq.11

    TODO Is this not acctually Z_eff, or some contribution to it? Should that not be a method of IonizationState?
    """
    ai = a_states[i]
    denominator = sum(
        state.number_density * state.particle.integer_charge ** 2 for state in a_states
    )
    return ai.number_density * ai.particle.integer_charge ** 2 / denominator
