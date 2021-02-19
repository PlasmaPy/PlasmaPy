import numpy as np

from astropy import constants
from scipy.special import erf
from typing import List

from plasmapy.formulary import thermal_speed
from plasmapy.formulary.mathematics import Chandrasekhar_G


def xab_ratio(a, b):
    return thermal_speed(b.T_e, b.base_particle) / thermal_speed(a.T_e, a.base_particle)


def M_matrix(species_a, species_b):
    a, b = species_a, species_b
    xab = xab_ratio(a, b)
    temperature_ratio = a.T_e / b.T_e
    mass_ratio = a._particle.mass / b._particle.mass
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
    temperature_ratio = a.T_e / b.T_e
    mass_ratio = a._particle.mass / b._particle.mass
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

CL = lambda a, b: Coulomb_logarithm(
    b.T_e,
    b.n_elem,
    (a.base_particle, b.base_particle),  # simplifying assumption after A4
)


def effective_momentum_relaxation_rate(charge_states_a, charge_states_b):
    def contributions():
        CL = lambda ai, bj: Coulomb_logarithm(
            charge_states_b.T_e,
            charge_states_b.n_elem,
            (ai._particle, bj._particle),  # simplifying assumption after A4
        )
        for i in range(
            len(charge_states_a.number_densities)
        ):  # These really need to implement len.
            ai = charge_states_a[i]
            if ai._particle.charge == 0:
                continue
            for j in range(len(charge_states_b.number_densities)):
                bj = charge_states_b[j]
                if bj._particle.charge == 0:
                    continue
                # Eq. A4, Houlberg_1997
                # Eq. A3, Houlberg_1997
                collision_frequency_ai_bj = (
                    4
                    / (3 * np.sqrt(np.pi))
                    * (
                        4
                        * np.pi
                        * ai._particle.charge ** 2
                        * bj._particle.charge ** 2
                        * bj.number_density
                        * CL(ai, bj)
                    )
                    / (
                        (4 * np.pi * constants.eps0) ** 2
                        * ai._particle.mass ** 2
                        * thermal_speed(charge_states_a.T_e, ai._particle) ** 3
                    )
                ).si

                yield (
                    ai.number_density * ai._particle.mass
                ) * collision_frequency_ai_bj

    return sum(contributions())


def charge_weighting_factor(i, a_states):
    """Charge weighting factor; Houlberg_1997, eq.11

    TODO Is this not acctually Z_eff, or some contribution to it? Should that not be a method of IonizationState?
    """
    ai = a_states[i]
    denominator = sum(
        state.number_density * state.integer_charge ** 2 for state in a_states
    )
    return ai.number_density * ai.integer_charge ** 2 / denominator


def N_script(species_a, species_b):
    N = N_matrix(species_a, species_b)
    # Equation A2b
    N_script = effective_momentum_relaxation_rate(species_a, species_b) * N
    return N_script


def M_script(species_a, all_species):
    # Equation A2a
    def gener():
        for species_b in all_species:
            if species_b is not species_a:  # direct comparison glitches out
                yield M_matrix(
                    species_a, species_b
                ) * effective_momentum_relaxation_rate(species_a, species_b)

    return sum(gener())


def L_friction_coefficient(species_a, i, species_b, j, all_species):
    parenthesis = N_script(species_a, species_b) * charge_weighting_factor(j, species_b)
    if species_a is species_b and i == j:
        parenthesis += M_script(species_a, all_species)
    return charge_weighting_factor(i, species_a) * parenthesis


def pitch_angle_diffusion_rate(x, index, a_states, all_species):
    # Houlberg_1997, equation B4b,
    ai = a_states[index]

    def sum_items():
        for b in all_species:
            xab = xab_ratio(
                a_states, b
            )  # TOZO should be over all other species, not ionization states
            numerator = erf(x / xab) - Chandrasekhar_G(x / xab)
            denominator = x ** 3
            fraction = numerator / denominator
            result = fraction * effective_momentum_relaxation_rate(a_states, b)
            yield result

    return (
        charge_weighting_factor(index, a_states)
        / (ai.number_density * ai._particle.mass)
        * 3
        * np.sqrt(np.pi)
        / 4
        * sum(sum_items())
    )
