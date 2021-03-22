import numpy as np

from astropy import constants
from astropy import units as u
from scipy import integrate
from scipy.special import erf

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
            (ai.ion, bj.ion),  # simplifying assumption after A4
        )
        for ai in charge_states_a:
            if ai.ion.charge == 0:
                continue
            for bj in charge_states_b:
                if bj.ion.charge == 0:
                    continue
                # Eq. A4, Houlberg_1997
                # Eq. A3, Houlberg_1997
                collision_frequency_ai_bj = (
                    4
                    / (3 * np.sqrt(np.pi))
                    * (
                        4
                        * np.pi
                        * ai.ion.charge ** 2
                        * bj.ion.charge ** 2
                        * bj.number_density
                        * CL(ai, bj)
                    )
                    / (
                        (4 * np.pi * constants.eps0) ** 2
                        * ai.ion.mass ** 2
                        * thermal_speed(charge_states_a.T_e, ai.ion) ** 3
                    )
                ).si

                yield (ai.number_density * ai.ion.mass) * collision_frequency_ai_bj

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


# TODO refactor for ionizationstatecollection
def L_friction_coefficient(species_a, i, species_b, j, all_species):
    # Houlberg_1997, equation 10
    parenthesis = N_script(species_a, species_b) * charge_weighting_factor(j, species_b)
    if species_a == species_b and i == j:
        parenthesis += M_script(species_a, all_species)
    return charge_weighting_factor(i, species_a) * parenthesis


def pitch_angle_diffusion_rate(x, index, a_states, all_species):
    # Houlberg_1997, equation B4b,
    ai = a_states[index]  # TODO I wouldn't need to carry the index around, if...

    def sum_items():
        for b in all_species:
            xab = xab_ratio(
                a_states, b
            )  # TODO should be over all other species, not ionization states
            numerator = erf(x / xab) - Chandrasekhar_G(x / xab)
            denominator = x ** 3
            fraction = numerator / denominator
            result = fraction * effective_momentum_relaxation_rate(a_states, b)
            yield result

    return (
        charge_weighting_factor(index, a_states)
        / (ai.number_density * ai.ion.mass)
        * 3
        * np.sqrt(np.pi)
        / 4
        * sum(sum_items())
    )


def K_B_ai(x, index, a_states, all_species, flux_surface, *, orbit_squeezing=False):
    # eq. B1-B4, Houlberg_1997
    f_t = flux_surface.trapped_fraction()
    f_c = 1 - f_t
    if orbit_squeezing:
        raise NotImplementedError(
            "TODO allow for non-zero, changing radial electric fields (orbit squeezing)"
        )
    else:
        S_ai = B2 = 1
    return (
        pitch_angle_diffusion_rate(x, index, a_states, all_species)
        * f_t
        / f_c
        / S_ai ** 1.5
    )


LaguerrePolynomials = [
    lambda x: 1,
    lambda x: 5 / 2 - x,
    lambda x: 35 / 8 - 7 / 2 * x + 1 / 2 * x ** 2,
]


def K(ai, x, *args, **kwargs):
    return np.ones_like(x)


def mu_hat(ai, return_with_unc: bool = False):
    orders = range(1, 4)
    mu_hat_ai = np.zeros((3, 3))
    if return_with_unc:
        dmu_hat_ai = np.zeros((3, 3))
    π = np.pi
    for α in orders:
        for β in orders:
            integrand = lambda x: (
                x ** 4
                * np.exp(-(x ** 2))
                * LaguerrePolynomials[α - 1](x ** 2)
                * LaguerrePolynomials[β - 1](x ** 2)
                * K(ai, x)
            )
            integral = integrate.quad(integrand, 0, np.inf)
            mass_density_probably = ai.number_density * ai.ion.mass
            value, stderr = integral
            mu_hat_ai[α - 1, β - 1] = value * (-1) ** (α + β)
            if return_with_unc:
                dmu_hat_ai[α - 1, β - 1] = stderr
    actual_units = lambda value: (
        8 / 3 / np.sqrt(π) * value / u.s * mass_density_probably
    )
    mu_hat_ai = actual_units(mu_hat_ai)
    if return_with_unc:
        dmu_hat_ai = actual_units(dmu_hat_ai)
        return mu_hat_ai, dmu_hat_ai
    else:
        return mu_hat_ai


def ξ(isotope):
    array = u.Quantity(
        [ai.number_density * ai.ion.integer_charge ** 2 for ai in isotope]
    )
    return array / array.sum()


def rbar(a, all_species, beta_coeffs=None) -> u.Quantity:
    if beta_coeffs is not None:
        # TODO should be a dict or sth
        raise NotImplementedError
    else:
        beta_cx = np.zeros(3)  # TODO
        beta_an = np.zeros(3)  # TODO
        beta_coeffs = np.diag(beta_cx + beta_an) * u.kg / u.m ** 3 / u.s

    def gen():
        for i, ai in enumerate(a):
            if ξ(a)[i] == 0:
                continue  # won't add anything to sum anyway, and matrix gets singular
            Aai = ξ(a)[i] * M_script(a, all_species) - mu_hat(ai) - beta_coeffs
            S_matrix = ξ(a)[i] * np.eye(3)
            rai_as_rows = np.linalg.solve(Aai, S_matrix)
            # TODO does not include r_pT, r_E, r_NBI yet
            rbar_ingredient = ξ(a)[i] * rai_as_rows
            yield rbar_ingredient

    return sum(gen())


def eq34matrix(all_species, beta_coeffs=None):
    output_matrix = u.Quantity(np.eye(3 * len(all_species)))

    for I, a in enumerate(all_species):
        i = 3 * I
        rarray = rbar(a, all_species, beta_coeffs)
        for J, b in enumerate(all_species):
            j = 3 * J
            narray = N_script(a, b).sum(axis=0, keepdims=True)
            result = narray * rarray.T
            output_matrix[i : i + 3, j : j + 3] += result

    return output_matrix
