from __future__ import annotations

import numpy as np

from astropy import constants
from astropy import units as u
from scipy.special import erf
from typing import Union

from plasmapy.formulary import thermal_speed
from plasmapy.formulary.collisions import Coulomb_logarithm
from plasmapy.formulary.mathematics import Chandrasekhar_G
from plasmapy.particles import IonizationState, IonizationStateCollection
from plasmapy.plasma.fluxsurface import FluxSurface

try:
    from scipy.integrate import trapz as trapezoid
except ImportError:
    from scipy.integrate import trapezoid


def xab_ratio(a: IonizationState, b: IonizationState):
    return thermal_speed(b.T_e, b.base_particle) / thermal_speed(a.T_e, a.base_particle)


def M_matrix(species_a: IonizationState, species_b: IonizationState):
    a, b = species_a, species_b
    xab = xab_ratio(a, b)
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


def N_matrix(species_a: IonizationState, species_b: IonizationState):
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


CL = lambda a, b: Coulomb_logarithm(
    b.T_e,
    b.n_elem,
    (a.base_particle, b.base_particle),  # simplifying assumption after A4
)


def effective_momentum_relaxation_rate(
    charge_states_a: IonizationState, charge_states_b: IonizationState
):
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


def ξ(isotope: IonizationState):
    array = u.Quantity(
        [ai.number_density * ai.ion.integer_charge ** 2 for ai in isotope]
    )
    return array / array.sum()


def N_script(species_a: IonizationState, species_b: IonizationState):
    N = N_matrix(species_a, species_b)
    # Equation A2b
    N_script = effective_momentum_relaxation_rate(species_a, species_b) * N
    return N_script


def M_script(species_a: IonizationState, all_species: IonizationStateCollection):
    # Equation A2a
    def gener():
        for species_b in all_species:
            if species_b is not species_a:  # direct comparison glitches out
                # TODO am I sure this if is necessary here?
                yield M_matrix(
                    species_a, species_b
                ) * effective_momentum_relaxation_rate(species_a, species_b)

    return sum(gener())


def pitch_angle_diffusion_rate(
    x: np.ndarray,
    index: int,
    a_states: IonizationState,
    all_species: IonizationStateCollection,
):
    # Houlberg_1997, equation B4b,
    ai = a_states[index]  # TODO I wouldn't need to carry the index around, if...
    xi = ξ(a_states)[index]

    def sum_items():
        for b in all_species:
            xab = xab_ratio(a_states, b)
            numerator = erf(x / xab) - Chandrasekhar_G(x / xab)
            denominator = x ** 3
            fraction = numerator / denominator
            result = fraction * effective_momentum_relaxation_rate(a_states, b)
            yield result

    result = (
        xi
        / (ai.number_density * ai.ion.mass)
        * 3
        * np.sqrt(np.pi)
        / 4
        * sum(sum_items())
    )
    return result


def K_B_ai(
    x: np.ndarray,
    index: int,
    a_states: IonizationState,
    all_species: IonizationStateCollection,
    flux_surface: FluxSurface,
    *,
    orbit_squeezing=False
):
    # eq. B1-B4, Houlberg_1997
    f_t = flux_surface.trapped_fraction()
    f_c = 1 - f_t
    if orbit_squeezing:
        raise NotImplementedError(
            "TODO allow for non-zero, changing radial electric fields (orbit squeezing)"
        )
    else:
        S_ai = 1  # Equation B2
    padr = pitch_angle_diffusion_rate(x, index, a_states, all_species)
    return padr * f_t / f_c / S_ai ** 1.5


LaguerrePolynomials = [
    lambda x: 1,
    lambda x: 5 / 2 - x,
    lambda x: 35 / 8 - 7 / 2 * x + 1 / 2 * x ** 2,
]


def F_m(m: Union[int, np.ndarray], flux_surface: FluxSurface, g=1):
    fs = flux_surface
    B20 = fs.Brvals * fs.Bprimervals + fs.Bzvals * fs.Bprimezvals
    if g == 0:
        return np.nan
    under_average_B16 = np.sin(m * fs.Theta) * B20
    under_average_B15 = under_average_B16 / fs.Bmag
    under_average_B16_cos = np.cos(m * fs.Theta) * B20
    under_average_B15_cos = under_average_B16_cos / fs.Bmag
    #     plt.plot(fs.lp, under_average_B15)
    #     plt.plot(fs.lp, under_average_B16)
    B15 = fs.flux_surface_average(under_average_B15)
    B16 = fs.gamma * fs.flux_surface_average(under_average_B16)
    B15_cos = fs.flux_surface_average(under_average_B15_cos)
    B16_cos = fs.gamma * fs.flux_surface_average(under_average_B16_cos)

    B2mean = fs.flux_surface_average(fs.B2)

    # equation B9
    F_m = 2 / B2mean / fs.BDotNablaThetaFSA * (B15 * B16 + B15_cos * B16_cos)
    return F_m


def ωm(x: np.ndarray, m: Union[int, np.ndarray], a: IonizationState, fs: FluxSurface):
    B11 = (
        x * thermal_speed(a.T_e, a._particle) * m * fs.gamma / u.m
    )  # TODO why the u.m?
    return B11


def ν_T_ai(
    x: np.ndarray, i: int, a: IonizationState, all_species: IonizationStateCollection
):
    ai = a[i]
    prefactor = 3 * np.pi ** 0.5 / 4 * ξ(a)[i] / ai.number_density / ai.ion.mass

    def gen():
        for b in all_species:
            if b.base_particle != a.base_particle:  # TODO is not should work
                x_over_xab = x / xab_ratio(a, b).value
                part1 = (erf(x_over_xab) - 3 * Chandrasekhar_G(x_over_xab)) / x ** 3
                part2 = 4 * (
                    a.T_e / b.T_e + xab_ratio(a, b) ** -2
                )  # TODO double check this ratio
                part2full = part2 * Chandrasekhar_G(x_over_xab) / x
                result = (part1 + part2full) * effective_momentum_relaxation_rate(a, b)
                # print(f"{b=} {result=}")
                yield result

    result = prefactor * sum(gen())
    return result


def K_ps_ai(
    x: np.ndarray,
    i: int,
    a: IonizationState,
    all_species: IonizationStateCollection,
    flux_surface: FluxSurface,
    *,
    m_max=100,
    g=1  # TODO get from m_max
):
    ν = ν_T_ai(x, i, a, all_species)

    m = np.arange(1, m_max + 1)
    F = F_m(m[:, np.newaxis], flux_surface, g=g)  # TODO replace
    ω = ωm(x, m[:, np.newaxis], a, flux_surface)
    B10 = (
        1.5 * (ν / ω) ** 2
        - 9 / 2 * (ν / ω) ** 4
        + (1 / 4 + (3 / 2 + 9 / 4 * (ν / ω) ** 2) * (ν / ω) ** 2)
        * (2 * ν / ω)
        * np.arctan(ω / ν).si.value
    )
    # print(F.shape, B10.shape)
    onepart = F[:, np.newaxis] * B10
    full_sum = np.sum(onepart / ν, axis=0)
    # print(f"{full_sum=}")

    return (
        3
        / 2
        * thermal_speed(a.T_e, a.base_particle) ** 2
        * x ** 2
        * full_sum
        / u.m ** 2
    )


def K(
    x: np.ndarray,
    i: int,
    a: IonizationState,
    all_species: IonizationStateCollection,
    flux_surface: FluxSurface,
    *,
    m_max=100,
    orbit_squeezing=False,
    g=1
):
    # Eq 16
    kb = K_B_ai(x, i, a, all_species, flux_surface, orbit_squeezing=orbit_squeezing)
    # print(f"got {kb=}")
    kps = K_ps_ai(x, i, a, all_species, flux_surface, m_max=m_max, g=g)
    # print(f"got {kps=}")
    return 1 / (1 / kb + 1 / kps)


def _integrand(
    x: np.ndarray,
    α: int,
    β: int,
    i: int,
    a: IonizationState,
    all_species: IonizationStateCollection,
    flux_surface: FluxSurface,
    **kwargs
):
    laguerreterm1 = LaguerrePolynomials[α - 1](x ** 2)
    laguerreterm2 = LaguerrePolynomials[β - 1](x ** 2)
    kterm = K(x, i, a, all_species, flux_surface, **kwargs)
    xterm = x ** 4 * np.exp(-(x ** 2))
    result = laguerreterm1 * laguerreterm2 * kterm * xterm
    return result


def mu_hat(
    i: int,
    a: IonizationState,
    all_species: IonizationStateCollection,
    flux_surface: FluxSurface,
    *,
    xmin=0.0015,
    xmax=10,
    N=1000,
    **kwargs
):
    # TODO need to rework how this works... needs to be calculated for all i, earlier, otherwise - plenty of duplication
    ai = a[i]
    orders = range(1, 4)
    mu_hat_ai = u.Quantity(np.zeros((3, 3)), 1 / u.s)

    π = np.pi
    x = np.logspace(np.log10(xmin), np.log10(xmax), N)
    for α in orders:
        for β in orders:
            y = _integrand(x, α, β, i, a, all_species, flux_surface)
            integral = trapezoid(y, x)
            mu_hat_ai[α - 1, β - 1] = integral * (-1) ** (α + β)
    mass_density_probably = ai.number_density * ai.ion.mass
    actual_units = 8 / 3 / np.sqrt(π) * mu_hat_ai * mass_density_probably
    return actual_units
