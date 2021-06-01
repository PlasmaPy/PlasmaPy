r"""This module implements helper functions from [1]_.

.. [1] Houlberg et al, Bootstrap current and neoclassical transport in tokamaks of arbitrary collisionality and aspect ratio, 1997,
   Physics of Plasmas 4, 3230 (1997); , JGR, 117, A12219, doi: `10.1063/1.872465
   <https://aip.scitation.org/doi/10.1063/1.872465>`_.

.. _houlberg-1997:
This work is referenced in docstrings as `Houlberg_1997`.
"""

# TODO rename to Houlberg1997.py
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

__all__ = [
    "contributing_states",
    "mu_hat",
    "K",
    "K_ps_ai",
    "ν_T_ai",
    "ωm",
    "F_m",
    "_B17",
    "K_B_ai",
    "pitch_angle_diffusion_rate",
    "M_script",
    "N_script",
    "ξ",
    "effective_momentum_relaxation_rate",
    "N_matrix",
    "M_matrix",
    "xab_ratio",
    "ionizationstate_mass_densities",
]


def ionizationstate_mass_densities(a: IonizationState) -> u.kg / u.m ** 3:
    r"""Mass densities for each ionic level in an |IonizationState|.

    Parameters
    ----------
    a : IonizationState
        a

    Returns
    -------
    u.kg/u.m**3

    """

    return a.number_densities * u.Quantity([ai.ion.mass for ai in a])


def xab_ratio(a: IonizationState, b: IonizationState) -> u.dimensionless_unscaled:
    """Ratio of thermal speeds as defined by |Houlberg_1997|

    Parameters
    ----------
    a : IonizationState
    b : IonizationState

    Returns
    -------
    u.dimensionless_unscaled

    """

    return thermal_speed(b.T_e, b.base_particle) / thermal_speed(a.T_e, a.base_particle)


def M_matrix(species_a: IonizationState, species_b: IonizationState):
    """Test particle matrix - equation A5a through A5f from |Houlberg_1997|

    Parameters
    ----------
    species_a : IonizationState
    species_b : IonizationState
    """
    a, b = species_a, species_b
    xab = xab_ratio(a, b)
    mass_ratio = a._particle.mass / b._particle.mass
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
    """Test particle matrix - equation A6a through A6f from |Houlberg_1997|

    Parameters
    ----------
    species_a : IonizationState
    species_b : IonizationState
    """
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
    species_a: IonizationState, species_b: IonizationState
):
    """Equations A3, A4 from |Houlberg_1997|

    Parameters
    ----------
    species_a : IonizationState
    species_b : IonizationState
    """
    # TODO could be refactored using numpy
    def contributions():
        CL = lambda ai, bj: Coulomb_logarithm(
            species_b.T_e,
            species_b.n_elem,
            (ai.ion, bj.ion),  # simplifying assumption after A4
        )
        for ai in species_a:
            if ai.ion.charge == 0:
                continue
            for bj in species_b:
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
                        * thermal_speed(species_a.T_e, ai.ion) ** 3
                    )
                ).si

                yield (ai.number_density * ai.ion.mass) * collision_frequency_ai_bj

    return sum(contributions())


def ξ(isotope: IonizationState):
    """The charge state weighting factor - equation 11 from |Houlberg_1997|

    Parameters
    ----------
    isotope : IonizationState
    """
    array = u.Quantity(
        [ai.number_density * ai.ion.charge_number ** 2 for ai in isotope]
    )
    return array / array.sum()


def N_script(species_a: IonizationState, species_b: IonizationState):
    """Weighted field particle matrix - equation A2b from |Houlberg_1997|

    Parameters
    ----------
    species_a : IonizationState
    species_b : IonizationState
    """
    N = N_matrix(species_a, species_b)
    # Equation A2b
    N_script = effective_momentum_relaxation_rate(species_a, species_b) * N
    return N_script


def M_script(species_a: IonizationState, all_species: IonizationStateCollection):
    """Weighted test particle matrix - equation A2a from |Houlberg_1997|

    Parameters
    ----------
    species_a : IonizationState
    all_species : IonizationStateCollection
    """
    # Equation A2a
    def gener():
        for species_b in all_species:
            if species_b is not species_a:
                # TODO am I sure this if is necessary here?
                yield M_matrix(
                    species_a, species_b
                ) * effective_momentum_relaxation_rate(species_a, species_b)

    return sum(gener())


def pitch_angle_diffusion_rate(
    x: np.ndarray,
    a: IonizationState,
    all_species: IonizationStateCollection,
):
    """The pitch angle diffusion rate, equation B4b from |Houlberg_1997|

    Parameters
    ----------
    x : np.ndarray
        distribution velocity relative to the thermal velocity.
    a : IonizationState
    all_species : IonizationStateCollection
    """
    # Houlberg_1997, equation B4b,
    xi = ξ(a)
    denominator = x ** 3

    def sum_items():
        """sum_items."""
        for b in all_species:
            xab = xab_ratio(a, b)
            x_over_xab = (x / xab).value
            numerator = erf(x_over_xab) - Chandrasekhar_G(x_over_xab)
            fraction = numerator / denominator
            result = fraction * effective_momentum_relaxation_rate(a, b)
            yield result

    mass_density = ionizationstate_mass_densities(a)
    result = (
        xi[:, np.newaxis]
        / mass_density[:, np.newaxis]
        * 3
        * np.sqrt(np.pi)
        / 4
        * sum(sum_items())[np.newaxis, :]
    )
    return result


def K_B_ai(
    x: np.ndarray,
    a_states: IonizationState,
    all_species: IonizationStateCollection,
    flux_surface: FluxSurface,
    *,
    orbit_squeezing: bool = False,
):
    """Banana regime contribution to effective viscosity - eq. B1 from |Houlberg_1997|

    Parameters
    ----------
    x : np.ndarray
        x
    a_states : IonizationState
        a_states
    all_species : IonizationStateCollection
        all_species
    flux_surface : FluxSurface
        flux_surface
    orbit_squeezing : bool (default `False`)
        orbit_squeezing
    """
    # eq. B1-B4, Houlberg_1997
    f_t = flux_surface.trapped_fraction()
    f_c = 1 - f_t
    if orbit_squeezing:
        raise NotImplementedError(
            "TODO allow for non-zero, changing radial electric fields (orbit squeezing)"
        )
    else:
        S_ai = 1  # Equation B2
    padr = pitch_angle_diffusion_rate(x, a_states, all_species)
    return padr * f_t / f_c / S_ai ** 1.5


LaguerrePolynomials = [
    lambda x: np.ones_like(x),
    lambda x: 5 / 2 - x,
    lambda x: 35 / 8 - 7 / 2 * x + 1 / 2 * x ** 2,
]


def _B17(flux_surface):
    """Equation B17 from |Houlberg_1997|. Likely bugged!

    Notes
    -----
    Eventually this should allow picking the right `m` in `K` below.

    Parameters
    ----------
    flux_surface :
        flux_surface
    """
    fs = flux_surface
    B20 = fs.Brvals * fs.Bprimervals + fs.Bzvals * fs.Bprimezvals
    under_average_B17 = (B20 / fs.Bmag) ** 2
    return fs.flux_surface_average(under_average_B17) / fs.flux_surface_average(fs.B2)


def F_m(m: Union[int, np.ndarray], flux_surface: FluxSurface):
    """Mode weights for the Pfirsch-Schlüter contribution.

    Equation B9 from |Houlberg_1997|.

    Parameters
    ----------
    m : Union[int, np.ndarray]
        m
    flux_surface : FluxSurface
        flux_surface
    """
    fs = flux_surface
    B20 = fs.Brvals * fs.Bprimervals + fs.Bzvals * fs.Bprimezvals
    under_average_B16 = np.sin(m * fs.Theta) * B20
    under_average_B15 = under_average_B16 / fs.Bmag
    under_average_B16_cos = np.cos(m * fs.Theta) * B20
    under_average_B15_cos = under_average_B16_cos / fs.Bmag
    B15 = fs.flux_surface_average(under_average_B15)
    B16 = fs.gamma * fs.flux_surface_average(under_average_B16)
    B15_cos = fs.flux_surface_average(under_average_B15_cos)
    B16_cos = fs.gamma * fs.flux_surface_average(under_average_B16_cos)

    B2mean = fs.flux_surface_average(fs.B2)

    F_m = 2 / B2mean / fs.BDotNablaThetaFSA * (B15 * B16 + B15_cos * B16_cos)
    return F_m


def ωm(x: np.ndarray, m: Union[int, np.ndarray], a: IonizationState, fs: FluxSurface):
    """Equation B11 from |Houlberg_1997|.

    Parameters
    ----------
    x : np.ndarray
        distribution velocity relative to the thermal velocity.
    m : Union[int, np.ndarray]
    a : IonizationState
    fs : FluxSurface
    """
    """Equation B11 of Houlberg_1997."""
    B11 = (
        x * thermal_speed(a.T_e, a._particle) * m * fs.gamma / u.m
    )  # TODO why the u.m?
    return B11


def ν_T_ai(x: np.ndarray, a: IonizationState, all_species: IonizationStateCollection):
    """Characteristic rate of anisotropy relaxation.

    Equation B12 from |Houlberg_1997|.

    Parameters
    ----------
    x : np.ndarray
        distribution velocity relative to the thermal velocity.
    a : IonizationState
        a
    all_species : IonizationStateCollection
        all_species
    """
    mass_density = ionizationstate_mass_densities(a)
    prefactor = 3 * np.pi ** 0.5 / 4 * ξ(a) / mass_density

    def gen():
        """gen."""
        for b in all_species:
            if b.base_particle != a.base_particle:  # TODO is not should work
                x_over_xab = (x / xab_ratio(a, b)).value
                part1 = (erf(x_over_xab) - 3 * Chandrasekhar_G(x_over_xab)) / x ** 3
                part2 = 4 * (
                    a.T_e / b.T_e + xab_ratio(a, b) ** -2
                )  # TODO double check this ratio
                part2full = part2 * Chandrasekhar_G(x_over_xab) / x
                result = (part1 + part2full) * effective_momentum_relaxation_rate(a, b)
                yield result

    result = prefactor[:, np.newaxis] * sum(gen())[np.newaxis, :]
    return result


def K_ps_ai(
    x: np.ndarray,
    a: IonizationState,
    all_species: IonizationStateCollection,
    flux_surface: FluxSurface,
    *,
    m_max=100,
):
    """Pfirsch-Schlüter regime contribution to effective viscosity - eq. B8 from |Houlberg_1997|

    Parameters
    ----------
    x : np.ndarray
        x
    a : IonizationState
        a
    all_species : IonizationStateCollection
        all_species
    flux_surface : FluxSurface
        flux_surface
    m_max :
        m_max
    """
    ν = ν_T_ai(x, a, all_species)[:, np.newaxis, :]

    m = np.arange(1, m_max + 1)
    F = F_m(m[:, np.newaxis], flux_surface)
    ω = ωm(x, m[:, np.newaxis], a, flux_surface)[np.newaxis, ...]
    B10 = (
        1.5 * (ν / ω) ** 2
        - 9 / 2 * (ν / ω) ** 4
        + (1 / 4 + (3 / 2 + 9 / 4 * (ν / ω) ** 2) * (ν / ω) ** 2)
        * (2 * ν / ω)
        * np.arctan(ω / ν).si.value
    )
    onepart = F[:, np.newaxis] * B10
    full_sum = np.sum(onepart / ν, axis=1)

    return (
        3
        / 2
        * thermal_speed(a.T_e, a.base_particle)
        ** 2  # TODO replace T_e with T_i in the fullness of time
        * x ** 2
        * full_sum
        / u.m ** 2
    )


def K(
    x: np.ndarray,
    a: IonizationState,
    all_species: IonizationStateCollection,
    flux_surface: FluxSurface,
    *,
    m_max: int = 100,
    orbit_squeezing: bool = False,
):
    """Total effective velocity-dependent viscosity with contributions from -
    eq. 16 from |Houlberg_1997|

    Notes
    -----
    This expression originally comes from
    K. C. Shaing, C. T. Hsu, M. Yokoyama, and M. Wakatani, Phys. Plasmas
    2, 349 (1995); and K. C. Shaing, M. Yokoyama, M. Wakatani, and C. T. Hsu,
    ibid. 3, 965 (1996).

    Parameters
    ----------
    x : np.ndarray
        distribution velocity relative to the thermal velocity.
    a : IonizationState
    all_species : IonizationStateCollection
    flux_surface : FluxSurface
    m_max : int
    orbit_squeezing : bool
        orbit_squeezing
    """
    # Eq 16
    kb = K_B_ai(x, a, all_species, flux_surface, orbit_squeezing=orbit_squeezing)
    # print(f"got {kb=}")
    kps = K_ps_ai(x, a, all_species, flux_surface, m_max=m_max)
    # print(f"got {kps=}")
    return 1 / (1 / kb + 1 / kps)


def mu_hat(
    a: IonizationState,
    all_species: IonizationStateCollection,
    flux_surface: FluxSurface,
    *,
    xmin: float = 0.0015,
    xmax: float = 10,
    N: int = 1000,
    **kwargs,
):
    """Viscosity coefficients - equation 15 from |Houlberg_1997|.

    Parameters
    ----------
    a : IonizationState
    all_species : IonizationStateCollection
    flux_surface : FluxSurface
    xmin :
        xmin
    xmax :
        xmax
    N :
        N
    kwargs :
        kwargs
    """
    if N is None:
        N = 1000
    orders = np.arange(1, 4)
    π = np.pi
    x = np.logspace(np.log10(xmin), np.log10(xmax), N)

    α = orders
    β = orders
    len_a = len(a.number_densities)
    signs = (-1) * (α[:, None] + β[None, :])
    laguerres = np.vstack([LaguerrePolynomials[o - 1](x ** 2) for o in orders])
    kterm = K(x, a, all_species, flux_surface, **kwargs)
    kterm = kterm.reshape(len_a, N, 1, 1)
    xterm = (x ** 4 * np.exp(-(x ** 2))).reshape(1, N, 1, 1)
    y = laguerres.reshape(1, N, 3, 1) * laguerres.reshape(1, N, 1, 3) * kterm * xterm
    integral = trapezoid(y, x, axis=1)
    mu_hat_ai = integral * signs[None, ...]
    mass_density = ionizationstate_mass_densities(a)
    actual_units = (8 / 3 / np.sqrt(π)) * mu_hat_ai * mass_density[:, None, None]
    return actual_units


def contributing_states(a: IonizationState):
    """Helper iterator over |IonizationState|s.

    Return a generator of tuples (`ξ(a)[i], a[i]`).

    Parameters
    ----------
    a : IonizationState

    """
    xi = ξ(a)
    for i, ai in enumerate(a):
        if xi[i] == 0:
            continue
        yield xi[i], ai
