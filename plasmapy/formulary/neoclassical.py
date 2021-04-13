from __future__ import annotations

__all__ = [
    "get_flows",
]


import numpy as np

from astropy import constants
from astropy import units as u
from scipy import integrate
from scipy.special import erf

from plasmapy.formulary import thermal_speed
from plasmapy.formulary.mathematics import Chandrasekhar_G

try:
    from scipy.integrate import trapz as trapezoid
except ImportError:
    from scipy.integrate import trapezoid


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


def ξ(isotope):
    array = u.Quantity(
        [ai.number_density * ai.ion.integer_charge ** 2 for ai in isotope]
    )
    return array / array.sum()


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
                # TODO am I sure this if is necessary here?
                yield M_matrix(
                    species_a, species_b
                ) * effective_momentum_relaxation_rate(species_a, species_b)

    return sum(gener())


def pitch_angle_diffusion_rate(x, index, a_states, all_species):
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
    padr = pitch_angle_diffusion_rate(x, index, a_states, all_species)
    return padr * f_t / f_c / S_ai ** 1.5


LaguerrePolynomials = [
    lambda x: 1,
    lambda x: 5 / 2 - x,
    lambda x: 35 / 8 - 7 / 2 * x + 1 / 2 * x ** 2,
]


def rbar(a, all_species, flux_surface, beta_coeffs=None) -> u.Quantity:
    if beta_coeffs is not None:
        raise NotImplementedError
    else:
        beta_cx = np.zeros(3)  # TODO
        beta_an = np.zeros(3)  # TODO
        beta_coeffs = np.diag(beta_cx + beta_an) * u.kg / u.m ** 3 / u.s

    def gen():
        for i, ai in enumerate(a):
            if ξ(a)[i] == 0:
                continue  # won't add anything to sum anyway, and matrix gets singular
            Aai = (
                ξ(a)[i] * M_script(a, all_species)
                - mu_hat(i, a, all_species, flux_surface)
                - beta_coeffs
            )
            S_matrix = ξ(a)[i] * np.eye(3)
            rai_as_rows = np.linalg.solve(Aai, S_matrix)
            # TODO does not include r_pT, r_E, r_NBI yet. Should it?
            rbar_ingredient = ξ(a)[i] * rai_as_rows
            yield rbar_ingredient

    return sum(gen())


def eq34matrix(all_species, flux_surface, beta_coeffs=None):
    output_matrix = u.Quantity(np.eye(3 * len(all_species)))

    for I, a in enumerate(all_species):
        i = 3 * I
        # this is probably how rbar should work!
        # original_rhs = np.concatenate([rbar(a, all_species, fs) for a in all_species])

        rarray = rbar(a, all_species, flux_surface, beta_coeffs)
        for J, b in enumerate(all_species):
            j = 3 * J
            narray = N_script(a, b).sum(axis=0, keepdims=True)
            result = narray * rarray.T
            output_matrix[i : i + 3, j : j + 3] += result

    return output_matrix


def F_m(m, flux_surface, g=1):
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

    jacobian = g ** 0.5
    BdotNablatheta = fs.BDotNablaThetaFSA
    B2mean = fs.flux_surface_average(fs.B2)

    # equation B9
    F_m = 2 / B2mean / BdotNablatheta * (B15 * B16 + B15_cos * B16_cos)
    return F_m


def ωm(x, m, a: isotopelike, fs):
    B11 = (
        x * thermal_speed(a.T_e, a._particle) * m * fs.gamma / u.m
    )  # TODO why the u.m?
    return B11


def ν_T_ai(x, i, a, all_species):
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
    x, i, a, all_species, flux_surface, *, m_max=100, g=1  # TODO get from m_max
):
    ai = a[i]
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


def K(x, i, a, all_species, flux_surface, *, m_max=100, orbit_squeezing=False, g=1):
    # Eq 16
    kb = K_B_ai(x, i, a, all_species, flux_surface)
    # print(f"got {kb=}")
    kps = K_ps_ai(x, i, a, all_species, flux_surface, m_max=m_max, g=g)
    # print(f"got {kps=}")
    return 1 / (1 / kb + 1 / kps)


def _integrand(x, α, β, i, a, all_species, flux_surface, **kwargs):
    laguerreterm1 = LaguerrePolynomials[α - 1](x ** 2)
    laguerreterm2 = LaguerrePolynomials[β - 1](x ** 2)
    kterm = K(x, i, a, all_species, flux_surface, **kwargs)
    xterm = x ** 4 * np.exp(-(x ** 2))
    result = laguerreterm1 * laguerreterm2 * kterm * xterm
    return result


def mu_hat(i, a, all_species, flux_surface, *, xmin=0.0015, xmax=10, N=1000, **kwargs):
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


if __name__ == "__main__":
    from plasmapy.particles import (
        IonizationState,
        IonizationStateCollection,
        Particle,
        proton,
    )

    all_species = IonizationStateCollection(
        {
            "H": [0, 1],
            #      "D": [0, 1],   raises ParticleError, why?
            "C": [0, 1 / 1.1, 0.1 / 1.1, 0, 0, 0, 0],
        },
        n0=1e20 * u.m ** -3,
        abundances={"H": 1, "C": 0.11},
        T_e=10 * u.eV,
    )
    hydrogen = all_species["H"]
    ν_T_ai(1.0495932305582267e-05, 1, hydrogen, all_species)


def rbar_sources(a, all_species, flux_surface, beta_coeffs=None) -> u.Quantity:
    fs = flux_surface
    if beta_coeffs is not None:
        # TODO should be a dict or sth
        raise NotImplementedError
    else:
        beta_cx = np.zeros(3)  # TODO
        beta_an = np.zeros(3)  # TODO
        beta_coeffs = np.diag(beta_cx + beta_an) * u.kg / u.m ** 3 / u.s
    num_charge_states = len(a.integer_charges)
    temperatures = u.Quantity([a.T_e] * num_charge_states)

    # TODO make these inputs
    density_gradient = u.Quantity(num_charge_states * [1e18 * u.m ** -3 / u.m])
    temperature_gradient = u.Quantity(num_charge_states * [10 * u.K / u.m])

    # but not this, this is dervied
    pressure_gradient = constants.k_B * (
        temperatures * density_gradient + a.number_densities * temperature_gradient
    )

    def gen():
        for i, ai in enumerate(a):
            if ξ(a)[i] == 0:
                continue  # won't add anything to sum anyway, and matrix gets singular
            Aai = (
                ξ(a)[i] * M_script(a, all_species)
                - mu_hat(i, a, all_species, flux_surface)
                - beta_coeffs
            )
            μ = mu_hat(i, a, all_species, fs)
            Spt = (
                fs.Fhat
                / ai.ion.charge
                / ai.number_density
                * u.Quantity(
                    [
                        pressure_gradient[i] * μ[0, 0]
                        + a.number_densities[i]
                        * constants.k_B
                        * temperature_gradient[i]
                        * μ[0, 1],
                        pressure_gradient[i] * μ[1, 0]
                        + a.number_densities[i]
                        * constants.k_B
                        * temperature_gradient[i]
                        * μ[1, 1],
                    ]
                )
            ).si
            Spt = np.append(Spt, 0)
            S_matrix = Spt
            rai_as_rows = np.linalg.solve(Aai, S_matrix)
            # TODO does not include r_pT, r_E, r_NBI yet
            rbar_ingredient = ξ(a)[i] * rai_as_rows
            yield rbar_ingredient

    return sum(gen())


def get_flows(
    # TODO make these inputs
    all_species,
    flux_surface,
    density_gradient,
    temperature_gradient,
    beta_coeffs=None,
):
    fs = flux_surface
    # rbar_sources ABSOLUTELY needs a bloody rework to do stuff simultaneously
    rhs = np.concatenate(
        [rbar_sources(a, all_species, fs, beta_coeffs=beta_coeffs) for a in all_species]
    ).si
    if beta_coeffs is not None:
        # TODO should be a dict or sth
        raise NotImplementedError
    else:
        beta_cx = np.zeros(3)  # TODO
        beta_an = np.zeros(3)  # TODO
        beta_coeffs = np.diag(beta_cx + beta_an) * u.kg / u.m ** 3 / u.s

    lhs = eq34matrix(all_species, fs)
    ubar = np.linalg.solve(lhs, rhs)

    outputs = {}
    for I, a in enumerate(all_species):
        # use Eq31 to get charge state flows from isotopic flows
        def gen():
            i = 3 * I
            for J, b in enumerate(all_species):
                j = 3 * J
                ubar_b = ubar[j : j + 3]
                yield (N_script(a, b) * ubar_b.reshape(1, -1)).sum(axis=1)

        Λ = -sum(gen())
        M = M_script(a, all_species)
        xi = ξ(a)
        # TODO need to rework how mu_hat works... needs to be calculated for all i, earlier, otherwise - plenty of duplication
        for i, ai in enumerate(a):
            if i == 0:
                continue
            μ = mu_hat(i, a, all_species, fs)
            Aai = xi[i] * M - μ - beta_coeffs
            S_ai = xi[i] * np.diag(Λ)  # TODO np.diag(Δ)?
            rai_as_rows = np.linalg.solve(Aai, S_ai)
            order_flow_sum = (
                (Λ.reshape(-1, 1) * rai_as_rows).sum(axis=0).si.value
            )  # TODO fix units

            temperature = ai.T_i
            density = a.number_densities[i]  # TODO cleanup
            ne_grad = density_gradient.get(ai.ionic_symbol, 0 * u.m ** -4)
            T_grad = temperature_gradient.get(ai.ionic_symbol, 0 * u.K / u.m)
            pressure_gradient = constants.k_B * (
                temperature * ne_grad + density * T_grad
            )
            Spt = (
                fs.Fhat
                / ai.ion.charge
                / ai.number_density
                * u.Quantity(
                    [
                        pressure_gradient * μ[0, 0]
                        + density * constants.k_B * T_grad * μ[0, 1],
                        pressure_gradient * μ[1, 0]
                        + density * constants.k_B * T_grad * μ[1, 1],
                    ]
                )
            ).si
            Spt = np.append(Spt, 0)
            rpt_row = np.linalg.solve(
                Aai, Spt
            ).si.value  # TODO units are wrong here too; but I think the mechanics should just about work
            flows = order_flow_sum + rpt_row  # Eq31
            outputs[ai.ionic_symbol] = flows
    return outputs
