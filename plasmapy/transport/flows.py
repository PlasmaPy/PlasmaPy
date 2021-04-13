from __future__ import annotations

__all__ = [
    "get_flows",
    "get_fluxes",
    "linear_transport_coefficients",
]


import numpy as np

from astropy import constants
from astropy import units as u

from plasmapy.formulary.neoclassical import M_script, mu_hat, N_script, ξ

try:
    from functools import cached_property
except ImportError:
    from cached_property import cached_property


def S_pt(ai, μ, fs, density_gradient, temperature_gradient):
    # TODO gradients should be attached to ai
    ne_grad = density_gradient.get(ai.ionic_symbol, 0 * u.m ** -4)
    T_grad = temperature_gradient.get(ai.ionic_symbol, 0 * u.K / u.m)
    pressure_gradient = constants.k_B * (ai.T_i * ne_grad + ai.number_density * T_grad)
    Spt = (
        fs.Fhat
        / ai.ion.charge
        / ai.number_density
        * u.Quantity(
            [
                pressure_gradient * μ[0, 0]
                + ai.number_density * constants.k_B * T_grad * μ[0, 1],
                pressure_gradient * μ[1, 0]
                + ai.number_density * constants.k_B * T_grad * μ[1, 1],
            ]
        )
    ).si
    Spt = np.append(Spt, 0)
    return Spt


class FlowCalculator:
    def __init__(
        self, all_species, flux_surface,
    ):
        self.all_species = all_species
        self.flux_surface = flux_surface

    @cached_property
    def μ(self):  # this would be better as a cached method... depending on ai
        results = {}
        for a in self.all_species:
            xi = ξ(a)
            for i, ai in enumerate(a):
                if i == 0 or xi[i] == 0:
                    continue
                results[ai.ionic_symbol] = mu_hat(
                    i, a, self.all_species, self.flux_surface
                )
        return results

    def rbar(self, a, beta_coeffs=None) -> u.Quantity:
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
                μ = self.μ[ai.ionic_symbol]
                Aai = ξ(a)[i] * M_script(a, self.all_species) - μ - beta_coeffs
                S_matrix = ξ(a)[i] * np.eye(3)
                rai_as_rows = np.linalg.solve(Aai, S_matrix)
                # TODO does not include r_pT, r_E, r_NBI yet. Should it?
                rbar_ingredient = ξ(a)[i] * rai_as_rows
                yield rbar_ingredient

        return sum(gen())

    def rbar_sources(
        self, density_gradient, temperature_gradient, beta_coeffs=None
    ) -> u.Quantity:
        fs = self.flux_surface
        if beta_coeffs is not None:
            # TODO should be a dict or sth
            raise NotImplementedError
        else:
            beta_cx = np.zeros(3)  # TODO
            beta_an = np.zeros(3)  # TODO
            beta_coeffs = np.diag(beta_cx + beta_an) * u.kg / u.m ** 3 / u.s

        results = []
        for a in self.all_species:

            def gen():
                for i, ai in enumerate(a):
                    if ξ(a)[i] == 0:
                        continue  # won't add anything to sum anyway, and matrix gets singular
                    μ = self.μ[ai.ionic_symbol]
                    Aai = ξ(a)[i] * M_script(a, self.all_species) - μ - beta_coeffs
                    Spt = S_pt(ai, μ, fs, density_gradient, temperature_gradient)
                    rai_as_rows = np.linalg.solve(Aai, Spt)
                    # TODO does not include r_pT, r_E, r_NBI yet
                    rbar_ingredient = ξ(a)[i] * rai_as_rows
                    yield rbar_ingredient

            results.append(sum(gen()))
        return np.concatenate(results).si

    def get_flows(
        self,
        density_gradient,
        temperature_gradient,
        # TBH could probably pass profile shapes here, instead...
        beta_coeffs=None,
    ):
        fs = self.flux_surface
        rhs = self.rbar_sources(
            density_gradient, temperature_gradient, beta_coeffs=beta_coeffs
        )
        if beta_coeffs is not None:
            # TODO should be a dict or sth
            raise NotImplementedError
        else:
            beta_cx = np.zeros(3)  # TODO
            beta_an = np.zeros(3)  # TODO
            beta_coeffs = np.diag(beta_cx + beta_an) * u.kg / u.m ** 3 / u.s

        lhs = self.eq34matrix()
        ubar = np.linalg.solve(lhs, rhs)

        outputs = {}
        for I, a in enumerate(self.all_species):
            # use Eq31 to get charge state flows from isotopic flows
            def gen():
                i = 3 * I
                for J, b in enumerate(self.all_species):
                    j = 3 * J
                    ubar_b = ubar[j : j + 3]
                    yield (N_script(a, b) * ubar_b.reshape(1, -1)).sum(axis=1)

            Λ = -sum(gen())
            M = M_script(a, self.all_species)
            xi = ξ(a)
            for i, ai in enumerate(a):
                if i == 0 or xi[i] == 0:
                    continue
                μ = self.μ[ai.ionic_symbol]
                Aai = xi[i] * M - μ - beta_coeffs
                S_ai = xi[i] * np.diag(Λ)
                rai_as_rows = np.linalg.solve(Aai, S_ai)
                order_flow_sum = (
                    (Λ.reshape(-1, 1) * rai_as_rows).sum(axis=0).si.value
                )  # TODO fix units

                Spt = S_pt(ai, μ, fs, density_gradient, temperature_gradient)
                rpt_row = np.linalg.solve(
                    Aai, Spt
                ).si.value  # TODO units are wrong here too; but I think the mechanics should just about work
                flows = order_flow_sum + rpt_row  # Eq31
                outputs[ai.ionic_symbol] = flows
        return outputs

    def eq34matrix(self, beta_coeffs=None):
        output_matrix = u.Quantity(np.eye(3 * len(self.all_species)))

        for I, a in enumerate(self.all_species):
            i = 3 * I
            rarray = self.rbar(a, beta_coeffs)
            for J, b in enumerate(self.all_species):
                j = 3 * J
                narray = N_script(a, b).sum(axis=0, keepdims=True)
                result = narray * rarray.T
                output_matrix[i : i + 3, j : j + 3] += result

        return output_matrix
