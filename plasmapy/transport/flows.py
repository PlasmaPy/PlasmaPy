from __future__ import annotations

__all__ = [
    "FlowCalculator",
]


import numpy as np

from astropy import constants
from astropy import units as u
from collections import defaultdict, namedtuple

from .neoclassical import M_script, mu_hat, N_script, ξ

try:
    from functools import cached_property
except ImportError:
    from cached_property import cached_property

Fluxes = namedtuple("Fluxes", ["particle_flux", "heat_flux"])



class FlowCalculator:
    """
    This does, in fact, do most things for my thesis.
    """

    def __init__(
        self,
        all_species,
        flux_surface,
        density_gradient,
        temperature_gradient,
    ):
        self.all_species = all_species
        self.flux_surface = fs = flux_surface
        self.density_gradient = density_gradient
        self.temperature_gradient = temperature_gradient
        self.M_script_matrices = {}
        self.N_script_matrices = {}

        self.S_pt = {}
        self.μ = {}
        self.Aai = {}
        self.thermodynamic_forces = {}
        self.pressure_gradient = {}
        self.ξ = {}

        r_flows_list = []
        r_sources_list = []
        rbar_flows_list  = []
        rbar_sources_list = []
        self.r_pt = {}
        S_pt_list = []
        for a in self.all_species:
            sym = a.base_particle
            charges = a.integer_charges * constants.e.si
            n_charge_states = len(a.integer_charges)
            xi = ξ(a)
            T_i = u.Quantity([ai.T_i for ai in a]) # TODO should be an attribute of a
            n_i = a.number_densities
            density_gradient = self.density_gradient.get(sym, np.zeros(n_charge_states) * (u.m**-4))
            temperature_gradient = self.temperature_gradient.get(sym, np.zeros(n_charge_states) * (u.K / u.m))
            pressure_gradient_over_n_i = constants.k_B * (T_i * density_gradient / n_i + temperature_gradient)
            # we divide by n_i, which can be zero, leadning to inf, so to correct that...
            pressure_gradient_over_n_i[np.isinf(pressure_gradient_over_n_i)] = 0
            μ = u.Quantity([mu_hat(i, a, self.all_species, self.flux_surface) for i in a.integer_charges]) # TODO rework to work on arrays
            
            Aai = xi[:, np.newaxis, np.newaxis] * self.M_script(a)[np.newaxis, ...] - μ
            # --- TD forces eq21
            thermodynamic_forces = (fs.Fhat / charges * u.Quantity([pressure_gradient_over_n_i,
                                                                    constants.k_B * temperature_gradient,
                                                                    np.zeros(n_charge_states) * (u.J / u.m)])).T
            CHARGE_STATE_AXIS = 0
            BETA_AXIS = 2
            S_pt = (thermodynamic_forces[:, np.newaxis, :] * μ).sum(axis=BETA_AXIS)
            S_pt_list.append(S_pt)
            r_pt = np.linalg.solve(Aai, S_pt)
            # TODO r_E itd
            r_sources = r_pt + 0
            r_sources_list.append(r_sources)
            rbar_sources = (xi[:, np.newaxis] * r_sources).nansum(axis=CHARGE_STATE_AXIS)
            rbar_sources_list.append(rbar_sources)
            S_flows = (xi[:, np.newaxis, np.newaxis] * np.eye(3)) * u.Unit("N T / m3")
            r_flows = np.linalg.solve(Aai, S_flows)
            r_flows_list.append(r_flows)
            rbar_flows = (xi[:, np.newaxis, np.newaxis] * r_flows).nansum(axis=CHARGE_STATE_AXIS)
            rbar_flows_list.append(rbar_flows)
            for i, ai in enumerate(a):
                sym = ai.ionic_symbol
                self.density_gradient[sym] = density_gradient[i]
                self.temperature_gradient[sym] = temperature_gradient[i]
                self.r_pt[sym] = r_pt[i]
                self.S_pt[sym] = S_pt[i]
                self.thermodynamic_forces[sym] = thermodynamic_forces[i]




        lhs = u.Quantity(np.eye(3 * len(all_species)), "J2 / (A m6)")  # TODO verify
        for i, a in enumerate(all_species):
            rarray = rbar_flows_list[i]
            for j, b in enumerate(self.all_species):
                narray = self.N_script(a, b).sum(axis=0, keepdims=True)
                result = narray * rarray.T
                lhs[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] += result
        rhs = u.Quantity(rbar_sources_list)
        ubar = np.linalg.solve(lhs, rhs.ravel())

        self._charge_state_flows = {}
        for r_flows, r_sources, a in zip(r_flows_list,
                                         r_sources_list,
                                         self.all_species):
            def gen():
                for j, b in enumerate(self.all_species):
                    ubar_b = ubar[3 * j : 3 * j + 3]
                    yield (self.N_script(a, b) * ubar_b.reshape(1, -1)).sum(axis=1)
            Λ = -sum(gen())

            # TODO units are wrong between the two here
            self_consistent_u = np.sum(Λ[np.newaxis, :, np.newaxis] * r_flows, axis=2)
            u_velocity = self_consistent_u + r_sources

            for i, ai in enumerate(a):
                if np.isfinite(u_velocity[i]).all():
                    self._charge_state_flows[ai.ionic_symbol] = u_velocity[i]

    @staticmethod
    def contributing_states(a):
        xi = ξ(a)
        for i, ai in enumerate(a):
            if xi[i] == 0:
                continue
            yield xi[i], ai

    def M_script(self, a):
        sym = a.base_particle
        if sym not in self.M_script_matrices:
            self.M_script_matrices[sym] = M_script(a, self.all_species)
        return self.M_script_matrices[sym]

    def N_script(self, a, b):
        sym_tuple = a.base_particle, b.base_particle
        if sym_tuple not in self.N_script_matrices:
            self.N_script_matrices[sym_tuple] = N_script(a, b)
        return self.N_script_matrices[sym_tuple]

    def funnymatrix(self, a_symbol):
        a = self.all_species[a_symbol]
        M = self.M_script(a)
        outputs = {}
        for _, ai in self.contributing_states(a):
            sym = ai.ionic_symbol
            S = self.S_pt[sym]
            output = S * M
            for b in self.all_species:
                N = N_script(a, b)
                for xj, bj in self.contributing_states(b):
                    output += xj * N * self.S_pt[bj.ionic_symbol]
            outputs[sym] = output
        return outputs

    @cached_property
    def _fluxes_BP(self):
        Fhat = self.flux_surface.Fhat
        fs = self.flux_surface
        B2fsav = fs.flux_surface_average(fs.B2) * u.T ** 2  # flux surface averaged B^2
        results = {}
        for a in self.all_species:
            for (
                ai
            ) in (
                a
            ):  # this could be rfactored out by iterating over self._charge_state_flows, instead, given a way to access ionizationstate back from ioniclevel
                sym = ai.ionic_symbol
                if sym not in self._charge_state_flows:
                    continue
                u_velocity = self._charge_state_flows[sym]

                u_θ = (u_velocity + self.thermodynamic_forces[sym]) / B2fsav
                μ = self.μ[sym]
                Γ_BP = -(Fhat / ai.ion.charge * (μ[0, :] * u_θ).sum()).si
                # TODO verify unit; does not look right
                q_BP = -(
                    fs.Fhat
                    * constants.k_B
                    * ai.T_i
                    / ai.ion.charge
                    * (μ[1, :] * u_θ).sum()
                ).si
                results[sym] = Fluxes(Γ_BP.si, q_BP.si)
        return results

    @cached_property
    def _fluxes_PS(self):
        fs = self.flux_surface
        B2fsav = fs.flux_surface_average(fs.B2) * u.T ** 2  # flux surface averaged B^2
        Binv2fsav = fs.flux_surface_average(1 / fs.B2) / u.T ** 2
        results = {}
        fs = self.flux_surface
        for a in self.all_species:
            silly = self.funnymatrix(a.base_particle)
            for xi, ai in self.contributing_states(a):
                sym = ai.ionic_symbol
                prefactor = (
                    -fs.Fhat / ai.ion.charge * xi / B2fsav * (1 - B2fsav * Binv2fsav)
                )
                Γ_PS = prefactor * silly[sym][0].sum()
                q_PS = prefactor * constants.k_B * ai.T_i * silly[sym][1].sum()
                results[sym] = Fluxes(Γ_PS.si, q_PS.si)
        return results

    @cached_property
    def _fluxes_CL(self):
        fs = self.flux_surface
        FSA = fs.flux_surface_average(fs.GradRho2 / fs.B2) / u.T ** 2 / u.m
        # TODO if FSA does not drop units, the above line is completely wrong
        Fhat = self.flux_surface.Fhat
        results = {}
        for a in self.all_species:
            silly = self.funnymatrix(a.base_particle)
            for xi, ai in self.contributing_states(a):
                sym = ai.ionic_symbol
                prefactor = FSA / Fhat * xi / ai.ion.charge
                Γ_CL = prefactor * silly[sym][0].sum()
                q_CL = prefactor * constants.k_B * ai.T_i * silly[sym][1].sum()
                results[sym] = Fluxes(Γ_CL.si, q_CL.si)
        return results

    @cached_property
    def fluxes(self):
        results = {}
        for a in self.all_species:
            for _, ai in self.contributing_states(a):
                sym = ai.ionic_symbol
                Γ_BP, q_BP = self._fluxes_BP[sym]
                Γ_PS, q_PS = self._fluxes_PS[sym]
                Γ_CL, q_CL = self._fluxes_CL[sym]
                results[sym] = Fluxes(Γ_BP + Γ_PS + Γ_CL, q_BP + q_PS + q_CL)
        return results

    @cached_property
    def diffusion_coefficient(self):
        results = {}
        for a in self.all_species:
            for _, ai in self.contributing_states(a):
                sym = ai.ionic_symbol
                flux = self.fluxes[sym].particle_flux
                results[sym] = (
                    -flux / self.density_gradient[sym]
                )  # Eq48 TODO this is a partial adaptation
        return results

    @cached_property
    def thermal_conductivity(self):
        results = {}
        for a in self.all_species:
            for _, ai in self.contributing_states(a):
                sym = ai.ionic_symbol
                flux = self.fluxes[ai].heat_flux
                results[sym] = -flux / self.temperature_gradient[sym]
        return results

    @cached_property
    def bootstrap_current(self):
        def gen():
            for a in self.all_species:
                for _, ai in self.contributing_states(a):
                    sym = ai.ionic_symbol
                    yield ai.ion.charge * ai.number_density * self.r_pt[sym][
                        0
                    ]  # eq 37, second term

        return sum(gen())
