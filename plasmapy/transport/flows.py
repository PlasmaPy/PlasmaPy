from __future__ import annotations

__all__ = [
    "FlowCalculator",
]


import numpy as np
import typing
import xarray

from astropy import constants
from astropy import units as u
from collections import defaultdict, namedtuple

from plasmapy.particles import IonizationStateCollection, Particle
from plasmapy.plasma.fluxsurface import FluxSurface

from .Houlberg1997 import ExtendedParticleList

try:
    from functools import cached_property
except ImportError:
    from cached_property import cached_property
particle_flux_unit = u.m ** -2 / u.s
heat_flux_unit = u.J * particle_flux_unit
from plasmapy.particles.exceptions import InvalidElementError
from plasmapy.utils.decorators import validate_quantities

Fluxes = namedtuple("Fluxes", ["particle_flux", "heat_flux"])


class FlowCalculator:
    r"""Interface to neoclassical transport calculations.

    Reimplements NCLASS from |Houlberg_1997|[1]_

    References
    ----------
    .. [1] Houlberg et al, Bootstrap current and neoclassical transport in tokamaks of arbitrary collisionality and aspect ratio, 1997,
       Physics of Plasmas 4, 3230 (1997); , JGR, 117, A12219, doi: `10.1063/1.872465
       <https://aip.scitation.org/doi/10.1063/1.872465>`_.

    Parameters
    ----------
    all_species : List[ExtendedParticleList]
        all_species
    flux_surface : FluxSurface
        flux_surface
    mu_N : int
        mu_N
    dataset_input : xarray.Dataset
        dataset_input
    """

    @classmethod
    def from_xarray_surface(cls, dataset: xarray.Dataset, flux_surface: FluxSurface):
        """
        Alternate constructor from an `xarray.Dataset` and an associated `~FluxSurface`.
        """

        arrays = {}
        lists = []
        data_vars = ["n", "gradn", "T", "gradT"]
        dataset["basic_elements"] = "particle", [
            Particle(i).element
            if Particle(i).element is not None
            else Particle(i).symbol
            for i in dataset.particle.values
        ]
        # TODO does not handle deuterium because H 1+.isotope = None, why?
        # TODO to wciaż jest przekombinowane i wszystko to trzeba spłaszczyć, bo bez sensu
        for basic_element, a in dataset.groupby("basic_elements"):
            lists.append(
                ExtendedParticleList(
                    [Particle(i) for i in a.particle.values],
                    u.Quantity(a.T, dataset.attrs["T unit"]),
                    u.Quantity(a.n, dataset.attrs["n unit"]),
                    u.Quantity(a.gradT, dataset.attrs["gradT unit"]),
                    u.Quantity(a.gradn, dataset.attrs["gradn unit"]),
                )
            )

        return cls(
            lists,
            flux_surface,
            dataset_input=dataset,
        )

    # profile
    def __init__(
        self,
        all_species: IonizationStateCollection,
        flux_surface: FluxSurface,
        *,
        mu_N: int = None,
        dataset_input: xarray.Dataset = None,
    ):
        self.all_species = all_species
        N = len(all_species)
        self.flux_surface = fs = flux_surface
        self._dataset_input = dataset_input

        charges = all_species.charge
        xi = all_species.ξ
        T_i = all_species.T.to(u.K, equivalencies=u.temperature_energy())
        n_i = all_species.n
        density_gradient = all_species.dn
        temperature_gradient = (all_species.dT * u.m).to(
            u.K, equivalencies=u.temperature_energy()
        ) / u.m
        pressure_gradient_over_n_i = constants.k_B * (
            T_i * density_gradient / n_i + temperature_gradient
        )
        μ = all_species.mu_hat(self.flux_surface, N=mu_N)

        Aai = (
            xi.reshape(1, 1, -1) * all_species.decompress(all_species.M_script, axis=2)
            - μ
        ).T
        # --- TD forces eq21
        # r"""$S_{\theta,\beta}^{ai}"""
        thermodynamic_forces = S_pt_θ = (
            fs.Fhat
            / charges
            * u.Quantity(
                [
                    pressure_gradient_over_n_i,
                    constants.k_B * temperature_gradient,
                    u.Quantity(np.zeros_like(temperature_gradient).value, u.J / u.m),
                ]
            )
        )

        BETA_AXIS = 1
        S_pt = (thermodynamic_forces[:, np.newaxis, :] * μ).sum(
            axis=BETA_AXIS
        )  # Equation 29
        r_pt = np.linalg.solve(Aai, S_pt.T).T
        r_sources = r_pt + 0  # TODO r_E itd
        rbar_sources_presum = xi[np.newaxis, :] * r_sources  # (3, N)
        rbar_sources = all_species.compress(rbar_sources_presum, axis=1)
        N_isotopes = all_species.num_isotopes
        assert rbar_sources.shape == (N_isotopes, 3)
        S_flows = (xi[:, np.newaxis, np.newaxis] * np.eye(3)) * u.Unit("N T / m3")
        r_flows = np.linalg.solve(Aai, S_flows)
        assert r_flows.shape == (N, 3, 3)
        np.testing.assert_allclose(
            np.swapaxes(r_flows, 1, 2), r_flows
        )  # symmetric matrix
        rbar_flows_presum = xi[:, np.newaxis, np.newaxis] * r_flows
        rbar_flows = all_species.compress(rbar_flows_presum, axis=0)
        assert rbar_flows.shape == (N_isotopes, 3, 3)
        self.r_pt = r_pt
        self.S_pt = S_pt
        self.thermodynamic_forces = thermodynamic_forces
        self.μ = μ

        # equation 34
        lhs = (
            u.Quantity(
                np.stack(N_isotopes * [np.eye(3)], axis=0),
                "J2 / (A m6)",
                dtype=np.float64,
            )
            + all_species.N_script.sum(axis=-1).T * rbar_flows
        )
        ubar = np.linalg.solve(lhs, rbar_sources).si

        Λ = -(all_species.N_script * ubar.reshape(1, 3, 1, -1)).sum(axis=(1, 3))

        self_consistent_u = (
            all_species.decompress(Λ, axis=1).reshape(-1, 3, 1) * r_flows
        ).sum(
            axis=1
        )  # TODO is axis=1 right?
        u_velocity = (self_consistent_u + r_sources.T).si

        self._charge_state_flows = {
            ai.symbol: u_velocity[i] for i, ai in enumerate(all_species)
        }

    def all_contributing_states_symbols(self) -> typing.Iterator[str]:
        """Helper iterator over all charge levels of all isotopes in the calculation."""
        for a in self.all_species:
            for _, ai in contributing_states(a):
                yield ai.symbol

    def M_script(self, a: IonizationState) -> np.ndarray:
        """Thin, cached wrapper on top of `~plasmapy.transport.Houlberg1997.M_script`."""
        sym = a[0].symbol
        if sym not in self.M_script_matrices:
            self.M_script_matrices[sym] = M_script(a, self.all_species)
        return self.M_script_matrices[sym]

    def N_script(self, a: IonizationState, b: IonizationState) -> np.ndarray:
        """Thin, cached wrapper on top of `~plasmapy.transport.Houlberg1997.N_script`."""
        sym_tuple = a[0].symbol, b[0].symbol
        if sym_tuple not in self.N_script_matrices:
            self.N_script_matrices[sym_tuple] = N_script(a, b)
        return self.N_script_matrices[sym_tuple]

    # profile
    def _funnymatrix(self, a_symbol):
        a = self._all_species_map[a_symbol]
        M = self.M_script(a)
        outputs = {}
        for _, ai in contributing_states(a):
            sym = ai.symbol
            output = self.thermodynamic_forces[sym] * M
            for b in self.all_species:
                N = self.N_script(a, b)
                for xj, bj in contributing_states(b):
                    output += xj * N * self.thermodynamic_forces[bj.symbol]
            outputs[sym] = output.sum(axis=1)
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
                sym = ai.symbol
                if sym not in self._charge_state_flows:
                    continue

                u_θ = (
                    self._charge_state_flows[sym] + self.thermodynamic_forces[sym]
                ) / B2fsav
                μ = self.μ[sym]
                Γ_BP = -(Fhat / ai.charge * (μ[0, :] * u_θ).sum()).si
                q_BP = -(
                    fs.Fhat * constants.k_B * a.T / ai.charge * (μ[1, :] * u_θ).sum()
                ).si
                results[sym] = Fluxes(
                    Γ_BP.to(particle_flux_unit), q_BP.to(heat_flux_unit)
                )
        return results

    @cached_property
    def _fluxes_PS(self):
        fs = self.flux_surface
        B2fsav = fs.flux_surface_average(fs.B2) * u.T ** 2  # flux surface averaged B^2
        Binv2fsav = fs.flux_surface_average(1 / fs.B2) / u.T ** 2
        results = {}
        fs = self.flux_surface
        for a in self.all_species:
            silly = self._funnymatrix(a[0].symbol)
            for xi, ai in contributing_states(a):
                sym = ai.symbol
                prefactor = (
                    -fs.Fhat / ai.charge * xi / B2fsav * (1 - B2fsav * Binv2fsav)
                )
                Γ_PS = prefactor * silly[sym][0]  # overlarge by s/m5
                q_PS = (
                    prefactor * constants.k_B * a.T * silly[sym][1]
                )  # overlarge by μ.unit
                results[sym] = Fluxes(
                    Γ_PS.to(particle_flux_unit), q_PS.to(heat_flux_unit)
                )
        return results

    @cached_property
    def _fluxes_CL(self):
        fs = self.flux_surface
        FSA = fs.flux_surface_average(fs.GradRho2 / fs.B2) / u.T ** 2
        # TODO fs.rho is [m]; fs.GradRho2 is actually [fs.gradRho]^2, gradRho is [1]
        # TODO FSA does not drop units; B2 and the others are unitless
        Fhat = self.flux_surface.Fhat
        results = {}
        for a in self.all_species:
            silly = self._funnymatrix(a[0].symbol)
            for xi, ai in contributing_states(a):
                sym = ai.symbol
                prefactor = FSA / Fhat * xi / ai.charge
                Γ_CL = prefactor * silly[sym][0]
                q_CL = prefactor * constants.k_B * a.T * silly[sym][1]
                results[sym] = Fluxes(
                    Γ_CL.to(particle_flux_unit), q_CL.to(heat_flux_unit)
                )
        return results

    @cached_property
    def fluxes(self):
        results = {}
        for a in self.all_species:
            for _, ai in contributing_states(a):
                sym = ai.symbol
                Γ_BP, q_BP = self._fluxes_BP[sym]
                Γ_PS, q_PS = self._fluxes_PS[sym]
                Γ_CL, q_CL = self._fluxes_CL[sym]
                results[sym] = Fluxes(Γ_BP + Γ_PS + Γ_CL, q_BP + q_PS + q_CL)
        return results

    def to_dataset(self, *, with_input=True) -> xarray.Dataset:
        r"""
        Converts the outputs to `~xarray.Dataset`.

        Parameters
        ----------
        with_input: bool (default: True)
            if True and `self` was initialized with an xarray object through
            `~self.from_xarray_surface`, return a dataset merged with the input
            object. Otherwise, return just the calculation results.

        """
        result = xarray.Dataset(
            {
                "total_particle_flux": (
                    "particle",
                    u.Quantity(
                        [flux.particle_flux for flux in self.fluxes.values()]
                    ).ravel(),
                ),
                "total_heat_flux": (
                    "particle",
                    u.Quantity(
                        [flux.heat_flux for flux in self.fluxes.values()]
                    ).ravel(),
                ),
                "BP_particle_flux": (
                    "particle",
                    u.Quantity(
                        [flux.particle_flux for flux in self._fluxes_BP.values()]
                    ).ravel(),
                ),
                "BP_heat_flux": (
                    "particle",
                    u.Quantity(
                        [flux.heat_flux for flux in self._fluxes_BP.values()]
                    ).ravel(),
                ),
                "CL_particle_flux": (
                    "particle",
                    u.Quantity(
                        [flux.particle_flux for flux in self._fluxes_CL.values()]
                    ).ravel(),
                ),
                "CL_heat_flux": (
                    "particle",
                    u.Quantity(
                        [flux.heat_flux for flux in self._fluxes_CL.values()]
                    ).ravel(),
                ),
                "PS_particle_flux": (
                    "particle",
                    u.Quantity(
                        [flux.particle_flux for flux in self._fluxes_PS.values()]
                    ).ravel(),
                ),
                "PS_heat_flux": (
                    "particle",
                    u.Quantity(
                        [flux.heat_flux for flux in self._fluxes_PS.values()]
                    ).ravel(),
                ),
                "diffusion_coefficient": (
                    "particle",
                    u.Quantity(list(self.diffusion_coefficient.values())).ravel(),
                ),
                "thermal_conductivity": (
                    "particle",
                    u.Quantity(list(self.thermal_conductivity.values())).ravel(),
                ),
                "bootstrap_current": self.bootstrap_current,
                # TODO this probably won't fit in the common array because of spatial dependence
                # "local_heat_flux": (
                #     ("particle", "directions", "lp", ),
                #     u.Quantity(list(self.local_heat_flux_components.values())),
                # ),
                # "local_particle_velocities": (
                #     ("particle", "directions", "lp", ),
                #     u.Quantity(list(self.local_flow_velocities.values())),
                # ),
            },
            {
                "particle": list(self.all_contributing_states_symbols()),
                "psi": self.flux_surface.psi,
                # "directions": ["poloidal", "toroidal", "parallel", "perpendicular"],
                # "lp": self.flux_surface.lp,
            },
        )

        if with_input and self._dataset_input is not None:
            return xarray.merge([result, self._dataset_input])
        else:
            return result

    @cached_property
    def diffusion_coefficient(self):
        results = {}
        for a in self.all_species:
            for (_, ai), dn in zip(contributing_states(a), a.dn):
                sym = ai.symbol
                flux = self.fluxes[sym].particle_flux
                results[sym] = -flux / dn  # Eq48 TODO this is a partial adaptation
        return results

    @cached_property
    def thermal_conductivity(self):
        results = {}
        for a in self.all_species:
            for (_, ai), dT in zip(contributing_states(a), a.dT):
                sym = ai.symbol
                flux = self.fluxes[sym].heat_flux
                results[sym] = -flux / dT
        return results

    @cached_property
    def bootstrap_current(self) -> u.J / u.m ** 2 / u.T:
        """
        Bootstrap current caused by the charge state flows.
        """

        def gen():
            for a in self.all_species:
                for (_, ai), n in zip(contributing_states(a), a.n):
                    sym = ai.symbol
                    u_velocity = self._charge_state_flows[sym][0]
                    yield ai.charge * n * u_velocity
                    # eq 37, second term

        return sum(gen())

    @cached_property
    def local_flow_velocities(self):
        fs = self.flux_surface
        B2fsav = fs.flux_surface_average(fs.B2) * u.T ** 2  # flux surface averaged B^2
        B_p = fs.Bp * u.T
        B_t = fs.Bphivals * u.T  # TODO needs renaming T_T
        B = fs.Bmag * u.T
        results = {}
        for a in self.all_species:
            for _, ai in contributing_states(a):
                sym = ai.symbol
                u_θ = (
                    self._charge_state_flows[sym] + self.thermodynamic_forces[sym]
                ) / B2fsav
                u_hat_theta_1_ai = u_θ[0]
                S_theta_1_ai = self.thermodynamic_forces[sym][0]
                u_p_ai = B_p * u_hat_theta_1_ai
                u_t_ai = B_t * u_hat_theta_1_ai - S_theta_1_ai / B_t
                u_parallel_ai = B * u_hat_theta_1_ai - S_theta_1_ai / B
                u_perp_ai = B_p / B / B_t * S_theta_1_ai
                results[sym] = u.Quantity([u_p_ai, u_t_ai, u_parallel_ai, u_perp_ai]).si
        return results

    @cached_property
    def local_heat_flux_components(self):
        fs = self.flux_surface
        B2fsav = fs.flux_surface_average(fs.B2) * u.T ** 2  # flux surface averaged B^2
        B_p = fs.Bp * u.T
        B_t = fs.Bphivals * u.T  # TODO needs renaming T_T
        B = fs.Bmag * u.T
        results = {}
        for a in self.all_species:
            n_i = a.n
            T_i = a.T
            p = n_i * T_i
            for _, ai in contributing_states(a):
                sym = ai.symbol
                u_θ = (
                    self._charge_state_flows[sym] + self.thermodynamic_forces[sym]
                ) / B2fsav
                u_hat_theta_2_ai = u_θ[1]
                S_theta_2_ai = self.thermodynamic_forces[sym][1]
                p_ai = p[ai.charge_number]
                q_p_ai = 5 / 2 * p_ai * B_p * u_hat_theta_2_ai
                q_t_ai = 5 / 2 * p_ai * (B_t * u_hat_theta_2_ai - S_theta_2_ai / B_t)
                q_parallel_ai = 5 / 2 * p_ai * (B * u_hat_theta_2_ai - S_theta_2_ai / B)
                q_perp_ai = 5 / 2 * p_ai * B_p / B / B_t * S_theta_2_ai
                results[sym] = u.Quantity([q_p_ai, q_t_ai, q_parallel_ai, q_perp_ai]).si
        return results
