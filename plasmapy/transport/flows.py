__all__ = [
    "FlowCalculator",
]


import numpy as np
import typing
import xarray

from astropy import constants
from astropy import units as u
from collections import defaultdict, namedtuple
import scipy.linalg

from plasmapy.particles import IonizationStateCollection, Particle
from plasmapy.plasma.fluxsurface import FluxSurface

from .Houlberg1997 import ExtendedParticleList

A_AXIS = 0
B_AXIS = 1
ALPHA_AXIS = 2
BETA_AXIS = 3

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

        # TODO does not handle deuterium because H 1+.isotope = None, why?
        lists = ExtendedParticleList(
            [Particle(i) for i in dataset.particle.values],
            u.Quantity(dataset.T, dataset.attrs["T unit"]),
            u.Quantity(dataset.n, dataset.attrs["n unit"]),
            u.Quantity(dataset.gradT, dataset.attrs["gradT unit"]),
            u.Quantity(dataset.gradn, dataset.attrs["gradn unit"]),
        )
        return cls(
            lists,
            flux_surface,
            dataset_input=dataset,
        )

    # profile
    def __init__(
        self,
        all_species: ExtendedParticleList,
        flux_surface: FluxSurface,
        *,
        mu_N: int = None,
        dataset_input: xarray.Dataset = None,
    ):
        self.all_species = all_species
        self.flux_surface = flux_surface
        self._dataset_input = dataset_input

        self.μ = all_species.mu_hat(self.flux_surface, N=mu_N)

    @cached_property
    @validate_quantities
    def thermodynamic_forces(self) -> u.Unit("V / m"):
        charges = self.all_species.charge
        xi = self.all_species.ξ
        T_i = self.all_species.T.to(u.K, equivalencies=u.temperature_energy())
        n_i = self.all_species.n
        density_gradient = self.all_species.dn
        temperature_gradient = self.all_species.dT
        pressure_gradient_over_n_i = self.all_species.dP / self.all_species.n

        # --- TD forces eq21
        # r"""$S_{\theta,\beta}^{ai}"""
        thermodynamic_forces = (   # S_pt_θ
            self.flux_surface.Fhat
                / charges[np.newaxis, :]
            * u.Quantity(
                [
                    pressure_gradient_over_n_i,
                    constants.k_B * temperature_gradient,
                    u.Quantity(np.zeros_like(temperature_gradient).value, u.J / u.m),
                ]
            )
        ).T    # (3, N) -> (N, 3)
        return thermodynamic_forces

    @cached_property
    @validate_quantities
    def _Aai(self) -> u.Unit("kg / (m3 s)"):
        all_species = self.all_species
        M_script_particlewise = all_species.decompress(all_species.M_script, axis=A_AXIS)
        xi = all_species.ξ
        Aai = xi.reshape(-1, 1, 1) * M_script_particlewise - self.μ
        return Aai

    @cached_property
    @validate_quantities
    def _charge_state_flows(self) -> u.Unit("V / m"):
        """Reconstructs charge state flows using equation 31.

            TODO

            There's a bug in this code.

            Find it.

            Remove it.
        """
        p_eb = 0 # TODO
        rhat, crhat = self._rhat_crhat
        Λ_ai = self.all_species.decompress(self.Λ, axis = 0)

        r_flows = rhat[:,:,:3]
        r_sources = rhat[:,:,3:]
        r_sources[:,:,2] *= p_eb
        velocity_responses = r_flows[:,:,0]
        heat_responses = r_flows[:,:,1]
        u2_responses = r_flows[:,:,2]
        u_hat = Λ_ai[:,0][:,np.newaxis] * velocity_responses +\
                Λ_ai[:,1][:,np.newaxis] * heat_responses +\
                Λ_ai[:,2][:,np.newaxis] * u2_responses +\
                r_sources.sum(axis=-1) # TODO check
        return u_hat * u.Unit("V / m")

    @cached_property
    @validate_quantities
    def Λ(self):
        """Equation 27 for Λ"""
        Nscript = self.all_species.N_script

        M = self.all_species.num_isotopes
        xab = self._xab
        isotopic_velocities = xab[:M]
        isotopic_heat_flows = xab[M:2*M]
        u2_flows = xab[2*M:]
        u_bar = np.stack([isotopic_velocities, isotopic_heat_flows, u2_flows], axis=1)
        Λ = -(Nscript * u_bar[np.newaxis, :, :, :]).sum(axis=(B_AXIS, BETA_AXIS))
        return Λ

    @cached_property
    @validate_quantities
    def _xab(self) -> u.Unit("m3 s / kg"):   #TODO
        r"""Solves equation 34 for $\bar{u}_k^a$."""
        rhat, crhat = self._rhat_crhat
        lhs = self._lhs
        rhs = np.vstack(crhat[:,:,3:])
        xab = np.linalg.solve(lhs, rhs).si
        return xab

    @cached_property
    @validate_quantities
    def _lhs(self) -> u.Unit("kg / (m3 s)"):
        r"""The left hand side of equation 34."""
        M = self.all_species.num_isotopes
        crhat = self._rhat_crhat[1]
        # The $\bar{u}_k^a$ term:
        ab = u.Quantity(np.eye(3 * M), "kg / (m3 s)")

        # The summation term:
        Nscript = self.all_species.N_script
        for im in range(M):
            for m in range(3):
                m1 = im + m * M
                for jm in range(M):
                    for l in range(3):
                        l1 = jm + l * M
                        for k in range(3):
                            value = Nscript[im, jm, k,l] * crhat[im, m,k]
                            ab[m1, l1] += value
        return ab
        
    @cached_property
    def _rhat_crhat(self):
        p_eb = 0
        N = len(self.all_species)
        M = self.all_species.num_isotopes
        A_matrices = self._Aai
        charges = self.all_species.charge
        temperatures = self.all_species.T
        densities = self.all_species.n
        orders = np.arange(3)
        xi = self.all_species.ξ.value
        rhat = np.zeros((N, 3, 6))
        crhat = np.zeros((M, 3, 6))
        rhatp = np.zeros((N, N, 3))
        rhatt = np.zeros((N, N, 3))
        srcth = self.thermodynamic_forces
        for i in range(N):
            isotope_index = self.all_species.isotope_indices[i]
            assert A_matrices.shape[0] == N
            lu = A_lu, A_pivots = scipy.linalg.lu_factor(A_matrices[i])
            response1 = scipy.linalg.lu_solve(lu, np.eye(3) * xi[i])
            rhat[i, :, :3] = response1
            crhat[isotope_index, :, :3] += xi[i] * response1

            response2 = scipy.linalg.lu_solve(lu, (srcth[i] * self.μ[i]).sum(axis=1))
            rhat[i, :, 3] = response2
            crhat[isotope_index, :, 3] += xi[i] * response2

            input4 = np.append(-charges[i] * densities[i], [0, 0])
            response4 = scipy.linalg.lu_solve(lu, input4)
            rhat[i, :, 4] = response4
            crhat[isotope_index, :, 4] += xi[i] * p_eb * response4

            # external forces are skipped
        return rhat, crhat

    # @cached_property
    # def _crhat(self):
    #     p_eb 
    #     rhat = self._rhat.copy()


    @cached_property
    @validate_quantities
    def _funnymatrix(self) -> u.Unit("J2 / (A m6)"):
        M = self.all_species.decompress(
            self.all_species.M_script,
            axis=0,
        )
        output = self.thermodynamic_forces[:, np.newaxis, :] * M
        N = self.all_species.decompress(
                self.all_species.N_script,
                axis=(0, 1),
        )
        ξ = self.all_species.ξ.reshape(1, -1, 1, 1)
        presum = N * ξ * self.thermodynamic_forces[np.newaxis, :, np.newaxis, :]
        sum_b = presum.sum(axis=B_AXIS)
        return (output + sum_b).sum(axis=-1)

    @cached_property
    def _fluxes_BP(self):
        fs = self.flux_surface
        u_θ = (self._charge_state_flows + self.thermodynamic_forces) / fs.fsa_B2
        μ = self.μ
        charge = self.all_species.charge
        Γ_BP = -(fs.Fhat / charge * (μ[:, 0] * u_θ).sum(axis=-1)).si
        q_BP = -(
            fs.Fhat
            * constants.k_B
            * self.all_species.T
            / charge
            * (μ[:, 1] * u_θ).sum(axis=1)
        ).si
        return Fluxes(Γ_BP.to(particle_flux_unit), q_BP.to(heat_flux_unit))

    @cached_property
    def _fluxes_PS(self):
        fs = self.flux_surface
            # TODO move these and others into fluxsurface.py
        B2fsav = fs.fsa_B2
        Binv2fsav = fs.fsa_invB2
        fs = self.flux_surface
        ξ = self.all_species.ξ
        prefactor = (
            -fs.Fhat / self.all_species.charge * ξ / B2fsav * (1 - B2fsav * Binv2fsav)
        )
        Γ_PS = prefactor * self._funnymatrix[:,0]
        q_PS = prefactor * constants.k_B * self.all_species.T * self._funnymatrix[:,1]
        return Fluxes(Γ_PS.to(particle_flux_unit), q_PS.to(heat_flux_unit))

    @cached_property
    def _fluxes_CL(self):
        fs = self.flux_surface
        FSA = fs.grbm2
        # TODO fs.rho is [m]; fs.GradRho2 is actually [fs.gradRho]^2, gradRho is [1]
        # TODO FSA does not drop units; B2 and the others are unitless
        Fhat = fs.Fhat
        ξ = self.all_species.ξ
        prefactor = FSA / Fhat * self.all_species.ξ / self.all_species.charge
        Γ_CL = prefactor * self._funnymatrix[:,0]
        q_CL = prefactor * constants.k_B * self.all_species.T * self._funnymatrix[:,1]
        return Fluxes(Γ_CL.to(particle_flux_unit), q_CL.to(heat_flux_unit))

    @cached_property
    def fluxes(self):
        Γ_BP, q_BP = self._fluxes_BP
        Γ_PS, q_PS = self._fluxes_PS
        Γ_CL, q_CL = self._fluxes_CL
        return Fluxes(Γ_BP + Γ_PS + Γ_CL, q_BP + q_PS + q_CL)

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
                "total_particle_flux": ("particle", self.fluxes.particle_flux),
                "total_heat_flux": ("particle", self.fluxes.heat_flux),
                "BP_particle_flux": ("particle", self._fluxes_BP.particle_flux),
                "BP_heat_flux": ("particle", self._fluxes_BP.heat_flux),
                "CL_particle_flux": ("particle", self._fluxes_CL.particle_flux),
                "CL_heat_flux": ("particle", self._fluxes_CL.heat_flux),
                "PS_particle_flux": ("particle", self._fluxes_PS.particle_flux),
                "PS_heat_flux": ("particle", self._fluxes_PS.heat_flux),
                "diffusion_coefficient": ("particle", self.diffusion_coefficient),
                "thermal_conductivity": ("particle", self.thermal_conductivity),
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
                "particle": self.all_species.symbols,
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
        flux = self.fluxes.particle_flux
        dn = self.all_species.dn
        return -flux / dn  # Eq48 TODO this is a partial adaptation

    @cached_property
    def thermal_conductivity(self):
        dT = self.all_species.dT
        flux = self.fluxes.heat_flux
        return -flux / dT

    @cached_property
    def bootstrap_current(self) -> u.J / u.m ** 2 / u.T:
        """
        Bootstrap current caused by the charge state flows.

        Needs normalizing by a magnetic field.
        """

        return (
            self._charge_state_flows[:, 0]
            * self.all_species.charge
            * self.all_species.n
        ).sum()

    @cached_property
    def local_flow_velocities(self):
        fs = self.flux_surface
        B2fsav = fs.fsa_B2
        B_p = fs.Bp 
        B_t = fs.Bphivals # TODO needs renaming T_T
        B = fs.Bmag 
        u_θ = (self._charge_state_flows + self.thermodynamic_forces) / B2fsav
        u_hat_theta_1_ai = u_θ[:, 0][:, np.newaxis]
        S_theta_1_ai = self.thermodynamic_forces[:, 0][:, np.newaxis]
        u_p_ai = B_p * u_hat_theta_1_ai
        u_t_ai = B_t * u_hat_theta_1_ai - S_theta_1_ai / B_t
        u_parallel_ai = B * u_hat_theta_1_ai - S_theta_1_ai / B
        u_perp_ai = B_p / B / B_t * S_theta_1_ai
        return u.Quantity([u_p_ai, u_t_ai, u_parallel_ai, u_perp_ai]).si

    @cached_property
    def local_heat_flux_components(self):
        fs = self.flux_surface
        B2fsav = fs.fsa_B2
        B_p = fs.Bp
        B_t = fs.Bphivals # TODO needs renaming T_T
        B = fs.Bmag
        n_i = self.all_species.n
        T_i = self.all_species.T
        p = n_i * T_i
        u_θ = (self._charge_state_flows + self.thermodynamic_forces) / B2fsav
        u_hat_theta_2_ai = u_θ[:, 1][:, np.newaxis]
        S_theta_2_ai = self.thermodynamic_forces[:, 1][:, np.newaxis]
        p_ai = p[:, np.newaxis]
        q_p_ai = 5 / 2 * p_ai * B_p * u_hat_theta_2_ai
        q_t_ai = 5 / 2 * p_ai * (B_t * u_hat_theta_2_ai - S_theta_2_ai / B_t)
        q_parallel_ai = 5 / 2 * p_ai * (B * u_hat_theta_2_ai - S_theta_2_ai / B)
        q_perp_ai = 5 / 2 * p_ai * B_p / B / B_t * S_theta_2_ai
        return u.Quantity([q_p_ai, q_t_ai, q_parallel_ai, q_perp_ai]).si
