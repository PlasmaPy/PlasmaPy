r"""This module implements helper functions from [1]_.

.. [1] Houlberg et al, Bootstrap current and neoclassical transport in tokamaks of arbitrary collisionality and aspect ratio, 1997,
   Physics of Plasmas 4, 3230 (1997); , JGR, 117, A12219, doi: `10.1063/1.872465
   <https://aip.scitation.org/doi/10.1063/1.872465>`_.

This work is referenced in docstrings as `Houlberg_1997`.
"""

import numpy as np
from numpy import pi
import warnings

from astropy import constants
from astropy import units as u
from functools import cached_property
from scipy.special import erf
from typing import Iterable, Union, Optional

from plasmapy.formulary import thermal_speed
from plasmapy.formulary.collisions import Coulomb_logarithm
from plasmapy.formulary.mathematics import Chandrasekhar_G
from plasmapy.particles import Particle, electron
from plasmapy.plasma.fluxsurface import FluxSurface
import xarray
from plasmapy.utils.decorators import validate_quantities

A_AXIS = 0
B_AXIS = 1
ALPHA_AXIS = 2
BETA_AXIS = 3

temperature_unit = u.K
try:
    from scipy.integrate import trapz as trapezoid
except ImportError:
    from scipy.integrate import trapezoid
1
__all__ = [
    "ExtendedParticleList",
    "ωm",
    "_B17",
]

from plasmapy.particles.particle_collections import ParticleList

LaguerrePolynomials = [
    lambda x: np.ones_like(x),
    lambda x: 5 / 2 - x,
    lambda x: 35 / 8 - 7 / 2 * x + 1 / 2 * x ** 2,
]


class ExtendedParticleList(ParticleList):
    def __init__(
        self,
        particles: Iterable,
        T: u.eV,
        n: u.m ** -3,
        dT: u.eV / u.m = None,
        dn: u.m ** -4 = None,
    ):
        symbols = np.array([particle.symbol for particle in particles])
        sorted_particles = np.sort(symbols)
        if not (sorted_particles == symbols).all():
            warnings.warn("Your particles array was not sorted. Sorting!")
            indices = np.argsort(symbols)
            particles = [Particle(s) for s in sorted_particles]
            if not T.size == 1:
                T = T[indices]
            n = n[indices]
            dT = dT[indices]
            dn = dn[indices]

        super().__init__(particles)
        T = u.Quantity(T).to(temperature_unit, equivalencies=u.temperature_energy())
        self.n = u.Quantity(n)
        if T.size == 1:
            T = np.broadcast_to(T, self.n.shape) * T.unit
        self.T = T

        assert len(self.n) == len(particles)
        # TODO replace instances of element by isotope when that's fixed
        self.basic_elements = [
            p.element if p.element is not None else p.symbol for p in self
        ]
        indices = np.argsort(self.basic_elements)
        is_sorted = (np.diff(indices) > 0).all()
        if not is_sorted:
            warnings.warn("Unsorted input array! Sorting.")
            self._data = np.array(self._data)[indices].tolist()
            self.n = self.n[indices]
            self.T = self.T[indices]

        if dT is not None:
            self.dT = (u.Quantity(dT[indices]) * u.m).to(temperature_unit, equivalencies = u.temperature_energy()) / u.m
            assert len(self.dT) == len(self.T)
        if dn is not None:
            self.dn = u.Quantity(dn[indices])
            assert len(self.dn) == len(self.n)

        uniq_all = np.unique(
            self.basic_elements, return_index=True, return_inverse=True
        )
        unique, self.split_index, self.inverse_index = uniq_all
        self.num_isotopes = unique.size
        index = 0
        isotope_indices = {}
        for elem in self.basic_elements:
            if elem not in isotope_indices:
                isotope_indices[elem] = index
                index += 1
        self.isotope_indices = np.array(list(map(isotope_indices.get, self.basic_elements)))  # TODO could definitely be refactored

    @cached_property
    def dP(self) -> u.Pa / u.m:
        return constants.k_B * (
            self.T * self.dn + self.dT * self.n
        )

    def split_isotopes(self, arr, axis):
        assert arr.shape[axis] == len(self)
        output = np.split(arr, self.split_index, axis=axis)
        return [a for a in output if a.size > 0]

    def compress(self, arr, axis, aggregator=np.sum):
        """Compresses an N-sized array into an M sized array.

        N is the number of particles (different charge states)
        in the system, while M is the number of isotopes.

        The input must be N-sized along the dimension specified
        by `axis`.

        If the size of axis is equal to M (already compressed) and
        N = M, the matrix is passed through without change.
        This helps preserve ordering.
        If N != M in this case, an exception is raised.

        
        """
        if isinstance(axis, tuple):
            for integer in axis:
                arr = self.compress(arr, integer, aggregator)
            return arr
        elif arr.shape[axis] == self.num_isotopes:
            if self.num_isotopes == len(self):
                return arr
            else:
                raise ValueError("Are you sure you should be compressing here?")
        else:
            split = self.split_isotopes(arr, axis=axis)
            return u.Quantity([aggregator(a, axis=axis) for a in split])
        # TODO needs a test - breaks ordering

    def decompress(self, arr, axis):
        # this is *not* quite invertible...
        if isinstance(axis, tuple):
            for integer in axis:
                arr = self.decompress(arr, integer)
            return arr
        else:
            return np.take(arr, self.inverse_index, axis=axis)

    @cached_property
    def mass_density(self) -> u.kg / u.m**3:
        r"""Mass densities for each particle.

        Returns
        -------
        u.kg/u.m**3

        """
        return self.n * self.mass

    @cached_property
    def thermal_speed(self) -> u.m / u.s:
        return u.Quantity([thermal_speed(T, p, mass=m) for T, p, m in zip(self.T, self, self.mass)])

    @cached_property
    def ξ(self) -> np.ndarray:
        # assert np.all(np.array(self.basic_elements) == np.sort(self.basic_elements))
        input_indices = np.unique(
            self.basic_elements,
            return_index=True,
        )[1]
        charge_numbers = np.split(self.charge_number, input_indices)
        ns = np.split(self.n, input_indices)
        outputs = []
        for charge_number, n in zip(charge_numbers, ns):
            array = charge_number ** 2 * n
            array /= array.sum()
            outputs.append(array)
        result = np.concatenate(outputs)
        assert result.size == len(self)
        return result

    @cached_property
    def isotopic_thermal_speed(self) -> u.m/u.s:
        return self.compress(self.thermal_speed, axis=0, aggregator = np.mean)

    @cached_property
    def mass(self) -> u.kg:
        return u.Quantity([p.mass_number * u.u for p in self])

    @cached_property
    def isotopic_mass(self) -> u.kg:
        return self.compress(self.mass, axis=0, aggregator = np.mean)

    @cached_property
    def isotopic_temperature(self) -> u.eV:
        return self.compress(self.T, axis=0, aggregator = np.mean)

    @cached_property
    def xab_ratio(self) -> np.ndarray:
        speed = self.isotopic_thermal_speed
        return (
            speed[np.newaxis, :]
            / speed[:, np.newaxis]
        )

    @cached_property
    def mass_ratio(self) -> np.ndarray:
        mass = self.isotopic_mass
        return mass[:, np.newaxis] / mass[np.newaxis, :]

    @cached_property
    def temperature_ratio(self) -> np.ndarray:
        temperature = self.isotopic_temperature
        return (
            temperature[:, np.newaxis] / temperature[np.newaxis, :]
        )

    @cached_property
    def M_matrix(self) -> np.ndarray:
        """(N, N, 3, 3)"""
        xab = self.xab_ratio
        mass_ratio = self.mass_ratio
        M11 = -(1 + mass_ratio) / (1 + xab ** 2) ** (3 / 2)
        M12 = 3 / 2 * (1 + mass_ratio) / (1 + xab ** 2) ** (5 / 2)
        M21 = M12
        M22 = -(13 / 4 + 4 * xab ** 2 + 15 / 2 * xab ** 4) / (1 + xab ** 2) ** (5 / 2)
        M13 = -15 / 8 * (1 + mass_ratio) / (1 + xab ** 2) ** (7 / 2)
        M31 = M13
        M23 = (69 / 16 + 6 * xab ** 2 + 63 / 4 * xab ** 4) / (1 + xab ** 2) ** (7 / 2)
        M32 = M23
        M33 = -(433 / 64 + 17 * xab ** 2 + (459 / 8) * xab ** 4 + 28 * xab**6 + (175 / 8) * xab ** 8) / (
            1 + xab ** 2
        ) ** (9 / 2)
        ordering = [[M11, M12, M13], [M21, M22, M23], [M31, M32, M33]]
        M = np.array(ordering)
        output = np.moveaxis(M, (2, 3), (0, 1))
        return output

    @cached_property
    def N_matrix(self) -> np.ndarray:
        """(N, N, 3, 3)"""
        # TODO np.where(all_species.N_matrix - np.swapaxes(all_species.N_matrix, 0, 1))  should maybe be symmetric, but doesn't seem to be
        # cross-ref with NCLASS site errata
        xab = self.xab_ratio
        temperature_ratio = self.temperature_ratio
        mass_ratio = self.mass_ratio
        N11 = (1 + mass_ratio) / (1 + xab ** 2) ** (3 / 2)
        N21 = -3 / 2 * (1 + mass_ratio) / (1 + xab ** 2) ** (5 / 2)
        N31 = 15 / 8 * (1 + mass_ratio) / (1 + xab ** 2) ** (7 / 2)
        M12 = 3 / 2 * (1 + mass_ratio) / (1 + xab ** 2) ** (5 / 2)
        N12 = -(xab ** 2) * M12
        N22 = (
            27
            / 4
            * (temperature_ratio) ** (1 / 2)
            * xab ** 2
            / (1 + xab ** 2) ** (5 / 2)
        )
        M13 = -15 / 8 * (1 + mass_ratio) / (1 + xab ** 2) ** (7 / 2)
        N13 = -(xab ** 4) * M13
        N23 = -225 / 16 * temperature_ratio * xab ** 4 / (1 + xab ** 2) ** (7 / 2)
        N32 = N23 / temperature_ratio
        N33 = (
            2625
            / 64
            * temperature_ratio ** (1 / 2)
            * xab ** 4
            / (1 + xab ** 2) ** (9 / 2)
        )
        N = np.array([[N11, N12, N13], [N21, N22, N23], [N31, N32, N33]])
        output = np.moveaxis(N, (2, 3), (0, 1))
        return output

    @cached_property
    def tau(self):
        """Equations A3, A4 from |Houlberg_1997|"""
        charge_a = self.charge[:, np.newaxis]
        charge_b = charge_a.T
        charge_number_a = self.charge[:, np.newaxis]
        charge_number_b = charge_number_a.T
        mass = self.isotopic_mass.to(u.u)
        mass_a = mass[:, np.newaxis]
        mass_b = mass_a.T
        T = self.isotopic_temperature.to(u.keV, equivalencies=u.temperature_energy())
        T_a = T[:, np.newaxis]
        T_b = T_a.T
        xn = self.compress(self.n, axis=0)
        xnz = self.compress(self.charge_number * self.n, axis=0)
        xz = xnz / xn
        xz_a = xz[:, np.newaxis]
        xz_b = xz_a.T
        n_a = self.n[:, np.newaxis]
        n_b = n_a.T
        xnz2 = self.compress(self.charge_number**2 * self.n, axis=0)
        xnz2_a = xnz2[:, np.newaxis]
        xnz2_b = xnz2_a.T
        t_ele = self.T.to(u.keV, equivalencies=u.temperature_energy())
        clnab = 37.8 - np.log((np.sqrt(self.n) / t_ele).value)
        xlnab = 40.3 - np.log((
            xz_a * xz_b * (mass_a + mass_b) /
            (mass_a * T_b + mass_b * T_a) *
            np.sqrt(xnz2_a / T_a + xnz2_b / T_b)
        ).value)
        for i, p1 in enumerate(self):
            for j, p2 in enumerate(self):
                index = None
                if p1 == electron:
                    index = i
                elif p2 == electron:
                    index = j

                if index is not None:
                    xlnab[self.isotope_indices[index], :] = clnab[index]
                    xlnab[:, self.isotope_indices[index]] = clnab[index]

        CL_matrix = self.decompress(xlnab, axis=(0, 1))
        c1=4.0/3.0/np.sqrt(pi)*4.0*np.pi/(4.0*np.pi*constants.eps0)**2
        c1_fortran_ver=4.0/3.0/np.sqrt(pi)*4.0*np.pi*(constants.e.si / 
            (4.0*np.pi*constants.eps0))**2*(constants.e.si/constants.m_p)**2
        c2 = self.thermal_speed**3 * self.mass **2 / c1
        c2_fortran_ver=(self.thermal_speed**3)*self.mass.to(u.u)**2/c1_fortran_ver
        tau = c2 / CL_matrix / charge_a ** 2 / n_b / charge_b ** 2
        tau_fortran_ver = c2_fortran_ver / CL_matrix / charge_number_a ** 2 / n_b / charge_number_b ** 2
        return tau.T

    @cached_property
    def effective_momentum_relaxation_rate(self):
        """Equations A3, A4 from |Houlberg_1997|"""
        mass_a = self.mass[:, np.newaxis]
        n_a = self.n[:, np.newaxis]
        amnt = self.compress(mass_a * n_a / self.tau, axis=(0,1))
        return amnt

    @cached_property
    @validate_quantities
    def N_script(self) -> u.kg / u.m**3 / u.s:
        """Weighted field particle matrix - equation A2b from |Houlberg_1997|"""
        N = self.N_matrix
        # Equation A2b
        N_script = self.effective_momentum_relaxation_rate[:, :, np.newaxis, np.newaxis] * N
        return N_script

    @cached_property
    @validate_quantities
    def M_script(self) -> u.kg / u.m**3 / u.s:
        """Weighted test particle matrix - equation A2a from |Houlberg_1997|"""
        # Equation A2a
        integrand = self.M_matrix * self.effective_momentum_relaxation_rate[:, :, np.newaxis, np.newaxis]
        return integrand.sum(axis=B_AXIS)

    def x_over_xab(self, x):
        xab = self.xab_ratio
        x_over_xab = (x[:, np.newaxis, np.newaxis] / xab[np.newaxis, ...]).value
        return x_over_xab

    @validate_quantities
    def pitch_angle_diffusion_rate(
        self,
        x: np.ndarray,
    ) -> u.s**-1:
        """The pitch angle diffusion rate, nu_{D,ai}, equation B4b from |Houlberg_1997|

        Parameters
        ----------
        x : np.ndarray
            distribution velocity relative to the thermal velocity.
        """
        # Houlberg_1997, equation B4b,
        denominator = x ** 3

        x_over_xab = self.x_over_xab(x)
        numerator = erf(x_over_xab) - Chandrasekhar_G(x_over_xab)
        fraction = numerator / denominator[:, np.newaxis, np.newaxis]
        sum_items = fraction / self.tau[np.newaxis, :, :]

        result = (
            (3 * np.sqrt(np.pi) / 4)
            * sum_items
        )
        return result.sum(axis=-1)

    @validate_quantities
    def K_B_ai(
        self,
        x: np.ndarray,
        trapped_fraction: float,
        *,
        orbit_squeezing: bool = False,
    ) -> u.s**-1:
        """Banana regime contribution to effective viscosity - eq. B1 from |Houlberg_1997|

        Parameters
        ----------
        x : np.ndarray
            x
        orbit_squeezing : bool (default `False`)
            orbit_squeezing
        """
        # eq. B1-B4, Houlberg_1997
        f_t = trapped_fraction
        f_c = 1 - f_t
        if orbit_squeezing:
            raise NotImplementedError(
                "TODO allow for non-zero, changing radial electric fields (orbit squeezing)"
            )
        else:
            S_ai = 1  # Equation B2
        padr = self.pitch_angle_diffusion_rate(x)
        return padr * f_t / f_c / S_ai ** 1.5

    @validate_quantities
    def ν_T_ai(self, x: np.ndarray) -> u.s**-1:
        """Characteristic rate of anisotropy relaxation.

        Equation B12 from |Houlberg_1997|.

        Parameters
        ----------
        x : np.ndarray
            distribution velocity relative to the thermal velocity ($v / v_th$).
        """
        prefactor = 3 * np.pi ** 0.5 / 4 * self.ξ / self.mass_density

        # TODO check      if b is not a:  # TODO is not should work
        x_over_xab = self.x_over_xab(x)
        x = x.reshape(-1, 1, 1)
        part1 = (erf(x_over_xab) - 3 * Chandrasekhar_G(x_over_xab)) / x ** 3
        part2 = 4 * (self.temperature_ratio + self.xab_ratio ** -2)
        part2full = part2 * Chandrasekhar_G(x_over_xab) / x
        summation = (part1 + part2full) * self.effective_momentum_relaxation_rate

        result = prefactor * self.decompress(summation.sum(axis=-1), axis=1)
        return result

    @validate_quantities
    def K_ps_ai(
        self,
        x: np.ndarray,
        flux_surface: FluxSurface,
        *,
        m_max=100,
        with_B10 = False,
    ):
        """Pfirsch-Schlüter regime contribution to effective viscosity - eq. B8 from |Houlberg_1997|

        Parameters
        ----------
        x : np.ndarray
            x
        flux_surface : FluxSurface
            flux_surface
        m_max :
            m_max
        """
        ν = self.ν_T_ai(x)[:, :, np.newaxis]

        F = flux_surface.F_m(m_max)
        m = np.arange(1, m_max + 1)
        ω = ωm(x, m, self.thermal_speed, flux_surface.gamma)
        B10 = (
            - 1.5 * (ν / ω) ** 2
            - 9 / 2 * (ν / ω) ** 4
            + (1 / 4 + (3 / 2 + 9 / 4 * (ν / ω) ** 2) * (ν / ω) ** 2)
            * (2 * ν / ω)
            * np.arctan(ω / ν).value
        )
        onepart = F * B10.si.value
        full_sum = np.sum(onepart / ν, axis=-1)

        output = (
            3
            / 2
            * self.thermal_speed.reshape(1, -1) ** 2
            * x.reshape(-1, 1) ** 2
            * full_sum
        )
        if with_B10:
            return output, B10
        else:
            return output

    @validate_quantities
    def K(
        self,
        x: np.ndarray,
        flux_surface: FluxSurface,
        *,
        m_max: int = 100,
        orbit_squeezing: bool = False,
    ) -> u.s**-1:
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
        flux_surface : FluxSurface
        m_max : int
        orbit_squeezing : bool
            orbit_squeezing
        """
        # Eq 16
        kb = self.K_B_ai(
            x, flux_surface.trapped_fraction, orbit_squeezing=orbit_squeezing
        )
        kps = self.K_ps_ai(x, flux_surface, m_max=m_max)
        return 1 / (1 / kb + 1 / kps)

    @validate_quantities
    def mu_hat(
        self,
        flux_surface: FluxSurface,
        *,
        xmin: float = 0.000015,
        xmax: float = 5,
        N: Optional[int] = None,
        **kwargs,
    ) -> u.Unit("kg / m3 / s"):
        """Viscosity coefficients - equation 15 from |Houlberg_1997|.

        Parameters
        ----------
        flux_surface : FluxSurface
        xmin :
            xmin
        xmax :
            xmax
        N :
            N 
        kwargs :
            kwargs

        Notes
        =====

        The array shapes in this function go as follows:
            (species, x_grid, alpha, beta)
        """
        # dirty hack because validate_quantities does not play well with **kwargs?
        if "kwargs" in kwargs:
            kwargs = kwargs["kwargs"]

        if N is None:
            N = 1000
        orders = np.arange(1, 4)
        π = np.pi
        x = np.linspace(xmin, xmax, N)

        # signs = np.array([
        #     [ 1, 1, -1],
        #     [ 1, 1,  1],
        #     [-1, 1,  1],
        # ])
        α = orders[:, None]
        β = orders[None, :]
        signs = (-1) ** (α + β)
        laguerres = np.vstack([LaguerrePolynomials[o - 1](x ** 2) for o in orders])
        kterm = self.K(x, flux_surface, **kwargs).T[..., np.newaxis, np.newaxis]
        xterm = (x ** 4 * np.exp(-(x ** 2)))[np.newaxis, :, np.newaxis, np.newaxis]
        # TODO try quadgk here, it should work decently!
        y = (
            laguerres.reshape(1, N, 3, 1)
            * laguerres.reshape(1, N, 1, 3)
            * kterm
            * xterm
        )
        integral = trapezoid(y, x, axis=1)
        mu_hat_ai = integral * signs[None, ...]
        actual_units = (8 / 3 / np.sqrt(π)) * mu_hat_ai * self.mass_density[:, None, None]
        return actual_units


def ωm(x: np.ndarray, m: Union[int, np.ndarray], thermal_speed: u.m / u.s, gamma):
    """Equation B11 from |Houlberg_1997|.

    Parameters
    ----------
    x : np.ndarray
        distribution velocity relative to the thermal velocity.
    m : Union[int, np.ndarray]
    thermal_speed : u.m/u.s
    """
    """Equation B11 of Houlberg_1997."""
    B11 = (
        x.reshape(-1, 1, 1)
        * thermal_speed.reshape(1, -1, 1)
        * m.reshape(1, 1, -1)
        * gamma
    )  # TODO why the u.m?
    return B11
