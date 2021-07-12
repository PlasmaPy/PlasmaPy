from functools import cached_property
import astropy.units as u
import numpy as np

rho = u.def_unit("rho")


class SimpleFluxSurface:
    def __init__(
        self,
        inverse_aspect_ratio: float,
        safety_factor: float,
        major_radius: u.m,
        minor_radius: u.m,
        axial_elongation: float,
        axial_toroidal_field: u.T,
        axial_safety_factor: float,
        q0: float,
        radial_electric_field: u.V / rho = None,
        radial_electric_field_gradient: u.V / rho ** 2 = None,
    ):
        self.p_eps = inverse_aspect_ratio
        self.p_grphi = radial_electric_field
        self.p_gr2phi = radial_electric_field_gradient
        # gr2phi is defined as Psi' (Phi'/Psi')'
        self.p_q = safety_factor
        self.r0 = major_radius  # think this is not the machine major radius but rather the flux surface major radius
        self.a0 = minor_radius
        self.e0 = axial_elongation
        self.bt0 = axial_toroidal_field
        self.q0 = axial_safety_factor

    @cached_property
    def fsa_B2(self):
        return self.bt0 ** 2 * (1 + 0.5 * self.p_eps ** 2)

    @cached_property
    def fsa_invB2(self):
        return (1.0 + 1.5 * self.p_eps ** 2) / self.bt0 ** 2

    def F_m(self, M=3):
        p_eps = self.p_eps
        i = np.arange(1, M + 1)
        C1 = np.sqrt(1 - p_eps ** 2)
        return (
            i
            * ((1 - C1) / p_eps) ** (2 * i)
            * (1 + i * C1)
            / (C1 ** 3 * (self.p_q * self.r0) ** 2)
        )

    @cached_property
    def trapped_fraction(self):
        return 1.46 * np.sqrt(self.p_eps)

    @cached_property
    def grbm2(self):
        # <grad(ρ)**2 / B **2>    #ρ**2 / m**2 / T**2
        return 1 / self.bt0 ** 2

    @cached_property
    def gamma(self):
        return 1.0 / (self.p_q * self.r0)

    @cached_property
    def F_m3(self):
        return self.F_m()
