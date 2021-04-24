__all__ = [
    "AbstractClassicalTransportCoefficients",
    "AbstractInterpolatedCoefficients",
]

import astropy.constants as const
import astropy.units as u
import numpy as np
import os

from abc import ABC
from collections import namedtuple
from scipy.interpolate import interp2d

from plasmapy import particles
from plasmapy.formulary.parameters import gyrofrequency
from plasmapy.particles import Particle


m_e = const.m_e.si

class AbstractClassicalTransportCoefficients(ABC):
    r"""
    An abstract class representing classical transport coefficients.

    Subclasses representing different classical transport models re-implement
    the transport cofficient methods of this abstract class.
    """

    # @particles.particle_input
    def __init__(
        self,
        chi,
        Z,
        B: u.T = None,
        ne: u.m ** -3 = None,
        ni: u.m ** -3 = None,
        Te: u.K = None,
        Ti: u.K = None,
        particle: Particle = "",
    ):

        self.chi = chi
        self.Z = Z

        self.B = B
        self.ne = ne
        self.ni = ni
        self.Te = Te
        self.Ti = Ti
        self.particle = particle

        try:
            self.ei_collision_period = chi / gyrofrequency(B, "e-")
        except ValueError:
            self.ei_collision_period = None

        if Z < 0:
            raise ValueError(f"Z > 0 required: supplied value was {Z}")

        # Normalizations are defined such that multiplying the dimensionless
        # quantity by the normalization constant will return the dimensional
        # quantity
        try:
            self.alpha_normalization = (
                self.particle.mass * self.n_e / self.ei_collision_period
            )
        except (TypeError, AttributeError):
            self.alpha_normalization = None

        self.beta_normalization = 1.0  # Beta is naturally dimensionless

        try:
            self.kappa_normalization = (
                self.n_e * self.T_e * self.ei_collision_period / m_e
            )
        except (TypeError, AttributeError):
            self.kappa_normalization = None

    # **********************************************************************
    # Resistivity (alpha)
    # **********************************************************************
    def norm_alpha_para(self):
        raise NotImplementedError

    def alpha_para(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to "
                "calculate alpha_para"
            )

        return self.norm_alpha_para() * self.alpha_normalization

    def norm_alpha_perp(self):
        raise NotImplementedError

    def alpha_perp(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to "
                "calculate alpha_perp"
            )

        return self.norm_alpha_perp() * self.alpha_normalization

    def norm_alpha_cross(self):
        raise NotImplementedError

    def alpha_cross(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to "
                "calculate alpha_cross"
            )
        return self.norm_alpha_cross() * self.alpha_normalization

    # **********************************************************************
    # Thermoelectric Coefficient (beta)
    # **********************************************************************
    def norm_beta_para(self):
        raise NotImplementedError

    def beta_para(self):
        return self.norm_beta_para() * self.beta_normalization

    def norm_beta_perp(self):
        raise NotImplementedError

    def beta_perp(self):
        return self.norm_beta_perp() * self.beta_normalization

    def norm_beta_cross(self):
        raise NotImplementedError

    def beta_cross(self):
        return self.norm_beta_cross() * self.beta_normalization

    # **********************************************************************
    # Electron Thermal Conductivity (kappa_e)
    # **********************************************************************
    def norm_kappa_e_para(self):
        raise NotImplementedError

    def kappa_e_para(self):
        if self.kappa_e_normalization is None:
            raise ValueError(
                "Keywords ne, Te, and B must be provided to "
                "calculate kappa_e_para"
            )

        return self.norm_kappa_e_para() * self.kappa_normalization

    def norm_kappa_perp(self):
        raise NotImplementedError

    def kappa_e_perp(self):
        if self.kappa_e_normalization is None:
            raise ValueError(
                "Keywords ne, Te, and B must be provided to "
                "calculate kappa_perp_e"
            )

        return self.norm_kappa_e_perp() * self.kappa_e_normalization

    def norm_kappa_e_cross(self):
        raise NotImplementedError

    def kappa_e_cross(self):
        if self.kappa_e_normalization is None:
            raise ValueError(
                "Keywords ne, Te, and B must be provided to "
                "calculate kappa_e_cross"
            )
        return self.norm_kappa_e_cross() * self.kappa_e_normalization
    
    
    # **********************************************************************
    # Ion Thermal Conductivity (kappa_e)
    # **********************************************************************
    
    # TODO: To implement these, need to figure out how to handle tau_i. 
    # Can't just calculate it from chi...


class AbstractInterpolatedCoefficients(AbstractClassicalTransportCoefficients):
    """
    Interpolates transport coefficients from arrays of calculated values
    """

    @property
    def _data_file(self):
        return None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        file = np.load(self._data_file)
        coefficients = [
            "alpha_para",
            "alpha_perp",
            "alpha_cross",
            "beta_para",
            "beta_perp",
            "beta_cross",
            "kappa_para",
            "kappa_perp",
            "kappa_cross",
        ]

        # Create an interpolator for each of the data tables
        # (All of the interpolators use the same chi and Z)
        self.interpolators = {}
        for coef in coefficients:
            self.interpolators[coef] = interp2d(
                file["Z"],
                file["chi"],
                file[coef],
                kind="cubic",
                bounds_error=False,
                fill_value=None,
            )

    def norm_alpha_para(self):
        return self.interpolators["alpha_para"](self.Z, self.chi)

    def norm_alpha_perp(self):
        return self.interpolators["alpha_perp"](self.Z, self.chi)

    def norm_alpha_cross(self):
        return self.interpolators["alpha_cross"](self.Z, self.chi)

    def norm_beta_para(self):
        return self.interpolators["beta_para"](self.Z, self.chi)

    def norm_beta_perp(self):
        return self.interpolators["beta_perp"](self.Z, self.chi)

    def norm_beta_cross(self):
        return self.interpolators["beta_cross"](self.Z, self.chi)

    def norm_kappa_para(self):
        return self.interpolators["kappa_para"](self.Z, self.chi)

    def norm_kappa_perp(self):
        return self.interpolators["kappa_perp"](self.Z, self.chi)

    def norm_kappa_cross(self):
        return self.interpolators["kappa_cross"](self.Z, self.chi)


if __name__ == "__main__":
    pass
