__all__ = [
    "AbstractClassicalTransportCoefficients",
    "AbstractInterpolatedCoefficients",
    "validate_object",
]

import astropy.constants as const
import astropy.units as u
import functools
import numpy as np

from abc import ABC, abstractmethod
from scipy.interpolate import interp2d

from plasmapy import particles
from plasmapy.formulary.collisions import (
    fundamental_electron_collision_freq,
    fundamental_ion_collision_freq,
)
from plasmapy.formulary.parameters import gyrofrequency
from plasmapy.particles import Particle

m_e = const.m_e.si


def validate_object(properties=[]):
    """
    This decorator for methods of a class validates that a list of
    attributes on that class are not None.

    Expand to a real separate decorator in the util package?
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            for p in properties:
                if getattr(self, p) is None:
                    raise ValueError(
                        f"{self.__class__.__name__}.{func.__name__}() "
                        f"cannot be calculated when the attribute {p} "
                        "is None."
                    )
            return func(self, *args, **kwargs)

        return wrapper

    return decorator


class AbstractClassicalTransportCoefficients(ABC):
    r"""
    An abstract class representing classical transport coefficients.

    Subclasses representing different classical transport models re-implement
    the transport cofficient methods of this abstract class.
    """

    @particles.particle_input(none_shall_pass=True)
    def __init__(
        self,
        chi_e=None,
        chi_i=None,
        Z=None,
        particle: Particle = None,
        B: u.T = None,
        ne: u.m ** -3 = None,
        ni: u.m ** -3 = None,
        Te: u.K = None,
        Ti: u.K = None,
        e_collision_freq=None,
        i_collision_freq=None,
    ):

        if Z is not None:
            self._constructor_normalized(chi_e, chi_i, Z)

        else:
            self._constructor_dimensional(
                particle=particle,
                B=B,
                ne=ne,
                ni=ni,
                Te=Te,
                Ti=Ti,
                e_collision_freq=e_collision_freq,
                i_collision_freq=i_collision_freq,
            )

    def _constructor_normalized(self, chi_e, chi_i, Z):
        self.chi_e = chi_e
        self.chi_i = chi_i
        self.Z = Z

        # Set the normalization coefficients to None, since the required
        # plasma parameters were not provided. This will cause the
        # dimensional coefficent functions to raise an error
        self.alpha_normalization = None
        self.beta_normalization = None
        self.kappa_e_normalization = None
        self.kappa_i_normalization = None

    def _constructor_dimensional(
        self,
        particle,
        B,
        ne=None,
        ni=None,
        Te=None,
        Ti=None,
        e_collision_freq=None,
        i_collision_freq=None,
    ):

        # Normalizations are defined such that multiplying the dimensionless
        # quantity by the normalization constant will return the dimensional
        # quantity

        self.Z = particle.integer_charge

        if all(v is not None for v in [ne, Te]):

            if e_collision_freq is None:
                e_collision_freq = fundamental_electron_collision_freq(Te, ne, particle)
            wce = gyrofrequency(B, "e-")

            self.chi_e = wce / e_collision_freq
            self.alpha_normalization = particle.mass * ne * e_collision_freq
            self.beta_normalization = 1.0
            self.kappa_e_normalization = ne * Te / e_collision_freq / m_e
        else:
            self.chi_e = None
            self.alpha_normalization = None
            self.beta_normalization = None
            self.kappa_e_normalization = None

        if all(v is not None for v in [ni, Ti]):

            if i_collision_freq is None:
                i_collision_freq = fundamental_ion_collision_freq(Ti, ni, particle)
            wci = gyrofrequency(B, particle)

            self.chi_i = wci / i_collision_freq
            self.kappa_i_normalization = ni * Ti / i_collision_freq / particle.mass
        else:
            self.chi_i = None
            self.kappa_i_normalization = None

    # **********************************************************************
    # Resistivity (alpha)
    # **********************************************************************
    def norm_alpha_para(self):
        raise NotImplementedError

    def alpha_para(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to " "calculate alpha_para"
            )

        return self.norm_alpha_para() * self.alpha_normalization

    def norm_alpha_perp(self):
        raise NotImplementedError

    def alpha_perp(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to " "calculate alpha_perp"
            )

        return self.norm_alpha_perp() * self.alpha_normalization

    def norm_alpha_cross(self):
        raise NotImplementedError

    def alpha_cross(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to " "calculate alpha_cross"
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
                "Keywords ne, Te, and B must be provided to " "calculate kappa_e_para"
            )

        return self.norm_kappa_e_para() * self.kappa_e_normalization

    def norm_kappa_e_perp(self):
        raise NotImplementedError

    def kappa_e_perp(self):
        if self.kappa_e_normalization is None:
            raise ValueError(
                "Keywords ne, Te, and B must be provided to " "calculate kappa_perp_e"
            )

        return self.norm_kappa_e_perp() * self.kappa_e_normalization

    def norm_kappa_e_cross(self):
        raise NotImplementedError

    def kappa_e_cross(self):
        if self.kappa_e_normalization is None:
            raise ValueError(
                "Keywords ne, Te, and B must be provided to " "calculate kappa_e_cross"
            )
        return self.norm_kappa_e_cross() * self.kappa_e_normalization

    # **********************************************************************
    # Ion Thermal Conductivity (kappa_i)
    # **********************************************************************

    def norm_kappa_i_para(self):
        raise NotImplementedError

    def kappa_i_para(self):
        if self.kappa_i_normalization is None:
            raise ValueError(
                "Keywords particle, B, ni, and Ti must be provided to "
                "calculate kappa_i_para"
            )

        return self.norm_kappa_i_para() * self.kappa_i_normalization

    def norm_kappa_i_perp(self):
        raise NotImplementedError

    def kappa_i_perp(self):
        if self.kappa_i_normalization is None:
            raise ValueError(
                "Keywords particle, B, ni, and Ti must be provided to "
                "calculate kappa_i_perp"
            )

        return self.norm_kappa_i_perp() * self.kappa_i_normalization

    def norm_kappa_i_cross(self):
        raise NotImplementedError

    def kappa_i_cross(self):
        if self.kappa_i_normalization is None:
            raise ValueError(
                "Keywords particle, B, ni, and Ti must be provided to "
                "calculate kappa_i_cross"
            )
        return self.norm_kappa_i_cross() * self.kappa_i_normalization


class AbstractInterpolatedCoefficients(AbstractClassicalTransportCoefficients):
    """
    Interpolates transport coefficients from arrays of calculated values
    """

    @property
    @abstractmethod
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
            "kappa_e_para",
            "kappa_e_perp",
            "kappa_e_cross",
        ]

        # Create an interpolator for each of the data tables
        # (All of the interpolators use the same chi and Z)
        self.interpolators = {}
        for coef in coefficients:
            self.interpolators[coef] = interp2d(
                file["Z"],
                file["chi_e"],
                file[coef],
                kind="cubic",
                bounds_error=False,
                fill_value=None,
            )

    def norm_alpha_para(self):
        return self.interpolators["alpha_para"](self.Z, self.chi_e)

    def norm_alpha_perp(self):
        return self.interpolators["alpha_perp"](self.Z, self.chi_e)

    def norm_alpha_cross(self):
        return self.interpolators["alpha_cross"](self.Z, self.chi_e)

    def norm_beta_para(self):
        return self.interpolators["beta_para"](self.Z, self.chi_e)

    def norm_beta_perp(self):
        return self.interpolators["beta_perp"](self.Z, self.chi_e)

    def norm_beta_cross(self):
        return self.interpolators["beta_cross"](self.Z, self.chi_e)

    def norm_kappa_para(self):
        return self.interpolators["kappa_e_para"](self.Z, self.chi_e)

    def norm_kappa_perp(self):
        return self.interpolators["kappa_e_perp"](self.Z, self.chi_e)

    def norm_kappa_cross(self):
        return self.interpolators["kappa_e_cross"](self.Z, self.chi_e)


if __name__ == "__main__":
    pass
