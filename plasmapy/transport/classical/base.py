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
    An abstract class representing classical transport coefficients. Subclasses
    representing different classical transport models re-implement
    the transport cofficient methods of this abstract class.

    The class can be initialized with one of two sets of keywords, allowing
    it to calcualte the dimensionless/normalized forms of the transport
    coefficients or the dimensional forms. If sets are provided, the dimensionless
    set takes precidence.

    Parameters (dimensionless)
    --------------------------

    chi_e : `numpy.ndarray` (N,)
        The electron Hall parameter

    chi_i : `numpy.ndarray` (N,)
        The ion Hall parameter

    Z : float
        The plasma mean ionization

    Parameters (dimensional)
    --------------------------

    particle : `~plasmapy.particles.Particle` instance or str
        The ion species

    B : `~astropy.units.Quantity`
        Magnetic field strength in units convertable to Tesla.

    ne : `~astropy.units.Quantity`
        Electron number density, in units convertable to cm^-3.

    ni : `~astropy.units.Quantity`
        Ion number density, in units convertable to cm^-3.

    Te : `~astropy.units.Quantity`
        Electron temperature, in units convertable to eV.

    Ti : `~astropy.units.Quantity`
        Ion temperature, in units convertable to eV.

    e_collision_freq : `~astropy.units.Quantity`
        The fundamental electron collision frequency, in units convertable
        to Hz. If not provided, this will be calculated using
        `~plasmapy.formulary.collisions.fundamental_electron_collision_freq()`

    i_collision_freq : `~astropy.units.Quantity`
        The fundamental ion collision frequency, in units convertable
        to Hz. If not provided, this will be calculated using
        `~plasmapy.formulary.collisions.fundamental_ion_collision_freq()`

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

        # The __init__ method calls one of two constructors based on which
        # keywords are provided.
        # The normalized/dimensionless keywords take precidence

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
    @property
    def norm_alpha_para(self):
        raise NotImplementedError

    @property
    def alpha_para(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to " "calculate alpha_para"
            )

        return self.norm_alpha_para * self.alpha_normalization

    @property
    def norm_alpha_perp(self):
        raise NotImplementedError

    @property
    def alpha_perp(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to " "calculate alpha_perp"
            )

        return self.norm_alpha_perp * self.alpha_normalization

    @property
    def norm_alpha_cross(self):
        raise NotImplementedError

    @property
    def alpha_cross(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to " "calculate alpha_cross"
            )
        return self.norm_alpha_cross * self.alpha_normalization

    # **********************************************************************
    # Thermoelectric Coefficient (beta)
    # **********************************************************************
    @property
    def norm_beta_para(self):
        raise NotImplementedError

    @property
    def beta_para(self):
        return self.norm_beta_para * self.beta_normalization

    @property
    def norm_beta_perp(self):
        raise NotImplementedError

    @property
    def beta_perp(self):
        return self.norm_beta_perp * self.beta_normalization

    @property
    def norm_beta_cross(self):
        raise NotImplementedError

    @property
    def beta_cross(self):
        return self.norm_beta_cross * self.beta_normalization

    # **********************************************************************
    # Electron Thermal Conductivity (kappa_e)
    # **********************************************************************
    @property
    def norm_kappa_e_para(self):
        raise NotImplementedError

    @property
    def kappa_e_para(self):
        if self.kappa_e_normalization is None:
            raise ValueError(
                "Keywords ne, Te, and B must be provided to " "calculate kappa_e_para"
            )

        return self.norm_kappa_e_para * self.kappa_e_normalization

    @property
    def norm_kappa_e_perp(self):
        raise NotImplementedError

    @property
    def kappa_e_perp(self):
        if self.kappa_e_normalization is None:
            raise ValueError(
                "Keywords ne, Te, and B must be provided to " "calculate kappa_perp_e"
            )

        return self.norm_kappa_e_perp * self.kappa_e_normalization

    @property
    def norm_kappa_e_cross(self):
        raise NotImplementedError

    @property
    def kappa_e_cross(self):
        if self.kappa_e_normalization is None:
            raise ValueError(
                "Keywords ne, Te, and B must be provided to " "calculate kappa_e_cross"
            )
        return self.norm_kappa_e_cross * self.kappa_e_normalization

    # **********************************************************************
    # Ion Thermal Conductivity (kappa_i)
    # **********************************************************************

    @property
    def norm_kappa_i_para(self):
        raise NotImplementedError

    @property
    def kappa_i_para(self):
        if self.kappa_i_normalization is None:
            raise ValueError(
                "Keywords particle, B, ni, and Ti must be provided to "
                "calculate kappa_i_para"
            )

        return self.norm_kappa_i_para * self.kappa_i_normalization

    @property
    def norm_kappa_i_perp(self):
        raise NotImplementedError

    @property
    def kappa_i_perp(self):
        if self.kappa_i_normalization is None:
            raise ValueError(
                "Keywords particle, B, ni, and Ti must be provided to "
                "calculate kappa_i_perp"
            )

        return self.norm_kappa_i_perp * self.kappa_i_normalization

    @property
    def norm_kappa_i_cross(self):
        raise NotImplementedError

    @property
    def kappa_i_cross(self):
        if self.kappa_i_normalization is None:
            raise ValueError(
                "Keywords particle, B, ni, and Ti must be provided to "
                "calculate kappa_i_cross"
            )
        return self.norm_kappa_i_cross * self.kappa_i_normalization

    # **********************************************************************
    # Electron Viscosity (pi_e)
    # **********************************************************************

    # TODO: implement electron viscosity

    # **********************************************************************
    # Electron Viscosity (pi_i)
    # **********************************************************************

    # TODO: implement ion viscosity

    # **********************************************************************
    # Resistive Velocity (delta_e)
    # "Symmetric" coefficent formulism of Sadler and Davies
    # **********************************************************************

    def norm_delta_e_perp(self):

        return self.norm_alpha_cross / self.chi_e

    def delta_e_perp(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to " "calculate delta_e_perp"
            )

        return self.norm_delta_e_perp * self.alpha_normalization

    def norm_delta_e_cross(self):
        return (self.norm_alpha_perp - self.norm_alpha_para) / self.chi_e

    def delta_e_cross(self):
        if self.alpha_normalization is None:
            raise ValueError(
                "Keywords ne and B must be provided to " "calculate delta_e_cross"
            )
        return self.norm_delta_e_cross * self.alpha_normalization

    # **********************************************************************
    # Nernst Coefficient (gamma_e)
    # "Symmetric" coefficent formulism of Sadler and Davies
    # **********************************************************************

    def norm_gamma_e_perp(self):
        return self.norm_beta_cross / self.chi_e

    def gamma_e_perp(self):
        return self.norm_gamma_e_perp * self.beta_normalization

    def norm_gamma_e_cross(self):
        return (self.norm_beta_para - self.norm_beta_perp) / self.chi_e

    def gamma_e_cross(self):
        return self.norm_gamma_e_cross * self.beta_normalization


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

        # Coefficients whose fits depend on chi_e
        e_coefficients = [
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

        # Coefficients whose fits depend on chi_i
        i_coefficients = [
            "kappa_i_para",
            "kappa_i_perp",
            "kappa_i_cross",
        ]

        # Create an interpolator for each of the data tables
        # using either the chi_e or chi_i tables as appropriate for that coefficient
        # (All of the interpolators use the same Z)
        self.interpolators = {}
        for parameter, coefficients in zip(
            ["chi_e", "chi_i"], [e_coefficients, i_coefficients]
        ):
            for coef in coefficients:
                if coef in list(file.keys()):
                    self.interpolators[coef] = interp2d(
                        file["Z"],
                        file[parameter],
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

    def norm_kappa_e_para(self):
        return self.interpolators["kappa_e_para"](self.Z, self.chi_e)

    def norm_kappa_e_perp(self):
        return self.interpolators["kappa_e_perp"](self.Z, self.chi_e)

    def norm_kappa_e_cross(self):
        return self.interpolators["kappa_e_cross"](self.Z, self.chi_e)

    def norm_kappa_i_para(self):
        return self.interpolators["kappa_i_para"](self.Z, self.chi_i)

    def norm_kappa_i_perp(self):
        return self.interpolators["kappa_i_perp"](self.Z, self.chi_i)

    def norm_kappa_i_cross(self):
        return self.interpolators["kappa_i_cross"](self.Z, self.chi_i)


if __name__ == "__main__":
    pass
