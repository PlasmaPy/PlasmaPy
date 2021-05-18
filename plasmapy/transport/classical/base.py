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
            missing = []
            
            for p in properties:
                if (getattr(self, p) is None):
                    missing.append(p)

            if len(missing) > 0:
                raise ValueError(
                    f"Keywords {properties} must be provided to calculate "
                    f"{self.__class__.__name__}.{func.__name__}(), but "
                    f"the following keywords were not provided {missing}."
                    )
            else:
                return func(self, *args, **kwargs)

        return wrapper

    return decorator



# This format string is used for cleaner NotImplemented errors for
# coefficients that are not included in a particular model
not_implemented = ("The {0} property is not implemented in {1}. "
                   "The {0} coefficient may not be defined as part of this model.")

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
        
        self.chi_e = chi_e
        self.chi_i = chi_i
        self.Z = Z
        self.particle = particle
        self.B = B
        self.ne = ne
        self.ni = ni
        self.Te = Te
        self.Ti = Ti
        self.e_collision_freq = e_collision_freq
        self.i_collision_freq = i_collision_freq

        

        
        
    
        # Ensure that only one set of keywords is being used
        # This code block checks to ensure that one of the sets of
        # keywords is all None, and raises an exception otherwise
        dim = [self.particle, self.B, self.ne, self.ni, 
                               self.Te, self.Ti, self.e_collision_freq, 
                               self.i_collision_freq]
        nodim = [self.chi_e, self.chi_i, self.Z]
        dim_none = all([v is None for v in dim])
        nodim_none = all([v is None for v in nodim])
        
        if not dim_none and not nodim_none:
            raise ValueError("The dimensionless keywords [chi_e, chi_i, Z] and "
                             "dimensional keywords (all others) are mutually "
                             "exlusive because, in the dimensional case, "
                             "chi and Z are calculated from the other parameters. "
                             "Use one or the other!")
        
    
        # The __init__ method calls one of two constructors based on which
        # keywords are provided.
        if Z is not None:
            self._constructor_normalized()

        else:
            self._constructor_dimensional()

    def _constructor_normalized(self):
        # Set the normalization coefficients to None, since the required
        # plasma parameters were not provided. This will cause the
        # dimensional coefficent functions to raise an error
        self.alpha_normalization = None
        self.beta_normalization = None
        self.kappa_e_normalization = None
        self.kappa_i_normalization = None
        self.eta_e_normalization = None
        self.eta_i_normalization = None

    def _constructor_dimensional(self):
        # Normalizations are defined such that multiplying the dimensionless
        # quantity by the normalization constant will return the dimensional
        # quantity

        self.Z = self.particle.integer_charge
    
        if all(v is not None for v in [self.ne, self.Te]):

            if self.e_collision_freq is None:
                self.e_collision_freq = fundamental_electron_collision_freq(self.Te, self.ne, self.particle)
            wce = gyrofrequency(self.B, "e-")

            self.chi_e = wce / self.e_collision_freq
            self.alpha_normalization = self.particle.mass * self.ne * self.e_collision_freq
            self.beta_normalization = 1.0
            self.kappa_e_normalization = self.ne * self.Te / self.e_collision_freq / m_e
            self.eta_e_normalization = self.ne * self.Te / self.e_collision_freq
        else:
            self.chi_e = None
            self.alpha_normalization = None
            self.beta_normalization = None
            self.kappa_e_normalization = None
            self.eta_e_normalization = None

        if all(v is not None for v in [self.ni, self.Ti]):

            if self.i_collision_freq is None:
                self.i_collision_freq = fundamental_ion_collision_freq(self.Ti, self.ni, self.particle)
            wci = gyrofrequency(self.B, self.particle)

            self.chi_i = wci / self.i_collision_freq
            self.kappa_i_normalization = self.ni * self.Ti / self.i_collision_freq / self.particle.mass
            self.eta_i_normalization = self.ni * self.Ti / self.i_collision_freq
        else:
            self.chi_i = None
            self.kappa_i_normalization = None
            self.eta_i_normalization = None

    # **********************************************************************
    # Resistivity (alpha)
    # **********************************************************************
    @property
    def norm_alpha_para(self):
        raise NotImplementedError(not_implemented.format("alpha_para",
                                                         self.__class__.__name__))

    @property
    @validate_object(properties=['ne', 'Te', 'B', 'particle'])
    def alpha_para(self):
        return self.norm_alpha_para * self.alpha_normalization

    @property
    def norm_alpha_perp(self):
        raise NotImplementedError(not_implemented.format("alpha_perp",
                                                         self.__class__.__name__))

    @property
    @validate_object(properties=['ne', 'Te', 'B', 'particle'])
    def alpha_perp(self):
        return self.norm_alpha_perp * self.alpha_normalization

    @property
    def norm_alpha_cross(self):
        raise NotImplementedError(not_implemented.format("alpha_cross",
                                                         self.__class__.__name__))

    @property
    @validate_object(properties=['ne', 'Te', 'B', 'particle'])
    def alpha_cross(self):
        return self.norm_alpha_cross * self.alpha_normalization

    # **********************************************************************
    # Thermoelectric Coefficient (beta)
    # **********************************************************************
    @property
    def norm_beta_para(self):
        raise NotImplementedError(not_implemented.format("beta_para",
                                                         self.__class__.__name__))

    @property
    def beta_para(self):
        return self.norm_beta_para * self.beta_normalization

    @property
    def norm_beta_perp(self):
        raise NotImplementedError(not_implemented.format("beta_perp",
                                                         self.__class__.__name__))

    @property
    def beta_perp(self):
        return self.norm_beta_perp * self.beta_normalization

    @property
    def norm_beta_cross(self):
        raise NotImplementedError(not_implemented.format("beta_cross",
                                                         self.__class__.__name__))

    @property
    def beta_cross(self):
        return self.norm_beta_cross * self.beta_normalization

    # **********************************************************************
    # Electron Thermal Conductivity (kappa_e)
    # **********************************************************************
    @property
    def norm_kappa_e_para(self):
        raise NotImplementedError(not_implemented.format("kappa_e_para",
                                                         self.__class__.__name__))

    @property
    @validate_object(properties=['ne', 'Te', 'B', 'particle'])
    def kappa_e_para(self):
        return self.norm_kappa_e_para * self.kappa_e_normalization

    @property
    def norm_kappa_e_perp(self):
        raise NotImplementedError(not_implemented.format("kappa_e_perp",
                                                         self.__class__.__name__))

    @property
    @validate_object(properties=['ne', 'Te', 'B', 'particle'])
    def kappa_e_perp(self):
        return self.norm_kappa_e_perp * self.kappa_e_normalization

    @property
    def norm_kappa_e_cross(self):
        raise NotImplementedError(not_implemented.format("kappa_e_cross",
                                                         self.__class__.__name__))

    @property
    @validate_object(properties=['ne', 'Te', 'B', 'particle'])
    def kappa_e_cross(self):
        return self.norm_kappa_e_cross * self.kappa_e_normalization

    # **********************************************************************
    # Ion Thermal Conductivity (kappa_i)
    # **********************************************************************

    @property
    def norm_kappa_i_para(self):
        raise NotImplementedError(not_implemented.format("kappa_i_para",
                                                         self.__class__.__name__))

    @property
    @validate_object(properties=['ni', 'Ti', 'B', 'particle'])
    def kappa_i_para(self):
        return self.norm_kappa_i_para * self.kappa_i_normalization

    @property
    def norm_kappa_i_perp(self):
        raise NotImplementedError(not_implemented.format("kappa_i_perp",
                                                         self.__class__.__name__))

    @property
    @validate_object(properties=['ni', 'Ti', 'B', 'particle'])
    def kappa_i_perp(self):
        return self.norm_kappa_i_perp * self.kappa_i_normalization

    @property
    def norm_kappa_i_cross(self):
        raise NotImplementedError(not_implemented.format("kappa_i_cross",
                                                         self.__class__.__name__))

    @property
    @validate_object(properties=['ni', 'Ti', 'B', 'particle'])
    def kappa_i_cross(self):   
        return self.norm_kappa_i_cross * self.kappa_i_normalization

    # **********************************************************************
    # Electron Viscosity (eta_e)
    # Note: returns an array [eta0, eta1, eta2, eta3, eta4]
    # **********************************************************************
    @property
    def norm_eta_e(self):
        raise NotImplementedError(not_implemented.format("eta_e",
                                                         self.__class__.__name__))
        
    @property
    @validate_object(properties=['ne', 'Te', 'B', 'particle'])
    def eta_e(self):
        return self.norm_eta_e * self.eta_e_normalization
    

    # **********************************************************************
    # Ion Viscosity (eta_i)
    # Note: returns an array [eta0, eta1, eta2, eta3, eta4]
    # **********************************************************************
    @property
    def norm_eta_i(self):
        raise NotImplementedError(not_implemented.format("eta_i",
                                                         self.__class__.__name__))
        
    @property
    @validate_object(properties=['ni', 'Ti', 'B', 'particle'])
    def eta_i(self):
        return self.norm_eta_i * self.eta_i_normalization

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
            "eta0_e",
            "eta1_e",
            "eta2_e",
            "eta3_e",
            "eta4_e",
        ]

        # Coefficients whose fits depend on chi_i
        i_coefficients = [
            "kappa_i_para",
            "kappa_i_perp",
            "kappa_i_cross",
            "eta0_i",
            "eta1_i",
            "eta2_i",
            "eta3_i",
            "eta4_i",
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
    @property
    def norm_alpha_para(self):
        return self.interpolators["alpha_para"](self.Z, self.chi_e)
    
    @property
    def norm_alpha_perp(self):
        return self.interpolators["alpha_perp"](self.Z, self.chi_e)
    
    @property
    def norm_alpha_cross(self):
        return self.interpolators["alpha_cross"](self.Z, self.chi_e)

    @property
    def norm_beta_para(self):
        return self.interpolators["beta_para"](self.Z, self.chi_e)
   
    @property
    def norm_beta_perp(self):
        return self.interpolators["beta_perp"](self.Z, self.chi_e)

    @property
    def norm_beta_cross(self):
        return self.interpolators["beta_cross"](self.Z, self.chi_e)
    
    @property
    def norm_kappa_e_para(self):
        return self.interpolators["kappa_e_para"](self.Z, self.chi_e)
    
    @property
    def norm_kappa_e_perp(self):
        return self.interpolators["kappa_e_perp"](self.Z, self.chi_e)

    @property
    def norm_kappa_e_cross(self):
        return self.interpolators["kappa_e_cross"](self.Z, self.chi_e)
    
    @property
    def norm_kappa_i_para(self):
        return self.interpolators["kappa_i_para"](self.Z, self.chi_i)

    @property
    def norm_kappa_i_perp(self):
        return self.interpolators["kappa_i_perp"](self.Z, self.chi_i)
    
    @property
    def norm_kappa_i_cross(self):
        return self.interpolators["kappa_i_cross"](self.Z, self.chi_i)
    
    @property
    def norm_eta_e(self):
        """
        
        Returns
        -------
        coef : `numpy.ndarray` (5,)
            The five electron viscosity coefficients:
                 eta0, eta1, eta2, eta3, eta4.
        """
        
        coef = []
        names = ['eta0_e', 'eta1_e', 'eta2_e', 'eta3_e', 'eta4_e']
        for c in names:
            coef.append(self.interpolators[c](self.Z, self.chi_e))
        return coef
    
    @property
    def norm_eta_i(self):
        """
        
        Returns
        -------
        coef : `numpy.ndarray` (5,)
            The five ion viscosity coefficients:
                 eta0, eta1, eta2, eta3, eta4.
        """
        
        coef = []
        names = ['eta0_i', 'eta1_i', 'eta2_i', 'eta3_i', 'eta4_i']
        for c in names:
            coef.append(self.interpolators[c](self.Z, self.chi_i))
        return coef


if __name__ == "__main__":
    chi = np.linspace(-2, 2, num=50)
    chi = 10 ** chi
    x = AbstractClassicalTransportCoefficients(chi_e = chi, Z=1)
    print(x.norm_alpha_para)

    
