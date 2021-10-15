__all__ = [
    "AbstractClassicalTransportCoefficients",
    "AbstractPolynomialCoefficients",
    "AbstractInterpolatedCoefficients",
    "validate_object",
]

import astropy.constants as const
import astropy.units as u
import functools
import numpy as np
import warnings

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
                if getattr(self, p) is None:
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
not_implemented = (
    "The {0} property is not implemented in {1}. "
    "The {0} coefficient may not be defined as part of this model."
)


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
        ne: u.cm ** -3 = None,
        ni: u.cm ** -3 = None,
        Te: u.K = None,
        Ti: u.K = None,
        e_collision_freq: u.Hz=None,
        i_collision_freq: u.Hz=None,
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
        self._dim_vars = {
            'particle':self.particle,
            'B':self.B,
            'ne':self.ne,
            'ni':self.ni,
            'Te':self.Te,
            'Ti':self.Ti,
            'e_collision_freq':self.e_collision_freq,
            'i_collision_freq':self.i_collision_freq,
        }
        self._nodim_vars = {'chi_e':self.chi_e, 
                            'chi_i':self.chi_i, 
                            'Z':self.Z}
        dim_none = all([v is None for v in self._dim_vars.values()])
        nodim_none = all([v is None for v in self._nodim_vars.values()])

        # Enforce that only one set of coefficents or the other is used
        if not dim_none and not nodim_none:
            raise ValueError(
                "The dimensionless keywords [chi_e, chi_i, Z] and "
                "dimensional keywords (all others) are mutually "
                "exlusive because, in the dimensional case, "
                "chi and Z are calculated from the other parameters. "
                "Use one or the other!"
            )

        # The __init__ method calls one of two constructors based on which
        # keywords are provided.
        # If the normalized constructor is used, only the normalized
        # coefficients can be used.
        # If the dimensional constructor is used, both the normalized
        # and dimensional forms of the coefficients can be used.
        if dim_none:
            self._dimensional = False
            self._constructor_dimensionless()
            

        else:
            self._constructor_dimensional()
            self._dimensional = True
            
    def _constructor_dimensionless(self):
        
        # Verify that Z is a number
        # TODO: Ensure support for arrays of Z
        if not isinstance(self.Z, (int,float, np.ndarray)):
            raise ValueError("Z must be a dimensionless number or "
                             "`~numpy.ndarray`.")
            
            
        # TODO: handle chi_e and chi_i mistakenly entered as a quantity array,
        # with or without dimensional units.
            
        # chi_e and chi_i must be a number or None, but both cannot be None
        if not (isinstance(self.chi_e, np.ndarray) or self.chi_e is None):
            raise ValueError("chi_e must be either a dimensionless "
                             "`~numpy.ndarray` or None.")
            
        if not (isinstance(self.chi_i, np.ndarray) or self.chi_i is None):
            raise ValueError("chi_i must be either a dimensionless "
                             "`~numpy.ndarray` or None.")
            
        if (self.chi_e is None) and (self.chi_i is None):
            raise ValueError("chi_e and chi_i cannot both be None when "
                             "using the dimensionless mode.")
            


    def _constructor_dimensional(self):
        # particle and B keyword arguments are always required
        if self.particle is None:
            raise ValueError("The 'particle' keyword is required in dimensional "
                             "mode")
        
        if self.B is None:
            raise ValueError("The 'B' keyword is required in dimensional "
                             "mode.")
        
        
        # Validate that input arrays have the same shape
        shape = self.B.shape
        array_quantities = ['ne', 'ni', 'Te', 'Ti', 
                            'e_collision_freq',
                            'i_collision_freq',
                            ]
        for q in array_quantities:
            if self._dim_vars[q] is not None and  \
                self._dim_vars[q].shape != shape:
                    raise ValueError("Shape of arrays must be equal, but "
                                     f"B: {self.B.shape} and "
                                     f"{q}:{self._dim_vars[q].shape} are not.")
                    
        
     
        self.Z = self.particle.integer_charge
        # If all the electron quantities are provided, calculate chi_e
        if all(v is not None for v in [self.ne, self.Te]):

            if self.e_collision_freq is None:
                self.e_collision_freq = fundamental_electron_collision_freq(
                    self.Te, self.ne, self.particle
                )
            wce = gyrofrequency(self.B, "e-")
            self.chi_e = (wce / self.e_collision_freq).to(u.rad).value

        else:
            self.chi_e = None

        # If all the ion quantities are provided, calculate chi_i
        if all(v is not None for v in [self.ni, self.Ti]):

            if self.i_collision_freq is None:
                self.i_collision_freq = fundamental_ion_collision_freq(
                    self.Ti, self.ni, self.particle
                )
            wci = gyrofrequency(self.B, self.particle)

            self.chi_i = (wci / self.i_collision_freq).to(u.rad).value
        else:
            self.chi_i = None
            
    # **********************************************************************
    # Normalization Constants
    # (defaults, overwritable by child classes)
    # Normalizations are defined such that multiplying the dimensionless
    # quantity by the normalization constant will return the dimensional
    # quantity
    # **********************************************************************    
    @property
    def alpha_normalization(self):
        """
        The normalization constant for alpha. 
        
        Defined such that multiplying the dimensionless quantity by the
        normalization constnat will return the dimensional quantity.
        """
        if self._dimensional:
            return self.particle.mass * self.ne * self.e_collision_freq
        else:
            return None
    
    @property
    def beta_normalization(self):
        """
        The normalization constant for beta. 
        
        Defined such that multiplying the dimensionless quantity by the
        normalization constnat will return the dimensional quantity.
        """
        if self._dimensional:
            return 1.0
        else:
            return None
        
    @property
    def kappa_e_normalization(self):
        """
        The normalization constant for kappa_e. 
        
        Defined such that multiplying the dimensionless quantity by the
        normalization constnat will return the dimensional quantity.
        """
        if self._dimensional:
            return self.ne * self.Te / self.e_collision_freq / m_e
        else:
            return None
        
    @property
    def eta_e_normalization(self):
        """
        The normalization constant for eta_e. 
        
        Defined such that multiplying the dimensionless quantity by the
        normalization constnat will return the dimensional quantity.
        """
        if self._dimensional:
            return self.ne * self.Te / self.e_collision_freq
        else:
            return None
        
    @property
    def kappa_i_normalization(self):
        """
        The normalization constant for kappa_i. 
        
        Defined such that multiplying the dimensionless quantity by the
        normalization constnat will return the dimensional quantity.
        """
        if self._dimensional:
            return self.ni * self.Ti / self.i_collision_freq / self.particle.mass
        else:
            return None
        
    @property
    def eta_i_normalization(self):
        """
        The normalization constant for eta_i. 
        
        Defined such that multiplying the dimensionless quantity by the
        normalization constnat will return the dimensional quantity.
        """
        if self._dimensional:
            return self.ni * self.Ti / self.i_collision_freq
        else:
            return None

    # **********************************************************************
    # Resistivity (alpha)
    # **********************************************************************
    @property
    def norm_alpha(self):
        raise NotImplementedError(
            not_implemented.format("alpha", self.__class__.__name__)
        )

    @property
    @validate_object(properties=["ne", "Te", "B", "particle"])
    def alpha(self):
        return self.norm_alpha * self.alpha_normalization

    # **********************************************************************
    # Thermoelectric Coefficient (beta)
    # **********************************************************************
    @property
    def norm_beta(self):
        raise NotImplementedError(
            not_implemented.format("beta", self.__class__.__name__)
        )

    @property
    @validate_object(properties=[])
    def beta(self):
        return self.norm_beta * self.beta_normalization

    # **********************************************************************
    # Electron Thermal Conductivity (kappa_e)
    # **********************************************************************
    @property
    def norm_kappa_e(self):
        raise NotImplementedError(
            not_implemented.format("kappa_e", self.__class__.__name__)
        )

    @property
    @validate_object(properties=["ne", "Te", "B", "particle"])
    def kappa_e(self):
        # newaxis is necessary if ne, B etc. are not scalars
        return self.norm_kappa_e * self.kappa_e_normalization[..., np.newaxis]

    # **********************************************************************
    # Ion Thermal Conductivity (kappa_i)
    # **********************************************************************

    @property
    def norm_kappa_i(self):
        raise NotImplementedError(
            not_implemented.format("kappa_i", self.__class__.__name__)
        )

    @property
    @validate_object(properties=["ni", "Ti", "B", "particle"])
    def kappa_i(self):
        return self.norm_kappa_i * self.kappa_i_normalization[..., np.newaxis]

    # **********************************************************************
    # Electron Viscosity (eta_e)
    # Note: returns an array [eta0, eta1, eta2, eta3, eta4]
    # **********************************************************************
    @property
    def norm_eta_e(self):
        raise NotImplementedError(
            not_implemented.format("eta_e", self.__class__.__name__)
        )

    @property
    @validate_object(properties=["ne", "Te", "B", "particle"])
    def eta_e(self):
        return self.norm_eta_e * self.eta_e_normalization

    # **********************************************************************
    # Ion Viscosity (eta_i)
    # Note: returns an array [eta0, eta1, eta2, eta3, eta4]
    # **********************************************************************
    @property
    def norm_eta_i(self):
        raise NotImplementedError(
            not_implemented.format("eta_i", self.__class__.__name__)
        )

    @property
    @validate_object(properties=["ni", "Ti", "B", "particle"])
    def eta_i(self):
        return self.norm_eta_i * self.eta_i_normalization

    # **********************************************************************
    # Resistive Velocity (delta)
    # "Symmetric" coefficent formulism of Sadler and Davies
    # **********************************************************************
    @property
    def norm_delta(self):
        """
        The normalized symmetric transport coefficient delta

        Returns
        -------
        coef : `~numpy.ndarray` (2,)
            An array of
            [norm_delta_perp, norm_delta_cross]
        """
        perp = self.norm_alpha[2] / self.chi_e
        cross = (self.norm_alpha[1] - self.norm_alpha[0]) / self.chi_e

        return np.array([perp, cross])

    @property
    @validate_object(properties=["ne", "Te", "B", "particle"])
    def delta(self):
        return self.norm_delta * self.alpha_normalization

    # **********************************************************************
    # Nernst Coefficient (gamma)
    # "Symmetric" coefficent formulism of Sadler and Davies
    # **********************************************************************

    @property
    def norm_gamma(self):
        """
        The normalized symmetric transport coefficient gamma

        Returns
        -------
        coef : `~numpy.ndarray` (2,)
            An array of
            [norm_gamma_perp, norm_gamma_cross]
        """

        perp = self.norm_beta_cross / self.chi_e
        cross = (self.norm_beta_para - self.norm_beta_perp) / self.chi_e

        return np.array([perp, cross])

    @property
    @validate_object(properties=[])
    def gamma(self):
        return self.norm_gamma * self.beta_normalization


# TODO: Should allow for each coefficient to have it's own Z and chi list, so they are
# independent (maybe a default if one isn't specifically listed)
# Braginskii only includes formulas for Z=1 for a few of the coefficients


class AbstractPolynomialCoefficients(AbstractClassicalTransportCoefficients):
    @property
    @abstractmethod
    def _c(self):
        """
        Dictionary of polynomial coefficients

        """
        ...

    def _find_nearest_Z(self, Z):
        """
        Finds the nearest Z-value to the given Z value in the coefficient tables.
        Prints a warning if the Z found is not equal to the Z requested.

        Parameters
        ----------
        Z : float
            An integer charge

        Returns
        -------
        i : int
            The index of the closest Z in the tables

        """
        if Z == np.inf:
            return -1

        i = np.argmin(np.abs(self._c["Z"] - Z))
        if self._c["Z"][i] != Z:
            warnings.warn(
                f"Value Z = {Z} is not in the coefficient table. "
                f"Using the nearest value, Z = {self._c['Z'][i]}. "
                f"The values in the table are {self._c['Z']}.",
                RuntimeWarning,
            )
        return i


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

                # If a coeffiecient table exists, create an interpolator
                # if it doesn't, set that interpolator in the dict to None
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
    def norm_alpha(self):
        coef = []
        names = ["alpha_para", "alpha_perp", "alpha_cross"]

        # Raise an exception of all of the required coefficients are not
        # existant in this class
        if not set(names).issubset(set(self.interpolators.keys())):
            raise ValueError(
                "alpha coefficient is not defined for class "
                f"{self.__class__.__name__}."
            )

        for c in names:
            coef.append(self.interpolators[c](self.Z, self.chi_e))
        return np.array(coef)

    @property
    def norm_beta(self):
        coef = []
        names = ["beta_para", "beta_perp", "beta_cross"]

        # Raise an exception of all of the required coefficients are not
        # existant in this class
        if not set(names).issubset(set(self.interpolators.keys())):
            raise ValueError(
                "beta coefficient is not defined for class "
                f"{self.__class__.__name__}."
            )

        for c in names:
            coef.append(self.interpolators[c](self.Z, self.chi_e))
        return np.array(coef)

    @property
    def norm_kappa_e(self):
        coef = []
        names = ["kappa_e_para", "kappa_e_perp", "kappa_e_cross"]

        # Raise an exception of all of the required coefficients are not
        # existant in this class
        if not set(names).issubset(set(self.interpolators.keys())):
            raise ValueError(
                "kappa_e coefficient is not defined for class "
                f"{self.__class__.__name__}."
            )

        for c in names:
            coef.append(self.interpolators[c](self.Z, self.chi_e))
        return np.array(coef)

    @property
    def norm_kappa_i(self):
        coef = []
        names = ["kappa_i_para", "kappa_i_perp", "kappa_i_cross"]

        # Raise an exception of all of the required coefficients are not
        # existant in this class
        if not set(names).issubset(set(self.interpolators.keys())):
            raise ValueError(
                "kappa_i coefficient is not defined for class "
                f"{self.__class__.__name__}."
            )

        for c in names:
            coef.append(self.interpolators[c](self.Z, self.chi_i))
        return np.array(coef)

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
        names = ["eta0_e", "eta1_e", "eta2_e", "eta3_e", "eta4_e"]

        # Raise an exception of all of the required coefficients are not
        # existant in this class
        if not set(names).issubset(set(self.interpolators.keys())):
            raise ValueError(
                "eta_e coefficient is not defined for class "
                f"{self.__class__.__name__}."
            )

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
        names = ["eta0_i", "eta1_i", "eta2_i", "eta3_i", "eta4_i"]

        # Raise an exception of all of the required coefficients are not
        # existant in this class
        if not set(names).issubset(set(self.interpolators.keys())):
            raise ValueError(
                "eta_i coefficient is not defined for class "
                f"{self.__class__.__name__}."
            )

        for c in names:
            coef.append(self.interpolators[c](self.Z, self.chi_i))
        return coef


if __name__ == "__main__":
    chi = np.linspace(-2, 2, num=50)
    chi = 10 ** chi
    x = AbstractClassicalTransportCoefficients(chi_e=chi, Z=1)
    print(x.beta)
