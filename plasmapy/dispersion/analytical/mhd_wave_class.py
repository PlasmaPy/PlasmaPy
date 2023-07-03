"""
Objects for storing magnetohydrodynamic wave parameters.
"""
__all__ = ["AlfvenWave", "FastMagnetosonicWave", "SlowMagnetosonicWave", "mhd_waves"]

import astropy.units as u
import numpy as np
import warnings

from abc import ABC, abstractmethod
from numbers import Integral, Real
from typing import Optional, Union

from plasmapy.formulary.speeds import Alfven_speed, ion_sound_speed
from plasmapy.particles import particle_input, ParticleLike
from plasmapy.utils.decorators import validate_quantities
from plasmapy.utils.exceptions import PhysicsWarning


class AbstractMHDWave(ABC):
    @particle_input
    @validate_quantities(
        B={"can_be_negative": False},
        n_i={"can_be_negative": False},
        T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
        T_i={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    )
    def __init__(
        self,
        B: u.T,
        ion: ParticleLike,
        n_i: u.m**-3,
        *,
        T_e: u.K = 0 * u.K,
        T_i: u.K = 0 * u.K,
        gamma_e: Union[float, int] = 1,
        gamma_i: Union[float, int] = 3,
        mass_numb: Optional[Integral] = None,
        Z: Optional[Real] = None,
    ):
        # validate arguments
        for arg_name in ("B", "n_i", "T_e", "T_i"):
            val = locals()[arg_name].squeeze()
            if val.shape != ():
                raise ValueError(
                    f"Argument '{arg_name}' must a single value and not an array of "
                    f"shape {val.shape}."
                )
            locals()[arg_name] = val

        # validate arguments
        for arg_name in ("gamma_e", "gamma_i"):
            if not isinstance(locals()[arg_name], Real):
                raise TypeError(
                    f"Expected int or float for argument '{arg_name}', but got "
                    f"{type(locals()[arg_name])}."
                )
        self._B = B
        self._ion = ion
        self._n_i = n_i
        self._T_e = T_e
        self._T_i = T_i
        self._gamma_e = gamma_e
        self._gamma_i = gamma_i
        self._Z = ion.charge_number

        self._v_A = Alfven_speed(self._B, self._n_i, ion=self._ion, Z=self._Z)

        self._n_e = self._Z * self._n_i
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=PhysicsWarning)
            self._c_s = ion_sound_speed(
                T_e=self._T_e,
                T_i=self._T_i,
                ion=self._ion,
                n_e=self._n_e,
                gamma_e=self._gamma_e,
                gamma_i=self._gamma_i,
                Z=self._Z,
            )

        # magnetosonic speed squared
        self._c_ms2 = self._v_A**2 + self._c_s**2

    @staticmethod
    def validate_k_theta(k: u.rad / u.m, theta: u.rad):
        # validate argument k
        k = k.squeeze()
        if k.ndim not in [0, 1]:
            raise ValueError(
                f"Argument 'k' needs to be a single valued or 1D array astropy Quantity,"
                f" got array of shape {k.shape}."
            )
        if np.any(k <= 0):
            raise ValueError("Argument 'k' can not be a or have negative values.")

        # validate argument theta
        theta = theta.squeeze()
        if theta.ndim not in [0, 1]:
            raise ValueError(
                f"Argument 'theta' needs to be a single valued or 1D array astropy "
                f"Quantity, got array of shape {k.shape}."
            )

        # return theta and k as coordinate arrays
        return np.meshgrid(theta, k)

    @abstractmethod
    def angular_frequency(self, k: u.rad / u.m, theta: u.rad):
        pass

    def phase_velocity(self, k: u.rad / u.m, theta: u.rad):
        return self.angular_frequency(k, theta) / k


class AlfvenWave(AbstractMHDWave):
    """
    Core class for MHD waves.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to T.
    ion : `str` or `~plasmapy.particles.particle_class.Particle`
        Representation of the ion species (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, ``'He-4 +1'`` for singly ionized
        helium-4, etc.). If no charge state information is provided,
        then the ions are assumed to be singly ionized.
    n_i : `~astropy.units.Quantity`
        Ion number density in units convertible to m\ :sup:`-3`.
    T_e : `~astropy.units.Quantity`, optional
        The electron temperature in units of K or eV, which defaults
        to zero.
    T_i : `~astropy.units.Quantity`, optional
        The ion temperature in units of K or eV, which defaults to
        zero.
    gamma_e : `float` or `int`, optional
        The adiabatic index for electrons, which defaults to 1.  This
        value assumes that the electrons are able to equalize their
        temperature rapidly enough that the electrons are effectively
        isothermal.
    gamma_i : `float` or `int`, optional
        The adiabatic index for ions, which defaults to 3. This value
        assumes that ion motion has only one degree of freedom, namely
        along magnetic field lines.
    Z : `float` or int, optional
        The average ionization state (arithmetic mean) of the ``ion``
        composing the plasma.  Will override any charge state defined
        by argument ``ion``.

    Raises
    ------
    TypeError
        If applicable arguments are not instances of
        `~astropy.units.Quantity` or cannot be converted into one.

    TypeError
        If ``ion`` is not of type or convertible to
        `~plasmapy.particles.particle_class.Particle`.

    TypeError
        If ``gamma_e``, ``gamma_i``, or ``Z`` are not of type `int`
        or `float`.

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the
        expected units.

    ValueError
        If any of ``B``, ``k``, ``n_i``, ``T_e``, or ``T_i`` is negative.

    ValueError
        If ``k`` is negative or zero.

    ValueError
        If ``ion`` is not of category ion or element.

    ValueError
        If ``B``, ``n_i``, ``T_e``, or ``T_I`` are not single valued
        `astropy.units.Quantity` (i.e. an array).
    """

    def angular_frequency(self, k: u.rad / u.m, theta: u.rad):
        """
        Calculate the analytical solution to a magnetohydrodynamic,
        low-frequency (:math:`ω/kc ≪ 1`) dispersion relation.

        Parameters
        ----------
        k : `~astropy.units.Quantity`, single valued or 1-D array
            Wavenumber in units convertible to rad/m`.  Either single
            valued or 1-D array of length :math:`N`.
        theta : `~astropy.units.Quantity`, single valued or 1-D array
            The angle of propagation of the wave with respect to the
            magnetic field, :math:`\\cos^{-1}(k_z / k)`, in units must be
            convertible to radians. Either single valued or 1-D array of
            size :math:`M`.

        Returns
        -------
        omega : `~astropy.units.Quantity`
            An :math:`N × M` array of computed wave frequencies in units
            rad/s.

        Raises
        ------
        ~astropy.units.UnitTypeError
            If applicable arguments do not have units convertible to the
            expected units.

        ValueError
            If ``k`` is negative or zero.

        ValueError
            If ``k`` or ``theta`` are not single valued or a 1-D array.
        """
        theta, k = super().validate_k_theta(k, theta)
        return np.squeeze(k * self._v_A * np.cos(theta))


class FastMagnetosonicWave(AbstractMHDWave):
    def angular_frequency(self, k: u.rad / u.m, theta: u.rad):
        theta, k = super().validate_k_theta(k, theta)
        return np.squeeze(
            k
            * np.sqrt(
                (
                    self._c_ms2
                    + np.sqrt(
                        self._c_ms2**2
                        - (2 * self._v_A * self._c_s * np.cos(theta)) ** 2
                    )
                )
                / 2
            )
        )


class SlowMagnetosonicWave(AbstractMHDWave):
    def angular_frequency(self, k: u.rad / u.m, theta: u.rad):
        theta, k = super().validate_k_theta(k, theta)
        return np.squeeze(
            k
            * np.sqrt(
                (
                    self._c_ms2
                    - np.sqrt(
                        self._c_ms2**2
                        - (2 * self._v_A * self._c_s * np.cos(theta)) ** 2
                    )
                )
                / 2
            )
        )


def mhd_waves(*args, **kwargs):
    return {
        "alfven": AlfvenWave(*args, **kwargs),
        "fast": FastMagnetosonicWave(*args, **kwargs),
        "slow": SlowMagnetosonicWave(*args, **kwargs),
    }
