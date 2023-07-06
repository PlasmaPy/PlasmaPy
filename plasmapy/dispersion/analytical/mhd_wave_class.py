"""
Objects for storing magnetohydrodynamic wave parameters.
"""
__all__ = ["AlfvenWave", "FastMagnetosonicWave", "SlowMagnetosonicWave", "mhd_waves"]

import astropy.units as u
import numpy as np

from abc import ABC, abstractmethod
from astropy.constants.si import k_B
from numbers import Integral, Real
from typing import Optional, Union

from plasmapy.formulary.speeds import Alfven_speed
from plasmapy.particles import electron, particle_input, ParticleLike
from plasmapy.utils.decorators import validate_quantities


class AbstractMHDWave(ABC):
    @particle_input
    @validate_quantities(
        B={"can_be_negative": False},
        density={"can_be_negative": False},
        T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    )
    def __init__(
        self,
        B: u.T,
        density: (u.m**-3, u.kg / u.m**3),
        ion: ParticleLike,
        *,
        T: u.K = 0 * u.K,
        gamma: Union[float, int] = 5 / 3,
        mass_numb: Optional[Integral] = None,
        Z: Optional[Real] = None,
    ):
        # validate arguments
        for arg_name in ("B", "density", "T"):
            val = locals()[arg_name].squeeze()
            if val.shape != ():
                raise ValueError(
                    f"Argument '{arg_name}' must a single value and not an array of "
                    f"shape {val.shape}."
                )
            locals()[arg_name] = val

        # validate gamma
        if not isinstance(gamma, Real):
            raise TypeError(
                f"Expected int or float for argument 'gamma', but got "
                f"{type(gamma)}."
            )

        if density.unit.physical_type == u.physical.mass_density:
            self._rho = density
        else:
            self._rho = (ion.mass + ion.charge_number * electron.mass) * density

        # Alfvén speed
        self._v_A = Alfven_speed(B, self._rho)
        # sound speed
        self._c_s = np.sqrt(gamma * k_B * T / ion.mass)
        # magnetosonic speed
        self._c_ms = np.sqrt(self._v_A**2 + self._c_s**2)

    @property
    def alfven_speed(self):
        """The Alfvén speed of the plasma."""
        return self._v_A

    @property
    def sound_speed(self):
        """The sound speed of the plasma."""
        return self._c_s

    @property
    def magnetosonic_speed(self):
        r"""
        The magnetosonic speed of the plasma.

        Defined as :math:`c_{ms} = \sqrt{v_A^2 + c_s^2}` where
        :math:`v_A` is the Alfvén speed and :math:`c_s` is the sound speed
        """
        return self._c_ms

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
    r"""
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
        r"""
        Calculate the analytical solution to a magnetohydrodynamic,
        low-frequency (:math:`ω/kc ≪ 1`) dispersion relation.

        Parameters
        ----------
        k : `~astropy.units.Quantity`, single valued or 1-D array
            Wavenumber in units convertible to rad/m`.  Either single
            valued or 1-D array of length :math:`N`.
        theta : `~astropy.units.Quantity`, single valued or 1-D array
            The angle of propagation of the wave with respect to the
            magnetic field, :math:`\cos^{-1}(k_z / k)`, in units must be
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
                    self._c_ms**2
                    + np.sqrt(
                        self._c_ms**4
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
                    self._c_ms**2
                    - np.sqrt(
                        self._c_ms**4
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
