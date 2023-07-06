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
    """
    Core class for magnetohydrodynamic waves.
    """

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
        self._v_a = Alfven_speed(B, self._rho)
        # sound speed
        self._c_s = np.sqrt(gamma * k_B * T / ion.mass)
        # magnetosonic speed
        self._c_ms = np.sqrt(self._v_a**2 + self._c_s**2)

    @property
    def alfven_speed(self):
        """The Alfvén speed of the plasma."""
        return self._v_a

    @property
    def sound_speed(self):
        """The sound speed of the plasma."""
        return self._c_s

    @property
    def magnetosonic_speed(self):
        r"""
        The magnetosonic speed of the plasma.

        Defined as :math:`c_{ms} = \sqrt{v_a^2 + c_s^2}` where
        :math:`v_a` is the Alfvén speed and :math:`c_s` is the sound speed
        """
        return self._c_ms

    @staticmethod
    def _validate_k_theta(k: u.rad / u.m, theta: u.rad):
        """Validate and return arguments."""
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
    A class to represent magnetohydrodynamic Alfvén waves.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to T.
    density : `~astropy.units.Quantity`
        Either the ion number density :math:`n_i` in units convertible
        to m\ :sup:`-3` or the total mass density :math:`ρ` in units
        convertible to kg m\ :sup:`-3`\ .
    ion : `str` or `~plasmapy.particles.particle_class.Particle`
        Representation of the ion species (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, ``'He-4 +1'`` for singly ionized
        helium-4, etc.). If no charge state information is provided,
        then the ions are assumed to be singly ionized.
    T : `~astropy.units.Quantity`, |keyword-only|, optional
        The plasma temperature in units of K or eV, which defaults
        to zero.
    gamma : `float` or `int`, |keyword-only|, optional
        The adiabatic index for the plasma, which defaults to 3/5.
    mass_numb : integer, |keyword-only|, optional
        The mass number corresponding to ``ion``.
    Z : `float` or int, optional, |keyword-only|
        The charge number corresponding to ``ion``.

    Raises
    ------
    TypeError
        If applicable arguments are not instances of
        `~astropy.units.Quantity` or cannot be converted into one.

    TypeError
        If ``ion`` is not |particle-like|.

    TypeError
        If ``gamma`` or ``Z`` are not of type `int` or `float`.

    TypeError
        If ``mass_numb`` is not of type `int`.

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the
        expected units.

    ValueError
        If any of ``B``, ``density``, or ``T`` is negative.

    ValueError
        If ``ion`` is not of category ion or element.

    ValueError
        If ``B``, ``density``, or ``T`` are not single valued
        `astropy.units.Quantity` (i.e. an array).
    """

    def angular_frequency(self, k: u.rad / u.m, theta: u.rad):
        r"""
        Calculate the angular frequency of magnetohydrodynamic
        Alfvén waves.

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
        theta, k = super()._validate_k_theta(k, theta)
        return np.squeeze(k * self._v_a * np.cos(theta))


class FastMagnetosonicWave(AbstractMHDWave):
    r"""
    A class to represent fast magnetosonic waves.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to T.
    density : `~astropy.units.Quantity`
        Either the ion number density :math:`n_i` in units convertible
        to m\ :sup:`-3` or the total mass density :math:`ρ` in units
        convertible to kg m\ :sup:`-3`\ .
    ion : `str` or `~plasmapy.particles.particle_class.Particle`
        Representation of the ion species (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, ``'He-4 +1'`` for singly ionized
        helium-4, etc.). If no charge state information is provided,
        then the ions are assumed to be singly ionized.
    T : `~astropy.units.Quantity`, |keyword-only|, optional
        The plasma temperature in units of K or eV, which defaults
        to zero.
    gamma : `float` or `int`, |keyword-only|, optional
        The adiabatic index for the plasma, which defaults to 3/5.
    mass_numb : integer, |keyword-only|, optional
        The mass number corresponding to ``ion``.
    Z : `float` or int, optional, |keyword-only|
        The charge number corresponding to ``ion``.

    Raises
    ------
    TypeError
        If applicable arguments are not instances of
        `~astropy.units.Quantity` or cannot be converted into one.

    TypeError
        If ``ion`` is not |particle-like|.

    TypeError
        If ``gamma`` or ``Z`` are not of type `int` or `float`.

    TypeError
        If ``mass_numb`` is not of type `int`.

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the
        expected units.

    ValueError
        If any of ``B``, ``density``, or ``T`` is negative.

    ValueError
        If ``ion`` is not of category ion or element.

    ValueError
        If ``B``, ``density``, or ``T`` are not single valued
        `astropy.units.Quantity` (i.e. an array).
    """

    def angular_frequency(self, k: u.rad / u.m, theta: u.rad):
        r"""
        Calculate the angular frequency of a fast magnetosonic waves.

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
        theta, k = super()._validate_k_theta(k, theta)
        return np.squeeze(
            k
            * np.sqrt(
                (
                    self._c_ms**2
                    + np.sqrt(
                        self._c_ms**4
                        - (2 * self._v_a * self._c_s * np.cos(theta)) ** 2
                    )
                )
                / 2
            )
        )


class SlowMagnetosonicWave(AbstractMHDWave):
    r"""
    A class to represent slow magnetosonic waves.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to T.
    density : `~astropy.units.Quantity`
        Either the ion number density :math:`n_i` in units convertible
        to m\ :sup:`-3` or the total mass density :math:`ρ` in units
        convertible to kg m\ :sup:`-3`\ .
    ion : `str` or `~plasmapy.particles.particle_class.Particle`
        Representation of the ion species (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, ``'He-4 +1'`` for singly ionized
        helium-4, etc.). If no charge state information is provided,
        then the ions are assumed to be singly ionized.
    T : `~astropy.units.Quantity`, |keyword-only|, optional
        The plasma temperature in units of K or eV, which defaults
        to zero.
    gamma : `float` or `int`, |keyword-only|, optional
        The adiabatic index for the plasma, which defaults to 3/5.
    mass_numb : integer, |keyword-only|, optional
        The mass number corresponding to ``ion``.
    Z : `float` or int, optional, |keyword-only|
        The charge number corresponding to ``ion``.

    Raises
    ------
    TypeError
        If applicable arguments are not instances of
        `~astropy.units.Quantity` or cannot be converted into one.

    TypeError
        If ``ion`` is not |particle-like|.

    TypeError
        If ``gamma`` or ``Z`` are not of type `int` or `float`.

    TypeError
        If ``mass_numb`` is not of type `int`.

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the
        expected units.

    ValueError
        If any of ``B``, ``density``, or ``T`` is negative.

    ValueError
        If ``ion`` is not of category ion or element.

    ValueError
        If ``B``, ``density``, or ``T`` are not single valued
        `astropy.units.Quantity` (i.e. an array).
    """

    def angular_frequency(self, k: u.rad / u.m, theta: u.rad):
        r"""
        Calculate the angular frequency of slow magnetosonic waves.

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
        theta, k = super()._validate_k_theta(k, theta)
        return np.squeeze(
            k
            * np.sqrt(
                (
                    self._c_ms**2
                    - np.sqrt(
                        self._c_ms**4
                        - (2 * self._v_a * self._c_s * np.cos(theta)) ** 2
                    )
                )
                / 2
            )
        )


def mhd_waves(*args, **kwargs):
    r"""
    Returns a dictionary containing objects of the three
    magnetohydrodynamic waves with identical parameters.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to T.
    density : `~astropy.units.Quantity`
        Either the ion number density :math:`n_i` in units convertible
        to m\ :sup:`-3` or the total mass density :math:`ρ` in units
        convertible to kg m\ :sup:`-3`\ .
    ion : `str` or `~plasmapy.particles.particle_class.Particle`
        Representation of the ion species (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, ``'He-4 +1'`` for singly ionized
        helium-4, etc.). If no charge state information is provided,
        then the ions are assumed to be singly ionized.
    T : `~astropy.units.Quantity`, |keyword-only|, optional
        The plasma temperature in units of K or eV, which defaults
        to zero.
    gamma : `float` or `int`, |keyword-only|, optional
        The adiabatic index for the plasma, which defaults to 3/5.
    mass_numb : integer, |keyword-only|, optional
        The mass number corresponding to ``ion``.
    Z : `float` or int, optional, |keyword-only|
        The charge number corresponding to ``ion``.

    Returns
    -------
    mhd_waves : Dict[str, `~plasmapy.dispersion.analytical.AlfvenWave` or `~plasmapy.dispersion.analytical.FastMagnetosonicWave` or `~plasmapy.dispersion.analytical.SlowMagnetosonicWave`]
        A dictionary of magnetohydrodynamic-wave objects. The
        dictionary contains three keys: ``'alfven'`` for the Alfvén
        mode, ``'fast'`` for the fast magnetosonic mode, and
        ``'slow'`` for the slow magnetosonic mode.  The value for
        each key will be of type

    Raises
    ------
    TypeError
        If applicable arguments are not instances of
        `~astropy.units.Quantity` or cannot be converted into one.

    TypeError
        If ``ion`` is not |particle-like|.

    TypeError
        If ``gamma`` or ``Z`` are not of type `int` or `float`.

    TypeError
        If ``mass_numb`` is not of type `int`.

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the
        expected units.

    ValueError
        If any of ``B``, ``density``, or ``T`` is negative.

    ValueError
        If ``ion`` is not of category ion or element.

    ValueError
        If ``B``, ``density``, or ``T`` are not single valued
        `astropy.units.Quantity` (i.e. an array).
    """
    return {
        "alfven": AlfvenWave(*args, **kwargs),
        "fast": FastMagnetosonicWave(*args, **kwargs),
        "slow": SlowMagnetosonicWave(*args, **kwargs),
    }
