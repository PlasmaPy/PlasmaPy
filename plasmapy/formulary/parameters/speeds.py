"""
Functions to calculated fundamental plasma speed/veolcity parameters.
"""
__all__ = [
    "Alfven_speed",
    "ion_sound_speed",
]
__aliases__ = [
    "cs_",
    "va_",
]

import astropy.units as u
import numbers
import numpy as np
import warnings

from astropy.constants.si import k_B, mu0
from typing import Optional

from plasmapy import particles
from plasmapy.formulary.parameters.lengths import Debye_length
from plasmapy.formulary.parameters.parameters_ import mass_density, _grab_charge
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import ChargeError
from plasmapy.utils.decorators import check_relativistic, validate_quantities
from plasmapy.utils.exceptions import PhysicsError, PhysicsWarning

__all__ += __aliases__


@check_relativistic
@validate_quantities(density={"can_be_negative": False})
def Alfven_speed(
    B: u.T,
    density: (u.m ** -3, u.kg / u.m ** 3),
    ion: Optional[Particle] = None,
    z_mean: Optional[numbers.Real] = None,
) -> u.m / u.s:
    r"""
    Calculate the Alfvén speed.

    The Alfvén speed :math:`V_A` is the typical propagation speed of magnetic
    disturbances in a plasma, and is given by:

    .. math::

        V_A = \frac{B}{\sqrt{μ_0 ρ}}

    where :math:`B` is the magnetic field and :math:`ρ = n_i m_i + n_e m_e`
    is the total mass density (:math:`n_i` is the ion number density,
    :math:`n_e` is the electron number density, :math:`m_i` is the ion mass,
    and :math:`m_e` is the electron mass) :cite:p:`alfven:1942`.

    **Aliases:** `va_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to tesla.

    density : `~astropy.units.Quantity`
        Either the ion number density :math:`n_i` in units convertible to
        m\ :sup:`-3` or the total mass density :math:`ρ` in units
        convertible to kg m\ :sup:`-3`\ .

    ion : `~plasmapy.particles.Particle`, optional
        Representation of the ion species (e.g., `'p'` for protons, `'D+'` for
        deuterium, `'He-4 +1'` for singly ionized helium-4, etc.). If no charge
        state information is provided, then the ions are assumed to be singly
        ionized. If the density is an ion number density, then this paramter
        is required in order to convert to mass density.

    z_mean : `~numbers.Real`, optional
        The average ionization state (arithmetic mean) of the ``ion`` composing
        the plasma.  This is used in calculating the mass density
        :math:`ρ = n_i (m_i + Z_{mean} m_e)`.  ``z_mean`` is ignored if
        ``density`` is passed as a mass density and overrides any charge state
        info provided by ``ion``.

    Returns
    -------
    V_A : `~astropy.units.Quantity`
        The Alfvén speed in units of m s\ :sup:`-1`.

    Raises
    ------
    `~plasmapy.utils.exceptions.RelativityError`
        If the Alfvén velocity is greater than or equal to the speed of light.

    `TypeError`
        If ``B`` and/or ``density`` are not of type `~astropy.units.Quantity`,
        or convertible.

    `TypeError`
        If ``ion`` is not of type or convertible to `~plasmapy.particles.Particle`.

    `TypeError`
        If ``z_mean`` is not of type `int` or `float`.

    `~astropy.units.UnitTypeError`
        If the magnetic field ``B`` does not have units equivalent to
        tesla.

    `~astropy.units.UnitTypeError`
        If the ``density`` does not have units equivalent to a number density
        or mass density.

    `ValueError`
        If ``density`` is negative.

    Warns
    -----
    : `~plasmapy.utils.exceptions.RelativityWarning`
        If the Alfvén velocity exceeds 5% of the speed of light.

    : `~astropy.units.UnitsWarning`
        If units are not provided for the magnetic field ``B``, units of
        tesla are assumed.

    Notes
    -----
    This expression does not account for relativistic effects, and
    loses validity when the resulting speed is a significant fraction
    of the speed of light.

    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy.constants.si import m_p, m_e
    >>> B = 0.014*u.T
    >>> n = 5e19*u.m**-3
    >>> rho = n*(m_p+m_e)
    >>> ion = 'p'
    >>> Alfven_speed(B, n, ion=ion)
    <Quantity 43173.870... m / s>
    >>> Alfven_speed(B, rho)
    <Quantity 43173.870... m / s>
    >>> Alfven_speed(B, rho).to(u.cm/u.us)
    <Quantity 4.317387 cm / us>
    >>> Alfven_speed(B, n, ion="He +2")
    <Quantity 21664.18... m / s>
    >>> Alfven_speed(B, n, ion="He++")
    <Quantity 21664.18... m / s>
    >>> Alfven_speed(B, n, ion="He", z_mean=1.8)
    <Quantity 21661.51... m / s>
    """
    if density.unit.is_equivalent(u.kg / u.m ** 3):
        rho = density
    else:
        if not isinstance(ion, Particle):
            try:
                ion = Particle(ion)
            except TypeError:
                raise TypeError(
                    f"If passing a number density, you must pass a plasmapy Particle "
                    f"(not type {type(ion)}) to calculate the mass density!"
                )
        if z_mean is None:
            try:
                z_mean = abs(ion.charge_number)
            except ChargeError:
                z_mean = 1

        z_mean = abs(z_mean)
        rho = mass_density(density, ion) + mass_density(density, "e", z_ratio=z_mean)

    V_A = np.abs(B) / np.sqrt(mu0 * rho)
    return V_A


va_ = Alfven_speed
"""Alias to `~plasmapy.formulary.parameters.parameters_.Alfven_speed`."""


@check_relativistic
@validate_quantities(
    T_i={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False, "none_shall_pass": True},
    k={"can_be_negative": False, "none_shall_pass": True},
)
def ion_sound_speed(
    T_e: u.K,
    T_i: u.K,
    ion: Particle,
    n_e: u.m ** -3 = None,
    k: u.m ** -1 = None,
    gamma_e=1,
    gamma_i=3,
    z_mean=None,
) -> u.m / u.s:
    r"""
    Return the ion sound speed for an electron-ion plasma.

    **Aliases:** `cs_`

    Parameters
    ----------
    T_e : `~astropy.units.Quantity`
        Electron temperature in units of temperature or energy per
        particle. If this is not given, then the electron temperature
        is assumed to be zero.

    T_i : `~astropy.units.Quantity`
        Ion temperature in units of temperature or energy per
        particle.  If this is not given, then the ion temperature is
        assumed to be zero.

    ion : `~plasmapy.particles.Particle`
        Representation of the ion species (e.g., `'p'` for protons,
        `'D+'` for deuterium, or 'He-4 +1' for singly ionized
        helium-4). If no charge state information is provided, then the
        ions are assumed to be singly charged.

    n_e : `~astropy.units.Quantity`
        Electron number density. If this is not given, then ion_sound_speed
        will be approximated in the non-dispersive limit
        (:math:`k^2 λ_{D}^2` will be assumed zero). If ``n_e`` is given,
        a value for ``k`` must also be given.

    k : `~astropy.units.Quantity`
        Wavenumber (in units of inverse length, e.g. m\ :sup:`-1`\ ). If this
        is not given, then ion_sound_speed will be approximated in the
        non-dispersive limit (:math:`k^2 λ_{D}^2` will be assumed zero).
        If ``k`` is given, a value for ``n_e`` must also be given.

    gamma_e : `float` or `int`
        The adiabatic index for electrons, which defaults to 1.  This
        value assumes that the electrons are able to equalize their
        temperature rapidly enough that the electrons are effectively
        isothermal.

    gamma_i : `float` or `int`
        The adiabatic index for ions, which defaults to 3.  This value
        assumes that ion motion has only one degree of freedom, namely
        along magnetic field lines.

    z_mean : `~astropy.units.Quantity`, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. If this quantity is not
        given then the charge number of the ion
        is used. This is effectively an average ion sound speed for the
        plasma where multiple charge states are present.

    Returns
    -------
    V_S : `~astropy.units.Quantity`
        The ion sound speed in units of meters per second.

    Raises
    ------
    `TypeError`
        If any of the arguments are not entered as keyword arguments
        or are of an incorrect type.

    `ValueError`
        If the ion mass, adiabatic index, or temperature are invalid.

    `~plasmapy.utils.exceptions.PhysicsError`
        If an adiabatic index is less than one.

    `~astropy.units.UnitConversionError`
        If the temperature, electron number density, or wavenumber
        is in incorrect units.

    Warns
    -----
    : `~plasmapy.utils.exceptions.RelativityWarning`
        If the ion sound speed exceeds 5% of the speed of light.

    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    : `~plasmapy.utils.exceptions.PhysicsWarning`
        If only one of ``k`` or ``n_e`` is given, the non-dispersive
        limit is assumed.

    Notes
    -----
    The ion sound speed :math:`V_S` is given by

    .. math::

        V_S = \sqrt{\frac{γ_e Z k_B T_e + γ_i k_B T_i}{m_i (1 + k^2 λ_{D}^2)}}

    where :math:`γ_e` and :math:`γ_i` are the electron and
    ion adiabatic indices, :math:`k_B` is the Boltzmann constant,
    :math:`T_e` and :math:`T_i` are the electron and ion temperatures,
    :math:`Z` is the charge state of the ion, :math:`m_i` is the
    ion mass, :math:`λ_D` is the Debye length, and :math:`k` is the
    wavenumber.

    In the non-dispersive limit (:math:`k^2 λ_D^2` is small) the
    equation for :math:`V_S` is approximated (the denominator reduces
    to :math:`m_i`).

    When the electron temperature is much greater than the ion
    temperature, the ion sound velocity reduces to
    :math:`\sqrt{γ_e k_B T_e / m_i}`. Ion acoustic waves can
    therefore occur even when the ion temperature is zero.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 5e19*u.m**-3
    >>> k_1 = 3e1*u.m**-1
    >>> k_2 = 3e7*u.m**-1
    >>> ion_sound_speed(T_e=5e6*u.K, T_i=0*u.K, ion='p', gamma_e=1, gamma_i=3)
    <Quantity 203155... m / s>
    >>> ion_sound_speed(T_e=5e6*u.K, T_i=0*u.K, n_e=n, k=k_1, ion='p', gamma_e=1, gamma_i=3)
    <Quantity 203155... m / s>
    >>> ion_sound_speed(T_e=5e6*u.K, T_i=0*u.K, n_e=n, k=k_2, ion='p', gamma_e=1, gamma_i=3)
    <Quantity 310.31... m / s>
    >>> ion_sound_speed(T_e=5e6*u.K, T_i=0*u.K, n_e=n, k=k_1, ion='p')
    <Quantity 203155... m / s>
    >>> ion_sound_speed(T_e=500*u.eV, T_i=200*u.eV, n_e=n, k=k_1, ion='D+')
    <Quantity 229585... m / s>

    """

    m_i = particles.particle_mass(ion)
    Z = _grab_charge(ion, z_mean)

    for gamma, species in zip([gamma_e, gamma_i], ["electrons", "ions"]):
        if not isinstance(gamma, (numbers.Real, numbers.Integral)):
            raise TypeError(
                f"The adiabatic index gamma for {species} must be a float or int"
            )
        if gamma < 1:
            raise PhysicsError(
                f"The adiabatic index for {species} must be between "
                f"one and infinity"
            )

    # Assume non-dispersive limit if values for n_e (or k) are not specified
    klD2 = 0.0
    if (n_e is None) ^ (k is None):
        warnings.warn(
            "The non-dispersive limit has been assumed for "
            "this calculation. To prevent this, values must "
            "be specified for both n_e and k.",
            PhysicsWarning,
        )
    elif n_e is not None and k is not None:
        lambda_D = Debye_length(T_e, n_e)
        klD2 = (k * lambda_D) ** 2

    try:
        V_S_squared = (gamma_e * Z * k_B * T_e + gamma_i * k_B * T_i) / (
            m_i * (1 + klD2)
        )
        V_S = np.sqrt(V_S_squared).to(u.m / u.s)
    except Exception:
        raise ValueError("Unable to find ion sound speed.")

    return V_S


cs_ = ion_sound_speed
"""Alias to `~plasmapy.formulary.parameters.parameters_.ion_sound_speed`."""
