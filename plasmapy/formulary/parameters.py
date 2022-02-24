"""Functions to calculate fundamental plasma parameters."""

__all__ = [
    "Alfven_speed",
    "Bohm_diffusion",
    "Debye_number",
    "Hall_parameter",
    "ion_sound_speed",
    "kappa_thermal_speed",
    "magnetic_energy_density",
    "magnetic_pressure",
    "mass_density",
    "plasma_frequency",
    "thermal_pressure",
    "thermal_speed",
    "thermal_speed_coefficients",
    "upper_hybrid_frequency",
]
__aliases__ = [
    "betaH_",
    "cs_",
    "DB_",
    "nD_",
    "pmag_",
    "pth_",
    "rho_",
    "ub_",
    "va_",
    "vth_",
    "vth_kappa_",
    "wp_",
    "wuh_",
]
__lite_funcs__ = ["thermal_speed_lite", "plasma_frequency_lite"]

import astropy.units as u
import numbers
import numpy as np
import warnings

from astropy.constants.si import e, eps0, k_B, mu0
from numba import njit
from typing import Optional, Union

from plasmapy import particles
from plasmapy.formulary import frequencies
from plasmapy.formulary import lengths
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import ChargeError
from plasmapy.utils import PhysicsError
from plasmapy.utils.decorators import (
    angular_freq_to_hz,
    bind_lite_func,
    check_relativistic,
    deprecated,
    preserve_signature,
    validate_quantities,
)
from plasmapy.utils.exceptions import (
    PhysicsWarning,
    PlasmaPyFutureWarning,
    RelativityWarning,
)

__aliases__ += frequencies.__aliases__ + lengths.__aliases__
__all__ += frequencies.__all__ + lengths.__all__ + __aliases__ + __lite_funcs__

e_si_unitless = e.value
eps0_si_unitless = eps0.value
k_B_si_unitless = k_B.value

funcs_to_deprecate_wrap = [  # (module_name, func_name)
    ("lengths", "Debye_length"),
    ("lengths", "lambdaD_"),
    ("lengths", "gyroradius"),
    ("lengths", "rc_"),
    ("lengths", "rhoc_"),
    ("lengths", "inertial_length"),
    ("lengths", "cwp_"),
    ("frequencies", "gyrofrequency"),
    ("frequencies", "oc_"),
    ("frequencies", "wc_"),
    ("frequencies", "lower_hybrid_frequency"),
    ("frequencies", "wlh_"),
]
for modname, name in funcs_to_deprecate_wrap:
    globals()[name] = deprecated(
        since="0.7.0",
        warning_type=PlasmaPyFutureWarning,
        message=(
            f"The {name}() function has been moved to "
            f"plasmapy.formulary.{modname}.  Update your import to get "
            f"rid of this warning."
        ),
    )(getattr(globals()[f"{modname}"], name))

del modname, name


def _grab_charge(ion: Particle, z_mean=None):
    """
    Merge two possible inputs for particle charge.

    Parameters
    ----------
    ion : `~plasmapy.particles.particle_class.Particle`
        a string representing a charged particle, or a Particle object.

    z_mean : `float`
        An optional float describing the average ionization of a particle
        species.

    Returns
    -------
    float
        if ``z_mean`` was passed, ``z_mean``, otherwise, the charge number
        of ``ion``.

    """
    if z_mean is None:
        # warnings.warn("No z_mean given, defaulting to atomic charge",
        #               PhysicsWarning)
        Z = particles.charge_number(ion)
    else:
        # using average ionization provided by user
        Z = z_mean
    return Z


@validate_quantities(
    density={"can_be_negative": False}, validations_on_return={"can_be_negative": False}
)
def mass_density(
    density: (u.m ** -3, u.kg / (u.m ** 3)),
    particle: Union[Particle, str],
    z_ratio: Optional[numbers.Real] = 1,
) -> u.kg / u.m ** 3:
    r"""
    Calculate the mass density from a number density.

    .. math::

        \rho = \left| \frac{Z_{s}}{Z_{particle}} \right| n_{s} m_{particle}
              = | Z_{ratio} | n_{s} m_{particle}

    where :math:`m_{particle}` is the particle mass, :math:`n_{s}` is a number
    density for plasma species :math:`s`, :math:`Z_{s}` is the charge number of
    species :math:`s`, and :math:`Z_{particle}` is the charge number of
    ``particle``.  For example, if the electron density is given for :math:`n_s`
    and ``particle`` is a doubly ionized atom, then :math:`Z_{ratio} = -1 / 2`\ .

    **Aliases:** `rho_`

    Parameters
    ----------
    density : `~astropy.units.Quantity`
        Either a particle number density (in units of m\ :sup:`-3` or
        equivalent) or a mass density (in units of kg / m\ :sup:`3` or
        equivalent).  If ``density`` is a mass density, then it will be passed
        through and returned without modification.

    particle : `~plasmapy.particles.particle_class.Particle`
        The particle for which the mass density is being calculated for.  Must
        be a `~plasmapy.particles.particle_class.Particle` or a value convertible to
        a `~plasmapy.particles.particle_class.Particle` (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, or ``'He-4 +1'`` for singly ionized helium-4).

    z_ratio : `int`, `float`, optional
        The ratio of the charge numbers corresponding to the plasma species
        represented by ``density`` and the ``particle``.  For example, if the
        given ``density`` is and electron density and ``particle`` is doubly
        ionized ``He``, then ``z_ratio = -0.5``.  Default is ``1``.

    Raises
    ------
    `~astropy.units.UnitTypeError`
        If the ``density`` does not have units equivalent to a number density
        or mass density.

    `TypeError`
        If ``density`` is not of type `~astropy.units.Quantity`, or convertible.

    `TypeError`
        If ``particle`` is not of type or convertible to
        `~plasmapy.particles.particle_class.Particle`.

    `TypeError`
        If ``z_ratio`` is not of type `int` or `float`.

    `ValueError`
        If ``density`` is negative.

    Returns
    -------
    `~astropy.units.Quantity`
        The mass density for the plasma species represented by ``particle``.

    Examples
    --------
    >>> import astropy.units as u
    >>> mass_density(1 * u.m ** -3, 'p')
    <Quantity 1.67262...e-27 kg / m3>
    >>> mass_density(4 * u.m ** -3, 'D+')
    <Quantity 1.33743...e-26 kg / m3>
    >>> mass_density(2.e12 * u.cm ** -3, 'He')
    <Quantity 1.32929...e-08 kg / m3>
    >>> mass_density(2.e12 * u.cm ** -3, 'He', z_ratio=0.5)
    <Quantity 6.64647...e-09 kg / m3>
    >>> mass_density(1.0 * u.g * u.m ** -3, "")
    <Quantity 0.001 kg / m3>
    """
    if density.unit.is_equivalent(u.kg / u.m ** 3):
        return density

    if not isinstance(particle, Particle):
        try:
            particle = Particle(particle)
        except TypeError:
            raise TypeError(
                f"If passing a number density, you must pass a plasmapy Particle "
                f"(not type {type(particle)}) to calculate the mass density!"
            )

    if not isinstance(z_ratio, (float, np.floating, int, np.integer)):
        raise TypeError(
            f"Expected type int or float for keyword z_ratio, got type {type(z_ratio)}."
        )

    return abs(z_ratio) * density * particle.mass


rho_ = mass_density
"""Alias to `~plasmapy.formulary.parameters.mass_density`."""


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

    ion : `~plasmapy.particles.particle_class.Particle`, optional
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
        If ``ion`` is not of type or convertible to `~plasmapy.particles.particle_class.Particle`.

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
"""Alias to `~plasmapy.formulary.parameters.Alfven_speed`."""


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

    ion : `~plasmapy.particles.particle_class.Particle`
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
        lambda_D = lengths.Debye_length(T_e, n_e)
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
"""Alias to `~plasmapy.formulary.parameters.ion_sound_speed`."""


def thermal_speed_coefficients(method: str, ndim: int) -> float:
    r"""
    Get the thermal speed coefficient corresponding to the desired
    thermal speed definition.

    See the `~plasmapy.formulary.parameters.thermal_speed`
    :ref:`Notes <thermal-speed-notes>` section for further details of
    the various thermal speed definitions.

    Parameters
    ----------
    method : `str`
        Method to be used for calculating the thermal speed. Valid
        values are ``"most_probable"``, ``"rms"``, ``"mean_magnitude"``,
        and ``"nrl"``.

    ndim : `int`
        Dimensionality (1D, 2D, 3D) of space in which to calculate
        thermal speed. Valid values are ``1``, ``2``, or ``3``.

    Raises
    ------
    `ValueError`
        If ``method`` or ``ndim`` are not a valid value.

    Notes
    -----
    For a detailed explanation of the different coefficients used to
    calculate the thermal speed, then look to the
    :ref:`Notes <thermal-speed-notes>` section for
    `~plasmapy.formulary.parameters.thermal_speed`.  The possible return
    values are listed the following table:

    .. table:: Thermal speed :math:`v_{th}` coefficients.
       :widths: 2 1 1 1 1
       :width: 100%

       +--------------+------------+---------------+---------------+---------------+
       | ↓ **method** | **ndim** → | ``1``         | ``2``         | ``3``         |
       +--------------+------------+---------------+---------------+---------------+
       | ``"most_probable"``       | .. math::     | .. math::     | .. math::     |
       |                           |    0          |    1          |    \sqrt{2}   |
       +--------------+------------+---------------+---------------+---------------+
       | ``"rms"``                 | .. math::     | .. math::     | .. math::     |
       |                           |    1          |    \sqrt{2}   |    \sqrt{3}   |
       +--------------+------------+---------------+---------------+---------------+
       | ``"mean_magnitude"``      | .. math::     | .. math::     | .. math::     |
       |                           |    \sqrt{2/π} |    \sqrt{π/2} |    \sqrt{8/π} |
       +--------------+------------+---------------+---------------+---------------+
       | ``"nrl"``                 | .. math::                                     |
       |                           |    1                                          |
       +--------------+------------+---------------+---------------+---------------+

    Examples
    --------
    >>> thermal_speed_coefficients(method="most_probable", ndim=3)
    1.414213...
    """
    _coefficients = {
        (1, "most_probable"): 0,
        (2, "most_probable"): 1,
        (3, "most_probable"): np.sqrt(2),
        (1, "rms"): 1,
        (2, "rms"): np.sqrt(2),
        (3, "rms"): np.sqrt(3),
        (1, "mean_magnitude"): np.sqrt(2 / np.pi),
        (2, "mean_magnitude"): np.sqrt(np.pi / 2),
        (3, "mean_magnitude"): np.sqrt(8 / np.pi),
        (1, "nrl"): 1,
        (2, "nrl"): 1,
        (3, "nrl"): 1,
    }

    try:
        coeff = _coefficients[(ndim, method)]
    except KeyError:
        raise ValueError(
            f"Value for (ndim, method) pair not valid, got '({ndim}, {method})'."
        )

    return coeff


@preserve_signature
@njit
def thermal_speed_lite(
    T: numbers.Real, mass: numbers.Real, coeff: numbers.Real
) -> numbers.Real:
    r"""
    The ":term:`lite-function`" version of
    `~plasmapy.formulary.parameters.thermal_speed`.  Performs the same
    thermal speed calculations as
    `~plasmapy.formulary.parameters.thermal_speed`, but is intended for
    computational use and, thus, has data conditioning safeguards
    removed.

    .. math::
        v_{th} = C_o \sqrt{\frac{k_B T}{m}}

    where :math:`T` is the temperature associated with the distribution,
    :math:`m` is the particle's mass, and :math:`C_o` is a constant of
    proportionality determined by the method in which :math:`v_{th}` is
    calculated and the dimensionality of the system (1D, 2D, 3D).  For
    further details see the :ref:`Notes <thermal-speed-notes>` section
    in the `~plasmapy.formulary.parameters.thermal_speed` documentation.

    Parameters
    ----------
    T : `~numbers.Real`
        The temperature of the particle distribution, in units of kelvin.

    mass : `~numbers.Real`
        Mass of the particle in kg.

    coeff : `~numbers.Real`
        The coefficient :math:`C_o` associated with the method used for
        calculating the thermal speed, see
        :ref:`Notes <thermal-speed-notes>` section in the
        `~plasmapy.formulary.parameters.thermal_speed` documentation.

    Returns
    -------
    vth : `~numbers.Real`
        Thermal speed of the Maxwellian distribution in units of m/s.

    Examples
    --------
    >>> from plasmapy.particles import Particle
    >>> mass = Particle("p").mass.value
    >>> coeff = thermal_speed_coefficients(method="most_probable", ndim=3)
    >>> thermal_speed_lite(T=1e6, mass=mass, coeff=coeff)
    128486...
    """
    return coeff * np.sqrt(k_B_si_unitless * T / mass)


@bind_lite_func(
    thermal_speed_lite,
    attrs={"coefficients": thermal_speed_coefficients},
)
@check_relativistic
@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    mass={"can_be_negative": False, "can_be_nan": True},
)
@particles.particle_input
def thermal_speed(
    T: u.K,
    particle: Particle,
    method="most_probable",
    mass: u.kg = None,
    ndim=3,
) -> u.m / u.s:
    r"""
    Calculate the speed of thermal motion for particles with a Maxwellian
    distribution.  (See the :ref:`Notes <thermal-speed-notes>` section for
    details.)

    **Aliases:** `~plasmapy.formulary.parameters.vth_`

    **Lite Version:** `~plasmapy.formulary.parameters.thermal_speed_lite`

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        The temperature of the particle distribution, in units of kelvin or
        energy.

    particle : `~plasmapy.particles.particle_class.Particle`
        Representation of the particle species (e.g., ``"p"`` for protons,
        ``"D+"`` for deuterium, or ``"He-4 +1"`` for singly ionized helium-4).
        If no charge state information is provided, then the particles are
        assumed to be singly charged.

    method : `str`, optional
        (Default ``"most_probable"``) Method to be used for calculating the
        thermal speed. Valid values are ``"most_probable"``, ``"rms"``,
        ``"mean_magnitude"``, and ``"nrl"``.

    mass : `~astropy.units.Quantity`
        Mass override in units convertible to kg.  If given, then ``mass`` will
        be used instead of the mass value associated with ``particle``.

    ndim : `int`
        (Default ``3``) Dimensionality (1D, 2D, 3D) of space in which to
        calculate thermal speed. Valid values are ``1``, ``2``, or ``3``.

    Returns
    -------
    vth : `~astropy.units.Quantity`
        Thermal speed of the Maxwellian distribution.

    Raises
    ------
    `TypeError`
        The particle temperature is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the particle temperature is not in units of temperature or
        energy per particle.

    `ValueError`
        The particle temperature is invalid or particle cannot be used to
        identify an isotope or particle.

    Warns
    -----
    : `~plasmapy.utils.exceptions.RelativityWarning`
        If the ion sound speed exceeds 5% of the speed of light.

    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.


    .. _thermal-speed-notes:

    Notes
    -----

    There are multiple methods (or definitions) for calculating the thermal
    speed, all of which give the expression

    .. math::
        v_{th} = C_o \sqrt{\frac{k_B T}{m}}

    where :math:`T` is the temperature associated with the distribution,
    :math:`m` is the particle's mass, and :math:`C_o` is a constant of
    proportionality determined by the method in which :math:`v_{th}` is
    calculated and the dimensionality of the system (1D, 2D, 3D).  The
    :math:`C_o` used for the ``thermal_speed`` calculation is determined from
    the input arguments ``method`` and ``ndim``, and the values can be seen in
    the table below:

    .. table:: Values for :math:`C_o`
       :widths: 2 1 1 1 1
       :width: 100%

       +--------------+------------+---------------+---------------+---------------+
       | ↓ **method** | **ndim** → | ``1``         | ``2``         | ``3``         |
       +--------------+------------+---------------+---------------+---------------+
       | ``"most_probable"``       | .. math::     | .. math::     | .. math::     |
       |                           |    0          |    1          |    \sqrt{2}   |
       +--------------+------------+---------------+---------------+---------------+
       | ``"rms"``                 | .. math::     | .. math::     | .. math::     |
       |                           |    1          |    \sqrt{2}   |    \sqrt{3}   |
       +--------------+------------+---------------+---------------+---------------+
       | ``"mean_magnitude"``      | .. math::     | .. math::     | .. math::     |
       |                           |    \sqrt{2/π} |    \sqrt{π/2} |    \sqrt{8/π} |
       +--------------+------------+---------------+---------------+---------------+
       | ``"nrl"``                 | .. math::                                     |
       |                           |    1                                          |
       +--------------+------------+---------------+---------------+---------------+

    The coefficents can be directly retrieved using
    `~plasmapy.formulary.parameters.thermal_speed_coefficients`.

        .. rubric:: The Methods

        In the following discussion the Maxwellian distribution
        :math:`f(\mathbf{v})` is assumed to be 3D, but similar expressions can
        be given for 1D and 2D.

        - **Most Probable** ``method = "most_probable"``

          This method expresses the thermal speed of the distribution by expressing
          it as the most probable speed a particle in the distribution may have.
          To do this we first define another function :math:`g(v)` given by

          .. math::
             \int_{0}^{\infty} g(v) dv
                = \int_{-\infty}^{\infty} f(\mathbf{v}) d^3\mathbf{v}
                \quad \rightarrow \quad
                g(v) = 4 \pi v^2 f(v)

          then

          .. math::
             g^{\prime}(v_{th}) = \left.\frac{dg}{dv}\right|_{v_{th}} = 0\\
             \implies v_{th} = \sqrt{\frac{2 k_B T}{m}}

        - **Root Mean Square** ``method = "rms"``

          This method uses the root mean square to calculate an expression for
          the thermal speed of the particle distribution, which is given by

          .. math::
             v_{th} = \left[\int v^2 f(\mathbf{v}) d^3 \mathbf{v}\right]^{1/2}
                        = \sqrt{\frac{3 k_B T}{m}}

        - **Mean Magnitude** ``method = "mean_magnitude"``

          This method uses the mean speed of the particle distribution to
          calculate an expression for the thermal speed, which is given by

          .. math::
             v_{th} = \int |\mathbf{v}| f(\mathbf{v}) d^3 \mathbf{v}
                         = \sqrt{\frac{8 k_B T}{\pi m}}


        - **NRL Formulary** ``method = "nrl"``

          The `NRL Plasma Formulary
          <https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary>`_
          uses the square root of the Normal distribution's variance
          as the expression for thermal speed.

          .. math::
             v_{th} = σ = \sqrt{\frac{k_B T}{m}} \quad
             \text{where} \quad f(v) \sim e^{v^2 / 2 σ^2}

    Examples
    --------
    >>> from astropy import units as u
    >>> thermal_speed(5*u.eV, 'p')
    <Quantity 30949.6... m / s>
    >>> thermal_speed(1e6*u.K, particle='p')
    <Quantity 128486... m / s>
    >>> thermal_speed(5*u.eV, particle='e-')
    <Quantity 132620... m / s>
    >>> thermal_speed(1e6*u.K, particle='e-')
    <Quantity 550569... m / s>
    >>> thermal_speed(1e6*u.K, "e-", method="rms")
    <Quantity 674307... m / s>
    >>> thermal_speed(1e6*u.K, "e-", method="mean_magnitude")
    <Quantity 621251... m / s>

    For user convenience `~plasmapy.formulary.parameters.thermal_speed_coefficients`
    and `~plasmapy.formulary.parameters.thermal_speed_lite` are bound to this function
    and can be used as follows.

    >>> from plasmapy.particles import Particle
    >>> mass = Particle("p").mass.value
    >>> coeff = thermal_speed.coefficients(method="most_probable", ndim=3)
    >>> thermal_speed.lite(T=1e6, mass=mass, coeff=coeff)
    128486...
    """
    if mass is None:
        mass = particles.particle_mass(particle)

    coeff = thermal_speed_coefficients(method=method, ndim=ndim)

    speed = thermal_speed_lite(T=T.value, mass=mass.value, coeff=coeff)
    return speed * u.m / u.s


vth_ = thermal_speed
"""Alias to :func:`~plasmapy.formulary.parameters.thermal_speed`."""


@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n={"can_be_negative": False},
)
def thermal_pressure(T: u.K, n: u.m ** -3) -> u.Pa:
    r"""
    Return the thermal pressure for a Maxwellian distribution.

    **Aliases:** `pth_`

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        The particle temperature in either kelvin or energy per particle.

    n : `~astropy.units.Quantity`
        The particle number density in units convertible to m\ :sup:`-3`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> thermal_pressure(1*u.eV, 1e20/u.m**3)
    <Quantity 16.021... Pa>
    >>> thermal_pressure(10*u.eV, 1e20/u.m**3)
    <Quantity 160.21... Pa>

    Returns
    -------
    p_th : `~astropy.units.Quantity`
        Thermal pressure.

    Raises
    ------
    `TypeError`
        The temperature or number density is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the particle temperature is not in units of temperature or
        energy per particle.

    Notes
    -----
    The thermal pressure is given by:

    .. math::
        T_{th} = n k_B T
    """
    return n * k_B * T


pth_ = thermal_pressure
"""Alias to `~plasmapy.formulary.parameters.thermal_pressure`."""


@check_relativistic
@validate_quantities(
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()}
)
def kappa_thermal_speed(
    T: u.K, kappa, particle: Particle, method="most_probable"
) -> u.m / u.s:
    r"""Return the most probable speed for a particle within a Kappa
    distribution.

    **Aliases:** `vth_kappa_`

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        The particle temperature in either kelvin or energy per particle

    kappa: `float`
        The ``kappa`` parameter is a dimensionless number which sets the slope
        of the energy spectrum of suprathermal particles forming the tail
        of the Kappa velocity distribution function. ``kappa`` must be greater
        than 3/2.

    particle : `~plasmapy.particles.particle_class.Particle`
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4). If no
        charge state information is provided, then the particles are
        assumed to be singly charged.

    method : `str`, optional
        Method to be used for calculating the thermal speed. Options are
        ``'most_probable'`` (default), ``'rms'``, and ``'mean_magnitude'``.

    Returns
    -------
    V : `~astropy.units.Quantity`
        Particle thermal speed.

    Raises
    ------
    `TypeError`
        The particle temperature is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the particle temperature is not in units of temperature or
        energy per particle.

    `ValueError`
        The particle temperature is invalid or particle cannot be used to
        identify an isotope or particle.

    Warns
    -----
    : `~plasmapy.utils.exceptions.RelativityWarning`
        If the particle thermal speed exceeds 5% of the speed of light.

    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The particle thermal speed is given by:

    .. math::
        V_{th,i} = \sqrt{(2 κ - 3)\frac{2 k_B T_i}{κ m_i}}

    For more discussion on the ``'mean_magnitude'`` calculation method,
    see `PlasmaPy issue #186
    <https://github.com/PlasmaPy/PlasmaPy/issues/186>`__.

    Examples
    --------
    >>> from astropy import units as u
    >>> kappa_thermal_speed(5*u.eV, 4, 'p') # defaults to most probable
    <Quantity 24467.87... m / s>
    >>> kappa_thermal_speed(5*u.eV, 4, 'p', 'rms')
    <Quantity 37905.47... m / s>
    >>> kappa_thermal_speed(5*u.eV, 4, 'p', 'mean_magnitude')
    <Quantity 34922.98... m / s>

    See Also
    --------
    ~plasmapy.formulary.kappa_thermal_speed
    ~plasmapy.formulary.kappa_velocity_1D
    """
    # Checking thermal units
    if kappa <= 3 / 2:
        raise ValueError(
            f"Must have kappa > 3/2, instead of {kappa}, for "
            "kappa distribution function to be valid."
        )
    # different methods, as per https://en.wikipedia.org/wiki/Thermal_velocity
    vTh = thermal_speed(T=T, particle=particle, method=method)

    if method == "most_probable":
        # thermal velocity of Kappa distribution function is just Maxwellian
        # thermal speed modulated by the following factor.
        # This is only true for "most probable" case. RMS and mean
        # magnitude velocities are same as Maxwellian.
        coeff = np.sqrt((kappa - 3 / 2) / kappa)
    else:
        coeff = 1

    return vTh * coeff


vth_kappa_ = kappa_thermal_speed
"""Alias to `~plasmapy.formulary.parameters.kappa_thermal_speed`."""


@validate_quantities(
    n={"can_be_negative": False},
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def Hall_parameter(
    n: u.m ** -3,
    T: u.K,
    B: u.T,
    ion: Particle,
    particle: Particle,
    coulomb_log=None,
    V=None,
    coulomb_log_method="classical",
):
    r"""
    Calculate the ``particle`` Hall parameter for a plasma.

    The Hall parameter for plasma species :math:`s` (``particle``) is given by:

    .. math::

        β_{s} = \frac{Ω_{c s}}{ν_{s s^{\prime}}}

    where :math:`Ω_{c s}` is the gyrofrequncy for plasma species :math:`s`
    (``particle``) and :math:`ν_{s s^{\prime}}` is the collision frequency
    between plasma species :math:`s` (``particle``) and species
    :math:`s^{\prime}` (``ion``).

    **Aliases:** `betaH_`

    Parameters
    ----------
    n : `~astropy.units.quantity.Quantity`
        The number density associated with ``particle``.

    T : `~astropy.units.quantity.Quantity`
        The temperature of associated with ``particle``.

    B : `~astropy.units.quantity.Quantity`
        The magnetic field.

    ion : `~plasmapy.particles.particle_class.Particle`
        The type of ion ``particle`` is colliding with.

    particle : `~plasmapy.particles.particle_class.Particle`
        The particle species for which the Hall parameter is calculated for.
        Representation of the particle species (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, or ``'He-4 +1'`` for singly ionized helium-4).
        If no charge state information is provided, then the particles are
        assumed to be singly charged.

    coulomb_log : `float`, optional
        Preset value for the Coulomb logarithm. Used mostly for testing purposes.

    V : `~astropy.units.quantity.Quantity`
        The relative velocity between ``particle`` and ``ion``.  If not provided,
        then the ``particle`` thermal velocity is assumed
        (`~plasmapy.formulary.parameters.thermal_speed`).

    coulomb_log_method : `str`, optional
        The method by which to compute the Coulomb logarithm.
        The default method is the classical straight-line Landau-Spitzer
        method (``"classical"`` or ``"ls"``). The other 6 supported methods
        are ``"ls_min_interp"``, ``"ls_full_interp"``, ``"ls_clamp_mininterp"``,
        ``"hls_min_interp"``, ``"hls_max_interp"``, and ``"hls_full_interp"``.
        Please refer to the docstring of
        `~plasmapy.formulary.collisions.Coulomb_logarithm` for more
        information about these methods.

    See Also
    --------
    ~plasmapy.formulary.parameters.gyrofrequency
    ~plasmapy.formulary.collisions.fundamental_electron_collision_freq
    ~plasmapy.formulary.collisions.fundamental_ion_collision_freq
    ~plasmapy.formulary.collisions.Coulomb_logarithm

    Returns
    -------
    `~astropy.units.quantity.Quantity`
        Hall parameter for ``particle``.

    Notes
    -----
    * For calculating the collision frequency
      `~plasmapy.formulary.collisions.fundamental_electron_collision_freq` is used
      when ``particle`` is an electron and
      `~plasmapy.formulary.collisions.fundamental_ion_collision_freq` when
      ``particle`` is an ion.
    * The collision frequencies are calculated assuming a slowly moving
      Maxwellian distribution.

    Examples
    --------
    >>> from astropy import units as u
    >>> import pytest
    >>> Hall_parameter(1e10 * u.m**-3, 2.8e2 * u.eV, 2.3 * u.T, 'He-4 +1', 'e-')
    <Quantity 2.500...e+15>
    >>> with pytest.warns(RelativityWarning): Hall_parameter(1e10 * u.m**-3, 5.8e3 * u.eV, 2.3 * u.T, 'He-4 +1', 'e-')
    <Quantity 2.11158...e+17>
    """
    from plasmapy.formulary.collisions import (
        fundamental_electron_collision_freq,
        fundamental_ion_collision_freq,
    )

    gyro_frequency = frequencies.gyrofrequency(B, particle)
    gyro_frequency = gyro_frequency / u.radian
    if particles.Particle(particle).symbol == "e-":
        coll_rate = fundamental_electron_collision_freq(
            T, n, ion, coulomb_log, V, coulomb_log_method=coulomb_log_method
        )
    else:
        coll_rate = fundamental_ion_collision_freq(T, n, ion, coulomb_log, V)
    return gyro_frequency / coll_rate


betaH_ = Hall_parameter
"""Alias to `~plasmapy.formulary.parameters.Hall_parameter`."""


@preserve_signature
@njit
def plasma_frequency_lite(
    n: numbers.Real, mass: numbers.Real, z_mean: numbers.Real, to_hz: bool = False
) -> numbers.Real:
    r"""
    The ":term:`lite-function`" version of
    `~plasmapy.formulary.parameters.plasma_frequency`.  Performs the
    same plasma frequency calculation as
    `~plasmapy.formulary.parameters.plasma_frequency`, but is intended
    for computational use and, thus, has all data conditioning
    safeguards removed.

    Parameters
    ----------
    n : `~numbers.Real`
        Particle number density, in units of m\ :sup:`-3`.

    mass : `~numbers.Real`
        Mass of the particle, in units of kg.

    z_mean : `~numbers.Real`
        The average ionization (arithmetic mean) for the particle
        species in the plasma.  For example, a proton would have a value
        of ``z_mean=1``.

    to_hz : `bool`
        (Default `False`) Set `True` to apply the factor of
        :math:`1/2π` and return a value in units of Hz.

    Returns
    -------
    wp : `~numbers.Real`
        The particle plasma frequency in radians per second.  Setting
        keyword ``to_hz=True`` will apply the factor of :math:`1/2π`
        and yield a value in Hz.

    Notes
    -----

    The particle plasma frequency is

    .. math::
        ω_{p} = Z |e| \sqrt{\frac{n}{\epsilon_0 m}}

    where :math:`m` is the mass of the particle, :math:`e` is the
    fundamental unit of charge, :math:`Z` is the average charge state
    ``z_mean`` of the particle species, :math:`n` is the particle number
    density.  This form of the plasma frequency has units of
    radians / s, but when using the ``to_hz`` keyword a factor of
    :math:`1/2π` will be applied to give a value in Hz.

    Examples
    --------

    >>> from plasmapy.particles import Particle
    >>> mass = Particle("p").mass.value
    >>> plasma_frequency_lite(n=1e19, mass=mass, z_mean=1)
    416329...
    >>> plasma_frequency_lite(n=1e19, mass=mass, z_mean=1, to_hz=True)
    662608...
    """
    omega_p = z_mean * e_si_unitless * np.sqrt(n / (eps0_si_unitless * mass))

    if to_hz:
        return omega_p / (2.0 * np.pi)

    return omega_p


@bind_lite_func(plasma_frequency_lite)
@validate_quantities(
    n={"can_be_negative": False},
    validations_on_return={
        "units": [u.rad / u.s, u.Hz],
        "equivalencies": [(u.cy / u.s, u.Hz)],
    },
)
@angular_freq_to_hz
def plasma_frequency(n: u.m ** -3, particle: Particle, z_mean=None) -> u.rad / u.s:
    r"""Calculate the particle plasma frequency.

    **Aliases:** `wp_`

    **Lite Version:** `~plasmapy.formulary.parameters.plasma_frequency_lite`

    Parameters
    ----------
    n : `~astropy.units.Quantity`
        Particle number density in units convertible to m\ :sup:`-3`.

    particle : `~plasmapy.particles.particle_class.Particle`
        Representation of the particle species (e.g., ``"p"`` for
        protons, ``"D+"`` for deuterium, or ``"He-4 +1"`` for singly
        ionized helium-4). If no charge state information is provided,
        then the particles are assumed to be singly charged.

    z_mean : `~numbers.Real`, optional
        The average ionization (arithmetic mean) for the particle
        species in the plasma.  Typically the charge state will be
        dervied from the ``particle`` argument, but this keyword will
        override that behavior.

    Returns
    -------
    omega_p : `~astropy.units.Quantity`
        The particle plasma frequency in radians per second.  Setting
        keyword ``to_hz=True`` will apply the factor of :math:`1/2π`
        and yield a value in Hz.

    Raises
    ------
    `TypeError`
        If ``n`` is not a `~astropy.units.Quantity` or particle is not
        of an appropriate type.

    `~astropy.units.UnitConversionError`
        If ``n`` is not in correct units.

    `ValueError`
        If ``n`` contains invalid values or particle cannot be used to
        identify an particle or isotope.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The particle plasma frequency is

    .. math::
        ω_{p} = Z |e| \sqrt{\frac{n}{\epsilon_0 m}}

    where :math:`m` is the mass of the particle, :math:`e` is the
    fundamental unit of charge, :math:`Z` is the average charge state
    ``z_mean`` of the particle species, :math:`n` is the particle number
    density.  This form of the plasma frequency has units of
    radians / s, but using the ``to_hz`` will apply the factor of
    :math:`1/2π` to give a value in Hz.

    Examples
    --------
    >>> from astropy import units as u
    >>> plasma_frequency(1e19*u.m**-3, particle='p')
    <Quantity 4.16329...e+09 rad / s>
    >>> plasma_frequency(1e19*u.m**-3, particle='p', to_hz=True)
    <Quantity 6.62608...e+08 Hz>
    >>> plasma_frequency(1e19*u.m**-3, particle='D+')
    <Quantity 2.94462...e+09 rad / s>
    >>> plasma_frequency(1e19*u.m**-3, 'e-')
    <Quantity 1.78398...e+11 rad / s>
    >>> plasma_frequency(1e19*u.m**-3, 'e-', to_hz=True)
    <Quantity 2.83930...e+10 Hz>

    For user convienence
    `~plasmapy.formulary.parameters.plasma_frequency_lite` is bound to
    this function and can be used as follows.

    >>> from plasmapy.particles import Particle
    >>> mass = Particle("p").mass.value
    >>> plasma_frequency.lite(n=1e19, mass=mass, z_mean=1)
    416329...
    >>> plasma_frequency.lite(n=1e19, mass=mass, z_mean=1, to_hz=True)
    662608...
    """

    try:
        m = particles.particle_mass(particle).value

        if z_mean is None:
            # warnings.warn("No z_mean given, defaulting to atomic charge",
            #               PhysicsWarning)
            try:
                Z = particles.charge_number(particle)
            except Exception:
                Z = 1
        else:
            # using user provided average ionization
            Z = z_mean
        Z = np.abs(Z)
        # TODO REPLACE WITH Z = np.abs(_grab_charge(particle, z_mean)), some bugs atm
    except Exception:
        raise ValueError(f"Invalid particle, {particle}, in plasma_frequency.")

    omega_p = plasma_frequency_lite(n=n, mass=m, z_mean=Z) * u.rad / u.s
    return omega_p


wp_ = plasma_frequency
"""Alias to `~plasmapy.formulary.parameters.plasma_frequency`."""


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def Debye_number(T_e: u.K, n_e: u.m ** -3) -> u.dimensionless_unscaled:
    r"""Return the number of electrons within a sphere with a radius
    of the Debye length.

    **Aliases:** `nD_`

    Parameters
    ----------
    T_e : `~astropy.units.Quantity`
        Electron temperature.

    n_e : `~astropy.units.Quantity`
        Electron number density.

    Raises
    ------
    `TypeError`
        If either argument is not a `~astropy.units.Quantity`.

    `astropy.units.UnitConversionError`
        If either argument is in incorrect units.

    `ValueError`
        If either argument contains invalid values.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Returns
    -------
    N_D : `~astropy.units.Quantity`
        Number of electrons within a sphere with a radius of the Debye
        length.

    Notes
    -----
    The Debye number is the number of electrons contained within a
    sphere with a radius of a Debye length and is given by

    .. math::
        N_D = \frac{4π}{3} n_e λ_D^3

    The Debye number is also known as the plasma parameter.

    Collective behavior requires :math:`N_D ≫ 1`\ .

    See Also
    --------
    ~plasmapy.formulary.lengths.Debye_length

    Examples
    --------
    >>> from astropy import units as u
    >>> Debye_number(5e6*u.K, 5e9*u.cm**-3)
    <Quantity 2.17658...e+08>

    """

    lambda_D = lengths.Debye_length(T_e, n_e)
    N_D = (4 / 3) * np.pi * n_e * lambda_D ** 3

    return N_D


nD_ = Debye_number
"""Alias to `~plasmapy.formulary.parameters.Debye_number`."""


@validate_quantities
def magnetic_pressure(B: u.T) -> u.Pa:
    r"""
    Calculate the magnetic pressure.

    **Aliases:** `pmag_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field in units convertible to tesla.

    Returns
    -------
    p_B : `~astropy.units.Quantity`
        The magnetic pressure in units in pascals (newtons per square meter).

    Raises
    ------
    `TypeError`
        If the input is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the input is not in units convertible to tesla.

    `ValueError`
        If the magnetic field strength is not a real number between
        :math:`±∞`\ .

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The magnetic pressure is given by:

    .. math::
        p_B = \frac{B^2}{2 μ_0}

    The motivation behind having two separate functions for magnetic
    pressure and magnetic energy density is that it allows greater
    insight into the physics that are being considered by the user and
    thus more readable code.

    See Also
    --------
    magnetic_energy_density : returns an equivalent `~astropy.units.Quantity`,
        except in units of joules per cubic meter.

    Examples
    --------
    >>> from astropy import units as u
    >>> magnetic_pressure(0.1*u.T).to(u.Pa)
    <Quantity 3978.87... Pa>

    """
    return (B ** 2) / (2 * mu0)


pmag_ = magnetic_pressure
"""Alias to `~plasmapy.formulary.parameters.magnetic_pressure`."""


@validate_quantities
def magnetic_energy_density(B: u.T) -> u.J / u.m ** 3:
    r"""
    Calculate the magnetic energy density.

    **Aliases:** `ub_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field in units convertible to tesla.

    Returns
    -------
    E_B : `~astropy.units.Quantity`
        The magnetic energy density in units of joules per cubic meter.

    Raises
    ------
    `TypeError`
        If the input is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the input is not in units convertible to tesla.

    `ValueError`
        If the magnetic field strength does not have an appropriate.
        value.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed

    Notes
    -----
    The magnetic energy density is given by:

    .. math::
        E_B = \frac{B^2}{2 μ_0}

    The motivation behind having two separate functions for magnetic
    pressure and magnetic energy density is that it allows greater
    insight into the physics that are being considered by the user and
    thus more readable code.

    See Also
    --------
    magnetic_pressure : Returns an equivalent `~astropy.units.Quantity`,
        except in units of pascals.

    Examples
    --------
    >>> from astropy import units as u
    >>> magnetic_energy_density(0.1*u.T)
    <Quantity 3978.87... J / m3>

    """
    return magnetic_pressure(B)


ub_ = magnetic_energy_density
"""Alias to `~plasmapy.formulary.parameters.magnetic_energy_density`."""


@validate_quantities(
    n_e={"can_be_negative": False},
    validations_on_return={
        "units": [u.rad / u.s, u.Hz],
        "equivalencies": [(u.cy / u.s, u.Hz)],
    },
)
@angular_freq_to_hz
def upper_hybrid_frequency(B: u.T, n_e: u.m ** -3) -> u.rad / u.s:
    r"""
    Return the upper hybrid frequency.

    **Aliases:** `wuh_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to tesla.

    n_e : `~astropy.units.Quantity`
        The electron number density.

    Returns
    -------
    omega_uh : `~astropy.units.Quantity`
        The upper hybrid frequency in radians per second.

    Raises
    ------
    `TypeError`
        If either of ``B`` or ``n_e`` is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If either of ``B`` or ``n_e`` is in incorrect units.

    `ValueError`
        If either of ``B`` or ``n_e`` contains invalid values or are of
        incompatible dimensions.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The upper hybrid frequency is given through the relation

    .. math::
        ω_{uh}^2 = ω_{ce}^2 + ω_{pe}^2

    where :math:`ω_{ce}` is the electron gyrofrequency and
    :math:`ω_{pe}` is the electron plasma frequency.

    The upper hybrid frequency is a resonance for electromagnetic
    waves in magnetized plasmas, namely for the X-mode. These are
    waves with their wave electric field being perpendicular to
    the background magnetic field. In the cold plasma model, i.e.
    without any finite temperature effects, the resonance acts
    merely as a resonance such that power can be deposited there.
    If finite temperature effects are considered, mode conversion
    can occur at the upper hybrid resonance, coupling to the
    electrostatic electron Bernstein wave.

    Examples
    --------
    >>> from astropy import units as u
    >>> upper_hybrid_frequency(0.2*u.T, n_e=5e19*u.m**-3)
    <Quantity 4.00459...e+11 rad / s>
    >>> upper_hybrid_frequency(0.2*u.T, n_e=5e19*u.m**-3, to_hz = True)
    <Quantity 6.37350...e+10 Hz>

    """
    omega_pe = plasma_frequency(n=n_e, particle="e-")
    omega_ce = frequencies.gyrofrequency(B, "e-")
    omega_uh = np.sqrt(omega_pe ** 2 + omega_ce ** 2)

    return omega_uh


wuh_ = upper_hybrid_frequency
"""Alias to `~plasmapy.formulary.parameters.upper_hybrid_frequency`."""


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    B={"can_be_negative": False},
)
def Bohm_diffusion(T_e: u.K, B: u.T) -> u.m ** 2 / u.s:
    r"""
    Return the Bohm diffusion coefficient.

    The Bohm diffusion coefficient was conjectured to follow Bohm model
    of the diffusion of plasma across a magnetic field and describe the
    diffusion of early fusion energy machines :cite:p:`bohm:1949`. The
    rate predicted by Bohm diffusion is much higher than classical
    diffusion, and if there were no exceptions, magnetically confined
    fusion would be impractical.

    .. math::

        D_B = \frac{1}{16} \frac{k_B T}{e B}

    where :math:`k_B` is the Boltzmann constant
    and :math:`e` is the fundamental charge.

    **Aliases:** `DB_`

    Parameters
    ----------
    T_e : `~astropy.units.Quantity`
        The electron temperature.

    B : `~astropy.units.Quantity`
        The magnitude of the magnetic field in the plasma.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Raises
    ------
    `TypeError`
        ``T_e`` is not a `~astropy.units.Quantity` and cannot be
        converted into one.

    `~astropy.units.UnitConversionError`
        If ``T_e`` is not in appropriate units.

    Examples
    --------
    >>> import astropy.units as u
    >>> T_e = 5000 * u.K
    >>> B = 10 * u.T
    >>> Bohm_diffusion(T_e, B)
    <Quantity 0.00269292 m2 / s>
    >>> T_e = 50 * u.eV
    >>> B = 10 * u.T
    >>> Bohm_diffusion(T_e, B)
    <Quantity 0.3125 m2 / s>

    Returns
    -------
    D_B : `~astropy.units.Quantity`
        The Bohm diffusion coefficient in meters squared per second.

    """
    D_B = k_B * T_e / (16 * e * B)
    return D_B


DB_ = Bohm_diffusion
"""Alias to `~plasmapy.formulary.parameters.Bohm_diffusion`."""
