"""Functions to calculate fundamental plasma length parameters."""

__all__ = ["Debye_length", "gyroradius", "inertial_length"]
__aliases__ = ["cwp_", "lambdaD_", "rc_", "rhoc_"]

import warnings
from numbers import Integral, Real
from typing import Optional

import astropy.units as u
import numpy as np
from astropy.constants.si import c, e, eps0, k_B

from plasmapy.formulary import frequencies, speeds
from plasmapy.formulary.relativity import RelativisticBody
from plasmapy.particles import ParticleLike, particle_input
from plasmapy.utils.decorators import validate_quantities

__all__ += __aliases__


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def Debye_length(T_e: u.Quantity[u.K], n_e: u.Quantity[u.m**-3]) -> u.Quantity[u.m]:
    r"""
    Calculate the exponential scale length for charge screening in an
    electron plasma with stationary ions.

    The Debye length is given by

    .. math::
        λ_D = \sqrt{\frac{ε_0 k_B T_e}{n_e q_e^2}},

    where :math:`n_e` is the electron number density, :math:`T_e` is the
    electron temperature, :math:`k_B` is the Boltzmann constant,
    :math:`q_e` is the elementary charge, and :math:`ε_0` is the vacuum
    permittivity.

    **Aliases:** `lambdaD_`

    Parameters
    ----------
    T_e : `~astropy.units.Quantity`
        Electron temperature.

    n_e : `~astropy.units.Quantity`
        Electron number density.

    Returns
    -------
    lambda_D : `~astropy.units.Quantity`
        The Debye length in meters.

    Raises
    ------
    `TypeError`
        If either argument is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If either argument is in incorrect units.

    `ValueError`
        If either argument contains invalid values.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    Plasmas will generally be quasineutral on length scales
    significantly longer than the Debye length.

    The electrical potential will drop by a factor of
    :math:`∼\frac{1}{e}` every Debye length away from the vicinity of a
    charged particle.

    See Also
    --------
    ~plasmapy.formulary.dimensionless.Debye_number

    Examples
    --------
    >>> import astropy.units as u
    >>> Debye_length(5e6 * u.K, 5e15 * u.m**-3)
    <Quantity 0.002182... m>
    """
    return np.sqrt(eps0 * k_B * T_e / (n_e * e**2))


lambdaD_ = Debye_length
"""Alias to `~plasmapy.formulary.lengths.Debye_length`."""


@validate_quantities(
    Vperp={"can_be_nan": True},  # none_shall_pass
    T={
        "can_be_nan": True,
        "equivalencies": u.temperature_energy(),
        "none_shall_pass": True,
    },
    validations_on_return={"equivalencies": u.dimensionless_angles()},
)
@particle_input(any_of={"charged", "uncharged"})
def gyroradius(  # noqa: C901
    B: u.Quantity[u.T],
    particle: ParticleLike,
    *,
    Vperp: u.Quantity[u.m / u.s] = np.nan * u.m / u.s,
    T: u.Quantity[u.K] = None,
    lorentzfactor=np.nan,
    relativistic: bool = True,
    mass_numb: Optional[Integral] = None,
    Z: Optional[Real] = None,
) -> u.Quantity[u.m]:
    r"""
    Calculate the radius of circular motion for a charged particle in a
    uniform magnetic field (including relativistic effects by default).

    **Aliases:** `rc_`, `rhoc_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to tesla.

    particle : |particle-like|
        Representation of the particle species (e.g., ``"p+"`` for
        protons, ``"D+"`` for a deuteron, or ``"He-4 1+"`` for singly
        ionized helium-4).

    Vperp : `~astropy.units.Quantity`, |keyword-only|, optional
        The component of particle velocity that is perpendicular to
        the magnetic field in units convertible to meters per second.

    T : `~astropy.units.Quantity`, |keyword-only|, optional
        The particle temperature in units convertible to kelvin or
        electron-volts. If provided, the perpendicular velocity gets set
        to the most probable *non-relativistic* thermal velocity for
        that particle at this temperature. Cannot be provided if
        ``Vperp`` is provided.

    lorentzfactor : `float` or `~numpy.ndarray`, |keyword-only|, optional
        The :wikipedia:`Lorentz factor` of the particle corresponding to
        the direction perpendicular to the magnetic field. Cannot be
        provided if ``Vperp`` or ``T`` is provided.

    relativistic : `bool`, |keyword-only|, default: `True`
        If `True`, the relativistic formula for gyroradius will be used.
        If `False`, the non-relativistic formula will be used.

    Returns
    -------
    r_L : `~astropy.units.Quantity`
        The particle gyroradius in units of meters. This
        `~astropy.units.Quantity` will be based on either the
        perpendicular component of particle velocity as inputted, or
        the most probable speed for a particle within a Maxwellian
        distribution for the particle temperature.

    Other Parameters
    ----------------
    mass_numb : integer, |keyword-only|, optional
        The mass number, if not provided in ``particle``.

    Z : real number, |keyword-only|, optional
        The |charge number|, if not provided in ``particle``.

    Raises
    ------
    `~astropy.units.UnitConversionError`
        If a |Quantity| argument has units of an incorrect physical
        type.

    `ValueError`
        If any argument contains invalid values.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        Issued if any of ``B``, ``Vperp``, or ``T`` do not have units,
        in which case SI units will be assumed.

    Warnings
    --------
    The Lorentz factor can be inferred from ``Vperp`` or ``T`` but near
    the speed of light, this can lead to rounding errors. For very high
    values of the Lo

    Notes
    -----
    The relativistic :wikipedia:`gyroradius` for a particle of species
    :math:`s` is given by

    .. math::

        r_{L,s} = \frac{γ m_s V_{⟂,s}}{ |q_s| B}

    where :math:`V_⟂` is the component of particle velocity that is
    perpendicular to the magnetic field, :math:`m_s` is the particle
    mass, :math:`q_s` is the particle charge, :math:`B` is the magnetic
    field magnitude, and :math:`γ` is the :wikipedia:`Lorentz factor`.

    In the non-relativistic limit, the gyroradius reduces to

    .. math::

        r_{Ls} = \frac{V_{⟂,s}{ω_{c,s}}

    where :math:`ω_{c,s}` is the particle gyrofrequency. To turn off
    relativistic effects, set the ``relativistic`` keyword to `False`.

    The gyroradius is sometimes called the Larmor radius, cyclotron
    radius, or radius of gyration.

    Examples
    --------
    >>> import astropy.units as u
    >>> from plasmapy.formulary import gyroradius
    >>> from astropy.constants import c

    Let's estimate the proton gyroradius in the solar corona.

    >>> gyroradius(B=0.2 * u.T, particle="p+", T=1e6 * u.K)
    <Quantity 0.0067... m>

    Let's estimate the gyroradius of a deuteron and a triton in ITER by
    providing the characteristic thermal energy per particle,
    :math:`k_B T`, to ``T``.

    >>> gyroradius(B=5 * u.T, particle=["D+", "T+"], T=13 * u.keV)
    <Quantity [0.00465..., 0.00570...] m>

    Relativistic effects are included by default, but can be turned
    off using the ``relativistic`` parameter. Let's use this in the
    calculation of the gyroradius of a cosmic ray in the interstellar
    medium (ISM). We will provide the magnetic field in units of
    microgauss (μG).

    >>> gyroradius(B=10 * u.uG, particle="p+", Vperp=0.99 * c)
    <Quantity 2.19642688e+10 m>
    >>> gyroradius(B=10 * u.uG, particle="p+", Vperp=0.99 * c, relativistic=False)
    <Quantity 3.09844141e+09 m>

    Let's calculate the gyroradius of a much higher energy cosmic ray
    in the ISM using ``lorentzfactor``.

    >>> gyroradius(B=10 * u.uG, particle="p+", lorentzfactor=3e6).to("pc")
    <Quantity 0.30428378 pc>
    """

    # Define helper functions for input processing and gyroradius calculation

    # check 1: ensure either Vperp or T invalid, keeping in mind that
    # the underlying values of the astropy quantity may be numpy arrays
    def _raise_error_if_both_vperp_and_t_are_given(isfinite_Vperp, isfinite_T):
        if np.any(np.logical_and(isfinite_Vperp, isfinite_T)):
            raise ValueError(
                "Must give Vperp or T, but not both, as arguments to gyroradius"
            )

    def _raise_error_if_lorentzfactor_not_scalar(lorentzfactor):
        if nans_in_both_T_and_Vperp and not np.isscalar(lorentzfactor):
            raise ValueError(
                "Inferring velocity(s) from more than one Lorentz "
                "factor is not currently supported"
            )

    def _calculate_vperp_from_lorentzfactor(
        isfinite_Vperp, Vperp, particle, lorentzfactor, relativistic
    ):
        if relativistic and nans_in_both_T_and_Vperp:
            Vperp = np.copy(Vperp)
            rbody = RelativisticBody(particle, lorentz_factor=lorentzfactor)
            Vperp[~isfinite_Vperp] = rbody.velocity
        return Vperp

    def _warn_if_lorentz_factor_and_relativistic(
        isfinite_lorentzfactor, relativistic
    ) -> None:
        if np.any(isfinite_lorentzfactor) and relativistic:
            warnings.warn(
                "lorentzfactor is given along with Vperp or T, will lead "
                "to inaccurate predictions unless they correspond"
            )

    # check 2: get Vperp as the thermal speed if is not already a valid input
    def get_Vperp(T, Vperp, particle, isfinite_T, isfinite_Vperp):
        is_scalar_Vperp = np.isscalar(Vperp.value)
        is_scalar_T = np.isscalar(T.value)

        if is_scalar_Vperp and is_scalar_T:
            # both T and Vperp are scalars
            # we know exactly one of them is nan from check 1
            return speeds.thermal_speed(T, particle=particle) if isfinite_T else Vperp
            # T is valid, so use it to determine Vperp
            # else: Vperp is already valid, do nothing

        elif is_scalar_Vperp:
            # this means either Vperp must be nan, or T must be an array of all nan,
            # or else we couldn't have gotten through check 1
            return (
                np.repeat(Vperp, len(T))
                if isfinite_Vperp
                # Vperp is valid, T is a vector that is all nan
                else speeds.thermal_speed(T, particle=particle)
                # normal case where Vperp is scalar nan and T is valid array
            )

        elif is_scalar_T:  # only Vperp is an array
            # this means either T must be nan, or V_perp must be an array of all nan,
            # or else we couldn't have gotten through check 1
            return (
                speeds.thermal_speed(np.repeat(T, len(Vperp)), particle=particle)
                if isfinite_T
                # T is valid, V_perp is an array of all nan
                else Vperp
                # else: normal case where T is scalar nan and Vperp is already a valid
                # array so, do nothing
            )
        else:  # both T and Vperp are arrays
            # we know all the elementwise combinations have one nan and one finite,
            # due to check 1 use the valid Vperps, and replace the others with those
            # calculated from T
            Vperp = Vperp.copy()
            Vperp[isfinite_T] = speeds.thermal_speed(T[isfinite_T], particle=particle)
            return Vperp

    def get_result(lorentzfactor, isfinite_lorentzfactor, particle, Vperp, omega_ci):
        if np.all(isfinite_lorentzfactor):
            return lorentzfactor * np.abs(Vperp) / omega_ci

        elif not np.any(isfinite_lorentzfactor):
            lorentzfactor = RelativisticBody(particle, V=Vperp).lorentz_factor
            return lorentzfactor * np.abs(Vperp) / omega_ci

        else:
            # the lorentzfactor is neither completely valid nor completely invalid,
            # so we have to correct the missing parts, note that we don't actually
            # have to check if it is a scalar since scalars cannot be partially valid
            rbody = RelativisticBody(particle, V=Vperp)
            lorentzfactor = np.copy(lorentzfactor)
            lorentzfactor[~isfinite_lorentzfactor] = rbody.lorentz_factor[
                ~isfinite_lorentzfactor
            ]
            return lorentzfactor * np.abs(Vperp) / omega_ci

    # Initial setup and input validation
    if not relativistic and not np.isnan(lorentzfactor):
        raise ValueError("Lorentz factor is provided but relativistic is set to false")
    lorentzfactor = lorentzfactor if relativistic else 1.0

    if T is None:
        T = np.nan * u.K

    isfinite_T = np.isfinite(T)
    isfinite_Vperp = np.isfinite(Vperp)
    isfinite_lorentzfactor = np.isfinite(lorentzfactor)
    nans_in_both_T_and_Vperp = np.any(np.logical_not(isfinite_T)) and np.any(
        np.logical_not(isfinite_Vperp)
    )

    # Call defined functions to process inputs
    _raise_error_if_both_vperp_and_t_are_given(isfinite_Vperp, isfinite_T)

    if nans_in_both_T_and_Vperp:
        _raise_error_if_lorentzfactor_not_scalar(lorentzfactor)
        Vperp = _calculate_vperp_from_lorentzfactor(
            isfinite_Vperp, Vperp, particle, lorentzfactor, relativistic
        )

    _warn_if_lorentz_factor_and_relativistic(isfinite_lorentzfactor, relativistic)

    Vperp = get_Vperp(T, Vperp, particle, isfinite_T, isfinite_Vperp)

    omega_ci = frequencies.gyrofrequency(B, particle)

    return get_result(lorentzfactor, isfinite_lorentzfactor, particle, Vperp, omega_ci)


rc_ = gyroradius
"""Alias to `~plasmapy.formulary.lengths.gyroradius`."""

rhoc_ = gyroradius
"""Alias to `~plasmapy.formulary.lengths.gyroradius`."""


@validate_quantities(
    n={"can_be_negative": False},
    validations_on_return={"equivalencies": u.dimensionless_angles()},
)
@particle_input(require="charged")
def inertial_length(
    n: u.Quantity[u.m**-3],
    particle: ParticleLike,
    *,
    mass_numb: Optional[Integral] = None,
    Z: Optional[Real] = None,
) -> u.Quantity[u.m]:
    r"""
    Calculate a charged particle's inertial length.

    The inertial length of a particle of species :math:`s` is given by

    .. math::

        d = \frac{c}{ω_{ps}}

    The inertial length is the characteristic length scale for a
    particle to be accelerated in a plasma. The Hall effect becomes
    important on length scales shorter than the ion inertial length.

    **Aliases:** `cwp_`

    Parameters
    ----------
    n : `~astropy.units.Quantity`
        Particle number density in units convertible to m\ :sup:`-3`\ .

    particle : |particle-like|
        Representation of the particle species (e.g., 'p+' for protons,
        'D+' for deuterium, or 'He-4 +1' for singly ionized helium-4).

    Returns
    -------
    d : `~astropy.units.Quantity`
        The particle's inertial length in meters.

    Other Parameters
    ----------------
    mass_numb : integer, |keyword-only|, optional
        The mass number, if not provided in ``particle``.

    Z : real number, |keyword-only|, optional
        The |charge number|, if not provided in ``particle``.

    Raises
    ------
    `TypeError`
        If ``n`` is not a `~astropy.units.Quantity` or ``particle`` is
        not a string.

    `~astropy.units.UnitConversionError`
        If ``n`` is not in units of a number density.

    `ValueError`
        The particle density does not have an appropriate value.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided and SI units are assumed.

    Notes
    -----
    The inertial length is also known as the skin depth.

    Examples
    --------
    >>> import astropy.units as u
    >>> inertial_length(5 * u.m**-3, "He+")
    <Quantity 2.02985...e+08 m>
    >>> inertial_length(5 * u.m**-3, "e-")
    <Quantity 2376534.75... m>
    """
    omega_p = frequencies.plasma_frequency(n, particle=particle)
    return c / omega_p


cwp_ = inertial_length
"""
Alias to `~plasmapy.formulary.lengths.inertial_length`.

* Name is shorthand for :math:`c / ω_p`.
"""
