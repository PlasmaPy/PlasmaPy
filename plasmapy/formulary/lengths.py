"""Functions to calculate fundamental plasma length parameters."""

__all__ = ["Debye_length", "gyroradius", "inertial_length"]
__aliases__ = ["cwp_", "lambdaD_", "rc_", "rhoc_"]

import astropy.units as u
import numpy as np
import warnings

from astropy.constants.si import c, e, eps0, k_B

from plasmapy.formulary import frequencies, speeds
from plasmapy.formulary.relativity import RelativisticBody
from plasmapy.particles import particle_input, ParticleLike
from plasmapy.utils.decorators import validate_quantities

__all__ += __aliases__


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def Debye_length(T_e: u.K, n_e: u.m**-3) -> u.m:
    r"""
    Calculate the characteristic decay length for electric fields due to
    charge screening.

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
    The Debye length is the exponential scale length for charge
    screening and is given by

    .. math::
        λ_D = \sqrt{\frac{ε_0 k_B T_e}{n_e e^2}}

    for an electron plasma with nearly stationary ions.

    The electrical potential will drop by a factor of :math:`1/e`
    every Debye length.

    Plasmas will generally be quasineutral on length scales
    significantly larger than the Debye length.

    See Also
    --------
    ~plasmapy.formulary.dimensionless.Debye_number

    Examples
    --------
    >>> from astropy import units as u
    >>> Debye_length(5e6*u.K, 5e15*u.m**-3)
    <Quantity 0.002182... m>
    """
    return np.sqrt(eps0 * k_B * T_e / (n_e * e**2))


lambdaD_ = Debye_length
"""Alias to `~plasmapy.formulary.lengths.Debye_length`."""


@validate_quantities(
    Vperp={"can_be_nan": True},
    T={
        "can_be_nan": True,
        "equivalencies": u.temperature_energy(),
        "none_shall_pass": True,
    },
    validations_on_return={"equivalencies": u.dimensionless_angles()},
)
@particle_input(any_of={"charged", "uncharged"})
def gyroradius(  # noqa: C901
    B: u.T,
    particle: ParticleLike,
    *,
    Vperp: u.m / u.s = np.nan * u.m / u.s,
    T: u.K = None,
    lorentzfactor=np.nan,
    relativistic: bool = True,
    mass_numb=None,
    Z=None,
) -> u.m:
    r"""
    Return the particle gyroradius.

    **Aliases:** `rc_`, `rhoc_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to tesla.

    particle : `~plasmapy.particles.particle_class.Particle`
        Representation of the particle species (e.g., ``'p'`` for
        protons, ``'D+'`` for deuterium, or ``'He-4 +1'`` for singly
        ionized helium-4). If no charge state information is
        provided, then the particles are assumed to be singly charged.

    Vperp : `~astropy.units.Quantity`, optional, |keyword-only|
        The component of particle velocity that is perpendicular to
        the magnetic field in units convertible to meters per second.

    T : `~astropy.units.Quantity`, optional, |keyword-only|
        The particle temperature in units convertible to kelvin.

    lorentzfactor : `float` or `~numpy.ndarray`, optional, |keyword-only|
        The Lorentz factor for the particles, use for high precision.

    relativistic : `bool`, optional, |keyword-only|
        Whether or not you want to use a relativistic approximation.
        `True` by default.

    mass_numb : integer, |keyword-only|, optional
        The mass number associated with ``particle``.

    Z : real number, |keyword-only|, optional
        The charge number associated with ``particle``.

    Returns
    -------
    r_Li : `~astropy.units.Quantity`
        The particle gyroradius in units of meters. This
        `~astropy.units.Quantity` will be based on either the
        perpendicular component of particle velocity as inputted, or
        the most probable speed for a particle within a Maxwellian
        distribution for the particle temperature. It is
        relativistically accurate.

    Raises
    ------
    `TypeError`
        The arguments are of an incorrect type.

    `~astropy.units.UnitConversionError`
        The arguments do not have appropriate units.

    `ValueError`
        If any argument contains invalid values.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    One but not both of ``Vperp`` and ``T`` must be inputted.

    ``lorentzfactor`` can be inferred from ``Vperp`` or ``T`` but near
    the speed of light, this can lead to rounding errors.

    If any of ``B``, ``Vperp``, or ``T`` is a number rather than a
    `~astropy.units.Quantity`, then SI units will be assumed and a
    warning will be raised.

    The particle gyroradius is also known as the particle Larmor
    radius and is given by

    .. math::
        r_{Li} = \frac{γ V_⟂}{ω_{ci}}

    where :math:`V_⟂` is the component of particle velocity that is
    perpendicular to the magnetic field, :math:`ω_{ci}` is the
    particle gyrofrequency, and :math:`γ` is the Lorentz factor.  If a
    temperature is provided, then :math:`V_⟂` will be the most
    probable thermal velocity of a particle at that temperature. The
    ``relativistic`` keyword can be set to `False` to avoid the
    relativistic correction.

    Examples
    --------
    >>> from astropy import units as u
    >>> gyroradius(0.2*u.T, particle='p+', T=1e5*u.K)
    <Quantity 0.002120... m>
    >>> gyroradius(0.2*u.T, particle='p+', T=1e5*u.K)
    <Quantity 0.002120... m>
    >>> gyroradius(5*u.uG, particle='alpha', T=1*u.eV)
    <Quantity 288002.38... m>
    >>> gyroradius(400*u.G, particle='Fe+++', Vperp=1e7*u.m/u.s)
    <Quantity 48.25815... m>
    >>> gyroradius(B=0.01*u.T, particle='e-', T=1e6*u.K)
    <Quantity 0.003130... m>
    >>> gyroradius(0.01*u.T, 'e-', Vperp=1e6*u.m/u.s)
    <Quantity 0.000568... m>
    >>> gyroradius(0.2*u.T, 'e-', T=1e5*u.K)
    <Quantity 4.94957...e-05 m>
    >>> gyroradius(5*u.uG, 'e-', T=1*u.eV)
    <Quantity 6744.27... m>
    >>> gyroradius(400*u.G, 'e-', Vperp=1e7*u.m/u.s)
    <Quantity 0.001422... m>
    >>> gyroradius(400*u.G, 'e-', Vperp=1e7*u.m/u.s, lorentzfactor=1.0)
    <Quantity 0.001421... m>
    >>> gyroradius(400*u.G, 'e-', Vperp=1e7*u.m/u.s, relativistic=False)
    <Quantity 0.001421... m>
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
                "Inferring velocity(s) from more than one Lorentz factor is not currently supported"
            )

    def _calculate_vperp_from_lorentzfactor(
        isfinite_Vperp, Vperp, particle, lorentzfactor, relativistic
    ):
        if relativistic and nans_in_both_T_and_Vperp:
            Vperp = np.copy(Vperp)
            rbody = RelativisticBody(particle, lorentz_factor=lorentzfactor)
            Vperp[~isfinite_Vperp] = rbody.velocity
        return Vperp

    def _warn_if_lorentz_factor_and_relativistic(isfinite_lorentzfactor, relativistic):
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
    n: u.m**-3, particle: ParticleLike, *, mass_numb=None, Z=None
) -> u.m:
    r"""
    Calculate a charged particle's inertial length.

    **Aliases:** `cwp_`

    Parameters
    ----------
    n : `~astropy.units.Quantity`
        Particle number density in units convertible to m\ :sup:`-3`\ .

    particle : `~plasmapy.particles.particle_class.Particle`
        Representation of the particle species (e.g., 'p+' for protons,
        'D+' for deuterium, or 'He-4 +1' for singly ionized helium-4).

    mass_numb : integer, |keyword-only|, optional
        The mass number associated with ``particle``.

    Z : real number, |keyword-only|, optional
        The charge number associated with ``particle``.

    Returns
    -------
    d : `~astropy.units.Quantity`
        The particle's inertial length in meters.

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
    The inertial length of a particle of species :math:`s` is given by

    .. math::

        d = \frac{c}{ω_{ps}}

    The inertial length is the characteristic length scale for a
    particle to be accelerated in a plasma. The Hall effect becomes
    important on length scales shorter than the ion inertial length.

    The inertial length is also known as the skin depth.

    Examples
    --------
    >>> from astropy import units as u
    >>> inertial_length(5 * u.m ** -3, 'He+')
    <Quantity 2.02985...e+08 m>
    >>> inertial_length(5 * u.m ** -3, 'e-')
    <Quantity 2376534.75... m>

    """
    omega_p = frequencies.plasma_frequency(n, particle=particle)

    return c / omega_p


cwp_ = inertial_length
"""
Alias to `~plasmapy.formulary.lengths.inertial_length`.

* Name is shorthand for :math:`c / ω_p`.
"""
