"""Functions to calculate fundamental plasma length parameters."""
__all__ = ["Debye_length", "gyroradius", "inertial_length"]
__aliases__ = ["cwp_", "lambdaD_", "rc_", "rhoc_"]

import astropy.units as u
import numpy as np
import warnings

from astropy.constants.si import c, e, eps0, k_B

from plasmapy.formulary import frequencies, speeds
from plasmapy.particles import particle_input, ParticleLike
from plasmapy.utils.decorators import validate_quantities
from plasmapy.utils.exceptions import PlasmaPyFutureWarning

__all__ += __aliases__


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def Debye_length(T_e: u.K, n_e: u.m**-3) -> u.m:
    r"""Calculate the characteristic decay length for electric fields,
     due to charge screening.

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
        λ_D = \sqrt{\frac{ε_0 k_b T_e}{n_e e^2}}

    for an electron plasma with nearly stationary ions.

    The electrical potential will drop by a factor of :math:`1/e` every Debye
    length.

    Plasmas will generally be quasineutral on length scales significantly
    larger than the Debye length.

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
    T_i={
        "can_be_nan": True,
        "equivalencies": u.temperature_energy(),
        "none_shall_pass": True,
    },
    T={
        "can_be_nan": True,
        "equivalencies": u.temperature_energy(),
        "none_shall_pass": True,
    },
    validations_on_return={"equivalencies": u.dimensionless_angles()},
)
@particle_input(any_of={"charged", "uncharged"})
def gyroradius(
    B: u.T,
    particle: ParticleLike,
    *,
    Vperp: u.m / u.s = np.nan * u.m / u.s,
    T_i: u.K = None,
    T: u.K = None,
) -> u.m:
    r"""Return the particle gyroradius.

    **Aliases:** `rc_`, `rhoc_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to tesla.

    particle : `~plasmapy.particles.particle_class.Particle`
        Representation of the particle species (e.g., ``'p'`` for protons, ``'D+'``
        for deuterium, or ``'He-4 +1'`` for singly ionized helium-4).  If no
        charge state information is provided, then the particles are assumed
        to be singly charged.

    Vperp : `~astropy.units.Quantity`, optional, |keyword-only|
        The component of particle velocity that is perpendicular to the
        magnetic field in units convertible to meters per second.

    T : `~astropy.units.Quantity`, optional, |keyword-only|
        The particle temperature in units convertible to kelvin.

    T_i : `~astropy.units.Quantity`, optional, |keyword-only|
        The particle temperature in units convertible to kelvin.
        Note: Deprecated. Use ``T`` instead.

    Returns
    -------
    r_Li : `~astropy.units.Quantity`
        The particle gyroradius in units of meters.  This
        `~astropy.units.Quantity` will be based on either the
        perpendicular component of particle velocity as inputted, or
        the most probable speed for a particle within a Maxwellian
        distribution for the particle temperature.

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

    If any of ``B``, ``Vperp``, or ``T`` is a number rather than a
    `~astropy.units.Quantity`, then SI units will be assumed and a
    warning will be raised.

    The particle gyroradius is also known as the particle Larmor
    radius and is given by

    .. math::
        r_{Li} = \frac{V_{\perp}}{ω_{ci}}

    where :math:`V_⟂` is the component of particle velocity that is
    perpendicular to the magnetic field and :math:`ω_{ci}` is the
    particle gyrofrequency.  If a temperature is provided, then
    :math:`V_⟂` will be the most probable thermal velocity of a
    particle at that temperature.

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
    <Quantity 48.23129... m>
    >>> gyroradius(B=0.01*u.T, particle='e-', T=1e6*u.K)
    <Quantity 0.003130... m>
    >>> gyroradius(0.01*u.T, 'e-', Vperp=1e6*u.m/u.s)
    <Quantity 0.000568... m>
    >>> gyroradius(0.2*u.T, 'e-', T=1e5*u.K)
    <Quantity 4.94949...e-05 m>
    >>> gyroradius(5*u.uG, 'e-', T=1*u.eV)
    <Quantity 6744.25... m>
    >>> gyroradius(400*u.G, 'e-', Vperp=1e7*u.m/u.s)
    <Quantity 0.001421... m>
    """

    # Backwards Compatibility and Deprecation check for keyword T_i
    if T_i is not None:
        warnings.warn(
            "Keyword T_i is deprecated, use T instead.",
            PlasmaPyFutureWarning,
        )
        if T is None:
            T = T_i
        else:
            raise ValueError(
                "Keywords T_i and T are both given.  T_i is deprecated, "
                "please use T only."
            )

    if T is None:
        T = np.nan * u.K

    isfinite_T = np.isfinite(T)
    isfinite_Vperp = np.isfinite(Vperp)

    # check 1: ensure either Vperp or T invalid, keeping in mind that
    # the underlying values of the astropy quantity may be numpy arrays
    if np.any(np.logical_and(isfinite_Vperp, isfinite_T)):
        raise ValueError(
            "Must give Vperp or T, but not both, as arguments to gyroradius"
        )

    # check 2: get Vperp as the thermal speed if is not already a valid input
    if np.isscalar(Vperp.value) and np.isscalar(
        T.value
    ):  # both T and Vperp are scalars
        # we know exactly one of them is nan from check 1
        if isfinite_T:
            # T is valid, so use it to determine Vperp
            Vperp = speeds.thermal_speed(T, particle=particle)
        # else: Vperp is already valid, do nothing
    elif np.isscalar(Vperp.value):  # only T is an array
        # this means either Vperp must be nan, or T must be an array of all nan,
        # or else we couldn't have gotten through check 1
        if isfinite_Vperp:
            # Vperp is valid, T is a vector that is all nan
            # uh...
            Vperp = np.repeat(Vperp, len(T))
        else:
            # normal case where Vperp is scalar nan and T is valid array
            Vperp = speeds.thermal_speed(T, particle=particle)
    elif np.isscalar(T.value):  # only Vperp is an array
        # this means either T must be nan, or V_perp must be an array of all nan,
        # or else we couldn't have gotten through check 1
        if isfinite_T:
            # T is valid, V_perp is an array of all nan
            # uh...
            Vperp = speeds.thermal_speed(np.repeat(T, len(Vperp)), particle=particle)
        # else: normal case where T is scalar nan and Vperp is already a valid
        # array so, do nothing
    else:  # both T and Vperp are arrays
        # we know all the elementwise combinations have one nan and one finite,
        # due to check 1 use the valid Vperps, and replace the others with those
        # calculated from T
        Vperp = Vperp.copy()  # avoid changing Vperp's value outside function
        Vperp[isfinite_T] = speeds.thermal_speed(T[isfinite_T], particle=particle)

    omega_ci = frequencies.gyrofrequency(B, particle)

    return np.abs(Vperp) / omega_ci


rc_ = gyroradius
"""Alias to `~plasmapy.formulary.lengths.gyroradius`."""

rhoc_ = gyroradius
"""Alias to `~plasmapy.formulary.lengths.gyroradius`."""


@validate_quantities(
    n={"can_be_negative": False},
    validations_on_return={"equivalencies": u.dimensionless_angles()},
)
@particle_input(require="charged")
def inertial_length(n: u.m**-3, particle: ParticleLike) -> u.m:
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
    particle to be accelerated in a plasma.  The Hall effect becomes
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
