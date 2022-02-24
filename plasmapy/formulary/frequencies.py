"""Functions to calculate fundamental plasma frequency parameters."""
__all__ = ["gyrofrequency", "lower_hybrid_frequency"]
__aliases__ = ["oc_", "wc_", "wlh_"]

import astropy.units as u
import numpy as np

from astropy.constants.si import e

from plasmapy import particles
from plasmapy.particles import Particle
from plasmapy.utils.decorators import (
    angular_freq_to_hz,
    validate_quantities,
)

__all__ += __aliases__


@validate_quantities(
    validations_on_return={
        "units": [u.rad / u.s, u.Hz],
        "equivalencies": [(u.cy / u.s, u.Hz)],
    }
)
@angular_freq_to_hz
def gyrofrequency(B: u.T, particle: Particle, signed=False, Z=None) -> u.rad / u.s:
    r"""
    Calculate the particle gyrofrequency in units of radians per second.

    **Aliases:** `oc_`, `wc_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to tesla.

    particle : `~plasmapy.particles.particle_class.Particle`
        Representation of the particle species (e.g., 'p' for protons, 'D+'
        for deuterium, or 'He-4 +1' for singly ionized helium-4). If no
        charge state information is provided, then the particles are assumed
        to be singly charged.

    signed : `bool`, optional
        The gyrofrequency can be defined as signed (negative for electron,
        positive for ion). Default is `False` (unsigned, i.e. always
        positive).

    Z : `float` or `~astropy.units.Quantity`, optional
        The average ionization (arithmetic mean) for a plasma where the
        a macroscopic description is valid. If this quantity is not
        given then the charge number of the ion
        is used. This is effectively an average gyrofrequency for the
        plasma where multiple charge states are present, and should
        not be interpreted as the gyrofrequency for any single particle.
        If not provided, it defaults to the charge number of the ``particle``.

    Returns
    -------
    omega_c : `~astropy.units.Quantity`
        The particle gyrofrequency in units of radians per second.

    Raises
    ------
    `TypeError`
        If the magnetic field is not a `~astropy.units.Quantity` or
        ``particle`` is not of an appropriate type.

    `ValueError`
        If the magnetic field contains invalid values or particle cannot
        be used to identify an particle or isotope.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The particle gyrofrequency is the angular frequency of particle
    gyration around magnetic field lines and is given by:

    .. math::
        ω_{c} = \frac{Z e B}{m}

    The particle gyrofrequency is also known as the particle cyclotron
    frequency or the particle Larmor frequency.

    The recommended way to convert from angular frequency to frequency
    is to use an equivalency between cycles per second and hertz, as
    Astropy's `~astropy.units.dimensionles_angles` equivalency does not
    account for the factor of :math:`2π` needed during this conversion.  The
    `~astropy.units.dimensionless_angles` equivalency is appropriate
    when dividing a velocity by an angular frequency to get a length scale.

    Examples
    --------
    >>> from astropy import units as u
    >>> gyrofrequency(0.1*u.T, 'e-')
    <Quantity 1.7588...e+10 rad / s>
    >>> gyrofrequency(0.1*u.T, 'e-', to_hz=True)
    <Quantity 2.79924...e+09 Hz>
    >>> gyrofrequency(0.1*u.T, 'e-', signed=True)
    <Quantity -1.75882...e+10 rad / s>
    >>> gyrofrequency(0.01*u.T, 'p')
    <Quantity 957883.32... rad / s>
    >>> gyrofrequency(0.01*u.T, 'p', signed=True)
    <Quantity 957883.32... rad / s>
    >>> gyrofrequency(0.01*u.T, particle='T+')
    <Quantity 319964.5... rad / s>
    >>> gyrofrequency(0.01*u.T, particle='T+', to_hz=True)
    <Quantity 50923.9... Hz>
    >>> omega_ce = gyrofrequency(0.1*u.T, 'e-')
    >>> print(omega_ce)
    1758820... rad / s
    >>> f_ce = omega_ce.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])
    >>> print(f_ce)
    279924... Hz

    """
    from plasmapy.formulary.parameters import _grab_charge

    m = particles.particle_mass(particle)
    Z = _grab_charge(particle, Z)
    if not signed:
        Z = abs(Z)

    omega_c = u.rad * (Z * e * np.abs(B) / m).to(1 / u.s)

    return omega_c


oc_ = gyrofrequency
"""Alias to `~plasmapy.formulary.parameters.gyrofrequency`."""

wc_ = gyrofrequency
"""Alias to `~plasmapy.formulary.parameters.gyrofrequency`."""


@validate_quantities(
    n_i={"can_be_negative": False},
    validations_on_return={
        "units": [u.rad / u.s, u.Hz],
        "equivalencies": [(u.cy / u.s, u.Hz)],
    },
)
@angular_freq_to_hz
def lower_hybrid_frequency(B: u.T, n_i: u.m ** -3, ion: Particle) -> u.rad / u.s:
    r"""
    Return the lower hybrid frequency.

    **Aliases:** `wlh_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to tesla.

    n_i : `~astropy.units.Quantity`
        Ion number density.

    ion : `~plasmapy.particles.particle_class.Particle`
        Representation of the ion species (e.g., ``'p'`` for protons, ``'D+'``
        for deuterium, or ``'He-4 +1'`` for singly ionized helium-4). If no
        charge state information is provided, then the ions are assumed to
        be singly charged.

    Returns
    -------
    omega_lh : `~astropy.units.Quantity`
        The lower hybrid frequency in radians per second.

    Raises
    ------
    `TypeError`
        If either of ``B`` or ``n_i`` is not a `~astropy.units.Quantity`,
        or ion is of an inappropriate type.

    `~astropy.units.UnitConversionError`
        If either of ``B`` or ``n_i`` is in incorrect units.

    `ValueError`
        If either of ``B`` or ``n_i`` contains invalid values or are of
        incompatible dimensions, or ion cannot be used to identify an
        ion or isotope.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The lower hybrid frequency is given through the relation

    .. math::
        \frac{1}{ω_{lh}^2} = \frac{1}{ω_{ci}^2 + ω_{pi}^2} +
        \frac{1}{ω_{ci}ω_{ce}}

    where :math:`ω_{ci}` is the ion gyrofrequency,
    :math:`ω_{ce}` is the electron gyrofrequency, and
    :math:`ω_{pi}` is the ion plasma frequency.

    The lower hybrid frequency consitutes a resonance for electromagnetic
    waves in magnetized plasmas, namely for the X-mode. These are waves
    with their wave electric field being perpendicular to the background
    magnetic field. For the lower hybrid frequency, ion and electron
    dynamics both play a role. As the name suggests, it has a lower frequency
    compared to the upper hybrid frequency. It can play an important role
    for heating and current drive in fusion plasmas.

    Examples
    --------
    >>> from astropy import units as u
    >>> lower_hybrid_frequency(0.2*u.T, n_i=5e19*u.m**-3, ion='D+')
    <Quantity 5.78372...e+08 rad / s>
    >>> lower_hybrid_frequency(0.2*u.T, n_i=5e19*u.m**-3, ion='D+', to_hz = True)
    <Quantity 92050879.3... Hz>

    """
    from plasmapy.formulary.parameters import plasma_frequency

    # We do not need a charge state here, so the sole intent is to
    # catch invalid ions.
    try:
        particles.charge_number(ion)
    except Exception:
        raise ValueError("Invalid ion in lower_hybrid_frequency.")

    omega_ci = gyrofrequency(B, particle=ion)
    omega_pi = plasma_frequency(n_i, particle=ion)
    omega_ce = gyrofrequency(B, particle="e-")
    omega_lh = ((omega_ci * omega_ce) ** -1 + omega_pi ** -2) ** -0.5
    # TODO possibly optimize the above line via np.sqrt
    omega_lh = omega_lh

    return omega_lh


wlh_ = lower_hybrid_frequency
"""Alias to `~plasmapy.formulary.parameters.lower_hybrid_frequency`."""
