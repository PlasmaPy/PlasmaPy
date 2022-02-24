"""Functions to calculate fundamental plasma frequency parameters."""
__all__ = ["gyrofrequency"]
__aliases__ = ["oc_", "wc_"]

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
