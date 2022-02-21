"""Functions to calculated fundamental plasma length parameters."""
__all__ = ["Debye_length", "inertial_length"]
__aliases__ = ["cwp_", "lambdaD_"]

import astropy.units as u
import numpy as np

from astropy.constants.si import c, e, eps0, k_B

from plasmapy import particles
from plasmapy.formulary.parameters import plasma_frequency
from plasmapy.particles import Particle
from plasmapy.utils.decorators import validate_quantities

__all__ += __aliases__


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    n_e={"can_be_negative": False},
)
def Debye_length(T_e: u.K, n_e: u.m ** -3) -> u.m:
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
    Debye_number

    Examples
    --------
    >>> from astropy import units as u
    >>> Debye_length(5e6*u.K, 5e15*u.m**-3)
    <Quantity 0.002182... m>

    """
    lambda_D = np.sqrt(eps0 * k_B * T_e / (n_e * e ** 2))
    return lambda_D


lambdaD_ = Debye_length
"""Alias to `~plasmapy.formulary.parameters.Debye_length`."""


@validate_quantities(
    n={"can_be_negative": False},
    validations_on_return={"equivalencies": u.dimensionless_angles()},
)
@particles.particle_input(require="charged")
def inertial_length(n: u.m ** -3, particle: Particle) -> u.m:
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
    omega_p = plasma_frequency(n, particle=particle)

    return c / omega_p


cwp_ = inertial_length
"""
Alias to `~plasmapy.formulary.parameters.inertial_length`.

* Name is shorthand for :math:`c / ω_p`.
"""
