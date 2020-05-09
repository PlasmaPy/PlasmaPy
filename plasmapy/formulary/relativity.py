r"""This module currently provides ample room for the Lorentz factor,
as it turned out we didn't really have much else of the relativistic
variety to add just yet! This is expected to change in the future.
"""
__all__ = ["Lorentz_factor"]

import numpy as np

from astropy import units as u
from astropy.constants import c
from plasmapy import utils
from plasmapy.utils.decorators import validate_quantities


@validate_quantities(V={"can_be_negative": True})
def Lorentz_factor(V: u.m / u.s):
    r"""
    Return the Lorentz factor.

    Parameters
    ----------
    V : ~astropy.units.Quantity
        The velocity in units convertible to meters per second.

    Returns
    -------
    gamma : float or ~numpy.ndarray
        The Lorentz factor associated with the inputted velocities.

    Raises
    ------
    TypeError
        The `V` is not a `~astropy.units.Quantity` and cannot be
        converted into a ~astropy.units.Quantity.

    ~astropy.units.UnitConversionError
        If the `V` is not in appropriate units.

    ValueError
        If the magnitude of `V` is faster than the speed of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.

    Notes
    -----
    The Lorentz factor is a dimensionless number given by

    .. math::
        \gamma = \frac{1}{\sqrt{1-\frac{V^2}{c^2}}}

    The Lorentz factor is approximately one for sub-relativistic
    velocities, and goes to infinity as the velocity approaches the
    speed of light.

    Examples
    --------
    >>> from astropy import units as u
    >>> velocity = 1.4e8 * u.m / u.s
    >>> Lorentz_factor(velocity)
    1.130885603948959
    >>> Lorentz_factor(299792458*u.m/u.s)
    inf
    """

    if not np.all(np.abs(V) <= c):
        raise utils.RelativityError(
            "The Lorentz factor cannot be calculated for "
            "speeds faster than the speed of light. "
        )

    if V.size > 1:

        gamma = np.zeros_like(V.value)

        equals_c = np.abs(V) == c
        is_slow = ~equals_c

        gamma[is_slow] = ((1 - (V[is_slow] / c) ** 2) ** -0.5).value
        gamma[equals_c] = np.inf

    else:
        if np.abs(V) == c:
            gamma = np.inf
        else:
            gamma = ((1 - (V / c) ** 2) ** -0.5).value

    return gamma


def relativistic_energy(V: u.m / u.s, m: u.kg):
    '''
    Return the relativistic energy of a particle
    
    Parameters
    ----------
    V : ~astropy.units.Quantity
        The velocity in units convertible to meters per second.
        
    m : ~astropy.units.Quantity
        The mass in units convertible to kilograms.

    Returns
    -------
    E : float or ~numpy.ndarray
        The Lorentz factor associated with the inputted velocities.
      
    Raises
    ------
    TypeError
        The `V` is not a `~astropy.units.Quantity` and cannot be
        converted into a ~astropy.units.Quantity.

    ~astropy.units.UnitConversionError
        If the `V` is not in appropriate units.

    ValueError
        If the magnitude of `V` is faster than the speed of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.
    
    Examples
    --------
    >>> from astropy import units as u
    >>> velocity = 1.4e8 * u.m / u.s
    >>> mass = 1 * u.kg
    >>> relativistic energy(velocity, mass)
    1.0163893e17 * u.kg * u.m ** 2 / u.s ** 2
    >>> relativistic(299792458*u.m/u.s)
    inf * u.kg * u.m ** 2 / u.s ** 2
    '''
    
    gamma = Lorentz_factor(V)
    E = gamma*m*c**2
    return E

