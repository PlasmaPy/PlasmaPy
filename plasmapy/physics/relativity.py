import numpy as np
from astropy import units

from ..constants import c
import plasmapy.atomic as atomic
from plasmapy.utils.checks import _check_quantity
from plasmapy.utils.exceptions import RelativityError


def Lorentz_factor(V):
    r"""Returns the Lorentz factor.

    Parameters
    ----------
    V : Quantity
        The velocity in units convertible to meters per second.

    Returns
    -------
    gamma : float or ndarray
        The Lorentz factor associated with the inputted velocities

    Raises
    ------
    TypeError
        The velocity is not a Quantity and cannot be converted into a
        Quantity.

    UnitConversionError
        If the velocity is not in appropriate units.

    ValueError
        If the magnitude of V is faster than the speed of light.

    UserWarning
        If V is not a Quantity, then a UserWarning will be raised and
        units of meters per second will be assumed.

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
    >>> velocity = 1.4e8*u.m/u.s
    >>> Lorentz_factor(velocity)
    1.130885603948959
    >>> Lorentz_factor(299792458*u.m/u.s)
    inf
    """

    _check_quantity(V, 'V', 'Lorentz_factor', units.m/units.s)

    if not np.all(np.abs(V) <= c):
        raise RelativityError("The Lorentz factor cannot be calculated for "
                              "speeds faster than the speed of light. ")

    if V.size > 1:

        gamma = np.zeros_like(V.value)

        equals_c = np.abs(V) == c
        is_slow = ~equals_c

        gamma[is_slow] = ((1 - (V[is_slow]/c)**2)**-0.5).value
        gamma[equals_c] = np.inf

    else:
        if np.abs(V) == c:
            gamma = np.inf
        else:
            gamma = ((1 - (V/c)**2)**-0.5).value

    return gamma
