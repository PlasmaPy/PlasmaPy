import numpy as np
from astropy import units
from ..constants import mu0, pi
from ..utils import _check_quantity

def infinite_wire(I, r):
    r"""Returns the magnetic field of an infinite,
    current-carrying wire at distance r.

    Parameters
    ----------
    I: Quantity
        Current, given in units convertible to Amperes.

    r: Quantity
        Distance from wire in units convertible to meters.

    Returns
    -------
    B: Quantity
        Magnetic field of the wire at distance r.

    Raises
    ------
    TypeError
        The parameter arguments are not Quantities and
        cannot be converted into Quantities.

    UnitConversionError
        If the parameters are not in appropriate units.

    Notes
    -----
    The magnetic field of an infinite, current-carrying wire is

    .. math::
    B = \frac{\mu_0 I}{2\pi r}

    where :math:`\mu_0` is the magnetic permeability of the vacuum,
    :math:`I` is the current passing through the wire,
    :math:`r` is the distance from the wire.

    Examples
    --------
    >>> from astropy import units as u
    >>> distance = 1.5*u.m
    >>> current = 2.45*u.ampere
    >>> infinite_wire(current, distance)
    <Quantity 3.2666666666666673e-07 T>
    """
    _check_quantity({"I": {"units" : units.ampere, "can_be_complex": False,
                          "can_be_inf": False},
                     "r": {"units" : units.m, "can_be_complex": False,
                          "can_be_inf": True, "can_be_negative": False}})

    if r > 0:
        B = np.divide(mu0 * I,2*pi*r) 
    else:
        B = np.inf * units.tesla

    return B.to(units.tesla)
