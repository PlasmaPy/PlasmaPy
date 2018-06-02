import astropy.units as u
from astropy.units.utils import np
from plasmapy import utils

@utils.check_quantity(E={'units':u.V/u.m}, B={'units':u.T})
def ExB_drift(E: u.V/u.m, B: u.T) -> u.m/u.s:
    """
    TODO

    Parameters
    ----------
    E : ~astropy.units.Quantity
        Electric field
    B : ~astropy.units.Quantity
        Magnetic field

    Returns
    -------
    v: ~astropy.units.Quantity
        Velocity, in m/s

    References
    ----------
    Bellan, Fundamentals of Plasma Physics, 3.57

    """

    # np.cross drops units right now, thus this hack
    cross = np.cross(E.si.value, B.si.value) * E.unit * B.unit
    return (cross / (B*B).sum(-1)).to(u.m/u.s)

@utils.check_quantity(F={'units':u.N}, B={'units':u.T}, q={'units':u.C})
def force_drift(F: u.N, B: u.T, q: u.C):
    """

    Parameters
    ----------
    F :
    B :
    q :

    Returns
    -------
    v: ~astropy.units.Quantity
        Velocity, in m/s

    References
    ----------
    Bellan, Fundamentals of Plasma Physics, 3.58

    """
    cross = np.cross(F.si.value, B.si.value) * F.unit * B.unit
    return (cross / (q * (B*B).sum(-1))).to(u.m/u.s)
