import astropy.units as u
from astropy.units.utils import np
from plasmapy import utils

@utils.check_quantity()
def ExB_drift(E: u.V/u.m, B: u.T) -> u.m/u.s:
    """
    Calculate the "electric cross magnetic" particle drift.

    Parameters
    ----------
    E : ~astropy.units.Quantity
        Electric field
    B : ~astropy.units.Quantity
        Magnetic field

    Returns
    -------
    v: ~astropy.units.Quantity
        Drift velocity, in m/s

    Notes
    -----
    The E cross B drift is given by

    .. math::

        \vec{v} = \frac{\vec{E} \times \vec{B}}{|B|^2}

    and is independent of particle charge.

    References
    ----------
    - PM Bellan, Fundamentals of Plasma Physics, 3.57

    """

    # np.cross drops units right now, thus this hack
    cross = np.cross(E.si.value, B.si.value) * E.unit * B.unit
    return (cross / (B*B).sum(-1)).to(u.m/u.s)

@utils.check_quantity()
def force_drift(F: u.N, B: u.T, q: u.C):
    """
    Calculate the general force drift for a particle in a magnetic field.

    Parameters
    ----------
    F : ~astropy.units.Quantity
        Force acting on particle
    B : ~astropy.units.Quantity
        Magnetic field
    q : ~astropy.units.Quantity
        Particle charge

    Returns
    -------
    v: ~astropy.units.Quantity
        Drift velocity, in m/s

    Notes
    -----
    The particle drift in a magnetic field and with a general
    force (e.g. gravity) applied to it is given by

    .. math::

        \vec{v} = \frac{\vec{F} \times \vec{B}}{q |B|^2}

    Note the charge dependency.

    References
    ----------
    - PM Bellan, Fundamentals of Plasma Physics, 3.58

    """
    cross = np.cross(F.si.value, B.si.value) * F.unit * B.unit
    return (cross / (q * (B*B).sum(-1))).to(u.m/u.s)
