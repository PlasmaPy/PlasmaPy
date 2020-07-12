"""
Formulas for calculating particle drifts.
"""
__all__ = [
    "ExB_drift",
    "force_drift",
    "veb_",
    "vfd_",
]

import astropy.units as u
import numpy as np

from plasmapy.utils.decorators import validate_quantities


@validate_quantities
def ExB_drift(E: u.V / u.m, B: u.T) -> u.m / u.s:
    r"""
    Calculate the "electric cross magnetic" particle drift.

    **Aliases:** `veb_`

    Parameters
    ----------
    E : ~astropy.units.Quantity
        Electric field vector
    B : ~astropy.units.Quantity
        Magnetic field vector

    Returns
    -------
    v: ~astropy.units.Quantity
        Drift velocity, in m/s

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.constants.si import g0, e, m_e
    >>> ex = np.array([1, 0, 0])
    >>> ez = np.array([0, 0, 1])
    >>> force_drift(-ez*g0*m_e, ex*0.01*u.T, e)
    <Quantity [ 0.0000000e+00, -5.5756984e-09,  0.0000000e+00] m / s>
    >>> force_drift(-ez*g0*m_e, ez*0.01*u.T, e)
    <Quantity [ 0., -0.,  0.] m / s>
    >>> force_drift(-ez*g0*m_e, ex*u.T, e)
    <Quantity [ 0.0000000e+00, -5.5756984e-11,  0.0000000e+00] m / s>

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

    # np.cross drops units right now, thus this hack: see
    # https://github.com/PlasmaPy/PlasmaPy/issues/59
    cross = np.cross(E.si.value, B.si.value) * E.unit * B.unit
    return cross / (B * B).sum(-1)


veb_ = ExB_drift
""" Alias to :func:`ExB_drift`. """


@validate_quantities
def force_drift(F: u.N, B: u.T, q: u.C) -> u.m / u.s:
    r"""
    Calculate the general force drift for a particle in a magnetic field.

    **Aliases:** `vfd_`

    Parameters
    ----------
    F : ~astropy.units.Quantity
        Force acting on particle
    B : ~astropy.units.Quantity
        Magnetic field
    q : ~astropy.units.Quantity
        Particle charge

    Examples
    --------
    >>> import astropy.units as u
    >>> ex = np.array([1, 0, 0])
    >>> ey = np.array([0, 1, 0])
    >>> ExB_drift(ex * u.V/u.m, ey * u.T)
    <Quantity [0., 0., 1.] m / s>
    >>> ExB_drift(ex * u.V/u.m, ex * u.T)
    <Quantity [0., 0., 0.] m / s>
    >>> ExB_drift(ex * u.V/u.m, 100 * ey * u.T)
    <Quantity [0.  , 0.  , 0.01] m / s>

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
    return cross / (q * (B * B).sum(-1))


vfd_ = force_drift
""" Alias to :func:`force_drift`. """
