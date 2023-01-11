"""Particle movement integrators, for particle simulations.

These do not have `astropy.units` support, choosing instead to
limit overhead and increase performance.

They act in-place on position and velocity arrays to reduce
memory allocation.
"""

__all__ = ["boris_push"]

import numpy as np


def boris_push(x, v, B, E, q, m, dt, inplace: bool = True):
    r"""
    The explicit Boris pusher.

    Parameters
    ----------
    x : `~numpy.ndarray`
        Particle position at full timestep, in SI (meter) units.

    v : `~numpy.ndarray`
        Particle velocity at half timestep, in SI (meter/second) units.

    B : `~numpy.ndarray`
        Magnetic field at full timestep, in SI (tesla) units.

    E : `float`
        Electric field at full timestep, in SI (V/m) units.

    q : `float`
        Particle charge, in SI (coulomb) units.

    m : `float`
        Particle mass, in SI (kg) units.

    dt : `float`
        Timestep, in SI (second) units.

    inplace : `bool`
        If set, data in the v and x arrays is changed during execution.
        Defaults to True for performance reasons.
        If False, this function returns `None`.

    Returns
    -------
    x : `~numpy.ndarray`
        Particle position x after Boris push algorithm, in SI (meter) units.

    v : `~numpy.ndarray`
        Particle velocity after Boris push algorithm, in SI (meter/second) units.

    Examples
    --------
    >>> B = np.array([[0.0, 0.0, 5.0]])
    >>> E = np.array([[0.0, 0.0, 0.0]])
    >>> x_t0 = np.array([[0.0, 0.0, 0.0]])
    >>> v_t0 = np.array([[5.0, 0.0, 0.0]])
    >>> x_t1, v_t1 = boris_push(x = x_t0, v = v_t0, B = B, E = E, q = 1.0, m = 1.0, dt = 0.01, inplace=False)
    >>> x_t1
    array([[ 0.04993754, -0.00249844,  0.        ]])
    >>> v_t1
    array([[ 4.9937539 , -0.24984385,  0.        ]])
    >>> boris_push(x = x_t0, v = v_t0, B = B, E = E, q = 1.0, m = 1.0, dt = 0.01, inplace = True)
    >>> x_t0
    array([[ 0.04993754, -0.00249844,  0.        ]])
    >>> v_t0
    array([[ 4.9937539 , -0.24984385,  0.        ]])
    >>> # B parallel to v
    >>> B = np.array([[0.0, 0.0, 5.0]])
    >>> v_t0 = np.array([[0.0, 0.0, 5.0]])
    >>> x_t0 = np.array([[0.0, 0.0, 0.0]])
    >>> boris_push(x = x_t0, v = v_t0, B = B, E = E, q = 1.0, m = 1.0, dt = 0.01, inplace = True)
    >>> # no rotation of vector v
    >>> v_t0
    array([[0., 0., 5.]])
    >>> x_t0
    array([[0.  , 0.  , 0.05]])
    >>> # B perpendicular to v
    >>> B = np.array([[5.0, 0.0, 0.0]])
    >>> v_t0 = np.array([[0.0, 5.0, 0.0]])
    >>> x_t0 = np.array([[0.0, 0.0, 0.0]])
    >>> boris_push(x = x_t0, v = v_t0, B = B, E = E, q = 1.0, m = 1.0, dt = 0.01, inplace = True)
    >>> # rotation of vector v
    >>> v_t0
    array([[ 0.        ,  4.9937539 , -0.24984385]])
    >>> x_t0
    array([[ 0.        ,  0.04993754, -0.00249844]])
    >>> # nonzero E and zero B
    >>> E = np.array([[1.0, 1.0, 1.0]])
    >>> B = np.array([[0.0, 0.0, 0.0]])
    >>> v_t0 = np.array([[0.0, 0.0, 5.0]])
    >>> x_t0 = np.array([[0.0, 0.0, 0.0]])
    >>> boris_push(x = x_t0, v = v_t0, B = B, E = E, q = 1.0, m = 1.0, dt = 0.01, inplace = True)
    >>> v_t0
    array([[0.01, 0.01, 5.01]])
    >>> x_t0
    array([[0.0001, 0.0001, 0.0501]])
    >>> boris_push(x = x_t0, v = v_t0, B = B, E = E, q = 1.0, m = 1.0, dt = 0.01, inplace = True)
    >>> v_t0
    array([[0.02, 0.02, 5.02]])
    >>> x_t0
    array([[0.0003, 0.0003, 0.1003]])

    Notes
    ----------
    The Boris algorithm :cite:p:`boris:1970` is the standard energy
    conserving algorithm for particle movement in plasma physics. See
    :cite:t:`birdsall:2004` for more details, and this `page on the
    Boris method <https://www.particleincell.com/2011/vxb-rotation>`__
    for a nice overview.

    Conceptually, the algorithm has three phases:

    1. Add half the impulse from electric field.
    2. Rotate the particle velocity about the direction of the magnetic
       field.
    3. Add the second half of the impulse from the electric field.

    This ends up causing the magnetic field action to be properly "centered" in
    time, and the algorithm, being a symplectic integrator, conserves energy.
    """
    hqmdt = 0.5 * dt * q / m
    vminus = v + hqmdt * E

    # rotate to add magnetic field
    t = B * hqmdt
    s = 2 * t / (1 + (t * t).sum(axis=1, keepdims=True))
    vprime = vminus + np.cross(vminus, t)
    vplus = vminus + np.cross(vprime, s)

    # add second half of electric impulse
    if inplace:
        v[...] = vplus + hqmdt * E
        x += v * dt
    else:
        v = vplus + hqmdt * E
        return x + v * dt, v
