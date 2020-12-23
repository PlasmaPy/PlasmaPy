"""Particle movement integrators, for particle simulations.

These do not have `astropy.units` support, choosing instead to
limit overhead and increase performance.

They act in-place on position and velocity arrays to reduce
memory allocation.
"""
import math
import numpy as np

from astropy import constants


def boris_push(x, v, b, e, q, m, dt):
    r"""
    The explicit Boris pusher.

    Parameters
    ----------
    x : np.ndarray
        particle position at full timestep, in SI (meter) units.
    v : np.ndarray
        particle velocity at half timestep, in SI (meter/second) units.
    B : np.ndarray
        magnetic field at full timestep, in SI (tesla) units.
    E : float
        electric field at full timestep, in SI (V/m) units.
    q : float
        particle charge, in SI (Coulomb) units.
    m : float
        particle mass, in SI (kg) units.
    dt : float
        timestep, in SI (second) units.

    Notes
    ----------
    The Boris algorithm is the standard energy conserving algorithm for
    particle movement in plasma physics. See [1]_ for more details, and
    [2]_ for a nice overview.

    Conceptually, the algorithm has three phases:

    1. Add half the impulse from electric field.
    2. Rotate the particle velocity about the direction of the magnetic
       field.
    3. Add the second half of the impulse from the electric field.

    This ends up causing the magnetic field action to be properly "centered" in
    time, and the algorithm, being a symplectic integrator, conserves energy.

    References
    ----------
    .. [1] C. K. Birdsall, A. B. Langdon, "Plasma Physics via Computer
           Simulation", 2004, p. 58-63
    .. [2] L. Brieda, "Particle Push in Magnetic Field (Boris Method)",
           https://www.particleincell.com/2011/vxb-rotation/
    """
    hqmdt = 0.5 * dt * q / m
    vminus = v + hqmdt * E

    # rotate to add magnetic field
    t = B * hqmdt
    s = 2 * t / (1 + (t * t).sum(axis=1, keepdims=True))
    vprime = vminus + np.cross(vminus, t)
    vplus = vminus + np.cross(vprime, s)

    # add second half of electric impulse
    v[...] = vplus + hqmdt * E

    x += v * dt
