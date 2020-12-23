import math
import numpy as np

from astropy import constants


def boris_push(x, v, b, e, q, m, dt):
    r"""
    Implement the explicit Boris pusher for moving and accelerating particles.

    Parameters
    ----------
    x : np.ndarray
        particle position at full timestep
    v : np.ndarray
        particle velocity at half timestep
    b : np.ndarray
        magnetic field at full timestep
    e : float
        electric field at full timestep
    q : float
        particle charge
    m : float
        particle mass
    dt : float
        timestep

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

    This ends up causing the magnetic field action to be properly
    "centered" in time, and the algorithm conserves energy.

    References
    ----------
    .. [1] C. K. Birdsall, A. B. Langdon, "Plasma Physics via Computer
           Simulation", 2004, p. 58-63
    .. [2] L. Brieda, "Particle Push in Magnetic Field (Boris Method)",
           https://www.particleincell.com/2011/vxb-rotation/
    """
    hqmdt = 0.5 * dt * q / m
    vminus = v + hqmdt * e

    # rotate to add magnetic field
    t = b * hqmdt
    s = 2 * t / (1 + (t * t).sum(axis=1, keepdims=True))
    vprime = vminus + np.cross(vminus, t)
    vplus = vminus + np.cross(vprime, s)

    # add second half of electric impulse
    v[...] = vplus + hqmdt * e

    x += v * dt
