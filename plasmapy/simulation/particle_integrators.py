import numba
import numpy as np


@numba.njit(parallel=True)
def _boris_push(x, v, b, e, hqmdt, dt):
    r"""
    Implement the explicit Boris pusher for moving and accelerating particles.

    Arguments
    ----------
    init : bool (optional)
        If `True`, does not change the particle positions and sets dt
        to -dt/2.

    Notes
    ----------
    The Boris algorithm is the standard energy conserving algorithm for
    particle movement in plasma physics. See [1]_ for more details.

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
    """
    for i in numba.prange(len(x)):
        # add first half of electric impulse
        vminus = v[i] + hqmdt * e[i]

        # rotate to add magnetic field
        t = -b[i] * hqmdt
        s = 2 * t / (1 + (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]))
        cross_result = np.cross(vminus, t)
        vprime = vminus + cross_result
        cross_result_2 = np.cross(vprime, s)
        vplus = vminus + cross_result_2

        # add second half of electric impulse
        v[i] = vplus + e[i] * hqmdt
        x[i] += v[i] * dt


@numba.njit(parallel=True)
def _boris_push_implicit(x, v, b, e, hqmdt, dt):
    r"""
    Implement the implicit Boris pusher for moving and accelerating particles.

    Arguments
    ----------
    init : bool (optional)
        If `True`, does not change the particle positions and sets dt
        to -dt/2.

    Notes
    ----------
    The Boris algorithm is the standard energy conserving algorithm for
    particle movement in plasma physics. See [1]_ for more details.

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
    """
    raise NotImplementedError
