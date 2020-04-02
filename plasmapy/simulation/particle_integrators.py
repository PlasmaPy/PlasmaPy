import math

import numba
import numpy as np
from astropy import constants


@numba.njit(parallel=True)
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
    hqmdt = 0.5 * dt * q / m
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


c = constants.c.si.value


@numba.njit()
def gamma_from_velocity(velocity):
    return np.sqrt(1 - ((np.linalg.norm(velocity) / c) ** 2))


@numba.njit()
def gamma_from_u(u):
    return np.sqrt(1 + ((np.linalg.norm(u) / c) ** 2))


@numba.njit(parallel=True)
def _zenitani(x, v, b, e, q, m, dt, B_numerical_threshold=1e-20):
    r"""
    Implement the Zenitani-Umeda pusher.

    This implementation is currently not complete - use at your own risk!

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
    TODO

    References
    ----------
    .. [1] Seiji Zenitani and Takayuki Umeda,
           On the Boris solver in particle-in-cell simulation
           Physics of Plasmas 25, 112110 (2018); https://doi.org/10.1063/1.5051077
    """
    C = q / m * dt
    for i in numba.prange(len(x)):
        # add first half of electric impulse
        epsilon = C / 2.0 * e[i]
        uminus = v[i] + epsilon
        magfield_norm = max((np.linalg.norm(b[i]), B_numerical_threshold))
        theta = C * magfield_norm / gamma_from_u(uminus)  # Eq. 6
        bnormed = b[i] / magfield_norm
        u_parallel_minus = np.dot(uminus, bnormed) * bnormed  # Eq. 11
        uplus = (
            u_parallel_minus
            + (uminus - u_parallel_minus) * math.cos(theta)
            + np.cross(uminus, bnormed) * math.sin(theta)
        )  # Eq. 12
        u_t_plus_half = uplus + epsilon
        v[i] = u_t_plus_half / gamma_from_u(u_t_plus_half)
        x[i] += v[i] * dt
