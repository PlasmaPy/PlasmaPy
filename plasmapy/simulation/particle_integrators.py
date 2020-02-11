import numba
import numpy as np


@numba.njit(parallel=True)
def _boris_push(x, v, b, e, q, m, dt):
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


@numba.njit(parallel=True)
def _boris_push_implicit(x, v, b, e, q, m, dt):
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
    C = q / m
    for i in numba.prange(len(x)):
        # add first half of electric impulse
        B_x, B_y, B_z = b[i]
        v_x, v_y, v_z = v[i]
        E_x, E_y, E_z = e[i]

        # calculated via sympy
        vx = (
            0.0625 * B_x ** 2 * C ** 3 * E_x * dt ** 3
            + 0.0625 * B_x ** 2 * C ** 2 * dt ** 2 * v_x
            + 0.0625 * B_x * B_y * C ** 3 * E_y * dt ** 3
            + 0.125 * B_x * B_y * C ** 2 * dt ** 2 * v_y
            + 0.0625 * B_x * B_z * C ** 3 * E_z * dt ** 3
            + 0.125 * B_x * B_z * C ** 2 * dt ** 2 * v_z
            - 0.0625 * B_y ** 2 * C ** 2 * dt ** 2 * v_x
            - 0.125 * B_y * C ** 2 * E_z * dt ** 2
            - 0.25 * B_y * C * dt * v_z
            - 0.0625 * B_z ** 2 * C ** 2 * dt ** 2 * v_x
            + 0.125 * B_z * C ** 2 * E_y * dt ** 2
            + 0.25 * B_z * C * dt * v_y
            + 0.25 * C * E_x * dt
            + 0.25 * v_x
        ) / (
            0.0625 * B_x ** 2 * C ** 2 * dt ** 2
            + 0.0625 * B_y ** 2 * C ** 2 * dt ** 2
            + 0.0625 * B_z ** 2 * C ** 2 * dt ** 2
            + 0.25
        )
        vy = (
            -0.0625 * B_x ** 2 * C ** 2 * dt ** 2 * v_y
            + 0.0625 * B_x * B_y * C ** 3 * E_x * dt ** 3
            + 0.125 * B_x * B_y * C ** 2 * dt ** 2 * v_x
            + 0.125 * B_x * C ** 2 * E_z * dt ** 2
            + 0.25 * B_x * C * dt * v_z
            + 0.0625 * B_y ** 2 * C ** 3 * E_y * dt ** 3
            + 0.0625 * B_y ** 2 * C ** 2 * dt ** 2 * v_y
            + 0.0625 * B_y * B_z * C ** 3 * E_z * dt ** 3
            + 0.125 * B_y * B_z * C ** 2 * dt ** 2 * v_z
            - 0.0625 * B_z ** 2 * C ** 2 * dt ** 2 * v_y
            - 0.125 * B_z * C ** 2 * E_x * dt ** 2
            - 0.25 * B_z * C * dt * v_x
            + 0.25 * C * E_y * dt
            + 0.25 * v_y
        ) / (
            0.0625 * B_x ** 2 * C ** 2 * dt ** 2
            + 0.0625 * B_y ** 2 * C ** 2 * dt ** 2
            + 0.0625 * B_z ** 2 * C ** 2 * dt ** 2
            + 0.25
        )
        vz = (
            -C
            * dt
            * (0.5 * B_x - 0.25 * B_y * B_z * C * dt)
            * (
                0.5 * B_x * C * dt * v_z
                - 0.5 * B_z * C * dt * v_x
                - 0.5
                * B_z
                * C
                * dt
                * (
                    -0.5 * B_y * C * dt * v_z
                    + 0.5 * B_z * C * dt * v_y
                    + C * E_x * dt
                    + v_x
                )
                + C * E_y * dt
                + v_y
            )
            + (0.25 * B_z ** 2 * C ** 2 * dt ** 2 + 1)
            * (
                -0.5 * B_x * C * dt * v_y
                + 0.5 * B_y * C * dt * v_x
                + 0.5
                * B_y
                * C
                * dt
                * (
                    -0.5 * B_y * C * dt * v_z
                    + 0.5 * B_z * C * dt * v_y
                    + C * E_x * dt
                    + v_x
                )
                + C * E_z * dt
                + v_z
            )
        ) / (
            C ** 2
            * dt ** 2
            * (0.5 * B_x - 0.25 * B_y * B_z * C * dt)
            * (0.5 * B_x + 0.25 * B_y * B_z * C * dt)
            + (0.25 * B_y ** 2 * C ** 2 * dt ** 2 + 1)
            * (0.25 * B_z ** 2 * C ** 2 * dt ** 2 + 1)
        )
        v[i] = (vx, vy, vz)
        x[i] += v[i] * dt


@numba.njit(parallel=True)
def _zenitani(x, v, b, e, q, m, dt):
    r"""
    Implement the Zenitani-Umeda pusher

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
    .. [1] Seiji Zenitani and Takayuki Umeda,
           On the Boris solver in particle-in-cell simulation
           Physics of Plasmas 25, 112110 (2018); https://doi.org/10.1063/1.5051077
    """
    C = q / m
    for i in numba.prange(len(x)):
        # add first half of electric impulse
        epsilon = q * dt / 2 / m * e[i]
        uminus = v[i] + epsilon
        magfield_norm = np.linalg.norm(b[i])
        if magfield_norm == 0.0:
            uplus = uminus
        else:
            theta = q * dt / m * magfield_norm  # Eq. 6
            bnormed = b[i] / magfield_norm
            u_parallel_minus = np.dot(uminus, bnormed) * bnormed  # Eq. 11
            uplus = (
                u_parallel_minus
                + (uminus - u_parallel_minus) * np.cos(theta)
                + np.cross(uminus, bnormed) * np.sin(theta)
            )  # Eq. 12
        v[i] = uplus + epsilon  # Eq. 5

        x[i] += v[i] * dt
