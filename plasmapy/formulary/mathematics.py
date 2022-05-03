"""Mathematical formulas relevant to plasma physics."""

__all__ = ["Fermi_integral", "rot_a_to_b"]

import numbers
import numpy as np

from mpmath import polylog
from typing import Union


def Fermi_integral(
    x: Union[float, int, complex, np.ndarray], j: Union[float, int, complex, np.ndarray]
) -> Union[float, complex, np.ndarray]:
    r"""
    Calculate the complete Fermi-Dirac integral.

    Parameters
    ----------
    x : `float`, `int`, `complex`, or `~numpy.ndarray`
        Argument of the Fermi-Dirac integral function.

    j : `float`, `int`, `complex`, or `~numpy.ndarray`
        Order/index of the Fermi-Dirac integral function.

    Returns
    -------
    integral : `float`, `complex`, or `~numpy.ndarray`
        Complete Fermi-Dirac integral for given argument and order.

    Raises
    ------
    `TypeError`
        If the argument is invalid.

    `~astropy.units.UnitsError`
        If the argument is a `~astropy.units.Quantity` but is not
        dimensionless.

    `ValueError`
        If the argument is not entirely finite.

    Notes
    -----
    The `complete Fermi-Dirac integral
    <https://en.wikipedia.org/wiki/Complete_Fermi-Dirac_integral>`_ is
    defined as:

    .. math::
        F_j (x) = \frac{1}{Γ(j+1)} \int_0^∞ \frac{t^j}{\exp{(t-x)} + 1} dt

    for :math:`j > 0`.

    This is equivalent to the following `polylogarithm
    <https://en.wikipedia.org/wiki/Polylogarithm>`_ function:

    .. math::
        F_j (x) = -Li_{j+1}\left(-e^{x}\right)

    Warnings
    --------
    At present this function is limited to relatively small arguments
    due to limitations in ``mpmath.polylog``.

    Examples
    --------
    >>> Fermi_integral(0, 0)
    (0.6931471805599453-0j)
    >>> Fermi_integral(1, 0)
    (1.3132616875182228-0j)
    >>> Fermi_integral(1, 1)
    (1.8062860704447743-0j)
    """
    if isinstance(x, (numbers.Integral, numbers.Real, numbers.Complex)):
        arg = -np.exp(x)
        return -1 * complex(polylog(j + 1, arg))
    elif isinstance(x, np.ndarray):
        integral_arr = np.zeros_like(x, dtype="complex")
        for idx, val in enumerate(x):
            integral_arr[idx] = -1 * complex(polylog(j + 1, -np.exp(val)))
        return integral_arr
    else:
        raise TypeError(f"Improper type {type(x)} given for argument x.")


def rot_a_to_b(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    r"""
    Calculates the 3D rotation matrix that will rotate vector ``a`` to be aligned
    with vector ``b``. The rotation matrix is calculated as follows. Let

    .. math::
        \vec v = \vec a \times \vec b

    and let :math:`\theta` be the angle between :math:`\vec a`
    and :math:`\vec b` such that the projection of :math:`\vec a` along
    :math:`\vec b` is

    .. math::
        c = \vec a \cdot \vec b \cos\theta

    Then the rotation matrix :math:`R` is

    .. math::
        R = I + v_x + v_x^2 \frac{1}{1 + c}


    where :math:`I` is the identity matrix and :math:`v_x` is the
    skew-symmetric cross-product matrix of :math:`v` defined as

    .. math::
        v_x = \begin{bmatrix}
                0 & -v_3 & v_2 \\
                v_3 & 0 & -v_1 \\
                -v_2 & v_1 & 0
            \end{bmatrix}

    Note that this algorithm fails when :math:`1+c=0`, which occurs when :math:`a` and
    :math:`b` are anti-parallel. However, since the correct rotation matrix
    in this case is simply :math:`R=-I`, this function just handles this
    special case explicitly.

    This algorithm is based on
    `this discussion
    <https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311>`_
    on StackExchange.

    Parameters
    ----------
    a : `~numpy.ndarray`, shape (3,)
        Vector to be rotated.  Should be a 1D, 3-element unit vector.  If ``a``
        is not normalize, then it will be normalized.

    b : `~numpy.ndarray`, shape (3,)
        Vector representing the desired orientation after rotation.  Should be
        a 1D, 3-element unit vector.  If ``b`` is not normalized, then it will
        be.

    Returns
    -------
    R : `~numpy.ndarray`, shape (3,3)
        The rotation matrix that will rotate vector ``a`` onto vector ``b``.
    """

    # Normalize and validate both vectors

    a = np.squeeze(a)
    if a.shape != (3,):
        raise ValueError(
            f"Argument 'a' must have shape (3,) but input has shape {a.shape}."
        )
    a = a / np.linalg.norm(a)

    b = np.squeeze(b)
    if b.shape != (3,):
        raise ValueError(
            f"Argument 'b' must have shape (3,) but input has shape {b.shape}."
        )
    b = b / np.linalg.norm(b)

    # Manually handle the case where a and b point in opposite directions
    if np.dot(a, b) == -1:
        return -np.identity(3)

    axb = np.cross(a, b)
    c = np.dot(a, b)
    vskew = np.array(
        [[0, -axb[2], axb[1]], [axb[2], 0, -axb[0]], [-axb[1], axb[0], 0]]
    ).T  # Transpose to get right orientation

    return np.identity(3) + vskew + np.dot(vskew, vskew) / (1 + c)
