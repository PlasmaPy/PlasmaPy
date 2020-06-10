"""
This module gathers highly theoretical mathematical formulas
relevant to plasma physics. Usually, those are used somewhere else in
the code but were deemed general enough for the mathematical apparatus
to be abstracted from the main function interface.
"""
__all__ = ["Fermi_integral"]

import numbers
from typing import Union

import numpy as np


def Fermi_integral(
    x: Union[float, int, complex, np.ndarray], j: Union[float, int, complex, np.ndarray]
) -> Union[float, complex, np.ndarray]:
    r"""
    Calculate the complete Fermi-Dirac integral.

    Parameters
    ----------
    x : float, int, complex, or ~numpy.ndarray
        Argument of the Fermi-Dirac integral function.

    j : float, int, complex, or ~numpy.ndarray
        Order/index of the Fermi-Dirac integral function.

    Returns
    -------
    integral : float, complex, or ~numpy.ndarray
        Complete Fermi-Dirac integral for given argument and order.

    Raises
    ------
    TypeError
        If the argument is invalid.

    ~astropy.units.UnitsError
        If the argument is a `~astropy.units.Quantity` but is not
        dimensionless.

    ValueError
        If the argument is not entirely finite.

    Notes
    -----
    The `complete Fermi-Dirac integral
    <https://en.wikipedia.org/wiki/Complete_Fermi-Dirac_integral>`_ is
    defined as:

    .. math::
        F_j (x) = \frac{1}{\Gamma (j+1)} \int_0^{\infty} \frac{t^j}{\exp{(t-x)} + 1} dt

    for :math:`j > 0`.

    This is equivalent to the following `polylogarithm
    <https://en.wikipedia.org/wiki/Polylogarithm>`_ function:

    .. math::
        F_j (x) = -Li_{j+1}\left(-e^{x}\right)

    Warning: at present this function is limited to relatively small
    arguments due to limitations in the `~mpmath` package's
    implementation of `~mpmath.polylog`.

    Examples
    --------
    >>> Fermi_integral(0, 0)
    (0.6931471805599453-0j)
    >>> Fermi_integral(1, 0)
    (1.3132616875182228-0j)
    >>> Fermi_integral(1, 1)
    (1.8062860704447743-0j)

    """
    try:
        from mpmath import polylog
    except (ImportError, ModuleNotFoundError) as e:
        from plasmapy.optional_deps import mpmath_import_error

        raise mpmath_import_error from e

    if isinstance(x, (numbers.Integral, numbers.Real, numbers.Complex)):
        arg = -np.exp(x)
        integral = -1 * complex(polylog(j + 1, arg))
        return integral
    elif isinstance(x, np.ndarray):
        integral_arr = np.zeros_like(x, dtype="complex")
        for idx, val in enumerate(x):
            integral_arr[idx] = -1 * complex(polylog(j + 1, -np.exp(val)))
        return integral_arr
    else:
        raise TypeError(f"Improper type {type(x)} given for argument x.")
