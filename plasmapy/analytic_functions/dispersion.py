"""Functions related to the plasma dispersion function"""

import numpy as np
import scipy as sp

def plasma_dispersion_func(zeta):
    """Calculate the plasma dispersion function

    Parameters
    ----------
    zeta : complex number or array
        Argument of plasma dispersion function.

    
    Returns
    -------
    Z : complex number or array
        Value of plasma dispersion function.

    Raises
    ------
    ValueError
        Invalid argument.

    """

    if not np.all(np.isfinite(zeta)):
        raise ValueError("Argument of plasma_dispersion_function " + \
                             "must be a complex number or array")

    Z = 1j * np.sqrt(np.pi) * np.exp(-zeta**2) * \
        (1.0 + sp.special.erf(1j*zeta))

    return Z

def plasma_dispersion_func_deriv(zeta):
    """Calculate the derivative of the plasma dispersion function
    
    Parameters
    ----------
    zeta : complex number or array
        Argument of plasma dispersion function.
 
    Returns
    -------
    Zprime : complex number or array
        Derivative of plasma dispersion function.

    Raises
    ------
    ValueError
        Invalid argument.

    """

    if not np.all(np.isfinite(zeta)):
        raise ValueError("Argument of plasma_dispersion_function " + \
                             "must be a complex number or array")
    
    Zprime = -2*(1+zeta*Z(zeta))

    return Zprime
