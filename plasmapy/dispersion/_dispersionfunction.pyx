import numpy as np

cimport numpy as np

from scipy.special import wofz

np.import_array()

def plasma_dispersion_func_lite(np.ndarray zeta):

    Z = 1j * np.sqrt(np.pi) * wofz(zeta)

    return Z

# def plasma_dispersion_func_deriv(
#     zeta: Union[complex, int, float, np.ndarray, u.Quantity]
# ) -> Union[complex, float, np.ndarray, u.Quantity]: