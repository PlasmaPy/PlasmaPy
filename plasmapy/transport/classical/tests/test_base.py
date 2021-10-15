import pytest
import numpy as np
import astropy.units as u

import plasmapy.transport.classical.base as base



chi_e = 10**(np.linspace(-2, 2, num=50))
chi_i = 10**(np.linspace(-3, 3, num=50))

test_input = [
                (chi_e, chi_i, 10, None), # Nominal input
                (chi_e, chi_i, None, ValueError), # Z cannot be None
                (None, None, 10, ValueError), # chi_e and chi_i cannot both be None
                (chi_e, None, 10, None), # Either chi_e or chi_i can be None
                (None, chi_i, 10, None), # Either chi_e or chi_i can be None
              ]

@pytest.mark.parametrize("chi_e,chi_i,Z,error", test_input)
def test_constructor_dimensionless(chi_e, chi_i, Z, error):
    
    if error is None:
        base.AbstractClassicalTransportCoefficients(chi_e, chi_i, Z)
        
    else:
        with pytest.raises(error):
            base.AbstractClassicalTransportCoefficients(chi_e, chi_i, Z)


