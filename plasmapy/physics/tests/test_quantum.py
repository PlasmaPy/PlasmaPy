if __name__ == "__main__":    
    import numpy as np
    import pytest
    import astropy.units as u
    from plasmapy.constants import c, h
    from plasmapy.physics.quantum import(deBroglie_wavelength, 
                                         thermal_deBroglie_wavelength, 
                                         Fermi_energy, 
                                         Thomas_Fermi_length)
else:
    import numpy as np
    import pytest
    import astropy.units as u
    from ...constants import c, h
    from ..quantum import (deBroglie_wavelength, 
                           thermal_deBroglie_wavelength, 
                           Fermi_energy, 
                           Thomas_Fermi_length)

def test_deBroglie_wavelength():

    dbwavelength1 = deBroglie_wavelength(2e7*u.cm/u.s, 'e')
    assert np.isclose(dbwavelength1.value, 3.628845222852886e-11)
    assert dbwavelength1.unit == u.m

    dbwavelength2 = deBroglie_wavelength(0*u.m/u.s, 'e')
    assert dbwavelength2 == np.inf*u.m

    V_array = np.array([2e5, 0])*u.m/u.s
    dbwavelength_arr = deBroglie_wavelength(V_array, 'e')

    assert np.isclose(dbwavelength_arr.value[0], 3.628845222852886e-11)
    assert dbwavelength_arr.value[1] == np.inf
    assert dbwavelength_arr.unit == u.m

    V_array = np.array([2e5, 2e5])*u.m/u.s
    dbwavelength_arr = deBroglie_wavelength(V_array, 'e')

    assert np.isclose(dbwavelength_arr.value[0], 3.628845222852886e-11)
    assert np.isclose(dbwavelength_arr.value[1], 3.628845222852886e-11)
    assert dbwavelength_arr.unit == u.m

    assert deBroglie_wavelength(-5e5*u.m/u.s, 'p') == \
        deBroglie_wavelength(5e5*u.m/u.s, 'p')

    assert deBroglie_wavelength(-5e5*u.m/u.s, 'e+') == \
        deBroglie_wavelength(5e5*u.m/u.s, 'e')

    assert deBroglie_wavelength(1*u.m/u.s, 5*u.kg) == \
        deBroglie_wavelength(100*u.cm/u.s, 5000*u.g)

    with pytest.raises(ValueError):
        deBroglie_wavelength(c*1.000000001, 'e')

    with pytest.raises(UserWarning):
        deBroglie_wavelength(0.79450719277, 'Be-7 1+')

    with pytest.raises(u.UnitConversionError):
        deBroglie_wavelength(8*u.m/u.s, 5*u.m)

    with pytest.raises(ValueError):
        deBroglie_wavelength(8*u.m/u.s, 'sddsf')

# defining some plasma parameters for tests
T_e = 1 * u.eV
n_e = 1e23 * u.cm**-3
# should probably change this to use unittest module
# add tests for numpy arrays as inputs
# add tests for different astropy units (random fuzzing method?)
def test_thermal_deBroglie_wavelength():
    r"""Test the thermal_deBroglie_wavelength function in quantum.py."""
    lambda_dbTh = thermal_deBroglie_wavelength(T_e)
    # test a simple case for expected value
    assert np.isclose(lambda_dbTh.value,
                      6.919367518364532e-10)
    # testing returned units
    assert lambda_dbTh.unit == u.m

def test_Fermi_energy():
    r"""Test the Fermi_energy function in quantum.py."""
    energy_F = Fermi_energy(n_e)
    # test a simple case for expected value
    assert np.isclose(energy_F.value,
                      1.2586761116196002e-18)
    # testing returned units
    assert energy_F.unit == u.J

def test_Thomas_Fermi_length():
    r"""Test the Thomas_Fermi_length function in quantum.py."""
    lambda_TF = Thomas_Fermi_length(n_e)
    # test a simple case for expected value
    assert np.isclose(lambda_TF.value,
                      5.379914085596706e-11)
    # testing returned units
    assert lambda_TF.unit == u.m
    