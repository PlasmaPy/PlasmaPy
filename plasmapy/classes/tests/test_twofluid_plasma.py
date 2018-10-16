import pytest
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
from plasmapy.classes.plasma import AnalyticTwoFluidPlasma

@pytest.fixture
def uniform():
    return AnalyticTwoFluidPlasma(magnetic_field = lambda x: 1*u.T,
                                  electron_density= lambda x: 1/u.m**3,
                                  ion_density = lambda x: 1/u.m**3,
                                  electron_temperature = lambda x: 1*u.eV,
                                  ion_temperature = lambda x: 1 * u.K)

def test_interface(uniform):
    assert_quantity_allclose(uniform.electron_density([0, 0, 0]), 1*u.m**-3)
    assert_quantity_allclose(uniform.ion_density([0, 0, 0]), 1*u.m**-3)
    assert_quantity_allclose(uniform.electron_temperature([0, 0, 0]), 1*u.eV)
    assert_quantity_allclose(uniform.ion_temperature([0, 0, 0]), 1*u.K)
    assert_quantity_allclose(uniform.magnetic_field([0, 0, 0]), 1*u.T)

