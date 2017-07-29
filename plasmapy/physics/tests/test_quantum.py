import numpy as np
import pytest
import astropy.units as u
from ...constants import c, h
from ..quantum import deBroglie_wavelength


def test_deBroglie_wavelength():

    dbwavelength1 = deBroglie_wavelength(2e7*u.cm/u.s, 'e')
    assert np.isclose(dbwavelength1.value, 3.628845222852886e-11)
    assert dbwavelength1.unit == u.m

    dbwavelength2 = deBroglie_wavelength(0*u.m/u.s, 'e')
    assert dbwavelength2 == np.inf*u.m

    assert deBroglie_wavelength(-5e5*u.m/u.s, 'p') == \
        deBroglie_wavelength(5e5*u.m/u.s, 'p')

    with pytest.raises(ValueError):
        deBroglie_wavelength(c*1.000000001, 'e')

    with pytest.raises(UserWarning):
        # Use the Kolakowski constant for no reason besides that it's cool
        deBroglie_wavelength(0.794507192779479276240362415636045646, 'Be-7 1+')
