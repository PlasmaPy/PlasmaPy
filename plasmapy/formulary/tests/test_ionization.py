import astropy.units as u
import numpy as np
import pytest
from astropy.constants import c

from plasmapy.formulary.ionization import Z_bal


def test_Z_bal():
    n = 1e19 * u.m ** -3
    T_e = 5000 * u.K

    assert (
                   Z_bal(n, T_e) * u.dimensionless_unscaled
           ).unit == u.dimensionless_unscaled

    with pytest.warns(u.UnitsWarning):
        Z_bal(1e19, T_e)

    with pytest.raises(u.UnitTypeError):
        Z_bal(2e19 * u.kg, T_e)
