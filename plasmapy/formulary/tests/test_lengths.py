"""Tests for functionality contained in `plasmapy.formulary.lengths`."""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.lengths import Debye_length, lambdaD_
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

Z = 1
n_e = Z * 5e19 * u.m ** -3
T_e = 1e6 * u.K


def test_Debye_length():
    r"""Test the Debye_length function in parameters.py."""

    assert Debye_length(T_e, n_e).unit.is_equivalent(u.m)

    assert np.isclose(Debye_length(1 * u.eV, 1 * u.cm ** -3).value, 7.43, atol=0.005)

    with pytest.warns(u.UnitsWarning):
        Debye_length(5, 5 * u.m ** -3)

    with pytest.raises(u.UnitTypeError):
        Debye_length(56 * u.kg, 5 * u.m ** -3)

    with pytest.raises(ValueError):
        Debye_length(5 * u.eV, -5 * u.m ** -3)

    with pytest.raises(ValueError):
        Debye_length(-45 * u.K, 5 * u.m ** -3)

    Tarr2 = np.array([1, 2]) * u.K
    narr3 = np.array([1, 2, 3]) * u.m ** -3
    with pytest.raises(ValueError):
        Debye_length(Tarr2, narr3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_length(2.0, 2.0) == Debye_length(2.0 * u.K, 2.0 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_length(2.0 * u.K, 2.0) == Debye_length(2.0, 2.0 * u.m ** -3)

    assert_can_handle_nparray(Debye_length)

@pytest.mark.parametrize(
    "alias, parent",
    [
        (lambdaD_, Debye_length),
    ],
)
def test_parameters_aliases(alias, parent):
    """Test all aliases defined in parameters.py"""
    assert alias is parent
