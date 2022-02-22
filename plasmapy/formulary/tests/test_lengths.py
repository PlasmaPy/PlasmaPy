"""Tests for functionality contained in `plasmapy.formulary.lengths`."""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.lengths import (
    cwp_,
    Debye_length,
    inertial_length,
    lambdaD_,
)
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


def test_inertial_length():
    r"""Test the inertial_length function in parameters.py."""

    assert inertial_length(n_i, particle="p").unit.is_equivalent(u.m)

    assert np.isclose(
        inertial_length(mu * u.cm ** -3, particle="p").cgs.value, 2.28e7, rtol=0.01
    )

    inertial_length_electron_plus = inertial_length(5.351 * u.m ** -3, particle="e+")
    assert inertial_length_electron_plus == inertial_length(
        5.351 * u.m ** -3, particle="e"
    )

    assert inertial_length(n_i, particle="p") == inertial_length(n_i, particle="p")

    with pytest.warns(u.UnitsWarning):
        inertial_length(4, particle="p")

    with pytest.raises(u.UnitTypeError):
        inertial_length(4 * u.m ** -2, particle="p")

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m ** -3, particle="p")

    with pytest.raises(InvalidParticleError):
        inertial_length(n_i, particle=-135)

    with pytest.warns(u.UnitsWarning):
        inertial_length_no_units = inertial_length(1e19, particle="p")
        assert inertial_length_no_units == inertial_length(
            1e19 * u.m ** -3, particle="p"
        )

    assert inertial_length(n_e, "e-").unit.is_equivalent(u.m)

    assert np.isclose(
        inertial_length(1 * u.cm ** -3, "e-").cgs.value, 5.31e5, rtol=1e-3
    )

    with pytest.warns(u.UnitsWarning):
        inertial_length(5, "e-")

    with pytest.raises(u.UnitTypeError):
        inertial_length(5 * u.m, "e-")

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m ** -3, "e-")

    with pytest.warns(u.UnitsWarning):
        assert inertial_length(1e19, "e-") == inertial_length(1e19 * u.m ** -3, "e-")

    assert_can_handle_nparray(inertial_length)


@pytest.mark.parametrize(
    "alias, parent",
    [
        (lambdaD_, Debye_length),
        (cwp_, inertial_length),
    ],
)
def test_parameters_aliases(alias, parent):
    """Test all aliases defined in parameters.py"""
    assert alias is parent
