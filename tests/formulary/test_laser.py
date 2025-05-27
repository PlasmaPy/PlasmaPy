"""Tests for functionality contained in `plasmapy.formulary.laser`."""


import astropy.units as u
import pytest
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.laser import (
    E0,
)


@pytest.mark.parametrize(
    ("Intensity", "expected"),
    [
        (1e-3 * u.watt / u.m**2, 0.8680211 * u.V / u.m),
        ([1, 0] * u.milliWatt, [0.8680211, 0] * u.V / u.m),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_E0(Intensity, expected) -> None:
    result = E0(Intensity = Intensity)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.V / u.m


@pytest.mark.parametrize(
    ("Intensity", "expected"),
    [
        (-5e4 * u.Watt / u.m**-2, ValueError),
        (1 * u.kg, u.UnitTypeError),
    ],
)
def test_E0_errors(Intensity, expected) -> None:
    with pytest.raises(expected):
        E0(Intensity = Intensity)


@pytest.mark.parametrize(
    ("Intensity", "expected_warning"),
    [
        (5, u.UnitsWarning),
    ],
)
def test_E0_warnings(Intensity, expected_warning) -> None:
    with pytest.warns(expected_warning):
        E0(Intensity = Intensity)


