"""Tests for functionality contained in `plasmapy.formulary.laser`."""

import astropy.units as u
import numpy as np
import pytest
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.laser import E0_, electric_field_amplitude


@pytest.mark.parametrize(
    ("intensity", "expected"),
    [
        (1e-3 * u.watt / u.m**2, 0.8680211 * u.V / u.m),
        ([1, 0] * u.milliWatt / u.m**2, [0.8680211, 0] * u.V / u.m),
        (np.nan * u.watt / u.m**2, np.nan * u.V / u.m),
    ],
)

@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_electric_field_amplitude(intensity, expected) -> None:
    result = electric_field_amplitude(intensity=intensity)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.V / u.m

@pytest.mark.parametrize(
    ("intensity", "expected"),
    [
        (-5e4 * u.Watt / u.m**2, ValueError),
        (1 * u.kg, u.UnitTypeError),
    ],
)
def test_electric_field_amplitude_errors(intensity, expected) -> None:
    with pytest.raises(expected):
        electric_field_amplitude(intensity=intensity)

@pytest.mark.parametrize(
    ("intensity", "expected_warning"),
    [
        (5, u.UnitsWarning),
    ],
)
def test_electric_field_amplitude_warnings(intensity, expected_warning) -> None:
    with pytest.warns(expected_warning):
        electric_field_amplitude(intensity=intensity)

@pytest.mark.parametrize(
    ("alias", "parent"),
    [
        (E0_, electric_field_amplitude),
    ],
)

def test_aliases(alias, parent) -> None:
    """Test all aliases defined in laser.py"""
    assert alias is parent