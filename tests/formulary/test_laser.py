"""Tests for functionality contained in `plasmapy.formulary.laser`."""

import astropy.units as u
import numpy as np
import pytest
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.laser import (
    E0,
)


@pytest.mark.parametrize(
    ("intensity", "expected"),
    [
        (1e-3 * u.watt / u.m**2, 0.8680211 * u.V / u.m),
        ([1, 0] * u.milliWatt / u.m**2, [0.8680211, 0] * u.V / u.m),
        (np.nan * u.watt / u.m**2, np.nan * u.V / u.m),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_E0(intensity, expected) -> None:
    result = E0(intensity=intensity)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.V / u.m


@pytest.mark.parametrize(
    ("intensity", "expected"),
    [
        (-5e4 * u.Watt / u.m**2, ValueError),
        (1 * u.kg, u.UnitTypeError),
    ],
)
def test_E0_errors(intensity, expected) -> None:
    with pytest.raises(expected):
        E0(intensity=intensity)


@pytest.mark.parametrize(
    ("intensity", "expected_warning"),
    [
        (5, u.UnitsWarning),
    ],
)
def test_E0_warnings(intensity, expected_warning) -> None:
    with pytest.warns(expected_warning):
        E0(intensity=intensity)
