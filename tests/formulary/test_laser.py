"""Tests for functionality contained in `plasmapy.formulary.laser`."""

import astropy.units as u
import numpy as np
import pytest
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.laser import (
    E0_,
    I_,
    Gaussian_beam_waist_radius,
    Gaussian_power,
    Gaussian_Rayleigh_length,
    Gaussian_spot_size_FWHM,
    a0_,
    electric_field_amplitude,
    em_angular_frequency,
    em_wavelength,
    intensity,
    normalized_vector_potential,
    omega_,
    w0_,
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
        (I_, intensity),
        (w0_, Gaussian_beam_waist_radius),
        (omega_, em_angular_frequency),
        (a0_, normalized_vector_potential),
    ],
)
def test_aliases(alias, parent) -> None:
    """Test all aliases defined in laser.py"""
    assert alias is parent


@pytest.mark.parametrize(
    ("electric_field_amplitude", "expected"),
    [
        (0.8680211 * u.V / u.m, 1e-3 * u.watt / u.m**2),
        ([0.8680211, 0] * u.V / u.m, [1, 0] * u.milliWatt / u.m**2),
        (np.nan * u.V / u.m, np.nan * u.watt / u.m**2),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_intensity(electric_field_amplitude, expected) -> None:
    result = intensity(electric_field_amplitude=electric_field_amplitude)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.Watt / u.m**2


@pytest.mark.parametrize(
    ("electric_field_amplitude", "expected"),
    [
        (-5e4 * u.V / u.m, ValueError),
        (1 * u.kg, u.UnitTypeError),
    ],
)
def test_intensity_errors(electric_field_amplitude, expected) -> None:
    with pytest.raises(expected):
        intensity(electric_field_amplitude=electric_field_amplitude)


@pytest.mark.parametrize(
    ("electric_field_amplitude", "expected_warning"),
    [
        (5, u.UnitsWarning),
    ],
)
def test_intensity_warnings(electric_field_amplitude, expected_warning) -> None:
    with pytest.warns(expected_warning):
        intensity(electric_field_amplitude=electric_field_amplitude)


@pytest.mark.parametrize(
    ("intensity", "beam_waist_radius", "expected"),
    [
        (1e-3 * u.watt / u.m**2, 1e-6 * u.m, 1.5707963267948967e-15 * u.Watt),
        (
            [1, 1e21] * u.milliWatt / u.m**2,
            [1, 1] * u.um,
            [1.5707963267948967e-15, 1570796.3267948967] * u.Watt,
        ),
        (
            [1, 1] * u.milliWatt / u.m**2,
            [1, 10] * u.um,
            [1.5707963267948967e-15, 1.5707963267948967e-13] * u.Watt,
        ),
        (
            1 * u.milliWatt / u.m**2,
            [1, 10] * u.um,
            [1.5707963267948967e-15, 1.5707963267948967e-13] * u.Watt,
        ),
        (
            [1, 1e21] * u.milliWatt / u.m**2,
            1 * u.um,
            [1.5707963267948967e-15, 1570796.3267948967] * u.Watt,
        ),
        (0 * u.watt / u.m**2, 1e-6 * u.m, 0 * u.Watt),
        (1e-3 * u.watt / u.m**2, 0 * u.m, 0 * u.Watt),
        (np.nan * u.watt / u.m**2, np.nan * u.m, np.nan * u.Watt),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_Gaussian_power(intensity, beam_waist_radius, expected) -> None:
    result = Gaussian_power(intensity=intensity, beam_waist_radius=beam_waist_radius)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.Watt


@pytest.mark.parametrize(
    ("intensity", "beam_waist_radius", "expected"),
    [
        (-5e4 * u.Watt / u.m**2, 2 * u.m, ValueError),
        (5e4 * u.Watt / u.m**2, -2 * u.m, ValueError),
        (1 * u.kg, 3 * u.s, u.UnitTypeError),
    ],
)
def test_Gaussian_power_errors(intensity, beam_waist_radius, expected) -> None:
    with pytest.raises(expected):
        Gaussian_power(intensity=intensity, beam_waist_radius=beam_waist_radius)


@pytest.mark.parametrize(
    ("intensity", "beam_waist_radius", "expected_warning"),
    [(5, 2 * u.um, u.UnitsWarning), (5 * u.Watt / u.m**2, 2, u.UnitsWarning)],
)
def test_Gaussian_power_errors_warnings(
    intensity, beam_waist_radius, expected_warning
) -> None:
    with pytest.warns(expected_warning):
        Gaussian_power(intensity=intensity, beam_waist_radius=beam_waist_radius)


@pytest.mark.parametrize(
    ("spot_size_FWHM", "expected"),
    [
        (8.242e-6 * u.m, 7.00011e-6 * u.m),
        ([8.242, 0] * u.um, [7.00011e-6, 0] * u.m),
        (np.nan * u.m, np.nan * u.m),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_Gaussian_beam_waist_radius(spot_size_FWHM, expected) -> None:
    result = Gaussian_beam_waist_radius(spot_size_FWHM=spot_size_FWHM)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.m


@pytest.mark.parametrize(
    ("spot_size_FWHM", "expected"),
    [
        (-5e4 * u.m, ValueError),
        (1 * u.kg, u.UnitTypeError),
    ],
)
def test_Gaussian_beam_waist_radius_errors(spot_size_FWHM, expected) -> None:
    with pytest.raises(expected):
        Gaussian_beam_waist_radius(spot_size_FWHM=spot_size_FWHM)


@pytest.mark.parametrize(
    ("spot_size_FWHM", "expected_warning"),
    [
        (5, u.UnitsWarning),
    ],
)
def test_Gaussian_beam_waist_radius_warnings(spot_size_FWHM, expected_warning) -> None:
    with pytest.warns(expected_warning):
        Gaussian_beam_waist_radius(spot_size_FWHM=spot_size_FWHM)


@pytest.mark.parametrize(
    ("beam_waist_radius", "expected"),
    [
        (7e-6 * u.m, 8.24187e-6 * u.m),
        ([7, 0] * u.um, [8.24187e-6, 0] * u.m),
        (np.nan * u.m, np.nan * u.m),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_Gaussian_spot_size_FWHM(beam_waist_radius, expected) -> None:
    result = Gaussian_spot_size_FWHM(beam_waist_radius=beam_waist_radius)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.m


@pytest.mark.parametrize(
    ("beam_waist_radius", "expected"),
    [
        (-5e4 * u.m, ValueError),
        (1 * u.kg, u.UnitTypeError),
    ],
)
def test_Gaussian_spot_size_FWHM_errors(beam_waist_radius, expected) -> None:
    with pytest.raises(expected):
        Gaussian_spot_size_FWHM(beam_waist_radius=beam_waist_radius)


@pytest.mark.parametrize(
    ("beam_waist_radius", "expected_warning"),
    [
        (5, u.UnitsWarning),
    ],
)
def test_Gaussian_spot_size_FWHM_warnings(beam_waist_radius, expected_warning) -> None:
    with pytest.warns(expected_warning):
        Gaussian_spot_size_FWHM(beam_waist_radius=beam_waist_radius)


@pytest.mark.parametrize(
    ("angular_frequency", "expected"),
    [
        (2.354307546e15 * u.rad / u.s, 800.0873e-9 * u.m),
        ([2.354307546, 0] * u.Prad / u.s, [800.0873e-9, np.inf] * u.m),
        (np.nan * u.rad / u.s, np.nan * u.m),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_em_wavelength(angular_frequency, expected) -> None:
    result = em_wavelength(angular_frequency=angular_frequency)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.m


@pytest.mark.parametrize(
    ("angular_frequency", "expected"),
    [
        (-5e4 * u.rad / u.s, ValueError),
        (1 * u.kg, u.UnitTypeError),
    ],
)
def test_em_wavelength_errors(angular_frequency, expected) -> None:
    with pytest.raises(expected):
        em_wavelength(angular_frequency=angular_frequency)


@pytest.mark.parametrize(
    ("angular_frequency", "expected_warning"),
    [
        (5, u.UnitsWarning),
    ],
)
def test_em_wavelength_warnings(angular_frequency, expected_warning) -> None:
    with pytest.warns(expected_warning):
        em_wavelength(angular_frequency=angular_frequency)


@pytest.mark.parametrize(
    ("wavelength", "expected"),
    [
        (800e-9 * u.m, 2.354564e15 * u.rad / u.s),
        ([800, 0] * u.nm, [2.354564e15, np.inf] * u.rad / u.s),
        (np.nan * u.m, np.nan * u.rad / u.s),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_em_angular_frequency(wavelength, expected) -> None:
    result = em_angular_frequency(wavelength=wavelength)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.rad / u.s


@pytest.mark.parametrize(
    ("wavelength", "expected"),
    [
        (-5e4 * u.m, ValueError),
        (1 * u.kg, u.UnitTypeError),
    ],
)
def test_em_angular_frequency_errors(wavelength, expected) -> None:
    with pytest.raises(expected):
        em_angular_frequency(wavelength=wavelength)


@pytest.mark.parametrize(
    ("wavelength", "expected_warning"),
    [
        (5, u.UnitsWarning),
    ],
)
def test_em_angular_frequency_warning(wavelength, expected_warning) -> None:
    with pytest.warns(expected_warning):
        em_angular_frequency(wavelength=wavelength)


@pytest.mark.parametrize(
    ("intensity", "wavelength", "expected"),
    [
        (1e-3 * u.watt / u.m**2, 800e-9 * u.m, 2.162820076644342e-13),
        (
            [1, 1e-3] * u.milliWatt / u.m**2,
            [800, 800] * u.nm,
            [2.162820076644342e-13, 6.839437611336064e-15],
        ),
        (
            [1, 1] * u.milliWatt / u.m**2,
            [800, 650] * u.nm,
            [2.162820076644342e-13, 1.7572913122735275e-13],
        ),
        (
            [1, 1e-3] * u.milliWatt / u.m**2,
            800 * u.nm,
            [2.162820076644342e-13, 6.839437611336064e-15],
        ),
        (
            1 * u.milliWatt / u.m**2,
            [800, 650] * u.nm,
            [2.162820076644342e-13, 1.7572913122735275e-13],
        ),
        (0 * u.watt / u.m**2, 800e-9 * u.m, 0),
        (1e-3 * u.watt / u.m**2, 0 * u.m, 0),
        (np.nan * u.Watt / u.m**2, np.nan * u.m, np.nan),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_normalized_vector_potential(intensity, wavelength, expected) -> None:
    result = normalized_vector_potential(intensity=intensity, wavelength=wavelength)
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)


#    assert result.unit == u.dimensionless_unscaled


@pytest.mark.parametrize(
    ("intensity", "wavelength", "expected"),
    [
        (-5e-3 * u.Watt / u.m**2, 5e4 * u.m, ValueError),
        (5e-3 * u.Watt / u.m**2, -5e4 * u.m, ValueError),
        (7 * u.s, 1 * u.kg, u.UnitTypeError),
    ],
)
def test_normalized_vector_potential_errors(intensity, wavelength, expected) -> None:
    with pytest.raises(expected):
        normalized_vector_potential(intensity=intensity, wavelength=wavelength)


@pytest.mark.parametrize(
    ("intensity", "wavelength", "expected_warning"),
    [
        (3, 5, u.UnitsWarning),
    ],
)
def test_normalized_vector_potential_warning(
    intensity, wavelength, expected_warning
) -> None:
    with pytest.warns(expected_warning):
        normalized_vector_potential(intensity=intensity, wavelength=wavelength)


@pytest.mark.parametrize(
    ("wavelength", "beam_waist_radius", "expected"),
    [
        (800e-9 * u.m, 1e-6 * u.m, 3.926990816987241e-06 * u.m),
        (
            [800, 650] * u.nm,
            [1, 1] * u.um,
            [3.926990816987241e-06, 4.8332194670612195e-06] * u.m,
        ),
        (
            [800, 800] * u.nm,
            [1, 20] * u.um,
            [3.926990816987241e-06, 0.0015707963267948962] * u.m,
        ),
        (
            [800, 650] * u.nm,
            1 * u.um,
            [3.926990816987241e-06, 4.8332194670612195e-06] * u.m,
        ),
        (
            800 * u.nm,
            [1, 20] * u.um,
            [3.926990816987241e-06, 0.0015707963267948962] * u.m,
        ),
        (800e-9 * u.m, 0 * u.m, 0 * u.m),
        (0 * u.m, 1e-6 * u.m, np.inf * u.m),
        (np.nan * u.m, np.nan * u.m, np.nan * u.m),
    ],
)
@pytest.mark.filterwarnings("ignore::astropy.units.UnitsWarning")
def test_Gaussian_Rayleigh_length(wavelength, beam_waist_radius, expected) -> None:
    result = Gaussian_Rayleigh_length(
        wavelength=wavelength, beam_waist_radius=beam_waist_radius
    )
    assert_quantity_allclose(result, expected, rtol=1e-6, equal_nan=True, verbose=True)
    assert result.unit == u.m


@pytest.mark.parametrize(
    ("wavelength", "beam_waist_radius", "expected"),
    [
        (-5e4 * u.m, 2 * u.m, ValueError),
        (5e4 * u.m, -2 * u.m, ValueError),
        (1 * u.kg, 3 * u.s, u.UnitTypeError),
    ],
)
def test_Gaussian_Rayleigh_length_errors(
    wavelength, beam_waist_radius, expected
) -> None:
    with pytest.raises(expected):
        Gaussian_Rayleigh_length(
            wavelength=wavelength, beam_waist_radius=beam_waist_radius
        )


@pytest.mark.parametrize(
    ("wavelength", "beam_waist_radius", "expected_warning"),
    [
        (5 * u.m, 2, u.UnitsWarning),
        (5, 2 * u.m, u.UnitsWarning),
    ],
)
def test_Rayleigh_length_warnings_warnings(
    wavelength, beam_waist_radius, expected_warning
) -> None:
    with pytest.warns(expected_warning):
        Gaussian_Rayleigh_length(
            wavelength=wavelength, beam_waist_radius=beam_waist_radius
        )
