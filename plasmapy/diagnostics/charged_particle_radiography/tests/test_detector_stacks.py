"""
Tests for proton radiography functions
"""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.diagnostics.charged_particle_radiography.detector_stacks import (
    Layer,
    Stack,
)
from plasmapy.utils.data.downloader import get_file


@pytest.fixture()
def hdv2_stack(tmp_path):
    # Fetch stopping power data files from data module
    tissue_path = get_file("NIST_PSTAR_tissue_equivalent.txt", directory=tmp_path)
    aluminum_path = get_file("NIST_PSTAR_aluminum.txt", directory=tmp_path)

    arr = np.loadtxt(tissue_path, skiprows=8)
    eaxis = arr[:, 0] * u.MeV
    tissue_density = 1.04 * u.g / u.cm**3
    tissue_equivalent = arr[:, 1] * u.MeV * u.cm**2 / u.g * tissue_density

    arr = np.loadtxt(aluminum_path, skiprows=8)
    aluminum_density = 2.7 * u.g / u.cm**3
    aluminum = arr[:, 1] * u.MeV * u.cm**2 / u.g * aluminum_density

    # Defines the order of layers that make up a single piece of
    # HDV2 film
    HDV2 = [
        Layer(12 * u.um, eaxis, tissue_equivalent, name="2-HDV2-active"),
        Layer(
            97 * u.um, eaxis, tissue_equivalent, name="2-HDV2-substrate", active=False
        ),
    ]

    # Define a film pack consisting of an aluminum filter followed by
    # 10 layers of HDV2
    layers = [*HDV2] * 10
    layers = [
        Layer(100 * u.um, eaxis, aluminum, name="1-aluminum filter", active=False),
        *layers,
    ]

    return Stack(layers)


def test_create_layer_with_different_stopping_powers(tmp_path) -> None:
    """
    Tests the input validation for creating a Layer with either the linear
    stopping power or the mass stopping power.
    """

    aluminum_path = get_file("NIST_PSTAR_aluminum.txt", directory=tmp_path)
    arr = np.loadtxt(aluminum_path, skiprows=8)
    eaxis = arr[:, 0] * u.MeV
    mass_stopping_power = arr[:, 1] * u.MeV * u.cm**2 / u.g
    mass_density = 2.7 * u.g / u.cm**3

    # No error should be raised initializing one of these ways
    Layer(100 * u.um, eaxis, mass_stopping_power * mass_density)
    Layer(100 * u.um, eaxis, mass_stopping_power, mass_density=mass_density)

    # Error should be raised if the wrong units are provided
    with pytest.raises(ValueError):
        Layer(100 * u.um, eaxis, mass_stopping_power.value * u.kg)
    with pytest.raises(ValueError):
        Layer(
            100 * u.um,
            eaxis,
            mass_stopping_power,
            mass_density=mass_density.value * u.m,
        )

    # Error should be raised if mass_density keyword is not provided when
    # the mass stopping power is given
    with pytest.raises(ValueError):
        Layer(100 * u.um, eaxis, mass_stopping_power, mass_density=None)


def test_film_stack_num_layers(hdv2_stack) -> None:
    # Test num_layers property
    assert hdv2_stack.num_layers == 21


def test_film_stack_num_active(hdv2_stack) -> None:
    # Test num_active property
    assert hdv2_stack.num_active == 10


def test_film_stack_thickness(hdv2_stack) -> None:
    # Test thickness property
    assert np.isclose(hdv2_stack.thickness.to(u.mm).value, 1.19)


def test_film_stack_deposition_curves(hdv2_stack) -> None:
    energies = np.arange(1, 60, 1) * u.MeV
    deposition_curves = hdv2_stack.deposition_curves(energies, return_only_active=False)

    integral = np.sum(deposition_curves, axis=0)
    assert np.allclose(
        integral, 1.0
    ), "The integral over all layers for each particle species is not unity."


def test_film_stack_energy_bands_active(hdv2_stack) -> None:
    # Test energy bands
    ebands = hdv2_stack.energy_bands([0.1, 60] * u.MeV, 0.1 * u.MeV, dx=1 * u.um)

    # Expected energy bands, in MeV (only in active layers)
    expected = np.array([[3.5, 3.8], [4.6, 4.9], [5.6, 5.7], [6.4, 6.5], [7.1, 7.2]])

    assert np.allclose(ebands.to(u.MeV).value[:5, :], expected, atol=0.15)


def test_film_stack_energy_bands_inum_active(hdv2_stack) -> None:
    # Test including inum_active layers
    ebands = hdv2_stack.energy_bands(
        [0.1, 60] * u.MeV, 0.1 * u.MeV, dx=1 * u.um, return_only_active=False
    )
    # Expected first 5 energy bands
    expected = np.array([[0.1, 4.2], [3.5, 3.8], [3.9, 5.1], [4.6, 4.9], [4.9, 6]])
    assert np.allclose(ebands.to(u.MeV).value[0:5, :], expected, atol=0.15)
