import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.radiation import thermal_bremsstrahlung
from plasmapy.utils.exceptions import PhysicsError


def test_thermal_bremsstrahlung():
    # Test correct spectrum created
    frequencies = (10 ** np.arange(15, 16, 0.01)) / u.s
    ne, Te = 1e22 * u.cm**-3, 1e2 * u.eV
    ion_species = "H+"

    spectrum = thermal_bremsstrahlung(frequencies, ne, Te, ion_species=ion_species)

    assert np.isclose(np.max(spectrum).value, 128.4, 1), (
        f"Spectrum maximum is {np.max(spectrum).value} "
        "instead of expected value 128.4"
    )

    # Test violates w > wpe limit
    small_frequencies = (10 ** np.arange(12, 16, 0.01)) / u.s
    with pytest.raises(PhysicsError):
        spectrum = thermal_bremsstrahlung(
            small_frequencies, ne, Te, ion_species=ion_species
        )

    # Test violates Rayleigh-Jeans limit
    small_Te = 1 * u.eV
    with pytest.raises(PhysicsError):
        spectrum = thermal_bremsstrahlung(
            frequencies, ne, small_Te, ion_species=ion_species
        )
