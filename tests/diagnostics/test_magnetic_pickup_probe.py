"""Unit test for the magnetic pickup probe (B-dot probe) diagnostics."""

import numpy as np
import astropy.units as u
import pytest

from plasmapy.diagnostics.magnetic_pickup_probe import compute_bfield


def test_compute_bfield_integrates_voltage_to_field():
    """
    Test compute_bfield by inputting a sine wave voltage and checking
    if the output magnetic field matches the integrated waveform.
    """
    # Create a time array in seconds
    time = np.linspace(0, 2 * np.pi, 1000) * u.s

    # Create a sine voltage signal: V(t) = sin(t)
    voltage = np.sin(time.value) * u.V

    # Set probe parameters
    loop_area = 1.0 * u.m**2
    num_loop = 1
    gain = 1.0

    # Expected magnetic field is the integral of V(t)/(N*A*gain)
    expected_field = -np.cos(time.value) + 1  # Integral of sin(t) is -cos(t), plus constant
    expected_field *= (1 / (loop_area * num_loop * gain)).value  # Apply scaling
    expected_field = expected_field * u.T  # Assign Tesla units

    # Compute B field using function
    computed_bfield = compute_bfield(voltage, time, loop_area, num_loop, gain)

    # Use a tolerance since numerical integration introduces small errors
    assert u.allclose(computed_bfield, expected_field, rtol=1e-3, atol=1e-5 * u.T)
