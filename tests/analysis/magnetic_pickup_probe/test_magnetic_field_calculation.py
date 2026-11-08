"""Tests for `plasmapy.analysis.magnetic_pickup_probe.magnetic_field_calculation`. 

Function integrates magnetic pickup probe (B-dot) input voltage 
into a magnetic-field function. Current unit tests check that sine 
and cosine voltages, which have known integral solutions, match
our functions expected results, and that constant (DC) and triangle (AC)
voltages, which the trapezoidal rule integrates exactly, match to machine
precision. Also implemented tests to ensure differnt input sizes raise `ValueError`.

From #3055 testing, found that atol=1e-5 was failing due to trapezoidal
integration error at 50 points exceeding this tolerance, (for sin func) therefore
by utilizing a new atol=2e-4 can allow tests to pass at 50 points.
"""

import numpy as np 
import pytest

from plasmapy.analysis.magnetic_pickup_probe import magnetic_field_calculation

# Set probe parameters
loop_area = 2.0
n_loop = 3
gain = 5

# Create time array in s
t = np.linspace(0, 2 * np.pi, 50)

@pytest.mark.parametrize(
    ("voltage", "integrated_voltage"),
    [
        (np.sin(t), 1 - np.cos(t)),
        (np.cos(t), np.sin(t)),
    ],
    ids=["sine", "cosine"],
)
def test_magnetic_field_calcultion_against_known_trig_integrals(
    voltage: np.ndarray, 
    integrated_voltage: np.ndarray
) -> None:
    """Check magnetic_field_calculation against sine & cosine voltages, whose
    integrals are known.
    """
    #expected magnetic field known from analytical solution, integral of -V(t)/(N*A*g)
    exp_field = -integrated_voltage / (loop_area * n_loop * gain)
    #call magnetic field calculation to find function's B field 
    comp_field = magnetic_field_calculation(voltage, t, loop_area, n_loop, gain)
    #ensure they're the same result with acceptable tolerances 
    np.testing.assert_allclose(comp_field, exp_field, rtol=1e-3, atol=2e-4)

def test_magnetic_field_calculation_constant_voltage() -> None:
    """Check a constant (DC) voltage, which the trapezoidal rule integrates
    exactly, so B(t) = -V*t/(N*A*g) holds to machine precision.
    """
    # constant 3 V input over the same time base
    voltage = 3.0 * np.ones_like(t)
    # exact integral of a constant is V*t, scaled by -1/(N*A*g)
    exp_field = -3.0 * t / (loop_area * n_loop * gain)
    comp_field = magnetic_field_calculation(voltage, t, loop_area, n_loop, gain)
    np.testing.assert_allclose(comp_field, exp_field, rtol=1e-14, atol=1e-15)


def test_magnetic_field_calculation_triangle_wave() -> None:
    """Check a triangle wave, which alternates (AC) but is piecewise linear, so
    the trapezoidal rule integrates it exactly to machine precision.
    """
    # one AC triangle cycle on an integer time base
    time = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    voltage = np.array([0.0, 1.0, 0.0, -1.0, 0.0])
    # exact area accumulated under the triangle at each sample
    integrated_voltage = np.array([0.0, 0.5, 1.0, 0.5, 0.0])
    exp_field = -integrated_voltage / (loop_area * n_loop * gain)
    comp_field = magnetic_field_calculation(voltage, time, loop_area, n_loop, gain)
    np.testing.assert_allclose(comp_field, exp_field, rtol=1e-14, atol=1e-15)

def test_magnetic_field_calculation_sizes() -> None:
    """Check that a size discrepancy between inputs raises `ValueError`."""
    with pytest.raises(ValueError, match="not equal"):
        magnetic_field_calculation(np.ones(5), np.linspace(0, 1, 4), loop_area=1.0)