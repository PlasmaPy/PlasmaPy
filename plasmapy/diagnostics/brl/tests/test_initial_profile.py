"""Test that the initial profile is correctly generated."""

import numpy as np
import pytest

from plasmapy.diagnostics.brl import initial_profile
from plasmapy.diagnostics.brl.net_spacing import get_s_points, get_x_and_dx_ds


@pytest.mark.parametrize(
    (
        "coefficients",
        "normalized_probe_potential",
        "normalized_probe_radius",
        "x_points",
        "dx_ds_points",
        "spherical",
        "expected_z",
        "expected_dxi_ds",
        "expected_eta",
    ),
    zip(
        [
            np.array([1., 1., 1., 1., 1., 1., 1., 1., 1.]),
            np.array([0., 1., 2., 3., 4., 5., 6., 7., 8.]),
        ],
        [
            1,
            10,
        ],
        [
            1,
            2**0.5,
        ],
        [
            np.array([1.0, 0.95652874, 0.91494723, 0.87517332, 0.83712843, 0.8007374, 0.76592834, 0.73263247, 0.70078401, 0.67032005]),
            np.array([1.0, 0.95122942, 0.90483742, 0.86070798, 0.81873075]),
        ],
        [
            np.array([-1.0, -0.95652874, -0.91494723, -0.87517332, -0.83712843, -0.8007374 , -0.76592834, -0.73263247, -0.70078401, -0.67032005]),
            np.array([-1.0, -0.95122942, -0.90483742, -0.86070798, -0.81873075])
        ],
        [
            True,
            False,
        ],
        [
            np.array([1.0, 0.80390446, 0.65135942, 0.53186217, 0.43758185, 0.36265505, 0.30267182, 0.25429964, 0.21500762, 0.1828636]),
            np.array([10.0, 7.68088972, 5.93381904, 4.61164307, 3.60619087])
        ],
        [
            np.array([-5.0, -3.87684065, -3.02671362, -2.37976806, -1.88464428, -1.50346003, -1.20818025, -0.97798354, -0.79734754, -0.65465412]),
            np.array([-53.33333333, -40.09340656, -30.2751221, -22.96902443, -17.51235483]),
        ],
        [
            np.array([-24.0, -16.48650496, -11.37863761, -7.89239534, -5.50292251, -3.85786865, -2.71998516, -1.92901722, -1.37635105, -0.98810852]),
            np.array([-153.33333333, -102.7242156, -69.05024154, -46.58060901, -31.54176446])
        ],
    )
)
def test_evaluate_polynomial(
    coefficients,
    normalized_probe_potential,
    normalized_probe_radius,
    x_points,
    dx_ds_points,
    spherical,
    expected_z,
    expected_dxi_ds,
    expected_eta,
):
    z, dxi_ds, eta = initial_profile.evaluate_polynomial(
        coefficients, normalized_probe_potential, normalized_probe_radius, x_points, dx_ds_points, spherical=spherical
    )

    assert np.isclose(z[0], normalized_probe_potential)
    assert np.allclose(z, expected_z)
    assert np.allclose(dxi_ds, expected_dxi_ds)
    assert np.allclose(eta, expected_eta)
    assert z.size == x_points.size
    assert dxi_ds.size == x_points.size
    assert eta.size == x_points.size
