"""Tests for calculating the net spacing."""

import numpy as np

from plasmapy.diagnostics.brl import net_spacing


def test_get_s_points():
    """Test the returned `s_points`."""
    num_points = 20
    s_end_point = 0.4

    expect_num_points = num_points
    expect_s_start_point = 0
    expect_s_end_point = s_end_point
    expect_s_spacing = s_end_point / (num_points - 1)

    s_points = net_spacing.get_s_points(num_points, s_end_point)

    assert s_points.size == expect_num_points
    assert s_points[0] == expect_s_start_point
    assert np.allclose(s_points[-1], expect_s_end_point)
    assert np.allclose(
        s_points[1:] - s_points[:-1], expect_s_spacing * np.ones(num_points - 1)
    )


class Test__private_x_and_dx_ds:
    r"""Test the different available functions for getting `x` and `dx_ds`."""

    s_points = net_spacing.get_s_points(1000, 0.8)
    ds = s_points[1] - s_points[0]

    def test_large_probe_x_and_dx_ds(self):
        r"""Test large probe x and dx/ds."""

        x_points, dx_ds_points = net_spacing._large_probe_x_and_dx_ds(self.s_points, 80)
        approximate_derivative = (x_points[2:] - x_points[:-2]) / (2 * self.ds)

        assert x_points[0] == 1
        assert np.allclose(dx_ds_points[1:-1], approximate_derivative)

    def test_medium_probe_x_and_dx_ds(self):
        r"""Test medium probe x and dx/ds."""

        x_points, dx_ds_points = net_spacing._medium_probe_x_and_dx_ds(self.s_points)
        approximate_derivative = (x_points[2:] - x_points[:-2]) / (2 * self.ds)

        assert x_points[0] == 1
        assert np.allclose(dx_ds_points[1:-1], approximate_derivative)

    def test_small_probe_x_and_dx_ds(self):
        r"""Test small probe x and dx/ds."""

        x_points, dx_ds_points = net_spacing._small_probe_x_and_dx_ds(self.s_points)
        approximate_derivative = (x_points[2:] - x_points[:-2]) / (2 * self.ds)

        assert x_points[0] == 1
        assert np.allclose(dx_ds_points[1:-1], approximate_derivative)
