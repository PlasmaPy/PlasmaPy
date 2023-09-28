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
