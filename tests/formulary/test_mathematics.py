"""Tests for the mathematics.py module."""

import numpy as np
import pytest

from plasmapy.formulary import mathematics
from plasmapy.formulary.mathematics import rot_a_to_b

# TODO: Move the Fermi integral tests over to this file?


@pytest.mark.parametrize(
    ("a", "b", "correct"),
    [
        # Test one rotation
        (
            np.array([1, 0, 0]),
            np.array([0, 1, 0]),
            np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
        ),
        # Test the special case where the two vectors are parallel
        (np.array([1, 0, 0]), np.array([1, 0, 0]), np.identity(3)),
        # Test the special case where the two vectors are anti-parallel
        (np.array([1, 0, 0]), np.array([-1, 0, 0]), -np.identity(3)),
    ],
)
def test_rot_a_to_b(a, b, correct) -> None:
    R = mathematics.rot_a_to_b(a, b)
    np.testing.assert_allclose(R, correct, atol=1e-6)


@pytest.mark.parametrize(
    ("a", "b", "_raises"),
    [
        # a is the wrong length
        (np.array([1, 0]), np.array([1, 0, 0]), ValueError),
        # b is the wrong length
        (np.array([1, 0, 0]), np.array([1, 0]), ValueError),
        # both are the wrong length
        (np.array([1, 0]), np.array([1, 0]), ValueError),
    ],
)
def test_rot_a_to_b_raises(a, b, _raises) -> None:
    with pytest.raises(_raises):
        mathematics.rot_a_to_b(a, b)


def rotation_angle_degrees(R):
    """
    Compute rotation angle from a 3x3 rotation matrix.

    angle = arccos((trace(R) - 1) / 2)
    """
    trace = np.trace(R)
    cos_theta = (trace - 1.0) / 2.0

    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    theta = np.arccos(cos_theta)
    return np.degrees(theta)


def test_rot_a_to_b_antiparallel_should_be_180_degrees():
    a = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    b = -a

    R = rot_a_to_b(a, b)

    angle = rotation_angle_degrees(R)

    assert angle > 179.9, f"Rotation angle too small: {angle} degrees"
