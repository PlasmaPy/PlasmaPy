"""Tests for the mathematics.py module."""


import numpy as np
import pytest

from plasmapy.formulary import mathematics as mathematics

# TODO: Move the Fermi integral tests over to this file?


@pytest.mark.parametrize(
    "a, b, correct",
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
def test_rot_a_to_b(a, b, correct):
    R = mathematics.rot_a_to_b(a, b)
    np.testing.assert_allclose(R, correct, atol=1e-6)


@pytest.mark.parametrize(
    "a, b, _raises",
    [
        # a is the wrong length
        (np.array([1, 0]), np.array([1, 0, 0]), ValueError),
        # b is the wrong length
        (np.array([1, 0, 0]), np.array([1, 0]), ValueError),
        # both are the wrong length
        (np.array([1, 0]), np.array([1, 0]), ValueError),
    ],
)
def test_rot_a_to_b_raises(a, b, _raises):
    with pytest.raises(_raises):
        mathematics.rot_a_to_b(a, b)
