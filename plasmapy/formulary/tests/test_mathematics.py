"""Tests for the mathematics.py module."""


import numpy as np
import pytest

from plasmapy.formulary import mathematics as mathematics


# TODO: Move the Fermi integral tests over to this file?


def test_rot_a_to_b():

    # Test one rotation
    a = np.array([1,0,0])
    b = np.array([0,1,0])
    R = mathematics.rot_a_to_b(a,b)
    correct = np.array([[0, 1, 0], [-1, 0, 0], [0,0,1]])
    assert np.allclose(R, correct, atol=1e-6)

    # Test the special case where the two vectors are parallel
    a = np.array([1,0,0])
    b = np.array([1,0,0])
    R = mathematics.rot_a_to_b(a,b)
    correct = np.identity(3)
    assert np.allclose(R, correct, atol=1e-6)

    # Test the special case where the two vectors are anti-parallel
    a = np.array([1,0,0])
    b = np.array([-1,0,0])
    R = mathematics.rot_a_to_b(a,b)
    correct = - np.identity(3)
    assert np.allclose(R, correct, atol=1e-6)




if __name__ == '__main__':
    """
    test_rot_a_to_b()
    """
