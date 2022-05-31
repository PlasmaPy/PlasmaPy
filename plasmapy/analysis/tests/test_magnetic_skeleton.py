"""
Tests for the null point finder class defined in `plasmapy.analysis.nullpoint`.
"""
import numpy as np
import pytest

from plasmapy.analysis.magnetic_skeleton import B_approx, spine_find
from plasmapy.analysis.nullpoint import _vector_space


def vspace_func_1(x, y, z):
    return [(y - 5.5), (z - 5.5), (x - 5.5)]


def test_1():
    vspace1_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_1,
    }
    vspace1 = _vector_space(**vspace1_args)
