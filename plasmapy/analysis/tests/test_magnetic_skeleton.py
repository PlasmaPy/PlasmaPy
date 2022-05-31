"""
Tests for the null point finder class defined in `plasmapy.analysis.nullpoint`.
"""
import numpy as np
import pytest

from plasmapy.analysis.magnetic_skeleton import (
    B_approx,
    magnetic_skeleton_find,
    spine_find,
)
from plasmapy.analysis.nullpoint import _vector_space


def vspace_func_1(x, y, z):
    return [(y - 5.5), (z - 5.5), (x - 5.5)]


def test_1():
    x_arr = np.linspace(0, 10, 46)
    y_arr = np.linspace(0, 10, 46)
    z_arr = np.linspace(0, 10, 46)
    x, y, z = np.meshgrid(
        x_arr,
        y_arr,
        z_arr,
        indexing="ij",
    )
    u, v, w = vspace_func_1(x, y, z)
    magnetic_skeleton_find(x_arr, y_arr, z_arr, u, v, w)
