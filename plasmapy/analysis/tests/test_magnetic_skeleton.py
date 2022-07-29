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
    return [(y - 5.5), (z - 4.5) * (z - 6.5), (x - 5.5)]


def vspace_func_2(x, y, z):
    return [(1 - z) * x + y + 0.1, -2 * x - (1 + z) * y, z**2 + 0.3 * y - 0.2]


def vspace_func_3(x, y, z):
    return [(2 - z) * x - y + 0.2, 3 * x - (2 + z) * y - 0.2, z**2 - 0.2]


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
    nullpoints, spines, fan, seperators = magnetic_skeleton_find(
        x_arr, y_arr, z_arr, u, v, w
    )
    # print(nullpoints)
    # print("##############")
    # print(spines)
    # print("##############")
    # print(fan)
    # print("##############")
    # print(seperators)


def test_2():
    x_arr = np.linspace(-5, 5, 60)
    y_arr = np.linspace(-5, 5, 60)
    z_arr = np.linspace(-5, 5, 60)
    x, y, z = np.meshgrid(
        x_arr,
        y_arr,
        z_arr,
        indexing="ij",
    )
    u, v, w = vspace_func_2(x, y, z)
    nullpoints, spines, fan, seperators = magnetic_skeleton_find(
        x_arr, y_arr, z_arr, u, v, w
    )
    # print(nullpoints)
    # print("##############")
    # print(spines)
    # print("##############")
    # print(fan)
    # print("##############")
    # print(seperators)


def test_3():
    x_arr = np.linspace(-5, 5, 60)
    y_arr = np.linspace(-5, 5, 60)
    z_arr = np.linspace(-5, 5, 60)
    x, y, z = np.meshgrid(
        x_arr,
        y_arr,
        z_arr,
        indexing="ij",
    )
    u, v, w = vspace_func_3(x, y, z)
    nullpoints, spines, fan, seperators = magnetic_skeleton_find(
        x_arr, y_arr, z_arr, u, v, w
    )
    # for elem in seperators:
    #     print("###########")
    #     for p in elem:
    #         print(p.loc)
