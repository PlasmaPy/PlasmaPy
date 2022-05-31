import numpy as np
import warnings

from plasmapy.analysis.nullpoint import (
    _bilinear_root,
    _locate_null_point,
    _reduction,
    _trilinear_analysis,
    _trilinear_coeff_cal,
    _trilinear_jacobian,
    _vector_space,
    _vspace_iterator,
    MultipleNullPointWarning,
    NonZeroDivergence,
    null_point_find,
    NullPointError,
    NullPointWarning,
    trilinear_approx,
    uniform_null_point_find,
)

_EQUALITY_ATOL = 1e-10


def cell_find(vspace, np):
    deltax, deltay, deltaz = vspace[2]
    location = np.loc
    i = 0
    j = 0
    k = 0
    xsum = 0
    ysum = 0
    zsum = 0
    # for dx in range(deltax-1):
    #     if location[0]

    return np.array([i, j, k])


def spine_find(
    x_arr=None,
    y_arr=None,
    z_arr=None,
    u_arr=None,
    v_arr=None,
    w_arr=None,
    nullpoint_list=np.array([]),
):
    vspace = _vector_space(
        x_arr,
        y_arr,
        z_arr,
        None,
        None,
        None,
        u_arr,
        v_arr,
        w_arr,
        None,
        None,
    )
    for np in nullpoint_list:
        jcb = _trilinear_jacobian(vspace, cell)
        M = jcb(np.loc[0], np.loc[1], np.loc[2])
        if not np.isclose(np.trace(M), 0, atol=_EQUALITY_ATOL):
            raise NonZeroDivergence()
        eigen_vals, eigen_vectors = np.linalg.eig(M)


def fan_find(
    x_arr=None,
    y_arr=None,
    z_arr=None,
    u_arr=None,
    v_arr=None,
    w_arr=None,
    nullpoint_list=np.array([]),
):
    vspace = _vector_space(
        x_arr,
        y_arr,
        z_arr,
        None,
        None,
        None,
        u_arr,
        v_arr,
        w_arr,
        None,
        None,
    )


def magnetic_skeleton_find(
    x_arr=None,
    y_arr=None,
    z_arr=None,
    u_arr=None,
    v_arr=None,
    w_arr=None,
):
    ring_layers = np.array([])
    nullpoints = null_point_find(x_arr, y_arr, z_arr, u_arr, v_arr, w_arr)
