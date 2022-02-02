"""
Tests for the null point finder class defined in `plasmapy.analysis.nullpoint`.
"""
import numpy as np
import pytest

from plasmapy.analysis.nullpoint import (
    _ATOL,
    bilinear_root,
    locate_null_point,
    null_point_find,
    reduction,
    trilinear_analysis,
    trilinear_approx,
    trilinear_coeff_cal,
    trilinear_jacobian,
    vector_space,
)


def vspace_func_1(x, y, z):
    return [(y - 5.5), (z - 5.5), (x - 5.5)]


def vspace_func_2(x, y, z):
    return [2 * y - z - 5.5, 3 * x + z - 22, x ** 2 - 11 * x + y + 24.75]


def test_trilinear_coeff_cal():
    r"""Test `~plasmapy.analysis.nullpoint.trilinear_coeff_cal`."""
    vspace1_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_1,
    }
    vspace1 = vector_space(**vspace1_args)
    vspace2_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_2,
    }
    vspace2 = vector_space(**vspace2_args)
    test_trilinear_coeff_cal_values = [
        (
            {"vspace": vspace2, "cell": [25, 25, 25]},
            [
                [-5.5, 0, 2, -1, 0, 0, 0, 0],
                [-22, 3, 0, 1, 0, 0, 0, 0],
                [-5.96833648, 0.08695652, 1, 0, 0, 0, 0],
            ],
        )
    ]

    @pytest.mark.parametrize("kwargs, expected", test_trilinear_coeff_cal_values)
    def test_trilinear_coeff_cal_vals(kwargs, expected):
        r"""Test expected values."""
        assert trilinear_coeff_cal(**kwargs) == expected


def test_trilinear_jacobian():
    r"""Test `~plasmapy.analysis.nullpoint.trilinear_jacobian`."""
    vspace_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [1, 1, 1],
        "func": vspace_func_1,
    }
    vspace = vector_space(**vspace_args)

    jcb = trilinear_jacobian(vspace, [0, 0, 0])
    mtrx = jcb(0.5, 0.5, 0.5)
    exact_mtrx = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    assert np.allclose(mtrx, exact_mtrx, atol=_ATOL)


def test_trilinear_approx():
    r"""Test `~plasmapy.analysis.nullpoint.trilinear_approx`."""
    vspace1_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_1,
    }
    vspace1 = vector_space(**vspace1_args)
    vspace2_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_2,
    }
    vspace2 = vector_space(**vspace2_args)
    dx, dy, dz = vspace2[2]
    dx = dx[0]
    dy = dy[0]
    dz = dz[0]
    f000 = [0, 0, 0]
    f001 = [0, 0, dz]
    f010 = [0, dy, 0]
    f011 = [0, dy, dz]
    f100 = [dx, 0, 0]
    f101 = [dx, 0, dz]
    f110 = [dx, dy, 0]
    f111 = [dx, dy, dz]
    mid = [dx / 2, dy / 2, dz / 2]
    corners = [f000, f001, f010, f011, f100, f101, f110, f111]
    tlApprox = trilinear_approx(vspace2, [0, 0, 0])
    # Testing trilinear approx function on the corners
    for p in corners:
        approx = tlApprox(p[0], p[1], p[2])
        exact = vspace_func_2(p[0], p[1], p[2])
        approx = approx.reshape(1, 3)
        assert np.allclose(approx, exact, atol=_ATOL)
    # Testing Trilinear Approx function on a midpoint
    approx = tlApprox(mid[0], mid[1], mid[2])
    approx = approx.reshape(1, 3)
    assert np.allclose(approx, [-5.39130435, -21.5652174, 23.68667299], atol=_ATOL)


class Test_reduction:
    r"""Test `~plasmapy.analysis.nullpoint.reduction`."""
    vspace_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_2,
    }
    vspace = vector_space(**vspace_args)

    test_reduction_values = [
        ({"vspace": vspace, "cell": [25, 25, 25]}, True),
        ({"vspace": vspace, "cell": [32, 14, 4]}, True),
        ({"vspace": vspace, "cell": [0, 0, 0]}, False),
        ({"vspace": vspace, "cell": [45, 45, 45]}, False),
        ({"vspace": vspace, "cell": [33, 12, 0]}, True),
        ({"vspace": vspace, "cell": [31, 16, 8]}, True),
        ({"vspace": vspace, "cell": [24, 25, 26]}, True),
    ]

    @pytest.mark.parametrize("kwargs, expected", test_reduction_values)
    def test_reduction_vals(self, kwargs, expected):
        r"""Test expected values."""
        assert reduction(**kwargs) == expected


class Test_trilinear_analysis:
    r"""Test `~plasmapy.analysis.nullpoint.trilinear_analysis`."""
    vspace_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_2,
    }
    vspace = vector_space(**vspace_args)
    test_trilinear_analysis_values = [
        ({"vspace": vspace, "cell": [25, 25, 25]}, True),
        ({"vspace": vspace, "cell": [32, 14, 4]}, True),
        ({"vspace": vspace, "cell": [33, 12, 0]}, False),
        ({"vspace": vspace, "cell": [31, 16, 8]}, False),
        ({"vspace": vspace, "cell": [24, 25, 26]}, False),
    ]

    @pytest.mark.parametrize("kwargs, expected", test_trilinear_analysis_values)
    def test_trilinear_analysis_vals(self, kwargs, expected):
        r"""Test expected values."""
        assert trilinear_analysis(**kwargs) == expected


class Test_bilinear_root:
    r"""Test `~plasmapy.analysis.nullpoint.bilinear_root`."""
    test_bilinear_root_values = [
        (
            {"a1": 1, "b1": 3, "c1": 5, "d1": 1, "a2": 2, "b2": 4, "c2": 6, "d2": 8},
            [0.358257569496, -0.387210334997, -0.558257569496, 0.15191621735],
        )
    ]

    @pytest.mark.parametrize("kwargs, expected", test_bilinear_root_values)
    def test_bilinear_root_vals(self, kwargs, expected):
        r"""Test expected values."""
        x1, y1 = bilinear_root(**kwargs)[0]
        x2, y2 = bilinear_root(**kwargs)[1]
        assert np.isclose(x1, expected[0], atol=_ATOL)
        assert np.isclose(y1, expected[1], atol=_ATOL)
        assert np.isclose(x2, expected[2], atol=_ATOL)
        assert np.isclose(y2, expected[3], atol=_ATOL)


class Test_locate_null_point:
    r"""Test `~plasmapy.analysis.nullpoint.locate_null_point`."""
    vspace_args = {
        "x_range": [5, 6],
        "y_range": [5, 6],
        "z_range": [5, 6],
        "precision": [1, 1, 1],
        "func": vspace_func_1,
    }
    vspace = vector_space(**vspace_args)

    test_locate_null_point_values = [
        (
            {"vspace": vspace, "cell": [0, 0, 0], "n": 500, "err": _ATOL},
            np.array([5.5, 5.5, 5.5]),
        )
    ]

    @pytest.mark.parametrize("kwargs, expected", test_locate_null_point_values)
    def test_locate_null_point_vals(self, kwargs, expected):
        r"""Test expected values."""
        assert np.isclose(
            locate_null_point(**kwargs).reshape(1, 3), expected, atol=_ATOL
        ).all()


def test_null_point_find1():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # Uniform grid
    nullpoint_args = {
        "x_range": [5, 6],
        "y_range": [5, 6],
        "z_range": [5, 6],
        "precision": [1, 1, 1],
        "func": vspace_func_1,
    }
    npoints = null_point_find(**nullpoint_args)
    loc = npoints[0].loc.reshape(1, 3)
    assert len(npoints) == 1
    assert np.isclose(loc, [5.5, 5.5, 5.5], atol=_ATOL).all()


def test_null_point_find2():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # Non-uniform grid
    nullpoint2_args = {
        "x_arr": [0, 1, 2, 3, 4, 5, 6],
        "y_arr": [0, 2, 4, 6, 8],
        "z_arr": [0, 2, 4, 6, 8],
        "func": vspace_func_1,
    }
    npoints2 = null_point_find(**nullpoint2_args)
    loc2 = npoints2[0].loc.reshape(1, 3)
    assert len(npoints2) == 1
    assert np.isclose(loc2, [5.5, 5.5, 5.5], atol=_ATOL).all()
