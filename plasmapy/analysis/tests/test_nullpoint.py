"""
Tests for the null point finder class defined in `plasmapy.analysis.nullpoint`.
"""
import numpy as np
import pytest

from plasmapy.analysis.nullpoint import (
    _bilinear_root,
    _EQUALITY_ATOL,
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

# Defining tolerance level for tests where the accuracy
# level can not rise too high due to run time issues.
_TESTING_ATOL = 1e-2


def vspace_func_1(x, y, z):
    return [(y - 5.5), (z - 5.5), (x - 5.5)]


def vspace_func_2(x, y, z):
    return [2 * y - z - 5.5, 3 * x + z - 22, x**2 - 11 * x + y + 24.75]


def vspace_func_3(x, y, z):
    return [(y - 5.3) * (y - 5.5), (z - 5.5), (x - 5.5)]


def vspace_func_4(x, y, z):
    return [y - y, (z - 5.5), (x - 5.5)]


def vspace_func_5(x, y, z):
    return [(y - y), (z - z), (x - x)]


def vspace_func_6(x, y, z):
    return [(-1 - x - z), (-x - z - x * z), (y)]


def vspace_func_7(x, y, z):
    return [(y - 5.5) * (y + 5.5), (z - 5.5), (x - 5.5)]


def test_trilinear_coeff_cal():
    r"""Test `~plasmapy.analysis.nullpoint.trilinear_coeff_cal`."""
    vspace1_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_1,
    }
    vspace1 = _vector_space(**vspace1_args)
    vspace2_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_2,
    }
    vspace2 = _vector_space(**vspace2_args)
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
        assert _trilinear_coeff_cal(**kwargs) == expected


def test_trilinear_jacobian():
    r"""Test `~plasmapy.analysis.nullpoint.trilinear_jacobian`."""
    vspace_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [1, 1, 1],
        "func": vspace_func_1,
    }
    vspace = _vector_space(**vspace_args)

    jcb = _trilinear_jacobian(vspace, [0, 0, 0])
    mtrx = jcb(0.5, 0.5, 0.5)
    exact_mtrx = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    assert np.allclose(mtrx, exact_mtrx, atol=_EQUALITY_ATOL)


def test_trilinear_approx():
    r"""Test `~plasmapy.analysis.nullpoint.trilinear_approx`."""
    vspace2_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_2,
    }
    vspace2 = _vector_space(**vspace2_args)
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
        assert np.allclose(approx, exact, atol=_EQUALITY_ATOL)
    # Testing Trilinear Approx function on a midpoint
    approx = tlApprox(mid[0], mid[1], mid[2])
    approx = approx.reshape(1, 3)
    assert np.allclose(
        approx, [-5.39130435, -21.5652174, 23.68667299], atol=_EQUALITY_ATOL
    )


class Test_reduction:
    r"""Test `~plasmapy.analysis.nullpoint.reduction`."""
    vspace_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_2,
    }
    vspace = _vector_space(**vspace_args)

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
        assert _reduction(**kwargs) == expected


class Test_trilinear_analysis:
    r"""Test `~plasmapy.analysis.nullpoint.trilinear_analysis`."""
    vspace_args = {
        "x_range": [0, 10],
        "y_range": [0, 10],
        "z_range": [0, 10],
        "precision": [10 / 46, 10 / 46, 10 / 46],
        "func": vspace_func_2,
    }
    vspace = _vector_space(**vspace_args)
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
        assert _trilinear_analysis(**kwargs) == expected


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
        x1, y1 = _bilinear_root(**kwargs)[0]
        x2, y2 = _bilinear_root(**kwargs)[1]
        assert np.isclose(x1, expected[0], atol=_EQUALITY_ATOL)
        assert np.isclose(y1, expected[1], atol=_EQUALITY_ATOL)
        assert np.isclose(x2, expected[2], atol=_EQUALITY_ATOL)
        assert np.isclose(y2, expected[3], atol=_EQUALITY_ATOL)


class Test_locate_null_point:
    r"""Test `~plasmapy.analysis.nullpoint.locate_null_point`."""
    vspace_args = {
        "x_range": [5, 6],
        "y_range": [5, 6],
        "z_range": [5, 6],
        "precision": [1, 1, 1],
        "func": vspace_func_1,
    }
    vspace = _vector_space(**vspace_args)

    test_locate_null_point_values = [
        (
            {"vspace": vspace, "cell": [0, 0, 0], "n": 500, "err": _EQUALITY_ATOL},
            np.array([5.5, 5.5, 5.5]),
        )
    ]

    @pytest.mark.parametrize("kwargs, expected", test_locate_null_point_values)
    def test_locate_null_point_vals(self, kwargs, expected):
        r"""Test expected values."""
        assert np.isclose(
            _locate_null_point(**kwargs).reshape(1, 3), expected, atol=_EQUALITY_ATOL
        ).all()


@pytest.mark.slow
def test_null_point_find1():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # Uniform grid
    nullpoint_args = {
        "x_range": [5, 6],
        "y_range": [5, 6],
        "z_range": [5, 6],
        "precision": [0.1, 0.1, 0.1],
        "func": vspace_func_1,
    }
    npoints = uniform_null_point_find(**nullpoint_args)
    loc = npoints[0].loc.reshape(1, 3)
    assert len(npoints) == 1
    assert np.isclose(loc, [5.5, 5.5, 5.5], atol=_EQUALITY_ATOL).all()


@pytest.mark.slow
def test_null_point_find2():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # Non-uniform grid
    vspace_args = {
        "x_arr": np.logspace(np.log10(5.48), np.log10(5.52), num=30),
        "y_arr": np.logspace(np.log10(5.48), np.log10(5.52), num=30),
        "z_arr": np.logspace(np.log10(5.48), np.log10(5.52), num=30),
        "func": vspace_func_1,
    }
    vspace = _vector_space(**vspace_args)
    npoints2 = _vspace_iterator(vspace)
    loc2 = npoints2[0].loc.reshape(1, 3)
    assert len(npoints2) == 1
    assert np.isclose(loc2, [5.5, 5.5, 5.5], atol=_EQUALITY_ATOL).all()


def test_null_point_find3():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # Vector values passed by hand
    nullpoint3_args = {
        "x_arr": [5, 6],
        "y_arr": [5, 6],
        "z_arr": [5, 6],
        "u_arr": np.array([[[-0.5, -0.5], [0.5, 0.5]], [[-0.5, -0.5], [0.5, 0.5]]]),
        "v_arr": np.array([[[-0.5, 0.5], [-0.5, 0.5]], [[-0.5, 0.5], [-0.5, 0.5]]]),
        "w_arr": np.array([[[-0.5, -0.5], [-0.5, -0.5]], [[0.5, 0.5], [0.5, 0.5]]]),
    }
    npoints3 = null_point_find(**nullpoint3_args)
    loc3 = npoints3[0].loc.reshape(1, 3)
    assert len(npoints3) == 1
    assert np.isclose(loc3, [5.5, 5.5, 5.5], atol=_EQUALITY_ATOL).all()


@pytest.mark.slow
def test_null_point_find4():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # Two null points
    nullpoint4_args = {
        "x_range": [5, 6],
        "y_range": [5, 6],
        "z_range": [5, 6],
        "precision": [0.07, 0.003, 0.07],
        "func": vspace_func_3,
    }
    npoints4 = uniform_null_point_find(**nullpoint4_args)
    first_loc4 = npoints4[0].loc.reshape(1, 3)
    second_loc4 = npoints4[1].loc.reshape(1, 3)
    assert len(npoints4) == 2
    assert np.isclose(first_loc4, [5.5, 5.3, 5.5], atol=_EQUALITY_ATOL).all()
    assert np.isclose(second_loc4, [5.5, 5.5, 5.5], atol=_EQUALITY_ATOL).all()


@pytest.mark.slow
def test_null_point_find5():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # Many null points because a y vector dimension is zero
    nullpoint5_args = {
        "x_range": [5.4, 5.6],
        "y_range": [5.4, 5.6],
        "z_range": [5.4, 5.6],
        "precision": [0.01, 0.01, 0.01],
        "func": vspace_func_4,
    }
    npoints5 = uniform_null_point_find(**nullpoint5_args)
    for p in npoints5:
        if np.allclose(p.loc, np.array([5.5, 5.5, 5.5]), _EQUALITY_ATOL):
            assert (
                p.classification
                == "Anti-parallel lines with null plane OR Planes of parabolae with null line"
            )
        if np.allclose(p.loc, np.array([5.47, 5.5, 5.5]), _EQUALITY_ATOL):
            assert (
                p.classification
                == "Anti-parallel lines with null plane OR Planes of parabolae with null line"
            )
        if np.allclose(p.loc, np.array([5.5, 5.5, 5.47]), _EQUALITY_ATOL):
            assert (
                p.classification
                == "Anti-parallel lines with null plane OR Planes of parabolae with null line"
            )


@pytest.mark.slow
def test_null_point_find6():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # Many null points; All vector dimensions zero
    nullpoint6_args = {
        "x_range": [5, 6],
        "y_range": [5, 6],
        "z_range": [5, 6],
        "precision": [0.08, 0.08, 0.08],
        "func": vspace_func_5,
    }
    npoints6 = uniform_null_point_find(**nullpoint6_args)
    assert len(npoints6) == 0


@pytest.mark.slow
def test_null_point_find7():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # No null points, discriminant less than zero
    nullpoint7_args = {
        "x_range": [-10, 10],
        "y_range": [-10, 10],
        "z_range": [-10, 10],
        "precision": [1, 1, 1],
        "func": vspace_func_6,
    }
    npoints7 = uniform_null_point_find(**nullpoint7_args)
    assert len(npoints7) == 0


@pytest.mark.slow
def test_null_point_find8():
    r"""Test `~plasmapy.analysis.nullpoint.null_point_find`."""
    # Non-linear field
    nullpoint8_args = {
        "x_range": [5, 6],
        "y_range": [-6, 6],
        "z_range": [5, 6],
        "precision": [0.3, 0.3, 0.3],
        "func": vspace_func_7,
    }
    npoints8 = uniform_null_point_find(**nullpoint8_args)
    assert len(npoints8) == 2
    loc1 = npoints8[0].loc.reshape(1, 3)
    loc2 = npoints8[1].loc.reshape(1, 3)
    assert np.allclose(loc1, [5.5, -5.5, 5.5], atol=_TESTING_ATOL)
    assert np.allclose(loc2, [5.5, 5.5, 5.5], atol=_TESTING_ATOL)


@pytest.mark.slow
class Test_classify_null_point:
    r"""Test `~plasmapy.analysis.nullpoint._classify_null_point`."""

    test_classify_null_point_values = [
        (
            {
                "x_range": [-0.1, 0.1],
                "y_range": [-0.1, 0.1],
                "z_range": [-0.1, 0.1],
                "precision": [0.03, 0.03, 0.03],
                "func": lambda x, y, z: [x, 2 * y, -3 * z],
            },
            "Improper radial null",
        ),
        (
            {
                "x_range": [-0.1, 0.1],
                "y_range": [-0.1, 0.1],
                "z_range": [-0.1, 0.1],
                "precision": [0.03, 0.03, 0.03],
                "func": lambda x, y, z: [
                    -1 * x + 2 * y - 4 * z,
                    2 * x + 2 * y + 2 * z,
                    -4 * x + 2 * y - 1 * z,
                ],
            },
            "Proper radial null",
        ),
        (
            {
                "x_range": [5, 6],
                "y_range": [-6, 6],
                "z_range": [5, 6],
                "precision": [0.3, 0.3, 0.3],
                "func": lambda x, y, z: [(y - 5.5) * (y + 5.5), (z - 5.5), (x - 5.5)],
            },
            "Spiral null",
        ),
        (
            {
                "x_range": [-0.1, 0.1],
                "y_range": [-0.1, 0.1],
                "z_range": [-0.1, 0.1],
                "precision": [0.03, 0.03, 0.03],
                "func": lambda x, y, z: [
                    0.5 * x - 2 * y + z,
                    x - y + z,
                    x + y + 0.5 * z,
                ],
            },
            "Critical spiral null",
        ),
        (
            {
                "x_range": [-0.1, 0.1],
                "y_range": [-0.1, 0.1],
                "z_range": [-0.1, 0.1],
                "precision": [0.03, 0.03, 0.03],
                "func": lambda x, y, z: [0.5 * x - y + z, x - y + z, x + y + 0.5 * z],
            },
            "Skewed improper null",
        ),
    ]

    @pytest.mark.parametrize("kwargs, expected", test_classify_null_point_values)
    def test_classify_null_point_vals(self, kwargs, expected):
        r"""Test expected values."""
        assert uniform_null_point_find(**kwargs)[0].classification == expected


def test_null_point_find9():
    """Testing a magnetic field that violates the divergence constraint"""
    nullpoint9_args = {
        "x_range": [-0.1, 0.1],
        "y_range": [-0.1, 0.1],
        "z_range": [-0.1, 0.1],
        "precision": [0.03, 0.03, 0.03],
        "func": lambda x, y, z: [x, y, z],
    }
    with pytest.raises(NonZeroDivergence):
        npoints = uniform_null_point_find(**nullpoint9_args)


# Tests that capture the degenerate nulls/2D nulls
@pytest.mark.slow
def test_null_point_find10():
    nullpoint10_args = {
        "x_range": [-0.1, 0.1],
        "y_range": [-0.1, 0.1],
        "z_range": [-0.1, 0.1],
        "precision": [0.01, 0.01, 0.01],
        "func": lambda x, y, z: [y * z, -x * z, x * y],
    }
    npoints = uniform_null_point_find(**nullpoint10_args)
    for p in npoints:
        if np.allclose(p.loc, np.array([0, 0, 0]), _EQUALITY_ATOL):
            assert p.classification == "Proper radial null"
        if np.allclose(p.loc, np.array([0.01, 0, 0]), _EQUALITY_ATOL):
            assert (
                p.classification
                == "Anti-parallel lines with null plane OR Planes of parabolae with null line"
            )
        if np.allclose(p.loc, np.array([0, 0, -0.01]), _EQUALITY_ATOL):
            assert p.classification == "Continuous concentric ellipses"


@pytest.mark.slow
def test_null_point_find11():
    nullpoint10_args = {
        "x_range": [-0.1, 0.1],
        "y_range": [-0.1, 0.1],
        "z_range": [-0.1, 0.1],
        "precision": [0.01, 0.01, 0.01],
        "func": lambda x, y, z: [1.01 * y * z, -x * z, x * y],
    }
    npoints = uniform_null_point_find(**nullpoint10_args)
    for p in npoints:
        if np.allclose(p.loc, np.array([0, 0, 0]), _EQUALITY_ATOL):
            assert p.classification == "Proper radial null"
        if np.allclose(p.loc, np.array([0.01, 0, 0]), _EQUALITY_ATOL):
            assert (
                p.classification
                == "Anti-parallel lines with null plane OR Planes of parabolae with null line"
            )
        if np.allclose(p.loc, np.array([0, 0, -0.01]), _EQUALITY_ATOL):
            assert p.classification == "Continuous concentric ellipses"
