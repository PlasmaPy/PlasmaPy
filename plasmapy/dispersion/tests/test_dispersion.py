"""Tests for the plasma dispersion function and its derivative"""

# This file contains experimental usage of unicode characters.

import numpy as np
import pytest

from astropy import units as u
from hypothesis import given
from hypothesis.strategies import complex_numbers
from numpy import pi as π
from scipy.special import gamma as Γ

from plasmapy.dispersion.dispersionfunction import (
    plasma_dispersion_func,
    plasma_dispersion_func_deriv,
    plasma_dispersion_func_deriv_lite,
    plasma_dispersion_func_lite,
)

# Expected errors table. Used for both plasma_dispersion_func
# and plasma_dispersion_func_deriv
# w, expected_error
plasma_disp_func_errors_table = [
    ("", TypeError),
    (7 * u.m, u.UnitsError),
    (np.inf, ValueError),
    (np.nan, ValueError),
]


# Array of expected values for plasma_dispersion_func
# (w, expected)
plasma_dispersion_func_table = [
    (0, 1j * np.sqrt(π)),
    (1, -1.076_159_013_825_536_8 + 0.652_049_332_173_292_2j),
    (1j, 0.757_872_156_141_311_87j),
    (1.2 + 4.4j, -0.054_246_157_069_223_27 + 0.207_960_584_359_855_62j),
    (9.2j, plasma_dispersion_func(9.2j * u.dimensionless_unscaled)),
    (5.4 - 3.1j, -0.139_224_873_051_713_11 - 0.082_067_822_640_155_802j),
    (9.9 - 10j, 2.013_835_257_947_027_6 - 25.901_274_737_989_727j),
    (4.5 - 10j, -1.367_495_046_340_094_7e35 - 6.853_923_234_842_270_6e34j),
]


class TestPlasmaDispersionFunction:
    """
    Test class for `plasmapy.dispersion.plasma_dispersion_func`.

    Note: Testing of `plasma_dispersion_func_lite` is done in a separate
    test class.
    """

    @pytest.mark.parametrize(
        "bound_name, bound_attr",
        [("lite", plasma_dispersion_func_lite)],
    )
    def test_lite_function_binding(self, bound_name, bound_attr):
        """Test expected attributes are bound correctly."""
        assert hasattr(plasma_dispersion_func, bound_name)
        assert getattr(plasma_dispersion_func, bound_name) is bound_attr

    def test_lite_function_marking(self):
        """
        Test function is marked as having a Lite-Function.
        """
        assert hasattr(plasma_dispersion_func, "__bound_lite_func__")
        assert isinstance(plasma_dispersion_func.__bound_lite_func__, dict)

        for (
            bound_name,
            bound_origin,
        ) in plasma_dispersion_func.__bound_lite_func__.items():
            assert hasattr(plasma_dispersion_func, bound_name)

            attr = getattr(plasma_dispersion_func, bound_name)
            origin = f"{attr.__module__}.{attr.__name__}"
            assert origin == bound_origin

    @pytest.mark.parametrize("w, expected", plasma_dispersion_func_table)
    def test_plasma_dispersion_func(self, w, expected):
        r"""Test plasma_dispersion_func against tabulated results and
        known symmetry properties."""

        # Many of the tabulated results originally came from the book
        # entitled "The Plasma Dispersion Function: The Hilbert Transform
        # of the Gaussian" by B. D. Fried and S. D. Conte (1961).

        Z_of_w = plasma_dispersion_func(w)

        assert np.isclose(Z_of_w, expected, atol=1e-12 * (1 + 1j), rtol=1e-12), (
            f"plasma_dispersion_func({w}) equals {Z_of_w} instead of the "
            f"expected approximate result of {expected}.  The difference between "
            f"the actual and expected results is {Z_of_w - expected}."
        )

    @given(complex_numbers(allow_infinity=False, allow_nan=False, max_magnitude=20))
    def test_plasma_dispersion_func_symmetry(self, w):
        r"""Test plasma_dispersion_func against its symmetry properties"""

        # The two symmetry properties of the plasma dispersion function
        # are taken from the bottom of page 30 of "NRL Plasma Formulary"
        # by A.S. Richardson (2019)

        Z_of_wconj = plasma_dispersion_func(w.conjugate())
        minusZ_of_minuswconj = -(plasma_dispersion_func(-w).conjugate())

        assert np.isclose(Z_of_wconj, minusZ_of_minuswconj, atol=0, rtol=1e-15), (
            "The symmetry property of the plasma dispersion function that "
            f"Z(w*) == -[Z(-w)]* is not met for w = {w}.  Instead, "
            f"plasma_dispersion_func({w.conjugate()}) = {Z_of_wconj} "
            f"whereas -plasma_dispersion_func({-w}).conjugate() = "
            f"{minusZ_of_minuswconj}.  "
            "The difference between Z(w*) and -[Z(-w)]* is "
            f"{Z_of_wconj - minusZ_of_minuswconj}."
        )

        if w.imag > 0:
            should_equal_Z_of_wconj = (
                plasma_dispersion_func(w)
            ).conjugate() + 2j * np.sqrt(π) * np.exp(-(w.conjugate() ** 2))

            assert np.isclose(Z_of_wconj, should_equal_Z_of_wconj, rtol=1e-13), (
                "The symmetry property of the plasma dispersion function that "
                "Z(w*) = Z(w)* + 2j * sqrt(pi) * exp[-(w*)**2] for Im(w) > 0 "
                f"is not met for w = {w}.  The value of "
                f"plasma_dispersion_func({w.conjugate()}) is {Z_of_wconj}, "
                f"which is different from {should_equal_Z_of_wconj}.  "
                "The difference between these two results is "
                f"{Z_of_wconj - should_equal_Z_of_wconj}."
            )

    def test_plasma_dispersion_func_power_series_expansion(self):
        """Test plasma_dispersion_func against a power series expansion of
        the plasma dispersion function."""

        w_array = np.array(
            [
                [0.1356 + 0.114j, -0.204 - 0.0012j],
                [-0.131 + 0.131j, 0.1313 - 0.125j],
                [-0.334 - 0.712j, 0.12411 + 0j],
                [0.1278 + 0.928j, 0 + 0j],
            ],
            dtype=np.complex128,
        )

        try:
            Z_of_w_array = plasma_dispersion_func(w_array)
        except Exception as exc:
            raise ValueError(
                "plasma_dispersion_func is unable to accept an "
                f"ndarray argument with values:\n{w_array}"
            ) from exc

        # The following power series expansion is given by equation (B.3)
        # on page 401 of Plasma Waves by D. G. Swanson (2003, 2nd
        # edition).  The range of convergence of this expansion is not
        # stated, but arguments are chosen to be close to the origin.

        Z_power_series = np.zeros_like(w_array)

        for n in range(200):
            Z_power_series += 1j * np.sqrt(π) * (1j * w_array) ** n / Γ(n / 2 + 1)

        assert np.allclose(
            Z_of_w_array, Z_power_series, atol=1e-15 * (1 + 1j), rtol=1e-15
        ), (
            "The values returned by plasma_dispersion_func are inconsistent "
            "with the power series expansion of the plasma dispersion function.  "
            f"The argument given to plasma_dispersion_func is:\n\n{w_array}\n\n"
            f"The results of plasma_dispersion_func are:\n\n{Z_of_w_array}\n\n"
            "The results from the power series expansion are:\n\n"
            f"{Z_power_series}\n\n"
            "The difference between these two results is:\n\n"
            f"{Z_of_w_array - Z_power_series}\n"
        )

    def test_plasma_dispersion_func_roots(self):
        """Test roots of the plasma dispersion function."""

        # The first five roots of the plasma dispersion function are given
        # on page 402 of Swanson (2003), with some roundoff or truncation
        # error in the final decimal point.  These roots were found to
        # higher precision using mpmath.findroot.

        roots = np.array(
            [
                1.991_466_842_833_879_6 - 1.354_810_128_112_006_2j,
                2.691_149_024_251_438_8 - 2.177_044_906_089_615_9j,
                3.235_330_868_352_816_5 - 2.784_387_613_230_428_2j,
                3.697_309_702_468_468_4 - 3.287_410_789_389_848_6j,
                4.106_107_284_682_632_1 - 3.725_948_719_445_790_4j,
            ],
            dtype=np.complex128,
        )

        for root in roots:
            Z_at_root = plasma_dispersion_func(root)
            assert np.isclose(Z_at_root, 0 + 0j, atol=1e-15 * (1 + 1j)), (
                "A root of the plasma dispersion function is expected at w = "
                f"{root}, but plasma_dispersion_func({root}) is equal to "
                f"{Z_at_root} instead of {0j}."
            )

    @pytest.mark.parametrize("w, expected_error", plasma_disp_func_errors_table)
    def test_plasma_dispersion_func_errors(self, w, expected_error):
        """Test errors that should be raised by plasma_dispersion_func."""

        with pytest.raises(expected_error):
            plasma_dispersion_func(w)
            pytest.fail(
                f"plasma_dispersion_func({w}) did not raise "
                f"{expected_error.__name__} as expected."
            )


class TestPlasmaDispersionFunctionLite:
    """Test class for `plasma_dispersion_func_lite`."""

    @pytest.mark.parametrize("w, expected", plasma_dispersion_func_table)
    def test_normal_vs_lite(self, w, expected):
        r"""Test that plasma_dispersion_func and plasma_dispersion_func_lite
        calculate the same values."""

        # Many of the tabulated results originally came from the book
        # entitled "The Plasma Dispersion Function: The Hilbert Transform
        # of the Gaussian" by B. D. Fried and S. D. Conte (1961).

        Z_of_w = plasma_dispersion_func(w)
        Z_of_w_lite = plasma_dispersion_func_lite(w)

        assert np.isclose(Z_of_w, Z_of_w_lite, atol=1e-12 * (1 + 1j), rtol=1e-12), (
            f"plasma_dispersion_func({w}) and plasma_dispersion_func_lite({w}) "
            "are not equal."
        )


# Array of expected values for plasma_dispersion_func_deriv
# w, expected
plasma_disp_deriv_table = [
    (0, -2),
    (1, 0.152_318 - 1.304_10j),
    (1j, -0.484_257),
    (1.2 + 4.4j, -0.397_561e-1 - 0.217_392e-1j),
    (9j, plasma_dispersion_func_deriv(9j * u.dimensionless_unscaled)),
    (5.4 - 3.1j, 0.012_449_1 + 0.023_138_3j),
    (9.9 - 10j, 476.153 + 553.121j),
    (5 + 7j, -4.591_20e-3 - 0.012_610_4j),
    (4.5 - 10j, 0.260_153e37 - 0.211_814e37j),
]


class TestPlasmaDispersionFunctionDeriv:
    """
    Test class for `plasmapy.dispersion.plasma_dispersion_func_deriv`.

    Note: Testing of `plasma_dispersion_func_deriv_lite` is done in a separate
    test class.
    """

    @pytest.mark.parametrize(
        "bound_name, bound_attr",
        [("lite", plasma_dispersion_func_deriv_lite)],
    )
    def test_lite_function_binding(self, bound_name, bound_attr):
        """Test expected attributes are bound correctly."""
        assert hasattr(plasma_dispersion_func_deriv, bound_name)
        assert getattr(plasma_dispersion_func_deriv, bound_name) is bound_attr

    def test_lite_function_marking(self):
        """
        Test function is marked as having a Lite-Function.
        """
        assert hasattr(plasma_dispersion_func_deriv, "__bound_lite_func__")
        assert isinstance(plasma_dispersion_func_deriv.__bound_lite_func__, dict)

        for (
            bound_name,
            bound_origin,
        ) in plasma_dispersion_func_deriv.__bound_lite_func__.items():
            assert hasattr(plasma_dispersion_func_deriv, bound_name)

            attr = getattr(plasma_dispersion_func_deriv, bound_name)
            origin = f"{attr.__module__}.{attr.__name__}"
            assert origin == bound_origin

    @pytest.mark.parametrize("w, expected", plasma_disp_deriv_table)
    def test_plasma_dispersion_func_deriv(self, w, expected):
        r"""Test plasma_dispersion_func_deriv against tabulated results"""

        # The tabulated results are taken from Fried & Conte (1961)

        Z_deriv = plasma_dispersion_func_deriv(w)

        assert np.isclose(Z_deriv, expected, atol=5e-5 * (1 + 1j), rtol=5e-6), (
            f"The derivative of the plasma dispersion function does not match "
            f"the expected value for w = {w}.  The value of "
            f"plasma_dispersion_func_deriv({w}) equals {Z_deriv} whereas the "
            f"expected value is {expected}.  The difference between the actual "
            f"and expected results is {Z_deriv - expected}."
        )

    @given(complex_numbers(allow_infinity=False, allow_nan=False, max_magnitude=20))
    def test_plasma_dispersion_func_deriv_characterization(self, w):
        r"""Test plasma_dispersion_func_deriv against an exact relationship."""

        # The exact analytical relationship comes from the bottom of
        # page 3 of Fried & Conte (1961).

        Z = plasma_dispersion_func(w)
        Z_deriv = plasma_dispersion_func_deriv(w)
        Z_deriv_characterization = -2 * (1 + w * Z)

        assert np.isclose(Z_deriv, Z_deriv_characterization, rtol=1e-15), (
            f"The relationship that Z'(w) = -2 * [1 + w * Z(w)] is not "
            f"met for w = {w}, where Z'(w) = {Z_deriv} and "
            f"-2 * [1 + w * Z(w)] = {Z_deriv_characterization}."
        )

    @pytest.mark.parametrize("w, expected_error", plasma_disp_func_errors_table)
    def test_plasma_dispersion_deriv_errors(self, w, expected_error):
        """Test errors that should be raised by plasma_dispersion_func_deriv."""

        with pytest.raises(expected_error):
            plasma_dispersion_func_deriv(w)
            pytest.fail(
                f"plasma_dispersion_func_deriv({w}) did not raise "
                f"{expected_error.__name__} as expected."
            )


class TestPlasmaDispersionFunctionDerivLite:
    """Test class for `plasma_dispersion_func_deriv_lite`."""

    @pytest.mark.parametrize("w, expected", plasma_disp_deriv_table)
    def test_normal_vs_lite(self, w, expected):
        r"""Test that plasma_dispersion_func_deriv and
        plasma_dispersion_func_deriv_lite
        calculate the same values."""

        # Many of the tabulated results originally came from the book
        # entitled "The Plasma Dispersion Function: The Hilbert Transform
        # of the Gaussian" by B. D. Fried and S. D. Conte (1961).

        Z_of_w = plasma_dispersion_func_deriv(w)
        Z_of_w_lite = plasma_dispersion_func_deriv_lite(w)

        assert np.isclose(Z_of_w, Z_of_w_lite, atol=1e-12 * (1 + 1j), rtol=1e-12), (
            f"plasma_dispersion_func_deriv({w}) and "
            "plasma_dispersion_func_deriv_lite({w}) "
            "are not equal."
        )
