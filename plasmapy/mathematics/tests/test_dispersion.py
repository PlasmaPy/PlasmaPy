"""Tests for the plasma dispersion function and its derivative"""

import numpy as np
from numpy import pi as π
import pytest
from astropy import units as u
from ..mathematics import plasma_dispersion_func, plasma_dispersion_func_deriv


# (ζ, expected)
plasma_dispersion_func_table = [
    (0, 1j * np.sqrt(π)),
    (1, -1.076_159_013_825_536_8 + 0.652_049_332_173_292_2j),
    (1j, 0.757_872_156_141_311_87j),
    (1.2 + 4.4j, -0.054_246_146_372_377_471 + 0.207_960_589_336_958_13j),
    (9.2j, plasma_dispersion_func(9.2j * u.dimensionless_unscaled)),
    (5.4 - 3.1j, -0.139_224_873_051_713_11 - 0.082_067_822_640_155_802j),
    (9.9 - 10j, 2.013_835_257_947_027_6 - 25.901_274_737_989_727j),
    (4.5 - 10j, -1.367_495_046_340_094_7e35 - 6.853_923_234_842_270_6e34j)]


@pytest.mark.parametrize('ζ, expected', plasma_dispersion_func_table)
def test_plasma_dispersion_func(ζ, expected):
    r"""Test the implementation of plasma_dispersion_func against
    exact results, quantities calculated by Fried & Conte (1961),
    symmetry properties, and analytic results."""

    Z0 = plasma_dispersion_func(ζ)

    assert np.isclose(Z0, expected, atol=1e-12*(1 + 1j), rtol=1e-12), \
        (f"plasma_dispersion_func({ζ}) equals {Z0} instead of the "
         f"expected approximate result of {expected}.  The difference between "
         f"the actual and expected results is {Z0 - expected}.")

    Z1 = plasma_dispersion_func(ζ.conjugate())
    Z2 = -(plasma_dispersion_func(-ζ).conjugate())

    assert np.isclose(Z1, Z2, atol=0, rtol=1e-15), \
        ("The symmetry property involving conjugates of the plasma dispersion "
         f"function and its arguments that Z(ζ*) == -[Z(-ζ)]* is not "
         f"met for ζ = {ζ}.  Instead, Z(ζ*) = {Z1} whereas "
         f"-[Z(-ζ)]* = {Z2}.  The difference between Z(ζ*) and "
         f"-[Z(-ζ)]* is {Z1 - Z2}.")

    if ζ.imag > 0:

        Z3 = plasma_dispersion_func(ζ.conjugate())

        Z4 = (plasma_dispersion_func(ζ)).conjugate() \
            + 2j * np.sqrt(π) * np.exp(-(ζ.conjugate()**2))

        assert np.isclose(Z3, Z4, atol=0, rtol=1e-15), \
            (f"A symmetry property of the plasma dispersion function "
             f"is not met for {ζ}.  The value of "
             f"plasma_dispersion_func({ζ.conjugate()}) is {Z3}, which "
             f"is not approximately equal to {Z4}.  The difference between "
             f"the two results is {Z3 - Z4}.")


def test_plasma_dispersion_func_power_series_expansion():
    """Test plasma_dispersion_func against a power series expansion of
    the plasma dispersion function."""

    ζ_array = np.array(
        [[0.1356 + 0.114j, -0.204 - 0.0012j],
         [-0.131 + 0.131j, 0.1313 - 0.125j],
         [-0.334 - 0.712j, 0.12411 + 0j],
         [0.1278 + 0.928j, 0 + 0j]],
        dtype=np.complex128)

    Z1 = plasma_dispersion_func(ζ_array)

    Z2 = np.zeros_like(ζ_array)

    for n in range(0, 200):
        Z2 += (1j * ζ_array)**n / np.math.gamma(n/2 + 1)

    Z2 = 1j * np.sqrt(π) * Z2

    assert np.allclose(Z1, Z2, atol=1e-10 * (1 + 1j), rtol=1e-10), \
        ("plasma_dispersion_func is returning values that are inconsistent "
         "with the power series expansion given by equation B.3 from Plasma "
         "Waves by D. G. Swanson (2nd edition, 2003).  The results from "
         f"plasma_dispersion_func are:\n\n{Z1}\n\n"
         f"whereas the results from the power series expansion are:\n\n{Z2}\n")


def test_plasma_dispersion_func_zeros():
    """Test the zeros of the plasma dispersion function."""
    zeros = np.array([
            1.991_466_842_833_879_6 - 1.354_810_128_112_006_2j,
            2.691_149_024_251_438_8 - 2.177_044_906_089_615_9j,
            3.235_330_868_352_816_5 - 2.784_387_613_230_428_2j,
            3.697_309_702_468_468_4 - 3.287_410_789_389_848_6j,
            4.106_107_284_682_632_1 - 3.725_948_719_445_790_4j],
            dtype=np.complex128)

    for zero in zeros:
        Z = plasma_dispersion_func(zero)
        assert np.isclose(Z, 0 + 0j, atol=1e-15 * (1 + 1j)), \
            ("A zero of the plasma dispersion function is expected at ζ = "
             f"{zero}, but plasma_dispersion_func({zero}) is equal to {Z}.")


# ζ, expected
plasma_disp_deriv_table = [
    (0, -2),
    (1, 0.152_318 - 1.304_10j),
    (1j, -0.484_257),
    (1.2 + 4.4j, -0.397_561e-1 - 0.217_392e-1j),
    (9j, plasma_dispersion_func_deriv(9j * u.dimensionless_unscaled)),
    (5.4 - 3.1j, 0.012_449_1 + 0.023_138_3j),
    (9.9 - 10j, 476.153 + 553.121j),
    (5 + 7j, -4.59120e-3 - 0.012_610_4j),
    (4.5 - 10j, 0.260_153e37 - 0.211_814e37j),
    ]


@pytest.mark.parametrize('ζ, expected', plasma_disp_deriv_table)
def test_plasma_dispersion_func_deriv(ζ, expected):
    r"""Test the implementation of plasma_dispersion_func_deriv
    against tabulated results and an analytical relationship from
    Fried & Conte (1961)."""

    Z_deriv = plasma_dispersion_func_deriv(ζ)

    assert np.isclose(Z_deriv, expected, atol=5e-5*(1+1j), rtol=5e-6), \
        (f"The derivative of the plasma dispersion function does not match "
         f"the expected value for ζ = {ζ}.  The value of "
         f"plasma_dispersion_func_deriv({ζ}) equals {Z_deriv} whereas the "
         f"expected value is {expected}.  The difference between the actual "
         f"and expected results is {Z_deriv - expected}.")

    Z = plasma_dispersion_func(ζ)
    Z_deriv_characterization = -2 * (1 + ζ*Z)

    assert np.isclose(Z_deriv, Z_deriv_characterization, rtol=1e-15), \
        (f"The relationship that Z'(ζ) = -2 * [1 + ζ * Z(ζ)] is not "
         f"met for ζ = {ζ}, where Z'(ζ) = {Z_deriv} and "
         f"-2 * [1 + ζ * Z(ζ)] = {Z_deriv_characterization}.")


# ζ, expected_error
plasma_disp_func_errors_table = [
    ('', TypeError),
    (7 * u.m, u.UnitsError),
    (np.inf, ValueError),
    (np.nan, ValueError),
    ]


@pytest.mark.parametrize('ζ, expected_error', plasma_disp_func_errors_table)
def test_plasma_dispersion_func_errors(ζ, expected_error):
    """Test errors that should be raised by plasma_dispersion_func."""

    err_msg = (f"plasma_dispersion_func({ζ}) did not raise "
               f"{expected_error.__name__}")

    with pytest.raises(expected_error, message=err_msg):
        plasma_dispersion_func(ζ)


@pytest.mark.parametrize('ζ, expected_error', plasma_disp_func_errors_table)
def test_plasma_dispersion_deriv_errors(ζ, expected_error):
    """Test errors that should be raised by plasma_dispersion_func_deriv."""

    err_msg = (f"plasma_dispersion_func_deriv({ζ}) did not raise "
               f"{expected_error.__name__}")

    with pytest.raises(expected_error, message=err_msg):
        plasma_dispersion_func_deriv(ζ)
