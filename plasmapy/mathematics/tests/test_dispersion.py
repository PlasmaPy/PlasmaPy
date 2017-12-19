"""Tests for the plasma dispersion function and its derivative"""

import numpy as np
import pytest
from astropy import units as u
from ..mathematics import plasma_dispersion_func, plasma_dispersion_func_deriv


# (zeta, expected)
plasma_dispersion_func_table = [
    (0, 1j * np.sqrt(np.pi)),
    (1, -1.076_159_01 + 0.652_049_33j),
    (1j, 0.757_872_156j),
    (1.2 + 4.4j, -0.054_246_146 + 0.207_960_589j),
    (9.2j, plasma_dispersion_func(9.2j * u.dimensionless_unscaled)),
    ]


@pytest.mark.parametrize('zeta, expected', plasma_dispersion_func_table)
def test_plasma_dispersion_func(zeta, expected):
    r"""Test the implementation of plasma_dispersion_func against
    exact results, quantities calculated by Fried & Conte (1961),
    symmetry properties, and analytic results."""

    Z0 = plasma_dispersion_func(zeta)

    assert np.isclose(Z0, expected, atol=1e-8*(1 + 1j), rtol=0), \
        (f"plasma_disperion_func({zeta}) equals {Z0} instead of the "
         f"expected approximate result of {expected}.  The difference between "
         f"the actual and expected results is {Z0 - expected}.")

    Z1 = plasma_dispersion_func(zeta.conjugate())
    Z2 = -(plasma_dispersion_func(-zeta).conjugate())

    assert np.isclose(Z1, Z2, atol=1e-8*(1 + 1j), rtol=0), \
        ("The symmetry property involving conjugates of the plasma dispersion "
         f"function and its arguments that Z(zeta*) == -[Z(-zeta)]* is not "
         f"met for zeta = {zeta}.  Instead, Z(zeta*) = {Z1} whereas "
         f"-[Z(-zeta)]* = {Z2}.  The difference between Z(zeta*) and "
         f"-[Z(-zeta)]* is {Z1 - Z2}.")

    if zeta.imag > 0:

        Z3 = plasma_dispersion_func(zeta.conjugate())

        Z4 = (plasma_dispersion_func(zeta)).conjugate() \
            + 2j * np.sqrt(np.pi) * np.exp(-(zeta.conjugate()**2))

        assert np.isclose(Z3, Z4, atol=3e-8 * (1 + 1j), rtol=0), \
            (f"A symmetry property of the plasma dispersion function "
             f"is not met for {zeta}.  The value of "
             f"plasma_dispersion_func({zeta.conjugate()}) is {Z3}, which "
             f"is not approximately equal to {Z4}.  The difference between "
             f"the two results is {Z3 - Z4}.")


# zeta, expected
plasma_disp_deriv_table = [
    (0, -2),
    (1, 0.152_318 - 1.304_10j),
    (1j, -0.484_257),
    (1.2 + 4.4j, -0.397_561e-1 - 0.217_392e-1j),
    (9j, plasma_dispersion_func_deriv(9j * u.dimensionless_unscaled)),
    ]


@pytest.mark.parametrize('zeta, expected', plasma_disp_deriv_table)
def test_plasma_dispersion_func_deriv(zeta, expected):
    r"""Test the implementation of plasma_dispersion_func_deriv against
    tabulated results from Fried & Conte (1961)."""

    Z_deriv = plasma_dispersion_func_deriv(zeta)

    assert np.isclose(Z_deriv, expected, atol=2e-6 * (1 + 1j), rtol=0), \
        (f"The derivative of the plasma dispersion function does not match "
         f"the expected value for zeta = {zeta}.  The value of "
         f"plasma_dispersion_func_deriv({zeta}) equals {Z_deriv} whereas the "
         f"expected value is {expected}.  The difference between the actual "
         f"and expected results is {Z_deriv - expected}.")


# zeta, expected_error
plasma_disp_func_errors_table = [
    ('', TypeError),
    (7 * u.m, u.UnitsError),
    (np.inf, ValueError),
    (np.nan, ValueError),
    ]


@pytest.mark.parametrize('zeta, expected_error', plasma_disp_func_errors_table)
def test_plasma_dispersion_func_errors(zeta, expected_error):
    """Test errors that should be raised by plasma_dispersion_func."""
    with pytest.raises(expected_error):
        plasma_dispersion_func(zeta)


@pytest.mark.parametrize('zeta, expected_error', plasma_disp_func_errors_table)
def test_plasma_dispersion_deriv_errors(zeta, expected_error):
    """Test errors that should be raised by plasma_dispersion_func_deriv."""
    with pytest.raises(expected_error):
        plasma_dispersion_func_deriv(expected_error)
