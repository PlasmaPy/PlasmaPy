"""Tests for the plasma dispersion function and its derivative"""

import numpy as np
import pytest
from astropy import units as u
from ..mathematics import plasma_dispersion_func, plasma_dispersion_func_deriv


# (argument, expected)
plasma_dispersion_func_table = [
    (0, 1j * np.sqrt(np.pi)),
    (1, -1.076_159_01 + 0.652_049_33j),
    (1j, 0.757_872_156j),
    (1.2 + 4.4j, -0.054_246_146 + 0.207_960_589j),
    (9.2j, plasma_dispersion_func(9.2j * u.dimensionless_unscaled))]


@pytest.mark.parametrize('argument, expected', plasma_dispersion_func_table)
def test_plasma_dispersion_func(argument, expected):

    atol = 1e-8*(1 + 1j)

    Z = plasma_dispersion_func(argument)

    assert np.isclose(Z, expected, atol=atol, rtol=0), \
        (f"plasma_disperion_func({argument}) equals {Z} instead of the "
         f"expected approximate result of {expected}")

    if argument.imag > 0:

        atol_conj = 4e-8 * (1 + 1j)

        Z1 = plasma_dispersion_func(argument.conjugate())

        Z2 = (plasma_dispersion_func(argument)).conjugate() \
            + 2j * np.sqrt(np.pi) * np.exp(-(argument.conjugate()**2))

        assert np.isclose(Z1, Z2, atol=atol_conj, rtol=0), \
            (f"The first symmetry property of the plasma dispersion function"
             f"is not met for an argument of {argument}.  The value of "
             f"plasma_dispersion_func({argument.conjugate()}) is {Z1}, which "
             f"is not approximately equal to {Z2}")


def test_plasma_dispersion_func_old():
    r"""Test the implementation of plasma_dispersion_func against exact
    results, quantities calculated by Fried & Conte (1961), and
    analytic results.
    """

    atol = 1e-8*(1 + 1j)

    assert np.isclose(plasma_dispersion_func(0), 1j*np.sqrt(np.pi),
                      atol=atol, rtol=0), \
        "Z(0) does not give accurate answer"

    assert np.isclose(plasma_dispersion_func(1), -1.07615901 + 0.65204933j,
                      atol=atol, rtol=0), \
        "Z(1) not consistent with tabulated results"

    assert np.isclose(plasma_dispersion_func(1j), 0.757872156j,
                      atol=atol, rtol=0), \
        "Z(1j) not consistent with tabulated results"

    assert np.isclose(plasma_dispersion_func(1.2 + 4.4j),
                      -0.054246146 + 0.207960589j, atol=atol, rtol=0), \
        "Z(1.2+4.4j) not consistent with tabulated results"

    zeta = -1.2 + 0.4j

    assert np.isclose(plasma_dispersion_func(zeta.conjugate()),
                      -(plasma_dispersion_func(-zeta).conjugate()),
                      atol=atol, rtol=0), \
        "Symmetry property of Z(zeta*) = -[Z(-zeta)]* not satisfied"

    if zeta.imag > 0:
        assert np.isclose(plasma_dispersion_func(zeta.conjugate()),
                          (plasma_dispersion_func(zeta)).conjugate() +
                          2j*np.sqrt(np.pi)*np.exp(-(zeta.conjugate()**2)),
                          atol=atol, rtol=0), \
            "Symmetry property of Z(zeta*) valid for zeta.imag>0 not satisfied"

    assert np.isclose(plasma_dispersion_func(9.2j),
                      plasma_dispersion_func(9.2j*u.dimensionless_unscaled),
                      atol=0, rtol=0)

    with pytest.raises(TypeError):
        plasma_dispersion_func('')

    with pytest.raises(u.UnitsError):
        plasma_dispersion_func(6*u.m)

    with pytest.raises(ValueError):
        plasma_dispersion_func(np.inf)


def test_plasma_dispersion_func_deriv():
    r"""Test the implementation of plasma_dispersion_func_deriv against
    tabulated results from Fried & Conte (1961)."""

    atol = 2e-6*(1 + 1j)

    assert np.isclose(plasma_dispersion_func_deriv(0), -2, atol=atol, rtol=0),\
        "Z'(0) not consistent with tabulated values"

    assert np.isclose(plasma_dispersion_func_deriv(1), 0.152318 - 1.30410j,
                      atol=atol, rtol=0), \
        "Z'(1) not consistent with tabulated values"

    assert np.isclose(plasma_dispersion_func_deriv(1j),
                      -0.484257, atol=atol, rtol=0), \
        "Z'(1j) not consistent with tabulated values"

    assert np.isclose(plasma_dispersion_func_deriv(1.2 + 4.4j),
                      -0.397561e-1 - 0.217392e-1j, atol=atol, rtol=0)

    assert np.isclose(plasma_dispersion_func_deriv(9j),
                      plasma_dispersion_func_deriv(9j*u.m/u.m),
                      atol=atol, rtol=0)

    with pytest.raises(TypeError):
        plasma_dispersion_func_deriv('')

    with pytest.raises(u.UnitsError):
        plasma_dispersion_func_deriv(6*u.m)

    with pytest.raises(ValueError):
        plasma_dispersion_func_deriv(np.inf)
