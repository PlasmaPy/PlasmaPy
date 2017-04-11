"""Tests for the plasma dispersion function and its derivative"""

import numpy as np

from ..dispersion import plasma_dispersion_func, \
    plasma_dispersion_func_deriv

def test_plasma_dispersion_func():
    """Test the implementation of plasma_dispersion_func against exact
    results, quantities calculated by Fried & Conte (1961), and
    analytic results.
    """

    atol = 1e-10*(1+1j)

    assert np.isclose(plasma_dispersion_func(0), 1j*np.sqrt(np.pi), \
                          atol=atol, rtol=0), \
        "Z(0) does not give accurate answer"

    assert np.isclose(plasma_dispersion_func(1), 
                      -1.076_159_014 + 0.652_049_332j, 
                      atol=atol, rtol=0), \
                      "Z(1) not consistent with tabulated results"

    assert np.isclose(plasma_dispersion_func(1j), 
                      0.757_872_156j,
                      atol=atol, rtol=0), \
                      "Z(1j) not consistent with tabulated results"

    assert np.isclose(plasma_dispersion_func(1.2+4.4j),
                      -0.054_246_146 + 0.207_960_589j,
                      atol=atol, rtol=0), \
                      "Z(1.2+4.4j) not consistent with tabulated results"

    zeta = -1.2 + 0.4j                      

    assert np.isclose(plasma_dispersion_func(zeta.conjugate()), \
                          -(plasma_dispersion_func(-zeta).conjugate()), \
                          atol=atol, rtol=0), \
        "Symmetry property of Z(zeta*) = -[Z(-zeta)]* not satisfied"

    if zeta.imag>0: 
        assert np.isclose(plasma_dispersion_func(zeta.conjugate()), \
            (plasma_dispersion_func(zeta)).conjugate() + \
            2j*np.sqrt(np.pi)*np.exp(-(zeta.conjugate()**2)), \
                              atol=atol, rtol=0), \
            "Symmetry property of Z(zeta*) valid for zeta.imag>0 not satisfied"

def test_plasma_dispersion_func_deriv():
    """Test the implementation of plasma_dispersion_func_deriv against
    tabulated results from Fried & Conte (1961)."""

    atol = 1e-6*(1+1j)

    assert np.isclose(plasma_dispersion_func_deriv(0),
                      -2.0, atol=atol, rtol=0), \
        "Z'(0) not consistent with tabulated values"

    assert np.isclose(plasma_dispersion_func_deriv(1),
                      0.152_318 - 0.130_410e1*1j, atol=atol, rtol=0), \
        "Z'(1) not consistent with tabulated values"

    assert np.isclose(plasma_dispersion_func_deriv(1j),
                      -0.484_257, atol=atol, rtol=0), \
        "Z'(1j) not consistent with tabulated values"

    assert np.isclose(plasma_dispersion_func_deriv(1.2+4.4j),
                      -0.397_561e-1 - 0.217_392e-1*1j, atol=atol, rtol=0)

    


    
