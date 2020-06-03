# -*- coding: utf-8 -*-
"""
Tests for Thomson scattering analysis functions
"""

import numpy as np
from astropy import units as u
import pytest
from astropy import units as u
from plasmapy.diagnostics import thomson


def width_at_value(x, y, val):
    """
    Calculates the width of a curve at a given value.
    """
    above = np.where(y > val, x, np.nan)
    return np.abs(np.nanmax(above) - np.nanmin(above))


def gen_collective_spectrum():
    """
    Generates an example Thomson scattering spectrum in the collective regime
    """
    wavelength = np.arange(520, 545, 0.01)*u.nm
    probe_wavelength = 532*u.nm
    ne = 5e17*u.cm**-3
    probe_n = np.array([1,0,0])
    scatter_n = np.array([0,1,0])
    fract = np.array([1.0])
    Te = 10*u.eV
    Ti = np.array([10])*u.eV
    ion_z = np.array([5])
    ion_mu =  np.array([12])
    fluid_vel = np.array([0,0,0])*u.km/u.s
    ion_vel = np.array([[0,0,0]])*u.km/u.s

    alpha, Skw = thomson.spectral_density(wavelength, probe_wavelength=probe_wavelength,
                     ne=ne,fract = fract, Te=Te, Ti=Ti,
                     ion_z = ion_z, ion_mu=ion_mu, ion_vel=ion_vel, fluid_vel=fluid_vel,
                     probe_n=probe_n, scatter_n=scatter_n)

    return alpha, wavelength, Skw



def gen_non_collective_spectrum():
    """
    Generates an example Thomson scattering spectrum in the non-collective
    regime
    """
    wavelength = np.arange(500, 570, 0.01)*u.nm
    probe_wavelength = 532*u.nm
    ne = 5e15*u.cm**-3
    probe_n = np.array([1,0,0])
    scatter_n = np.array([0,1,0])
    fract = np.array([1.0])
    Te = 100*u.eV
    Ti = np.array([10])*u.eV
    ion_z = np.array([1])
    ion_mu =  np.array([1])
    fluid_vel = np.array([0,0,0])*u.km/u.s
    ion_vel = np.array([[0,0,0]])*u.km/u.s

    alpha, Skw = thomson.spectral_density(wavelength, probe_wavelength=probe_wavelength,
                     ne=ne,fract = fract, Te=Te, Ti=Ti,
                     ion_z = ion_z, ion_mu=ion_mu, ion_vel=ion_vel, fluid_vel=fluid_vel,
                     probe_n=probe_n, scatter_n=scatter_n)

    return alpha, wavelength, Skw


def test_collective_spectrum():
    """
    Compares the generated spectrum to previously determined values
    """
    alpha, wavelength, Skw = gen_collective_spectrum()

    #Check that alpha is correct
    assert np.isclose(alpha, 1.27, atol=.01), f"Collective case alpha returns {alpha} instead of expected 1.27"

    i_width = width_at_value(wavelength.value, Skw.value, 2e-13)
    e_width = width_at_value(wavelength.value, Skw.value, .2e-13)

    #Check that the widths of the ion and electron features match expectations
    assert np.isclose(i_width, 0.14, 1e-3), f"Collective case ion feature width is {i_width} instead of expected 0.14"
    assert np.isclose(e_width, 16.32, 1e-3), f"Collective case electron feature width is {e_width} instead of expected 16.32"


def test_non_collective_spectrum():
    """
    Compares the generated spectrum to previously determined values
    """
    alpha, wavelength, Skw = gen_non_collective_spectrum()

    #Check that alpha is correct
    assert np.isclose(alpha, 0.0403, atol=.01), f"Non-collective case alpha returns {alpha} instead of expected 0.0403"

    e_width = width_at_value(wavelength.value, Skw.value, .2e-13)

    #Check that the widts of the electron feature matchs expectations
    assert np.isclose(e_width, 28.68, 1e-3), f"Non-collective case electron feature width is {e_width} instead of expected 28.68"
