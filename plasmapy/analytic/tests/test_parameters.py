"""Tests for functions that calculate plasma parameters."""

import numpy as np
import astropy
from astropy import units as u
from ...constants import c, m_p, m_e, e, mu0
import pytest
from ..parameters import (Alfven_speed,
                          electron_gyrofrequency, ion_gyrofrequency,
                          electron_gyroradius, ion_gyroradius,
                          electron_thermal_speed, ion_thermal_speed,
                          electron_plasma_frequency, ion_plasma_frequency,
                          Debye_length, Debye_number,
                          electron_inertial_length, ion_inertial_length,
                          ion_sound_speed,
                          magnetic_energy_density, magnetic_pressure)

B = 1.0*u.T
Z = 1
ion = 'p'
m_i = m_p
n_i = 5e19*u.m**-3
n_e = Z*5e19*u.m**-3
rho = n_i*m_i + n_e*m_e
T_e = 1e6*u.K
T_i = 1e6*u.K

B_arr = np.array([0.001, 0.002])*u.T
rho_arr = np.array([5e-10, 2e-10])*u.kg/u.m**3

B_nanarr = np.array([0.001, np.nan])*u.T
rho_infarr = np.array([np.inf, 5e19])*u.m**-3
rho_negarr = np.array([-5e19, 6e19])*u.m**-3
T_nanarr = np.array([1e6, np.nan])*u.K
T_negarr = np.array([1e6, -5151.])*u.K


def test_Alfven_speed():
    """Test the Alfven_speed function in parameters.py."""

    V_A = Alfven_speed(B, n_i)
    assert np.isclose(V_A.value, (B/np.sqrt(mu0*n_i*(m_p+m_e))).si.value)
    assert Alfven_speed(B, rho) == Alfven_speed(B, n_i)
    assert Alfven_speed(B, rho) == Alfven_speed(rho, B)
    assert Alfven_speed(B, rho).unit == 'm / s'
    assert Alfven_speed(B, rho) == Alfven_speed(-B, rho)
    assert Alfven_speed(B, 4*rho) == 0.5*Alfven_speed(B, rho)
    assert Alfven_speed(2*B, rho) == 2*Alfven_speed(B, rho)

    # Case where magnetic field and density are Quantity arrays
    V_A_arr = Alfven_speed(B_arr, rho_arr)
    V_A_arr0 = Alfven_speed(B_arr[0], rho_arr[0])
    V_A_arr1 = Alfven_speed(B_arr[1], rho_arr[1])
    assert np.isclose(V_A_arr0.value, V_A_arr[0].value)
    assert np.isclose(V_A_arr1.value, V_A_arr[1].value)

    # Case where magnetic field is an array but density is a scalar Quantity
    V_A_arr = Alfven_speed(B_arr, rho)
    V_A_arr0 = Alfven_speed(B_arr[0], rho)
    V_A_arr1 = Alfven_speed(B_arr[1], rho)
    assert np.isclose(V_A_arr0.value, V_A_arr[0].value)
    assert np.isclose(V_A_arr1.value, V_A_arr[1].value)

    with pytest.raises(ValueError):
        Alfven_speed(B_nanarr, rho_arr)

    with pytest.raises(ValueError):
        Alfven_speed(B_arr, rho_infarr)

    with pytest.raises(ValueError):
        Alfven_speed(B_arr, rho_negarr)

    with pytest.raises(u.UnitConversionError):
        Alfven_speed(5*u.A, n_i)

    with pytest.raises(TypeError):
        Alfven_speed(B, 5)

    with pytest.raises(u.UnitsError):
        Alfven_speed(B, 5*u.m**-2)

    with pytest.raises(ValueError):
        Alfven_speed(B, n_i, ion='spacecats')

    with pytest.raises(UserWarning):
        Alfven_speed(5e29*u.T, 5e19*u.m**-3)

    with pytest.raises(ValueError):
        Alfven_speed(0.001*u.T, -5e19*u.m**-3)

    with pytest.raises(ValueError):
        Alfven_speed(np.nan*u.T, 1*u.m**-3)

    with pytest.raises(ValueError):
        Alfven_speed(1*u.T, np.nan*u.m**-3)


def test_ion_sound_speed():
    """Test the ion_sound_velocity function in parameters.py."""

    assert ion_sound_speed(T_i, ion='p').unit == 'm / s'    
    assert ion_sound_speed(4*T_i) == 2*ion_sound_speed(T_i)
    assert ion_sound_speed(T_i, gamma=np.inf) == np.inf*u.m/u.s
    assert ion_sound_speed(T_i, gamma=5/3) == ion_sound_speed(T_i)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i, gamma=0.9999)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i, ion='cupcakes')

    with pytest.raises(ValueError):
        ion_sound_speed(-np.abs(T_i))

    with pytest.raises(UserWarning):
        ion_sound_speed(5e12*u.K)

    with pytest.raises(UserWarning):
        ion_sound_speed(5e19*u.K)

    with pytest.raises(u.UnitConversionError):
        ion_sound_speed(5*u.A)

    with pytest.raises(ValueError):
        ion_sound_speed(T_nanarr)

    with pytest.raises(ValueError):
        ion_sound_speed(T_negarr)


def test_electron_thermal_speed():
    """Test the electron_thermal_speed function in parameters.py."""
    assert electron_thermal_speed(T_e).unit == 'm / s'
    assert electron_thermal_speed(T_e) > ion_thermal_speed(T_e)

    with pytest.raises(u.UnitConversionError):
        electron_thermal_speed(5*u.m)

    with pytest.raises(ValueError):
        electron_thermal_speed(-T_e)

    with pytest.raises(UserWarning):
        electron_thermal_speed(1e9*u.K)

    with pytest.raises(UserWarning):
        electron_thermal_speed(5e19*u.K)


def test_ion_thermal_speed():
    """Test the ion_thermal_speed function in parameters.py"""
    assert ion_thermal_speed(T_i).unit == 'm / s'


def test_electron_gyrofrequency():
    """Test the electron_gyrofrequency function in parameters.py."""
    assert electron_gyrofrequency(B).unit == 'rad / s'


def test_ion_gyrofrequency():
    """Test the ion_gyrofrequency function in parameters.py."""
    assert ion_gyrofrequency(B, ion=ion).unit == 'rad / s'


def test_electron_gyroradius():
    """Test the electron_gyroradius function in parameters.py."""
    assert electron_gyroradius(B, T_e).unit == 'm'
    assert electron_gyroradius(B, 25*u.m/u.s).unit == 'm'


def test_ion_gyroradius():
    """Test the ion_gyroradius function in parameters.py."""
    return None


def test_electron_plasma_frequency():
    """Test the electron_plasma_frequency function in parameters.py."""
    assert electron_plasma_frequency(n_e).unit == 'rad / s'


def test_ion_plasma_frequency():
    """Test the ion_plasma_frequency function in parameters.py."""
    assert ion_plasma_frequency(n_i, Z=None, ion='p').unit == 'rad / s'
    assert ion_plasma_frequency(n_i, Z=1, ion='H-1') == \
        ion_plasma_frequency(n_i, Z=None, ion='p')


def test_Debye_length():
    """Test the Debye_length function in parameters.py."""
    assert Debye_length(T_e, n_e).unit == 'm'


def test_Debye_number():
    """Test the Debye_number function in parameters.py."""
    assert Debye_number(T_e, n_e).unit == u.dimensionless_unscaled


def test_ion_inertial_length():
    """Test the ion_inertial_length function in parameters.py."""
    assert ion_inertial_length(n_i, ion='p').unit == 'm'


def test_electron_inertial_length():
    """Test the electron_inertial_length function in parameters.py."""
    assert electron_inertial_length(n_e).unit == 'm'


def test_magnetic_pressure():
    """Test the magnetic_pressure function in parameters.py."""
    assert magnetic_pressure(B).unit == 'Pa'
    assert magnetic_pressure(B).value == magnetic_energy_density(B).value


def test_magnetic_energy_density():
    """Test the magnetic_energy_density function in parameters.py."""
    assert magnetic_energy_density(B).unit == 'J / m3'
