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

mu = m_p.to(u.u).value

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
    assert Alfven_speed(5*u.T, 5e19*u.m**-3, ion='H-1') == \
        Alfven_speed(5*u.T, 5e19*u.m**-3, ion='p')

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
        Alfven_speed(B_arr, rho_negarr)

    with pytest.raises(u.UnitConversionError):
        Alfven_speed(5*u.A, n_i, ion='p')

    with pytest.raises(TypeError):
        Alfven_speed(B, 5, ion='p')

    with pytest.raises(u.UnitsError):
        Alfven_speed(B, 5*u.m**-2, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(B, n_i, ion='spacecats')

    with pytest.raises(UserWarning):  # relativistic
        Alfven_speed(5e1*u.T, 5e19*u.m**-3, ion='p')

    with pytest.raises(UserWarning):  # superrelativistic
        Alfven_speed(5e8*u.T, 5e19*u.m**-3, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(0.001*u.T, -5e19*u.m**-3, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(np.nan*u.T, 1*u.m**-3, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(1*u.T, np.nan*u.m**-3, ion='p')

    with pytest.raises(UserWarning):
        assert Alfven_speed(np.inf*u.T, 1*u.m**-3, ion='p') == np.inf*u.m/u.s

    with pytest.raises(UserWarning):
        assert Alfven_speed(-np.inf*u.T, 1*u.m**-3, ion='p') == np.inf*u.m/u.s


def test_ion_sound_speed():
    """Test the ion_sound_speed function in parameters.py."""

    assert ion_sound_speed(T_i, ion='p').unit == 'm / s'
    assert ion_sound_speed(4*T_i) == 2*ion_sound_speed(T_i)
    assert ion_sound_speed(T_i, gamma_i=3) == ion_sound_speed(T_i)

    with pytest.raises(UserWarning):
        assert ion_sound_speed(T_i, gamma_i=np.inf) == np.inf*u.m/u.s

    # Add numerical tests!

    # Add an array test!

    with pytest.raises(ValueError):
        ion_sound_speed(T_i, gamma_i=0.9999)

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

    # Add a numerical test!

    # Add an array test!

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

    # Add a numerical test!

    # Add an array test!

    # Add tests of exceptions!


def test_electron_gyrofrequency():
    """Test the electron_gyrofrequency function in parameters.py."""
    assert electron_gyrofrequency(B).unit == 'rad / s'

    # Add a numerical test!

    # Add an array test!

    # Add tests of exceptions!


def test_ion_gyrofrequency():
    """Test the ion_gyrofrequency function in parameters.py."""
    assert ion_gyrofrequency(B, ion=ion).unit == 'rad / s'

    # Add a numerical test!

    # Add an array test!

    # Add tests of exceptions!


def test_electron_gyroradius():
    """Test the electron_gyroradius function in parameters.py."""
    assert electron_gyroradius(B, T_e).unit == 'm'
    assert electron_gyroradius(B, 25*u.m/u.s).unit == 'm'

    # Add a numerical test!

    # Add an array test!

    # Add tests of exceptions!


def test_ion_gyroradius():
    """Test the ion_gyroradius function in parameters.py."""
    assert ion_gyroradius(B, T_i).unit == 'm'
    assert ion_gyroradius(B, 25*u.m/u.s).unit == 'm'

    # Add a numerical test!

    # Add an array test!

    # Add tests of exceptions!


def test_electron_plasma_frequency():
    """Test the electron_plasma_frequency function in parameters.py."""
    assert electron_plasma_frequency(n_e).unit == 'rad / s'

    assert np.isclose(electron_plasma_frequency(1*u.cm**-3).value,
                      5.64e4, rtol=1e-2)

    # Add an array test!

    # Add tests of exceptions!


def test_ion_plasma_frequency():
    """Test the ion_plasma_frequency function in parameters.py."""
    assert ion_plasma_frequency(n_i, ion='p').unit == 'rad / s'
    assert (ion_plasma_frequency(n_i, ion='H-1 1+') ==
            ion_plasma_frequency(n_i, ion='p'))
    assert np.isclose(ion_plasma_frequency(mu*u.cm**-3, ion='p').value,
                      1.32e3, rtol=1e-2)

    # Add an array test!

    # Add tests of exceptions!


def test_Debye_length():
    """Test the Debye_length function in parameters.py."""

    assert Debye_length(T_e, n_e).unit == 'm'
    assert np.isclose(Debye_length(1*u.eV, 1*u.cm**-3).value, 7.43, atol=0.005)

    with pytest.raises(TypeError):
        Debye_length(5, 5*u.m**-3)

    with pytest.raises(u.UnitConversionError):
        Debye_length(56*u.kg, 5*u.m**-3)

    with pytest.raises(ValueError):
        Debye_length(5*u.eV, -5*u.m**-3)

    with pytest.raises(ValueError):
        Debye_length(-45*u.K, 5*u.m**-3)


def test_Debye_number():
    """Test the Debye_number function in parameters.py."""

    assert Debye_number(T_e, n_e).unit == u.dimensionless_unscaled

    T_e_eV = T_e.to(u.eV, equivalencies=u.temperature_energy())
    assert np.isclose(Debye_number(T_e, n_e).value,
                      Debye_number(T_e_eV, n_e).value)

    assert np.isclose(Debye_number(1*u.eV, 1*u.cm**-3).value, 1720862385.43342)

    with pytest.raises(TypeError):
        Debye_number(T_e, 4)

    with pytest.raises(TypeError):
        Debye_number(None, n_e)

    with pytest.raises(u.UnitConversionError):
        Debye_number(5*u.m, 5*u.m**-3)

    with pytest.raises(u.UnitConversionError):
        Debye_number(5*u.K, 5*u.m**3)

    with pytest.raises(ValueError):
        Debye_number(5j*u.K, 5*u.cm**-3)


def test_ion_inertial_length():
    """Test the ion_inertial_length function in parameters.py."""
    assert ion_inertial_length(n_i, ion='p').unit == 'm'

    # Add a numerical test!
    
    assert np.isclose(ion_inertial_length(mu*u.cm**-3, ion='p').cgs.value,
                      2.28e7, rtol=0.01)

    # Add an array test!

    with pytest.raises(TypeError):
        ion_inertial_length(4)

    with pytest.raises(u.UnitConversionError):
        ion_inertial_length(4*u.m**-2)

    with pytest.raises(ValueError):
        ion_inertial_length(-5*u.m**-3)


def test_electron_inertial_length():
    """Test the electron_inertial_length function in parameters.py."""
    assert electron_inertial_length(n_e).unit == 'm'

    # Add a numerical test!

    # Add an array test!

    with pytest.raises(TypeError):
        electron_inertial_length(5)

    with pytest.raises(u.UnitConversionError):
        electron_inertial_length(5*u.m)

    with pytest.raises(ValueError):
        electron_inertial_length(-5*u.m**-3)


def test_magnetic_pressure():
    """Test the magnetic_pressure function in parameters.py."""
    assert magnetic_pressure(B_arr).unit == 'J / m3'
    assert magnetic_pressure(B).unit == 'Pa'
    assert magnetic_pressure(B).value == magnetic_energy_density(B).value
    assert magnetic_pressure(B) == magnetic_energy_density(B.to(u.G))
    assert np.isclose(magnetic_pressure(B).value, 397887.35772973835)

    with pytest.raises(TypeError):
        magnetic_pressure(5)

    with pytest.raises(u.UnitConversionError):
        magnetic_pressure(5*u.m)

    with pytest.raises(ValueError):
        magnetic_pressure(np.nan*u.T)

    with pytest.raises(ValueError):
        magnetic_pressure(5j*u.T)

    with pytest.raises(ValueError):
        magnetic_pressure(B_nanarr)


def test_magnetic_energy_density():
    """Test the magnetic_energy_density function in parameters.py."""

    assert magnetic_energy_density(B_arr).unit == 'J / m3'
    assert magnetic_energy_density(B).unit == 'J / m3'
    assert magnetic_energy_density(B).value == magnetic_pressure(B).value
    assert magnetic_energy_density(2*B) == 4*magnetic_energy_density(B)
    assert np.isclose(magnetic_energy_density(B).value, 397887.35772973835)
    assert magnetic_energy_density(B) == magnetic_energy_density(B.to(u.G))

    # Add an array test!

    assert magnetic_energy_density(B_arr)

    with pytest.raises(TypeError):
        magnetic_energy_density(5)

    with pytest.raises(u.UnitConversionError):
        magnetic_energy_density(5*u.m)

    with pytest.raises(ValueError):
        magnetic_energy_density(np.nan*u.T)

    with pytest.raises(ValueError):
        magnetic_energy_density(5j*u.T)

    with pytest.raises(ValueError):
        magnetic_energy_density(B_nanarr)
