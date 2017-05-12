"""Tests for functions that calculate plasma parameters."""

from astropy import units as u, constants as const

from astropy.constants import c, m_p, m_e, e

from ..parameters import (Alfven_speed,
                          electron_gyrofrequency, ion_gyrofrequency,
                          electron_plasma_frequency, ion_plasma_frequency,
                          Debye_length, Debye_number,
                          electron_inertial_length, ion_inertial_length,
                          ion_sound_speed, magnetic_pressure)

B = 1.0*u.T
Z = 1
ion = 'p'
m_i = const.m_p
n_i = 5e19*u.m**-3
n_e = Z*5e19*u.m**-3
rho = n_i*m_i + n_e*m_e
T_e = 1e6*u.K
T_i = 1e6*u.K

def test_Alfven_speed():
    """Test the Alfven_speed function in parameters.py."""
    assert Alfven_speed(B, rho) == Alfven_speed(B, n_i)
    assert Alfven_speed(B, rho).unit == 'm / s'


def test_ion_sound_speed():
    """Test the ion_sound_velocity function in parameters.py."""
    assert ion_sound_speed(T_i, ion='p').unit == 'm / s'


def test_electron_thermal_speed():
    """Test the electron_thermal_speed function in parameters.py."""
    return None


def test_ion_thermal_speed():
    """Test the ion_thermal_speed function in parameters.py"""
    return None

    
def test_electron_gyrofrequency():
    """Test the electron_gyrofrequency function in parameters.py."""
    assert electron_gyrofrequency(B).unit == 'rad / s'


def test_ion_gyrofrequency():
    """Test the ion_gyrofrequency function in parameters.py."""
    assert ion_gyrofrequency(B, ion).unit == 'rad / s'


def test_electron_gyroradius():
    """Test the electron_gyroradius function in parameters.py."""
    return None


def test_ion_gyroradius():
    """Test the ion_gyroradius function in parameters.py."""
    return None


def test_electron_plasma_frequency():
    """Test the electron_plasma_frequency function in parameters.py."""
    assert electron_plasma_frequency(n_e).unit == 'rad / s'


def test_ion_plasma_frequency():
    """Test the ion_plasma_frequency function in parameters.py."""
    assert ion_plasma_frequency(n_i, Z, ion).unit == 'rad / s'


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


def test_magnetic_energy_density():
    """Test the magnetic_energy_density function in parameters.py."""
    assert magnetic_pressure(B).unit == 'J / m3'
