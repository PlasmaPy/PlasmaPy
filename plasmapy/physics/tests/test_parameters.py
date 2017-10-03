"""Tests for functions that calculate plasma parameters."""

import numpy as np
import pytest
from astropy import units as u


from ...constants import c, m_p, m_e, e, mu0

from ..parameters import (Alfven_speed,
                          electron_gyrofrequency,
                          ion_gyrofrequency,
                          electron_gyroradius,
                          ion_gyroradius,
                          thermal_speed,
                          electron_plasma_frequency,
                          ion_plasma_frequency,
                          Debye_length,
                          Debye_number,
                          electron_inertial_length,
                          ion_inertial_length,
                          ion_sound_speed,
                          magnetic_energy_density,
                          magnetic_pressure,
                          upper_hybrid_frequency,
                          lower_hybrid_frequency)

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

V = 25.2*u.m/u.s

# Assertions below that are in CGS units with 2-3 significant digits
# are generally from the NRL Plasma Formulary.


def test_Alfven_speed():
    """Test the Alfven_speed function in parameters.py."""

    assert np.isclose(Alfven_speed(1*u.T, 1e-8*u.kg*u.m**-3).value,
                      8920620.580763856)

    V_A = Alfven_speed(B, n_i)
    assert np.isclose(V_A.value, (B/np.sqrt(mu0*n_i*(m_p+m_e))).si.value)

    assert Alfven_speed(B, rho) == Alfven_speed(B, n_i)

    assert Alfven_speed(B, rho) == Alfven_speed(rho, B)

    assert Alfven_speed(B, rho).unit == u.m/u.s

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
        Alfven_speed(np.array([5, 6, 7])*u.T,
                     np.array([5, 6])*u.m**-3)

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

    with pytest.raises(UserWarning):  # super-relativistic
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

    with pytest.raises(UserWarning):
        assert Alfven_speed(1.0, n_i) == Alfven_speed(1.0*u.T, n_i)


def test_ion_sound_speed():
    """Test the ion_sound_speed function in parameters.py."""

    assert np.isclose(ion_sound_speed(T_i=1.3232*u.MK, T_e=1.831*u.MK,
                                      ion='p', gamma_e=1, gamma_i=3).value,
                      218816.06086407552)

    assert np.isclose(ion_sound_speed(T_i=0.88*u.MK, T_e=1.28*u.MK, ion='p',
                                      gamma_e=1.2, gamma_i=3.4).value,
                      193328.52857788358)

    assert ion_sound_speed(T_i=T_i, T_e=T_e, ion='p') == \
        ion_sound_speed(T_i=T_i, T_e=T_e, ion='H-1')

    assert ion_sound_speed(T_i=T_i, ion='p').unit == u.m/u.s

    assert ion_sound_speed(T_i=T_i, gamma_i=3) == ion_sound_speed(T_i=T_i)

    assert ion_sound_speed(T_e=T_e, gamma_e=1) == ion_sound_speed(T_e=T_e)

    with pytest.raises(UserWarning):
        assert ion_sound_speed(T_i=T_i, gamma_i=np.inf) == np.inf*u.m/u.s

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=np.array([5, 6, 5])*u.K, T_e=np.array([3, 4])*u.K)

    with pytest.raises(TypeError):    # Is this test right??????
        ion_sound_speed(5*u.T)

    with pytest.raises(TypeError):
        ion_sound_speed('p')

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_i, gamma_i=0.9999)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_i, gamma_e=0.9999)

    with pytest.raises(TypeError):
        ion_sound_speed(T_i=T_i, gamma_e='sdjklsf')

    with pytest.raises(TypeError):
        ion_sound_speed(T_i=T_i, gamma_i='fsdfas')

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_i, ion='cupcakes')

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=-np.abs(T_i))

    with pytest.raises(UserWarning):
        ion_sound_speed(T_i=5e12*u.K)

    with pytest.raises(UserWarning):
        ion_sound_speed(T_i=5e19*u.K)

    with pytest.raises(u.UnitConversionError):
        ion_sound_speed(T_i=5*u.A)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_nanarr)

    with pytest.raises(ValueError):
        ion_sound_speed(T_e=T_nanarr)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_negarr)

    with pytest.raises(ValueError):
        ion_sound_speed(T_e=T_negarr)

    with pytest.raises(UserWarning):
        assert ion_sound_speed(T_e=1.2e6) == ion_sound_speed(T_e=1.2e6*u.K)

    with pytest.raises(UserWarning):
        assert ion_sound_speed(T_i=1.3e6) == ion_sound_speed(T_i=1.3e6*u.K)



def test_thermal_speed():
    """Test the thermal_speed function in parameters.py"""
    assert thermal_speed(T_e).unit == u.m/u.s

    assert thermal_speed(T_e) > thermal_speed(T_e, 'p')

    # The NRL Plasma Formulary uses a definition of the electron
    # thermal speed that differs by a factor of sqrt(2).
    assert np.isclose(thermal_speed(1*u.MK).value,
                      5505694.743141063)

    with pytest.raises(u.UnitConversionError):
        thermal_speed(5*u.m)

    with pytest.raises(ValueError):
        thermal_speed(-T_e)

    with pytest.raises(UserWarning):
        thermal_speed(1e9*u.K)

    with pytest.raises(UserWarning):
        thermal_speed(5e19*u.K)

    with pytest.raises(UserWarning):
        assert thermal_speed(1e5) == thermal_speed(1e5*u.K)

    assert thermal_speed(T_i, ion='p').unit == u.m/u.s

    # The NRL Plasma Formulary uses a definition of the ion thermal
    # speed that differs by a factor of sqrt(2).
    assert np.isclose(thermal_speed(1*u.MK, ion='p').si.value,
                      128486.56960876315)

    assert thermal_speed(T_i, ion='p') == \
        thermal_speed(T_i, ion='H-1')

    assert thermal_speed(1*u.MK, ion='e+') == \
        thermal_speed(1*u.MK)

    with pytest.raises(u.UnitConversionError):
        thermal_speed(5*u.m, ion='p')

    with pytest.raises(ValueError):
        thermal_speed(-T_e, ion='p')

    with pytest.raises(UserWarning):
        thermal_speed(1e11*u.K, ion='p')

    with pytest.raises(UserWarning):
        thermal_speed(1e14*u.K, ion='p')

    with pytest.raises(ValueError):
        thermal_speed(T_i, ion='asdfasd')

    with pytest.raises(UserWarning):
        assert thermal_speed(1e6, ion='p') ==\
            thermal_speed(1e6*u.K, ion='p')


def test_electron_gyrofrequency():
    """Test the electron_gyrofrequency function in parameters.py."""

    assert electron_gyrofrequency(B).unit == u.rad/u.s

    assert np.isclose(electron_gyrofrequency(1*u.T).value, 175882008784.72018)

    assert np.isclose(electron_gyrofrequency(2.4*u.T).value,
                      422116821083.3284)

    assert np.isclose(electron_gyrofrequency(1*u.G).cgs.value,
                      1.76e7, rtol=1e-3)

    with pytest.raises(TypeError):
        electron_gyrofrequency(u.m)

    with pytest.raises(u.UnitConversionError):
        electron_gyrofrequency(u.m*1)

    with pytest.raises(ValueError):
        electron_gyrofrequency(B_nanarr)

    # The following is a test to check that equivalencies from astropy
    # are working.
    omega_ce = electron_gyrofrequency(2.2*u.T)
    f_ce = (omega_ce/(2*np.pi))/u.rad
    f_ce_use_equiv = omega_ce.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])
    assert np.isclose(f_ce.value, f_ce_use_equiv.value)

    with pytest.raises(UserWarning):
        assert electron_gyrofrequency(5.0) == electron_gyrofrequency(5.0*u.T)


def test_ion_gyrofrequency():
    """Test the ion_gyrofrequency function in parameters.py."""

    assert ion_gyrofrequency(B, ion=ion).unit == u.rad/u.s

    assert np.isclose(ion_gyrofrequency(1*u.T, ion='p').value,
                      95788335.834874)

    assert np.isclose(ion_gyrofrequency(2.4*u.T, ion='p').value,
                      229892006.00369796)

    assert np.isclose(ion_gyrofrequency(1*u.G, ion='p').cgs.value,
                      9.58e3, rtol=2e-3)

    assert ion_gyrofrequency(-5*u.T) == ion_gyrofrequency(5*u.T)

    assert ion_gyrofrequency(B, ion='p') == \
        ion_gyrofrequency(B, ion='H-1')

    assert ion_gyrofrequency(B, ion='e+') == \
        electron_gyrofrequency(B)

    with pytest.raises(UserWarning):
        ion_gyrofrequency(8)

    with pytest.raises(u.UnitConversionError):
        ion_gyrofrequency(5*u.m)

    with pytest.raises(ValueError):
        ion_gyrofrequency(8*u.T, ion='asdfasd')

    with pytest.raises(UserWarning):
        assert ion_gyrofrequency(5.0) == ion_gyrofrequency(5.0*u.T)


def test_electron_gyroradius():
    """Test the electron_gyroradius function in parameters.py."""

    assert electron_gyroradius(B, T_e).unit == u.m

    assert electron_gyroradius(B, 25*u.m/u.s).unit == u.m

    assert electron_gyroradius(T_e, B) == electron_gyroradius(B, T_e)

    assert electron_gyroradius(V, B) == electron_gyroradius(B, V)

    assert electron_gyroradius(B, V) == electron_gyroradius(B, -V)

    Vperp = 1e6*u.m/u.s
    Bmag = 1*u.T
    omega_ce = electron_gyrofrequency(Bmag)
    assert electron_gyroradius(Bmag, Vperp) == \
        (Vperp/omega_ce).to(u.m, equivalencies=u.dimensionless_angles())

    with pytest.raises(TypeError):
        electron_gyroradius(u.T, 8*u.m/u.s)

    with pytest.raises(u.UnitConversionError):
        electron_gyroradius(5*u.A, 8*u.m/u.s)

    with pytest.raises(u.UnitConversionError):
        electron_gyroradius(5*u.T, 8*u.m)

    with pytest.raises(ValueError):
        electron_gyroradius(np.array([5, 6])*u.T, np.array([5, 6, 7])*u.m/u.s)

    with pytest.raises(ValueError):
        electron_gyroradius(np.nan*u.T, 1*u.m/u.s)

    with pytest.raises(ValueError):
        electron_gyroradius(3.14159*u.T, -1*u.K)

    with pytest.raises(UserWarning):
        assert electron_gyroradius(1.0, Vperp=1.0) == \
            electron_gyroradius(1.0*u.T, Vperp=1.0*u.m/u.s)

    with pytest.raises(UserWarning):
        assert electron_gyroradius(1.1, T_e=1.2) == \
            electron_gyroradius(1.1*u.T, T_e=1.2*u.K)

    with pytest.raises(ValueError):
        electron_gyroradius(1.1*u.T, T_e=1.2*u.K, Vperp=1*u.m/u.s)

    with pytest.raises(ValueError):
        electron_gyroradius(1.1*u.T, 1.2*u.K, 1.1*u.m)


def test_ion_gyroradius():
    """Test the ion_gyroradius function in parameters.py."""

    assert ion_gyroradius(B, T_i).unit == u.m

    assert ion_gyroradius(B, 25*u.m/u.s).unit == u.m

    assert ion_gyroradius(B, T_i, ion='p') == \
        ion_gyroradius(B, T_i, ion='H-1')

    assert ion_gyroradius(T_i, B) == ion_gyroradius(B, T_i)

    assert ion_gyroradius(V, B) == ion_gyroradius(B, V)

    assert ion_gyroradius(B, V) == ion_gyroradius(B, -V)

    Vperp = 1e6*u.m/u.s
    Bmag = 1*u.T
    omega_ci = ion_gyrofrequency(Bmag, ion='p')
    assert ion_gyroradius(Bmag, Vperp) == \
        (Vperp/omega_ci).to(u.m, equivalencies=u.dimensionless_angles())

    T2 = 1.2*u.MK
    B2 = 123*u.G
    ion2 = 'alpha'
    Vperp2 = thermal_speed(T2, ion=ion2)
    assert ion_gyroradius(B2, Vperp=Vperp2, ion='alpha') == \
        ion_gyroradius(B2, T_i=T2, ion='alpha')

    assert ion_gyroradius(1*u.T, 1*u.MK, ion='positron') == \
        electron_gyroradius(1*u.T, 1*u.MK)

    with pytest.raises(TypeError):
        ion_gyroradius(u.T, 8*u.m/u.s)

    with pytest.raises(ValueError):
        ion_gyroradius(B, T_i, ion='asfdas')

    with pytest.raises(ValueError):
        ion_gyroradius(B, -1*u.K, ion='p')

    with pytest.raises(UserWarning):
        assert ion_gyroradius(1.0, Vperp=1.0) == \
            ion_gyroradius(1.0*u.T, Vperp=1.0*u.m/u.s)

    with pytest.raises(UserWarning):
        assert ion_gyroradius(1.1, T_i=1.2) == \
            ion_gyroradius(1.1*u.T, T_i=1.2*u.K)

    with pytest.raises(ValueError):
        ion_gyroradius(1.1*u.T, T_i=1.2*u.K, Vperp=1*u.m/u.s)

    with pytest.raises(ValueError):
        ion_gyroradius(1.1*u.T, 1.2*u.K, 1.1*u.m)

    with pytest.raises(ValueError):
        ion_gyroradius(1.1*u.T, 1.2*u.m, 1.1*u.K)


def test_electron_plasma_frequency():
    """Test the electron_plasma_frequency function in parameters.py."""

    assert electron_plasma_frequency(n_e).unit == u.rad/u.s

    assert np.isclose(electron_plasma_frequency(1*u.cm**-3).value,
                      5.64e4, rtol=1e-2)

    with pytest.raises(TypeError):
        electron_plasma_frequency(u.m**-3)

    with pytest.raises(u.UnitConversionError):
        electron_plasma_frequency(5*u.m**-2)

    with pytest.raises(ValueError):
        electron_plasma_frequency(np.nan*u.m**-3)

    with pytest.raises(UserWarning):
        assert electron_plasma_frequency(1e19) == \
            electron_plasma_frequency(1e19*u.m**-3)


def test_ion_plasma_frequency():
    """Test the ion_plasma_frequency function in parameters.py."""

    assert ion_plasma_frequency(n_i, ion='p').unit == u.rad/u.s

    assert ion_plasma_frequency(n_i, ion='H-1') == \
        ion_plasma_frequency(n_i, ion='p')

    assert np.isclose(ion_plasma_frequency(mu*u.cm**-3, ion='p').value,
                      1.32e3, rtol=1e-2)

    assert ion_plasma_frequency(n_i, ion='e+') == \
        electron_plasma_frequency(n_i)

    with pytest.raises(ValueError):
        ion_plasma_frequency(n_i=5*u.m**-3, ion='sdfas')

    with pytest.raises(UserWarning):
        assert ion_plasma_frequency(1e19) == ion_plasma_frequency(1e19*u.m**-3)


def test_Debye_length():
    """Test the Debye_length function in parameters.py."""

    assert Debye_length(T_e, n_e).unit == u.m

    assert np.isclose(Debye_length(1*u.eV, 1*u.cm**-3).value, 7.43, atol=0.005)

    with pytest.raises(UserWarning):
        Debye_length(5, 5*u.m**-3)

    with pytest.raises(u.UnitConversionError):
        Debye_length(56*u.kg, 5*u.m**-3)

    with pytest.raises(ValueError):
        Debye_length(5*u.eV, -5*u.m**-3)

    with pytest.raises(ValueError):
        Debye_length(-45*u.K, 5*u.m**-3)

    Tarr2 = np.array([1, 2])*u.K
    narr3 = np.array([1, 2, 3])*u.m**-3
    with pytest.raises(ValueError):
        Debye_length(Tarr2, narr3)

    with pytest.raises(UserWarning):
        assert Debye_length(2.0, 2.0) == Debye_length(2.0*u.K, 2.0*u.m**-3)

    with pytest.raises(UserWarning):
        assert Debye_length(2.0*u.K, 2.0) == Debye_length(2.0, 2.0*u.m**-3)


def test_Debye_number():
    """Test the Debye_number function in parameters.py."""

    assert Debye_number(T_e, n_e).unit == u.dimensionless_unscaled

    T_e_eV = T_e.to(u.eV, equivalencies=u.temperature_energy())
    assert np.isclose(Debye_number(T_e, n_e).value,
                      Debye_number(T_e_eV, n_e).value)

    assert np.isclose(Debye_number(1*u.eV, 1*u.cm**-3).value, 1720862385.43342)

    with pytest.raises(UserWarning):
        Debye_number(T_e, 4)

    with pytest.raises(TypeError):
        Debye_number(None, n_e)

    with pytest.raises(u.UnitConversionError):
        Debye_number(5*u.m, 5*u.m**-3)

    with pytest.raises(u.UnitConversionError):
        Debye_number(5*u.K, 5*u.m**3)

    with pytest.raises(ValueError):
        Debye_number(5j*u.K, 5*u.cm**-3)

    Tarr2 = np.array([1, 2])*u.K
    narr3 = np.array([1, 2, 3])*u.m**-3
    with pytest.raises(ValueError):
        Debye_number(Tarr2, narr3)

    with pytest.raises(UserWarning):
        assert Debye_number(1.1, 1.1) == Debye_number(1.1*u.K, 1.1*u.m**-3)

    with pytest.raises(UserWarning):
        assert Debye_number(1.1*u.K, 1.1) == Debye_number(1.1, 1.1*u.m**-3)


def test_ion_inertial_length():
    """Test the ion_inertial_length function in parameters.py."""

    assert ion_inertial_length(n_i, ion='p').unit == u.m

    assert np.isclose(ion_inertial_length(mu*u.cm**-3, ion='p').cgs.value,
                      2.28e7, rtol=0.01)

    assert ion_inertial_length(5.351*u.m**-3, ion='e+') == \
        electron_inertial_length(5.351*u.m**-3)

    assert ion_inertial_length(n_i, ion='p') == \
        ion_inertial_length(n_i, ion='H-1')

    with pytest.raises(UserWarning):
        ion_inertial_length(4)

    with pytest.raises(u.UnitConversionError):
        ion_inertial_length(4*u.m**-2)

    with pytest.raises(ValueError):
        ion_inertial_length(-5*u.m**-3)

    with pytest.raises(ValueError):
        ion_inertial_length(n_i, ion=-135)

    with pytest.raises(UserWarning):
        assert ion_inertial_length(1e19) == ion_inertial_length(1e19*u.m**-3)


def test_electron_inertial_length():
    """Test the electron_inertial_length function in parameters.py."""

    assert electron_inertial_length(n_e).unit == u.m

    assert np.isclose(electron_inertial_length(1*u.cm**-3).cgs.value,
                      5.31e5, rtol=1e-3)

    with pytest.raises(UserWarning):
        electron_inertial_length(5)

    with pytest.raises(u.UnitConversionError):
        electron_inertial_length(5*u.m)

    with pytest.raises(ValueError):
        electron_inertial_length(-5*u.m**-3)

    with pytest.raises(UserWarning):
        assert electron_inertial_length(1e19) == \
            electron_inertial_length(1e19*u.m**-3)


def test_magnetic_pressure():
    """Test the magnetic_pressure function in parameters.py."""

    assert magnetic_pressure(B_arr).unit == u.Pa

    assert magnetic_pressure(B).unit == u.Pa

    assert magnetic_pressure(B).unit.name == 'Pa'

    assert magnetic_pressure(B).value == magnetic_energy_density(B).value

    assert magnetic_pressure(B) == magnetic_energy_density(B.to(u.G))

    assert np.isclose(magnetic_pressure(B).value, 397887.35772973835)

    with pytest.raises(UserWarning):
        magnetic_pressure(5)

    with pytest.raises(u.UnitConversionError):
        magnetic_pressure(5*u.m)

    with pytest.raises(ValueError):
        magnetic_pressure(np.nan*u.T)

    with pytest.raises(ValueError):
        magnetic_pressure(5j*u.T)

    with pytest.raises(ValueError):
        magnetic_pressure(B_nanarr)

    with pytest.raises(UserWarning):
        assert magnetic_pressure(22.2) == magnetic_pressure(22.2*u.T)


def test_magnetic_energy_density():
    """Test the magnetic_energy_density function in parameters.py."""

    assert magnetic_energy_density(B_arr).unit == u.J/u.m**3

    assert str(magnetic_energy_density(B).unit) == 'J / m3'

    assert magnetic_energy_density(B).value == magnetic_pressure(B).value

    assert magnetic_energy_density(2*B) == 4*magnetic_energy_density(B)

    assert np.isclose(magnetic_energy_density(B).value, 397887.35772973835)

    assert magnetic_energy_density(B) == magnetic_energy_density(B.to(u.G))

    # Add an array test!

    assert magnetic_energy_density(B_arr)

    with pytest.raises(UserWarning):
        magnetic_energy_density(5)

    with pytest.raises(u.UnitConversionError):
        magnetic_energy_density(5*u.m)

    with pytest.raises(ValueError):
        magnetic_energy_density(np.nan*u.T)

    with pytest.raises(ValueError):
        magnetic_energy_density(5j*u.T)

    with pytest.raises(ValueError):
        magnetic_energy_density(B_nanarr)

    with pytest.raises(UserWarning):
        assert magnetic_energy_density(22.2) == \
            magnetic_energy_density(22.2*u.T)


def test_upper_hybrid_frequency():
    """Test the upper_hybrid_frequency function in parameters.py."""

    omega_uh = upper_hybrid_frequency(B, n_e=n_e)
    omega_ce = electron_gyrofrequency(B)
    omega_pe = electron_plasma_frequency(n_e=n_e)
    assert omega_ce.unit == u.rad/u.s
    assert omega_pe.unit == u.rad/u.s
    assert omega_uh.unit == u.rad/u.s
    LHS = omega_uh**2
    RHS = omega_ce**2 + omega_pe**2
    assert np.isclose(LHS.value, RHS.value)

    with pytest.raises(ValueError):
        upper_hybrid_frequency(5*u.T, n_e=-1*u.m**-3)

    with pytest.raises(UserWarning):
        assert upper_hybrid_frequency(1.2, 1.3) == \
            upper_hybrid_frequency(1.2*u.T, 1.3*u.m**-3)

    with pytest.raises(UserWarning):
        assert upper_hybrid_frequency(1.4*u.T, 1.3) == \
            upper_hybrid_frequency(1.4, 1.3*u.m**-3)


def test_lower_hybrid_frequency():
    """Test the lower_hybrid_frequency function in parameters.py."""

    ion = 'He-4 1+'
    omega_ci = ion_gyrofrequency(B, ion=ion)
    omega_pi = ion_plasma_frequency(n_i=n_i, ion=ion)
    omega_ce = electron_gyrofrequency(B)
    omega_lh = lower_hybrid_frequency(B, n_i=n_i, ion=ion)
    assert omega_ci.unit == u.rad/u.s
    assert omega_pi.unit == u.rad/u.s
    assert omega_ce.unit == u.rad/u.s
    assert omega_lh.unit == u.rad/u.s
    LHS = omega_lh**-2
    RHS = 1/(omega_ci**2 + omega_pi**2) + omega_ci**-1*omega_ce**-1
    assert np.isclose(LHS.value, RHS.value)

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2*u.T, n_i=5e19*u.m**-3, ion='asdfasd')

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2*u.T, n_i=-5e19*u.m**-3, ion='asdfasd')

    with pytest.raises(ValueError):
        lower_hybrid_frequency(np.nan*u.T, n_i=-5e19*u.m**-3, ion='asdfasd')

    with pytest.raises(UserWarning):
        assert lower_hybrid_frequency(1.3, 1e19) == \
            lower_hybrid_frequency(1.3*u.T, 1e19*u.m**-3)
