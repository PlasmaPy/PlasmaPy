"""Tests for functions that calculate plasma parameters."""

import numpy as np
import pytest
from astropy import units as u
from warnings import simplefilter

from ...utils.exceptions import RelativityWarning, RelativityError
from ...utils.exceptions import PhysicsError
from ...constants import c, m_p, m_e, e, mu0

from ..parameters import (Alfven_speed,
                          gyrofrequency,
                          gyroradius,
                          thermal_speed,
                          kappa_thermal_speed,
                          plasma_frequency,
                          Debye_length,
                          Debye_number,
                          inertial_length,
                          ion_sound_speed,
                          magnetic_energy_density,
                          magnetic_pressure,
                          upper_hybrid_frequency,
                          lower_hybrid_frequency)

B = 1.0 * u.T
Z = 1
ion = 'p'
m_i = m_p
n_i = 5e19 * u.m**-3
n_e = Z * 5e19 * u.m**-3
rho = n_i * m_i + n_e * m_e
T_e = 1e6 * u.K
T_i = 1e6 * u.K

B_arr = np.array([0.001, 0.002]) * u.T
rho_arr = np.array([5e-10, 2e-10]) * u.kg / u.m**3

B_nanarr = np.array([0.001, np.nan]) * u.T
rho_infarr = np.array([np.inf, 5e19]) * u.m**-3
rho_negarr = np.array([-5e19, 6e19]) * u.m**-3
T_nanarr = np.array([1e6, np.nan]) * u.K
T_negarr = np.array([1e6, -5151.]) * u.K

mu = m_p.to(u.u).value

V = 25.2 * u.m / u.s

# Assertions below that are in CGS units with 2-3 significant digits
# are generally from the NRL Plasma Formulary.


def test_Alfven_speed():
    r"""Test the Alfven_speed function in parameters.py."""

    assert np.isclose(Alfven_speed(1 * u.T, 1e-8 * u.kg * u.m**-3).value,
                      8920620.580763856)

    V_A = Alfven_speed(B, n_i)
    assert np.isclose(
        V_A.value, (B / np.sqrt(mu0 * n_i * (m_p + m_e))).si.value)

    assert Alfven_speed(B, rho) == Alfven_speed(B, n_i)

    assert Alfven_speed(B, rho).unit == u.m / u.s

    assert Alfven_speed(B, rho) == Alfven_speed(-B, rho)

    assert Alfven_speed(B, 4 * rho) == 0.5 * Alfven_speed(B, rho)

    assert Alfven_speed(2 * B, rho) == 2 * Alfven_speed(B, rho)

    # Case when Z=1 is assumed
    assert Alfven_speed(5 * u.T, 5e19 * u.m**-3, ion='H-1') == \
        Alfven_speed(5 * u.T, 5e19 * u.m**-3, ion='p')

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
        Alfven_speed(np.array([5, 6, 7]) * u.T,
                     np.array([5, 6]) * u.m**-3)

    with pytest.raises(ValueError):
        Alfven_speed(B_nanarr, rho_arr)

    with pytest.raises(ValueError):
        Alfven_speed(B_arr, rho_negarr)

    with pytest.raises(u.UnitConversionError):
        Alfven_speed(5 * u.A, n_i, ion='p')

    with pytest.raises(TypeError):
        Alfven_speed(B, 5, ion='p')

    with pytest.raises(u.UnitsError):
        Alfven_speed(B, 5 * u.m**-2, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(B, n_i, ion='spacecats')

    with pytest.warns(RelativityWarning):  # relativistic
        Alfven_speed(5e1 * u.T, 5e19 * u.m**-3, ion='p')

    with pytest.raises(RelativityError):  # super-relativistic
        Alfven_speed(5e8 * u.T, 5e19 * u.m**-3, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(0.001 * u.T, -5e19 * u.m**-3, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(np.nan * u.T, 1 * u.m**-3, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(1 * u.T, np.nan * u.m**-3, ion='p')

    with pytest.raises(RelativityError):
        assert Alfven_speed(np.inf * u.T, 1 * u.m**-3,
                            ion='p') == np.inf * u.m / u.s

    with pytest.raises(RelativityError):
        assert Alfven_speed(-np.inf * u.T, 1 * u.m**-3,
                            ion='p') == np.inf * u.m / u.s

    with pytest.raises(UserWarning):
        assert Alfven_speed(1.0, n_i) == Alfven_speed(1.0 * u.T, n_i)


def test_ion_sound_speed():
    r"""Test the ion_sound_speed function in parameters.py."""

    assert np.isclose(ion_sound_speed(T_i=1.3232 * u.MK, T_e=1.831 * u.MK,
                                      ion='p', gamma_e=1, gamma_i=3).value,
                      218816.06086407552)

    assert np.isclose(ion_sound_speed(T_i=0.88 * u.MK, T_e=1.28 * u.MK, ion='p',
                                      gamma_e=1.2, gamma_i=3.4).value,
                      193328.52857788358)

    # case when Z=1 is assumed
    assert ion_sound_speed(T_i=T_i, T_e=T_e, ion='p') == \
        ion_sound_speed(T_i=T_i, T_e=T_e, ion='H-1')

    assert ion_sound_speed(T_i=T_i, ion='p').unit == u.m / u.s

    assert ion_sound_speed(T_i=T_i, gamma_i=3) == ion_sound_speed(T_i=T_i)

    assert ion_sound_speed(T_e=T_e, gamma_e=1) == ion_sound_speed(T_e=T_e)

    with pytest.raises(RelativityError):
        ion_sound_speed(T_i=T_i, gamma_i=np.inf)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=np.array([5, 6, 5])
                        * u.K, T_e=np.array([3, 4]) * u.K)

    with pytest.raises(TypeError):    # Is this test right??????
        ion_sound_speed(5 * u.T)

    with pytest.raises(TypeError):
        ion_sound_speed('p')

    with pytest.raises(PhysicsError):
        ion_sound_speed(T_i=T_i, gamma_i=0.9999)

    with pytest.raises(PhysicsError):
        ion_sound_speed(T_i=T_i, gamma_e=0.9999)

    with pytest.raises(TypeError):
        ion_sound_speed(T_i=T_i, gamma_e='sdjklsf')

    with pytest.raises(TypeError):
        ion_sound_speed(T_i=T_i, gamma_i='fsdfas')

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_i, ion='cupcakes')

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=-np.abs(T_i))

    with pytest.warns(RelativityWarning):
        ion_sound_speed(T_i=5e11 * u.K)

    with pytest.raises(RelativityError):
        ion_sound_speed(T_i=5e19 * u.K)

    with pytest.raises(u.UnitConversionError):
        ion_sound_speed(T_i=5 * u.A)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_nanarr)

    with pytest.raises(ValueError):
        ion_sound_speed(T_e=T_nanarr)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_negarr)

    with pytest.raises(ValueError):
        ion_sound_speed(T_e=T_negarr)

    with pytest.raises(UserWarning):
        assert ion_sound_speed(T_e=1.2e6) == ion_sound_speed(T_e=1.2e6 * u.K)

    with pytest.raises(UserWarning):
        assert ion_sound_speed(T_i=1.3e6) == ion_sound_speed(T_i=1.3e6 * u.K)


def test_thermal_speed():
    r"""Test the thermal_speed function in parameters.py"""
    assert thermal_speed(T_e).unit == u.m / u.s

    assert thermal_speed(T_e) > thermal_speed(T_e, 'p')

    # The NRL Plasma Formulary uses a definition of the electron
    # thermal speed that differs by a factor of sqrt(2).
    assert np.isclose(thermal_speed(1 * u.MK).value,
                      5505694.743141063)

    with pytest.raises(u.UnitConversionError):
        thermal_speed(5 * u.m)

    with pytest.raises(ValueError):
        thermal_speed(-T_e)

    with pytest.warns(RelativityWarning):
        thermal_speed(1e9 * u.K)

    with pytest.raises(RelativityError):
        thermal_speed(5e19 * u.K)

    with pytest.raises(UserWarning):
        assert thermal_speed(1e5) == thermal_speed(1e5 * u.K)

    assert thermal_speed(T_i, particle='p').unit == u.m / u.s

    # The NRL Plasma Formulary uses a definition of the particle thermal
    # speed that differs by a factor of sqrt(2).
    assert np.isclose(thermal_speed(1 * u.MK, particle='p').si.value,
                      128486.56960876315)

    # Case when Z=1 is assumed
    assert thermal_speed(T_i, particle='p') == \
        thermal_speed(T_i, particle='H-1')

    assert thermal_speed(1 * u.MK, particle='e+') == \
        thermal_speed(1 * u.MK)

    with pytest.raises(u.UnitConversionError):
        thermal_speed(5 * u.m, particle='p')

    with pytest.raises(ValueError):
        thermal_speed(-T_e, particle='p')

    with pytest.warns(RelativityWarning):
        thermal_speed(1e11 * u.K, particle='p')

    with pytest.raises(RelativityError):
        thermal_speed(1e14 * u.K, particle='p')

    with pytest.raises(ValueError):
        thermal_speed(T_i, particle='asdfasd')

    with pytest.raises(UserWarning):
        assert thermal_speed(1e6, particle='p') ==\
            thermal_speed(1e6 * u.K, particle='p')

    assert np.isclose(thermal_speed(1e6 * u.K,
                                    method="mean_magnitude").si.value,
                      19517177.023383822)

    assert np.isclose(thermal_speed(1e6 * u.K, method="rms").si.value,
                      6743070.475775486)

    with pytest.raises(ValueError):
        thermal_speed(T_i, method="sadks")

# test class for kappa_thermal_speed() function:


class Test_kappa_thermal_speed(object):
    def setup_method(self):
        """initializing parameters for tests """
        self.T_e = 5 * u.eV
        self.kappaInvalid = 3 / 2
        self.kappa = 4
        self.particle = "p"
        self.probable1True = 24467.878463594963
        self.rms1True = 29966.908662120648
        self.mean1True = 86736.37081257407
    def test_invalid_kappa(self):
        """
        Checks if function raises error when kappa <= 3/2 is passed as an
        argument.
        """
        with pytest.raises(ValueError):
            kappa_thermal_speed(self.T_e,
                                self.kappaInvalid,
                                particle=self.particle)
        return
    def test_invalid_method(self):
        """
        Checks if function raises error when invalid method is passed as an
        argument.
        """
        with pytest.raises(ValueError):
            kappa_thermal_speed(self.T_e,
                                self.kappa,
                                particle=self.particle,
                                method="invalid")
        return
    def test_probable1(self):
        """
        Tests if expected value is returned for a set of regular inputs.
        """
        known1 = kappa_thermal_speed(self.T_e,
                                     self.kappa,
                                     particle=self.particle,
                                     method="most_probable")
        errStr = (f"Kappa thermal velocity should be {self.probable1True} "
                  f"and not {known1.si.value}.")
        assert np.isclose(known1.value,
                          self.probable1True,
                          rtol=1e-8,
                          atol=0.0), errStr
        return
    def test_rms1(self):
        """
        Tests if expected value is returned for a set of regular inputs.
        """
        known1 = kappa_thermal_speed(self.T_e,
                                     self.kappa,
                                     particle=self.particle,
                                     method="rms")
        errStr = (f"Kappa thermal velocity should be {self.rms1True} "
                  f"and not {known1.si.value}.")
        assert np.isclose(known1.value,
                          self.rms1True,
                          rtol=1e-8,
                          atol=0.0), errStr
        return
    def test_mean1(self):
        """
        Tests if expected value is returned for a set of regular inputs.
        """
        known1 = kappa_thermal_speed(self.T_e,
                                     self.kappa,
                                     particle=self.particle,
                                     method="mean_magnitude")
        errStr = (f"Kappa thermal velocity should be {self.mean1True} "
                  f"and not {known1.si.value}.")
        assert np.isclose(known1.value,
                          self.mean1True,
                          rtol=1e-8,
                          atol=0.0), errStr
        return
    
    


def test_gyrofrequency():
    r"""Test the gyrofrequency function in parameters.py."""

    assert gyrofrequency(B).unit == u.rad / u.s

    assert np.isclose(gyrofrequency(1 * u.T).value, 175882008784.72018)

    assert np.isclose(gyrofrequency(2.4 * u.T).value,
                      422116821083.3284)

    assert np.isclose(gyrofrequency(1 * u.G).cgs.value,
                      1.76e7, rtol=1e-3)

    with pytest.raises(TypeError):
        gyrofrequency(u.m)

    with pytest.raises(u.UnitConversionError):
        gyrofrequency(u.m * 1)

    with pytest.raises(ValueError):
        gyrofrequency(B_nanarr)

    # The following is a test to check that equivalencies from astropy
    # are working.
    omega_ce = gyrofrequency(2.2 * u.T)
    f_ce = (omega_ce / (2 * np.pi)) / u.rad
    f_ce_use_equiv = omega_ce.to(u.Hz, equivalencies=[(u.cy / u.s, u.Hz)])
    assert np.isclose(f_ce.value, f_ce_use_equiv.value)

    with pytest.raises(UserWarning):
        assert gyrofrequency(5.0) == gyrofrequency(5.0 * u.T)

    assert gyrofrequency(B, particle=ion).unit == u.rad / u.s

    assert np.isclose(gyrofrequency(1 * u.T, particle='p').value,
                      95788335.834874)

    assert np.isclose(gyrofrequency(2.4 * u.T, particle='p').value,
                      229892006.00369796)

    assert np.isclose(gyrofrequency(1 * u.G, particle='p').cgs.value,
                      9.58e3, rtol=2e-3)

    assert gyrofrequency(-5 * u.T, 'p') == gyrofrequency(5 * u.T, 'p')

    # Case when Z=1 is assumed
    assert gyrofrequency(B, particle='p') == \
        gyrofrequency(B, particle='H-1')

    assert gyrofrequency(B, particle='e+') == \
        gyrofrequency(B)

    with pytest.raises(UserWarning):
        gyrofrequency(8, 'p')

    with pytest.raises(u.UnitConversionError):
        gyrofrequency(5 * u.m, 'p')

    with pytest.raises(ValueError):
        gyrofrequency(8 * u.T, particle='asdfasd')

    with pytest.raises(UserWarning):
        assert gyrofrequency(5.0, 'p') == gyrofrequency(5.0 * u.T, 'p')


def test_gyroradius():
    r"""Test the gyroradius function in parameters.py."""

    assert gyroradius(B, T_e).unit == u.m

    assert gyroradius(B, 25 * u.m / u.s).unit == u.m

    assert gyroradius(T_e, B) == gyroradius(B, T_e)

    assert gyroradius(V, B) == gyroradius(B, V)

    assert gyroradius(B, V) == gyroradius(B, -V)

    Vperp = 1e6 * u.m / u.s
    Bmag = 1 * u.T
    omega_ce = gyrofrequency(Bmag)
    assert gyroradius(Bmag, Vperp) == \
        (Vperp / omega_ce).to(u.m, equivalencies=u.dimensionless_angles())

    with pytest.raises(TypeError):
        gyroradius(u.T, 8 * u.m / u.s)

    with pytest.raises(u.UnitConversionError):
        gyroradius(5 * u.A, 8 * u.m / u.s)

    with pytest.raises(u.UnitConversionError):
        gyroradius(5 * u.T, 8 * u.m)

    with pytest.raises(ValueError):
        gyroradius(np.array([5, 6]) * u.T, np.array([5, 6, 7]) * u.m / u.s)

    with pytest.raises(ValueError):
        gyroradius(np.nan * u.T, 1 * u.m / u.s)

    with pytest.raises(ValueError):
        gyroradius(3.14159 * u.T, -1 * u.K)

    with pytest.raises(UserWarning):
        assert gyroradius(1.0, Vperp=1.0) == \
            gyroradius(1.0 * u.T, Vperp=1.0 * u.m / u.s)

    with pytest.raises(UserWarning):
        assert gyroradius(1.1, T_i=1.2) == \
            gyroradius(1.1 * u.T, T_i=1.2 * u.K)

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, T_i=1.2 * u.K, Vperp=1 * u.m / u.s)

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, 1.2 * u.K, 1.1 * u.m)

    assert gyroradius(B, T_i, particle="p").unit == u.m

    assert gyroradius(B, 25 * u.m / u.s, particle="p").unit == u.m

    # Case when Z=1 is assumed
    assert gyroradius(B, T_i, particle='p') == \
        gyroradius(B, T_i, particle='H-1')

    assert gyroradius(T_i, B, particle="p") == gyroradius(B, T_i, particle="p")

    assert gyroradius(V, B, particle="p") == gyroradius(B, V, particle="p")

    assert gyroradius(B, V, particle="p") == gyroradius(B, -V, particle="p")

    Vperp = 1e6 * u.m / u.s
    Bmag = 1 * u.T
    omega_ci = gyrofrequency(Bmag, particle='p')
    assert gyroradius(Bmag, Vperp, particle="p") == \
        (Vperp / omega_ci).to(u.m, equivalencies=u.dimensionless_angles())

    T2 = 1.2 * u.MK
    B2 = 123 * u.G
    particle2 = 'alpha'
    Vperp2 = thermal_speed(T2, particle=particle2)
    assert gyroradius(B2, Vperp=Vperp2, particle='alpha') == \
        gyroradius(B2, T_i=T2, particle='alpha')

    assert gyroradius(1 * u.T, 1 * u.MK, particle='positron') == \
        gyroradius(1 * u.T, 1 * u.MK)

    with pytest.raises(TypeError):
        gyroradius(u.T, 8 * u.m / u.s, particle="p")

    with pytest.raises(ValueError):
        gyroradius(B, T_i, particle='asfdas')

    with pytest.raises(ValueError):
        gyroradius(B, -1 * u.K, particle='p')

    with pytest.raises(UserWarning):
        assert gyroradius(1.0, Vperp=1.0, particle="p") == \
            gyroradius(1.0 * u.T, Vperp=1.0 * u.m / u.s, particle="p")

    with pytest.raises(UserWarning):
        assert gyroradius(1.1, T_i=1.2, particle="p") == \
            gyroradius(1.1 * u.T, T_i=1.2 * u.K, particle="p")

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, T_i=1.2 * u.K, Vperp=1 * u.m / u.s, particle="p")

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, 1.2 * u.K, 1.1 * u.m, particle="p")

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, 1.2 * u.m, 1.1 * u.K, particle="p")


def test_plasma_frequency():
    r"""Test the plasma_frequency function in parameters.py."""

    assert plasma_frequency(n_e).unit == u.rad / u.s

    assert np.isclose(plasma_frequency(1 * u.cm**-3).value,
                      5.64e4, rtol=1e-2)

    with pytest.raises(TypeError):
        plasma_frequency(u.m**-3)

    with pytest.raises(u.UnitConversionError):
        plasma_frequency(5 * u.m**-2)

    with pytest.raises(ValueError):
        plasma_frequency(np.nan * u.m**-3)

    with pytest.raises(UserWarning):
        assert plasma_frequency(1e19) == \
            plasma_frequency(1e19 * u.m**-3)

        assert plasma_frequency(n_i, particle='p').unit == u.rad / u.s

    # Case where Z=1 is assumed
    assert plasma_frequency(n_i, particle='H-1') == \
        plasma_frequency(n_i, particle='p')

    assert np.isclose(plasma_frequency(mu * u.cm**-3, particle='p').value,
                      1.32e3, rtol=1e-2)

    with pytest.raises(ValueError):
        plasma_frequency(n=5 * u.m**-3, particle='sdfas')

    with pytest.raises(UserWarning):
        assert plasma_frequency(1e19, particle='p') ==\
            plasma_frequency(1e19 * u.m**-3, particle='p')


def test_Debye_length():
    r"""Test the Debye_length function in parameters.py."""

    assert Debye_length(T_e, n_e).unit == u.m

    assert np.isclose(Debye_length(
        1 * u.eV, 1 * u.cm**-3).value, 7.43, atol=0.005)

    with pytest.raises(UserWarning):
        Debye_length(5, 5 * u.m**-3)

    with pytest.raises(u.UnitConversionError):
        Debye_length(56 * u.kg, 5 * u.m**-3)

    with pytest.raises(ValueError):
        Debye_length(5 * u.eV, -5 * u.m**-3)

    with pytest.raises(ValueError):
        Debye_length(-45 * u.K, 5 * u.m**-3)

    Tarr2 = np.array([1, 2]) * u.K
    narr3 = np.array([1, 2, 3]) * u.m**-3
    with pytest.raises(ValueError):
        Debye_length(Tarr2, narr3)

    with pytest.raises(UserWarning):
        assert Debye_length(2.0, 2.0) == Debye_length(2.0 * u.K, 2.0 * u.m**-3)

    with pytest.raises(UserWarning):
        assert Debye_length(2.0 * u.K, 2.0) == Debye_length(2.0, 2.0 * u.m**-3)


def test_Debye_number():
    r"""Test the Debye_number function in parameters.py."""

    assert Debye_number(T_e, n_e).unit == u.dimensionless_unscaled

    T_e_eV = T_e.to(u.eV, equivalencies=u.temperature_energy())
    assert np.isclose(Debye_number(T_e, n_e).value,
                      Debye_number(T_e_eV, n_e).value)

    assert np.isclose(Debye_number(
        1 * u.eV, 1 * u.cm**-3).value, 1720862385.43342)

    with pytest.raises(UserWarning):
        Debye_number(T_e, 4)

    with pytest.raises(TypeError):
        Debye_number(None, n_e)

    with pytest.raises(u.UnitConversionError):
        Debye_number(5 * u.m, 5 * u.m**-3)

    with pytest.raises(u.UnitConversionError):
        Debye_number(5 * u.K, 5 * u.m**3)

    with pytest.raises(ValueError):
        Debye_number(5j * u.K, 5 * u.cm**-3)

    Tarr2 = np.array([1, 2]) * u.K
    narr3 = np.array([1, 2, 3]) * u.m**-3
    with pytest.raises(ValueError):
        Debye_number(Tarr2, narr3)

    with pytest.raises(UserWarning):
        assert Debye_number(1.1, 1.1) == Debye_number(1.1 * u.K, 1.1 * u.m**-3)

    with pytest.raises(UserWarning):
        assert Debye_number(1.1 * u.K, 1.1) == Debye_number(1.1, 1.1 * u.m**-3)


def test_inertial_length():
    r"""Test the inertial_length function in parameters.py."""

    assert inertial_length(n_i, particle='p').unit == u.m

    assert np.isclose(inertial_length(mu * u.cm**-3, particle='p').cgs.value,
                      2.28e7, rtol=0.01)

    assert inertial_length(5.351 * u.m**-3, particle='e+') == \
        inertial_length(5.351 * u.m**-3, particle='e')

    assert inertial_length(n_i, particle='p') == \
        inertial_length(n_i, particle='H-1')

    with pytest.raises(UserWarning):
        inertial_length(4, particle='p')

    with pytest.raises(u.UnitConversionError):
        inertial_length(4 * u.m**-2, particle='p')

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m**-3, particle='p')

    with pytest.raises(ValueError):
        inertial_length(n_i, particle=-135)

    with pytest.raises(UserWarning):
        assert inertial_length(1e19, particle='p') ==\
            inertial_length(1e19 * u.m**-3, particle='p')

    assert inertial_length(n_e).unit == u.m

    assert np.isclose(inertial_length(1 * u.cm**-3).cgs.value,
                      5.31e5, rtol=1e-3)

    with pytest.raises(UserWarning):
        inertial_length(5)

    with pytest.raises(u.UnitConversionError):
        inertial_length(5 * u.m)

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m**-3)

    with pytest.raises(UserWarning):
        assert inertial_length(1e19) == \
            inertial_length(1e19 * u.m**-3)


def test_magnetic_pressure():
    r"""Test the magnetic_pressure function in parameters.py."""

    assert magnetic_pressure(B_arr).unit == u.Pa

    assert magnetic_pressure(B).unit == u.Pa

    assert magnetic_pressure(B).unit.name == 'Pa'

    assert magnetic_pressure(B).value == magnetic_energy_density(B).value

    assert magnetic_pressure(B) == magnetic_energy_density(B.to(u.G))

    assert np.isclose(magnetic_pressure(B).value, 397887.35772973835)

    with pytest.raises(UserWarning):
        magnetic_pressure(5)

    with pytest.raises(u.UnitConversionError):
        magnetic_pressure(5 * u.m)

    with pytest.raises(ValueError):
        magnetic_pressure(np.nan * u.T)

    with pytest.raises(ValueError):
        magnetic_pressure(5j * u.T)

    with pytest.raises(ValueError):
        magnetic_pressure(B_nanarr)

    with pytest.raises(UserWarning):
        assert magnetic_pressure(22.2) == magnetic_pressure(22.2 * u.T)


def test_magnetic_energy_density():
    r"""Test the magnetic_energy_density function in parameters.py."""

    assert magnetic_energy_density(B_arr).unit == u.J / u.m**3

    assert str(magnetic_energy_density(B).unit) == 'J / m3'

    assert magnetic_energy_density(B).value == magnetic_pressure(B).value

    assert magnetic_energy_density(2 * B) == 4 * magnetic_energy_density(B)

    assert np.isclose(magnetic_energy_density(B).value, 397887.35772973835)

    assert magnetic_energy_density(B) == magnetic_energy_density(B.to(u.G))

    # Add an array test!

    assert magnetic_energy_density(B_arr)

    with pytest.raises(UserWarning):
        magnetic_energy_density(5)

    with pytest.raises(u.UnitConversionError):
        magnetic_energy_density(5 * u.m)

    with pytest.raises(ValueError):
        magnetic_energy_density(np.nan * u.T)

    with pytest.raises(ValueError):
        magnetic_energy_density(5j * u.T)

    with pytest.raises(ValueError):
        magnetic_energy_density(B_nanarr)

    with pytest.raises(UserWarning):
        assert magnetic_energy_density(22.2) == \
            magnetic_energy_density(22.2 * u.T)


def test_upper_hybrid_frequency():
    r"""Test the upper_hybrid_frequency function in parameters.py."""

    omega_uh = upper_hybrid_frequency(B, n_e=n_e)
    omega_ce = gyrofrequency(B)
    omega_pe = plasma_frequency(n=n_e)
    assert omega_ce.unit == u.rad / u.s
    assert omega_pe.unit == u.rad / u.s
    assert omega_uh.unit == u.rad / u.s
    LHS = omega_uh**2
    RHS = omega_ce**2 + omega_pe**2
    assert np.isclose(LHS.value, RHS.value)

    with pytest.raises(ValueError):
        upper_hybrid_frequency(5 * u.T, n_e=-1 * u.m**-3)

    with pytest.raises(UserWarning):
        assert upper_hybrid_frequency(1.2, 1.3) == \
            upper_hybrid_frequency(1.2 * u.T, 1.3 * u.m**-3)

    with pytest.raises(UserWarning):
        assert upper_hybrid_frequency(1.4 * u.T, 1.3) == \
            upper_hybrid_frequency(1.4, 1.3 * u.m**-3)


def test_lower_hybrid_frequency():
    r"""Test the lower_hybrid_frequency function in parameters.py."""

    ion = 'He-4 1+'
    omega_ci = gyrofrequency(B, particle=ion)
    omega_pi = plasma_frequency(n=n_i, particle=ion)
    omega_ce = gyrofrequency(B)
    omega_lh = lower_hybrid_frequency(B, n_i=n_i, ion=ion)
    assert omega_ci.unit == u.rad / u.s
    assert omega_pi.unit == u.rad / u.s
    assert omega_ce.unit == u.rad / u.s
    assert omega_lh.unit == u.rad / u.s
    LHS = omega_lh**-2
    RHS = 1 / (omega_ci**2 + omega_pi**2) + omega_ci**-1 * omega_ce**-1
    assert np.isclose(LHS.value, RHS.value)

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=5e19 * u.m**-3, ion='asdfasd')

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=-5e19 * u.m**-3, ion='asdfasd')

    with pytest.raises(ValueError):
        lower_hybrid_frequency(
            np.nan * u.T, n_i=-5e19 * u.m**-3, ion='asdfasd')

    with pytest.raises(UserWarning):
        assert lower_hybrid_frequency(1.3, 1e19) == \
            lower_hybrid_frequency(1.3 * u.T, 1e19 * u.m**-3)
