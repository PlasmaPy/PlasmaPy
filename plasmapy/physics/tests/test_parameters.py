"""Tests for functions that calculate plasma parameters."""

import numpy as np
import pytest
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from warnings import simplefilter

from plasmapy.utils.exceptions import RelativityWarning, RelativityError
from plasmapy.utils.exceptions import PhysicsError, InvalidParticleError
from astropy.constants import c, m_p, m_e, e, mu0

from plasmapy.physics.parameters import (mass_density,
                                         Alfven_speed,
                                         gyrofrequency,
                                         gyroradius,
                                         thermal_speed,
                                         thermal_pressure,
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
n_i = 5e19 * u.m ** -3
n_e = Z * 5e19 * u.m ** -3
rho = n_i * m_i + n_e * m_e
T_e = 1e6 * u.K
T_i = 1e6 * u.K

B_arr = np.array([0.001, 0.002]) * u.T
rho_arr = np.array([5e-10, 2e-10]) * u.kg / u.m ** 3

B_nanarr = np.array([0.001, np.nan]) * u.T
rho_infarr = np.array([np.inf, 5e19]) * u.m ** -3
rho_negarr = np.array([-5e19, 6e19]) * u.m ** -3
T_nanarr = np.array([1e6, np.nan]) * u.K
T_negarr = np.array([1e6, -5151.]) * u.K

mu = m_p.to(u.u).value

V = 25.2 * u.m / u.s


class Test_mass_density:
    def test_particleless(self):
        with pytest.raises(ValueError):
            mass_density(1 * u.m ** -3)
    def test_wrong_units(self):
        with pytest.raises(ValueError):
            mass_density(1 * u.J)


# Assertions below that are in CGS units with 2-3 significant digits
# are generally from the NRL Plasma Formulary.


def test_Alfven_speed():
    r"""Test the Alfven_speed function in parameters.py."""

    assert np.isclose(Alfven_speed(1 * u.T, 1e-8 * u.kg * u.m ** -3).value,
                      8920620.580763856)

    V_A = Alfven_speed(B, n_i)
    assert np.isclose(
        V_A.value, (B / np.sqrt(mu0 * n_i * (m_p + m_e))).si.value)

    assert Alfven_speed(B, rho) == Alfven_speed(B, n_i)

    assert Alfven_speed(B, rho).unit.is_equivalent(u.m / u.s)

    assert Alfven_speed(B, rho) == Alfven_speed(-B, rho)

    assert Alfven_speed(B, 4 * rho) == 0.5 * Alfven_speed(B, rho)

    assert Alfven_speed(2 * B, rho) == 2 * Alfven_speed(B, rho)

    # Case when Z=1 is assumed
    with pytest.warns(RelativityWarning):
        assert np.isclose(Alfven_speed(5 * u.T, 5e19 * u.m ** -3, ion='H+'),
                          Alfven_speed(5 * u.T, 5e19 * u.m ** -3, ion='p'),
                          atol=0 * u.m / u.s,
                          rtol=1e-3)

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
                     np.array([5, 6]) * u.m ** -3)

    with pytest.raises(ValueError):
        Alfven_speed(B_nanarr, rho_arr)

    with pytest.raises(ValueError):
        Alfven_speed(B_arr, rho_negarr)

    with pytest.raises(u.UnitConversionError):
        Alfven_speed(5 * u.A, n_i, ion='p')

    with pytest.raises(TypeError):
        Alfven_speed(B, 5, ion='p')

    with pytest.raises(u.UnitsError):
        Alfven_speed(B, 5 * u.m ** -2, ion='p')

    with pytest.raises(InvalidParticleError):
        Alfven_speed(B, n_i, ion='spacecats')

    with pytest.warns(RelativityWarning):  # relativistic
        Alfven_speed(5e1 * u.T, 5e19 * u.m ** -3, ion='p')

    with pytest.raises(RelativityError):  # super-relativistic
        Alfven_speed(5e8 * u.T, 5e19 * u.m ** -3, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(0.001 * u.T, -5e19 * u.m ** -3, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(np.nan * u.T, 1 * u.m ** -3, ion='p')

    with pytest.raises(ValueError):
        Alfven_speed(1 * u.T, np.nan * u.m ** -3, ion='p')

    with pytest.raises(RelativityError):
        assert Alfven_speed(np.inf * u.T, 1 * u.m ** -3,
                            ion='p') == np.inf * u.m / u.s

    with pytest.raises(RelativityError):
        assert Alfven_speed(-np.inf * u.T, 1 * u.m ** -3,
                            ion='p') == np.inf * u.m / u.s

    with pytest.warns(u.UnitsWarning):
        assert Alfven_speed(1.0, n_i) == Alfven_speed(1.0 * u.T, n_i)

    Alfven_speed(1 * u.T, 5e19 * u.m ** -3, ion='p')
    # testing for user input z_mean
    testMeth1 = Alfven_speed(1 * u.T,
                             5e19 * u.m ** -3,
                             ion='p',
                             z_mean=0.8).si.value
    testTrue1 = 3084015.75214846
    errStr = f"Alfven_speed() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1,
                      testTrue1,
                      atol=0.0,
                      rtol=1e-15), errStr


def test_ion_sound_speed():
    r"""Test the ion_sound_speed function in parameters.py."""

    assert np.isclose(ion_sound_speed(T_i=1.3232 * u.MK, T_e=1.831 * u.MK,
                                      ion='p', gamma_e=1, gamma_i=3).value,
                      218816.06086407552)

    assert np.isclose(ion_sound_speed(
        T_i=0.88 * u.MK, T_e=1.28 * u.MK, ion='p', gamma_e=1.2,
        gamma_i=3.4).value, 193328.52857788358)

    # case when Z=1 is assumed
    # assert ion_sound_speed(T_i=T_i, T_e=T_e, ion='p+') == ion_sound_speed(T_i=T_i, T_e=T_e,
    # ion='H-1')

    assert ion_sound_speed(T_i=T_i, T_e=0 * u.K, ion='p+').unit.is_equivalent(u.m / u.s)

    with pytest.raises(RelativityError):
        ion_sound_speed(T_i=T_i, T_e=T_e, gamma_i=np.inf)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=np.array([5, 6, 5]) * u.K,
                        T_e=np.array([3, 4]) * u.K)

    with pytest.raises(TypeError):  # Is this test right??????
        ion_sound_speed(5 * u.T)

    with pytest.raises(TypeError):
        ion_sound_speed('p')

    with pytest.raises(PhysicsError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, gamma_i=0.9999)

    with pytest.raises(PhysicsError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, gamma_e=0.9999)

    with pytest.raises(TypeError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, gamma_e='sdjklsf')

    with pytest.raises(TypeError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, gamma_i='fsdfas')

    with pytest.raises(InvalidParticleError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, ion='cupcakes')

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=-np.abs(T_i), T_e=0 * u.K, )

    with pytest.warns(RelativityWarning):
        ion_sound_speed(T_i=5e11 * u.K, T_e=0 * u.K)

    with pytest.raises(RelativityError):
        ion_sound_speed(T_i=5e19 * u.K, T_e=0 * u.K)

    with pytest.raises(u.UnitConversionError):
        ion_sound_speed(T_i=5 * u.A, T_e=0 * u.K)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_nanarr, T_e=0 * u.K)

    with pytest.raises(ValueError):
        ion_sound_speed(T_e=T_nanarr, T_i=0 * u.K)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_negarr, T_e=0 * u.K)

    with pytest.raises(ValueError):
        ion_sound_speed(T_e=T_negarr, T_i=0 * u.K)

    with pytest.warns(u.UnitsWarning):
        assert ion_sound_speed(T_e=1.2e6, T_i=0 * u.K) == ion_sound_speed(T_e=1.2e6 * u.K,
                                                                          T_i=0 * u.K)

    with pytest.warns(u.UnitsWarning):
        assert ion_sound_speed(T_i=1.3e6, T_e=0 * u.K) == ion_sound_speed(T_i=1.3e6 * u.K,
                                                                          T_e=0 * u.K)

    ion_sound_speed(T_e=1.2e6 * u.K, T_i=0 * u.K)
    # testing for user input z_mean
    testMeth1 = ion_sound_speed(T_e=1.2e6 * u.K, T_i=0 * u.K, z_mean=0.8).si.value
    testTrue1 = 89018.0944146141
    errStr = f"ion_sound_speed() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1,
                      testTrue1,
                      atol=0.0,
                      rtol=1e-15), errStr


def test_thermal_pressure():
    assert thermal_pressure(T_e, n_i).unit.is_equivalent(u.Pa)


def test_thermal_speed():
    r"""Test the thermal_speed function in parameters.py"""
    assert thermal_speed(T_e).unit.is_equivalent(u.m / u.s)

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

    with pytest.warns(u.UnitsWarning):
        assert thermal_speed(1e5) == thermal_speed(1e5 * u.K)

    assert thermal_speed(T_i, particle='p').unit.is_equivalent(u.m / u.s)

    # The NRL Plasma Formulary uses a definition of the particle thermal
    # speed that differs by a factor of sqrt(2).
    assert np.isclose(thermal_speed(1 * u.MK, particle='p').si.value,
                      128486.56960876315)

    # Case when Z=1 is assumed
    assert thermal_speed(T_i, particle='p') == thermal_speed(T_i, particle='H-1+')

    assert thermal_speed(1 * u.MK, particle='e+') == thermal_speed(1 * u.MK)

    with pytest.raises(u.UnitConversionError):
        thermal_speed(5 * u.m, particle='p')

    with pytest.raises(ValueError):
        thermal_speed(-T_e, particle='p')

    with pytest.warns(RelativityWarning):
        thermal_speed(1e11 * u.K, particle='p')

    with pytest.raises(RelativityError):
        thermal_speed(1e14 * u.K, particle='p')

    with pytest.raises(InvalidParticleError):
        thermal_speed(T_i, particle='asdfasd')

    with pytest.warns(u.UnitsWarning):
        assert thermal_speed(1e6, particle='p') == thermal_speed(1e6 * u.K, particle='p')

    assert np.isclose(thermal_speed(1e6 * u.K,
                                    method="mean_magnitude").si.value,
                      6212510.3969422)

    assert np.isclose(thermal_speed(1e6 * u.K, method="rms").si.value,
                      6743070.475775486)

    with pytest.raises(ValueError):
        thermal_speed(T_i, method="sadks")


# test class for kappa_thermal_speed() function:


class Test_kappa_thermal_speed(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T_e = 5 * u.eV
        self.kappaInvalid = 3 / 2
        self.kappa = 4
        self.particle = "p"
        self.probable1True = 24467.878463594963
        self.rms1True = 37905.474322612165
        self.mean1True = 34922.98563039583

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

    assert gyrofrequency(B).unit.is_equivalent(u.rad / u.s)

    assert np.isclose(gyrofrequency(1 * u.T).value, 175882008784.72018)

    assert np.isclose(gyrofrequency(2.4 * u.T).value,
                      422116821083.3284)

    assert np.isclose(gyrofrequency(2.4 * u.T, signed=True).value,
                      -422116821083.3284)

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

    with pytest.warns(u.UnitsWarning):
        assert gyrofrequency(5.0) == gyrofrequency(5.0 * u.T)

    assert gyrofrequency(B, particle=ion).unit.is_equivalent(u.rad / u.s)

    assert np.isclose(gyrofrequency(1 * u.T, particle='p').value,
                      95788335.834874)

    assert np.isclose(gyrofrequency(2.4 * u.T, particle='p').value,
                      229892006.00369796)

    assert np.isclose(gyrofrequency(1 * u.G, particle='p').cgs.value,
                      9.58e3, rtol=2e-3)

    assert gyrofrequency(-5 * u.T, 'p') == gyrofrequency(5 * u.T, 'p')

    # Case when Z=1 is assumed
    # assert gyrofrequency(B, particle='p+') == gyrofrequency(B, particle='H-1')

    assert gyrofrequency(B, particle='e+') == gyrofrequency(B)

    with pytest.warns(u.UnitsWarning):
        gyrofrequency(8, 'p')

    with pytest.raises(u.UnitConversionError):
        gyrofrequency(5 * u.m, 'p')

    with pytest.raises(InvalidParticleError):
        gyrofrequency(8 * u.T, particle='asdfasd')

    with pytest.warns(u.UnitsWarning):
        # TODO this should be WARNS, not RAISES. and it's probably still raised
        assert gyrofrequency(5.0, 'p') == gyrofrequency(5.0 * u.T, 'p')

    gyrofrequency(1 * u.T, particle='p')
    # testing for user input Z
    testMeth1 = gyrofrequency(1 * u.T, particle='p', Z=0.8).si.value
    testTrue1 = 76630665.79318453
    errStr = f"gyrofrequency() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1,
                      testTrue1,
                      atol=0.0,
                      rtol=1e-15), errStr


def test_gyroradius():
    r"""Test the gyroradius function in parameters.py."""

    assert gyroradius(B, T_i=T_e).unit.is_equivalent(u.m)

    assert gyroradius(B, Vperp=25 * u.m / u.s).unit.is_equivalent(u.m)

    Vperp = 1e6 * u.m / u.s
    Bmag = 1 * u.T
    omega_ce = gyrofrequency(Bmag)
    analytical_result = (Vperp / omega_ce).to(u.m, equivalencies=u.dimensionless_angles())
    assert gyroradius(Bmag, Vperp=Vperp) == analytical_result

    with pytest.raises(TypeError):
        gyroradius(u.T)

    with pytest.raises(u.UnitConversionError):
        gyroradius(5 * u.A, Vperp=8 * u.m / u.s)

    with pytest.raises(u.UnitConversionError):
        gyroradius(5 * u.T, Vperp=8 * u.m)

    with pytest.raises(ValueError):
        gyroradius(np.array([5, 6]) * u.T, Vperp=np.array([5, 6, 7]) * u.m / u.s)

    with pytest.raises(ValueError):
        gyroradius(np.nan * u.T, Vperp=1 * u.m / u.s)

    with pytest.raises(ValueError):
        gyroradius(3.14159 * u.T, T_i=-1 * u.K)

    with pytest.warns(u.UnitsWarning):
        assert gyroradius(1.0, Vperp=1.0) == gyroradius(1.0 * u.T, Vperp=1.0 * u.m / u.s)

    with pytest.warns(u.UnitsWarning):
        assert gyroradius(1.1, T_i=1.2) == gyroradius(1.1 * u.T, T_i=1.2 * u.K)

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, Vperp=1 * u.m / u.s, T_i=1.2 * u.K)

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, Vperp=1.1 * u.m, T_i=1.2 * u.K)

    assert gyroradius(B, particle="p", T_i=T_i).unit.is_equivalent(u.m)

    assert gyroradius(B, particle="p", Vperp=25 * u.m / u.s).unit.is_equivalent(u.m)

    # Case when Z=1 is assumed
    assert np.isclose(gyroradius(B, particle='p', T_i=T_i),
                      gyroradius(B, particle='H+', T_i=T_i),
                      atol=1e-6 * u.m)

    gyroPos = gyroradius(B, particle="p", Vperp=V)
    gyroNeg = gyroradius(B, particle="p", Vperp=-V)
    assert gyroPos == gyroNeg

    Vperp = 1e6 * u.m / u.s
    Bmag = 1 * u.T
    omega_ci = gyrofrequency(Bmag, particle='p')
    analytical_result = (Vperp / omega_ci).to(u.m, equivalencies=u.dimensionless_angles())
    assert gyroradius(Bmag, particle="p", Vperp=Vperp) == analytical_result

    T2 = 1.2 * u.MK
    B2 = 123 * u.G
    particle2 = 'alpha'
    Vperp2 = thermal_speed(T2, particle=particle2)
    gyro_by_vperp = gyroradius(B2, particle='alpha', Vperp=Vperp2)
    assert gyro_by_vperp == gyroradius(B2, particle='alpha', T_i=T2)

    explicit_positron_gyro = gyroradius(1 * u.T, particle='positron', T_i=1 * u.MK)
    assert explicit_positron_gyro == gyroradius(1 * u.T, T_i=1 * u.MK)

    with pytest.raises(TypeError):
        gyroradius(u.T, particle="p", Vperp=8 * u.m / u.s)

    with pytest.raises(ValueError):
        gyroradius(B, particle='p', T_i=-1 * u.K)

    with pytest.warns(u.UnitsWarning):
        gyro_without_units = gyroradius(1.0, particle="p", Vperp=1.0)
        gyro_with_units = gyroradius(1.0 * u.T, particle="p", Vperp=1.0 * u.m / u.s)
        assert gyro_without_units == gyro_with_units

    with pytest.warns(u.UnitsWarning):
        gyro_t_without_units = gyroradius(1.1, particle="p", T_i=1.2)
        gyro_t_with_units = gyroradius(1.1 * u.T, particle="p", T_i=1.2 * u.K)
        assert gyro_t_with_units == gyro_t_without_units

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, particle="p", Vperp=1 * u.m / u.s, T_i=1.2 * u.K)

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, particle="p", Vperp=1.1 * u.m, T_i=1.2 * u.K)

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, particle="p", Vperp=1.2 * u.m, T_i=1.1 * u.K)


def test_plasma_frequency():
    r"""Test the plasma_frequency function in parameters.py."""

    assert plasma_frequency(n_e).unit.is_equivalent(u.rad / u.s)

    assert np.isclose(plasma_frequency(1 * u.cm ** -3).value, 5.64e4, rtol=1e-2)

    with pytest.raises(TypeError):
        plasma_frequency(u.m ** -3)

    with pytest.raises(u.UnitConversionError):
        plasma_frequency(5 * u.m ** -2)

    with pytest.raises(ValueError):
        plasma_frequency(np.nan * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert plasma_frequency(1e19) == plasma_frequency(1e19 * u.m ** -3)

        assert plasma_frequency(n_i, particle='p').unit.is_equivalent(u.rad / u.s)

    # Case where Z=1 is assumed
    assert plasma_frequency(n_i, particle='H-1+') == plasma_frequency(n_i, particle='p')

    assert np.isclose(plasma_frequency(mu * u.cm ** -3, particle='p').value,
                      1.32e3, rtol=1e-2)

    with pytest.raises(ValueError):
        plasma_frequency(n=5 * u.m ** -3, particle='sdfas')

    with pytest.warns(u.UnitsWarning):
        plasma_freq_no_units = plasma_frequency(1e19, particle='p')
        assert plasma_freq_no_units == plasma_frequency(1e19 * u.m ** -3, particle='p')

    plasma_frequency(1e17 * u.cm ** -3, particle='p')
    # testing for user input z_mean
    testMeth1 = plasma_frequency(1e17 * u.cm ** -3,
                                 particle='p',
                                 z_mean=0.8).si.value
    testTrue1 = 333063562455.4028
    errStr = f"plasma_frequency() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1,
                      testTrue1,
                      atol=0.0,
                      rtol=1e-15), errStr


def test_Debye_length():
    r"""Test the Debye_length function in parameters.py."""

    assert Debye_length(T_e, n_e).unit.is_equivalent(u.m)

    assert np.isclose(Debye_length(
        1 * u.eV, 1 * u.cm ** -3).value, 7.43, atol=0.005)

    with pytest.warns(u.UnitsWarning):
        Debye_length(5, 5 * u.m ** -3)

    with pytest.raises(u.UnitConversionError):
        Debye_length(56 * u.kg, 5 * u.m ** -3)

    with pytest.raises(ValueError):
        Debye_length(5 * u.eV, -5 * u.m ** -3)

    with pytest.raises(ValueError):
        Debye_length(-45 * u.K, 5 * u.m ** -3)

    Tarr2 = np.array([1, 2]) * u.K
    narr3 = np.array([1, 2, 3]) * u.m ** -3
    with pytest.raises(ValueError):
        Debye_length(Tarr2, narr3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_length(2.0, 2.0) == Debye_length(2.0 * u.K, 2.0 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_length(2.0 * u.K, 2.0) == Debye_length(2.0, 2.0 * u.m ** -3)


def test_Debye_number():
    r"""Test the Debye_number function in parameters.py."""

    assert Debye_number(T_e, n_e).unit.is_equivalent(u.dimensionless_unscaled)

    T_e_eV = T_e.to(u.eV, equivalencies=u.temperature_energy())
    assert np.isclose(Debye_number(T_e, n_e).value,
                      Debye_number(T_e_eV, n_e).value)

    assert np.isclose(Debye_number(1 * u.eV, 1 * u.cm ** -3).value, 1720862385.43342)

    with pytest.warns(u.UnitsWarning):
        Debye_number(T_e, 4)

    with pytest.raises(ValueError):
        Debye_number(None, n_e)

    with pytest.raises(u.UnitConversionError):
        Debye_number(5 * u.m, 5 * u.m ** -3)

    with pytest.raises(u.UnitConversionError):
        Debye_number(5 * u.K, 5 * u.m ** 3)

    with pytest.raises(ValueError):
        Debye_number(5j * u.K, 5 * u.cm ** -3)

    Tarr2 = np.array([1, 2]) * u.K
    narr3 = np.array([1, 2, 3]) * u.m ** -3
    with pytest.raises(ValueError):
        Debye_number(Tarr2, narr3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_number(1.1, 1.1) == Debye_number(1.1 * u.K, 1.1 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert Debye_number(1.1 * u.K, 1.1) == Debye_number(1.1, 1.1 * u.m ** -3)


def test_inertial_length():
    r"""Test the inertial_length function in parameters.py."""

    assert inertial_length(n_i, particle='p').unit.is_equivalent(u.m)

    assert np.isclose(inertial_length(mu * u.cm ** -3, particle='p').cgs.value,
                      2.28e7, rtol=0.01)

    inertial_length_electron_plus = inertial_length(5.351 * u.m ** -3, particle='e+')
    assert inertial_length_electron_plus == inertial_length(5.351 * u.m ** -3, particle='e')

    assert inertial_length(n_i, particle='p') == inertial_length(n_i, particle='p')

    with pytest.warns(u.UnitsWarning):
        inertial_length(4, particle='p')

    with pytest.raises(u.UnitConversionError):
        inertial_length(4 * u.m ** -2, particle='p')

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m ** -3, particle='p')

    with pytest.raises(ValueError):
        inertial_length(n_i, particle=-135)

    with pytest.warns(u.UnitsWarning):
        inertial_length_no_units = inertial_length(1e19, particle='p')
        assert inertial_length_no_units == inertial_length(1e19 * u.m ** -3, particle='p')

    assert inertial_length(n_e).unit.is_equivalent(u.m)

    assert np.isclose(inertial_length(1 * u.cm ** -3).cgs.value, 5.31e5, rtol=1e-3)

    with pytest.warns(u.UnitsWarning):
        inertial_length(5)

    with pytest.raises(u.UnitConversionError):
        inertial_length(5 * u.m)

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert inertial_length(1e19) == inertial_length(1e19 * u.m ** -3)


def test_magnetic_pressure():
    r"""Test the magnetic_pressure function in parameters.py."""

    assert magnetic_pressure(B_arr).unit.is_equivalent(u.Pa)

    assert magnetic_pressure(B).unit.is_equivalent(u.Pa)

    assert magnetic_pressure(B).unit.name == 'Pa'

    assert magnetic_pressure(B).value == magnetic_energy_density(B).value

    assert magnetic_pressure(B) == magnetic_energy_density(B.to(u.G))

    assert np.isclose(magnetic_pressure(B).value, 397887.35772973835)

    with pytest.warns(u.UnitsWarning):
        magnetic_pressure(5)

    with pytest.raises(u.UnitConversionError):
        magnetic_pressure(5 * u.m)

    with pytest.raises(ValueError):
        magnetic_pressure(np.nan * u.T)

    with pytest.raises(ValueError):
        magnetic_pressure(5j * u.T)

    with pytest.raises(ValueError):
        magnetic_pressure(B_nanarr)

    with pytest.warns(u.UnitsWarning):
        assert magnetic_pressure(22.2) == magnetic_pressure(22.2 * u.T)


def test_magnetic_energy_density():
    r"""Test the magnetic_energy_density function in parameters.py."""

    assert magnetic_energy_density(B_arr).unit.is_equivalent(u.J / u.m ** 3)

    assert magnetic_energy_density(B).unit.is_equivalent('J / m3')

    assert magnetic_energy_density(B).value == magnetic_pressure(B).value

    assert_quantity_allclose(magnetic_energy_density(2 * B), 4 * magnetic_energy_density(B))

    assert_quantity_allclose(magnetic_energy_density(B).value, 397887.35772973835)

    assert_quantity_allclose(magnetic_energy_density(B), magnetic_energy_density(B.to(u.G)))

    # TODO Add an array test!

    assert magnetic_energy_density(B_arr)

    with pytest.warns(u.UnitsWarning):
        magnetic_energy_density(5)

    with pytest.raises(u.UnitConversionError):
        magnetic_energy_density(5 * u.m)

    with pytest.raises(ValueError):
        magnetic_energy_density(np.nan * u.T)

    with pytest.raises(ValueError):
        magnetic_energy_density(5j * u.T)

    with pytest.raises(ValueError):
        magnetic_energy_density(B_nanarr)

    with pytest.warns(u.UnitsWarning):
        assert magnetic_energy_density(22.2) == magnetic_energy_density(22.2 * u.T)


def test_upper_hybrid_frequency():
    r"""Test the upper_hybrid_frequency function in parameters.py."""

    omega_uh = upper_hybrid_frequency(B, n_e=n_e)
    omega_ce = gyrofrequency(B)
    omega_pe = plasma_frequency(n=n_e)
    assert omega_ce.unit.is_equivalent(u.rad / u.s)
    assert omega_pe.unit.is_equivalent(u.rad / u.s)
    assert omega_uh.unit.is_equivalent(u.rad / u.s)
    left_hand_side = omega_uh ** 2
    right_hand_side = omega_ce ** 2 + omega_pe ** 2
    assert np.isclose(left_hand_side.value, right_hand_side.value)

    with pytest.raises(ValueError):
        upper_hybrid_frequency(5 * u.T, n_e=-1 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert upper_hybrid_frequency(1.2, 1.3) == upper_hybrid_frequency(1.2 * u.T,
                                                                          1.3 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert upper_hybrid_frequency(1.4 * u.T, 1.3) == upper_hybrid_frequency(1.4,
                                                                                1.3 * u.m ** -3)


def test_lower_hybrid_frequency():
    r"""Test the lower_hybrid_frequency function in parameters.py."""

    ion = 'He-4 1+'
    omega_ci = gyrofrequency(B, particle=ion)
    omega_pi = plasma_frequency(n=n_i, particle=ion)
    omega_ce = gyrofrequency(B)
    omega_lh = lower_hybrid_frequency(B, n_i=n_i, ion=ion)
    assert omega_ci.unit.is_equivalent(u.rad / u.s)
    assert omega_pi.unit.is_equivalent(u.rad / u.s)
    assert omega_ce.unit.is_equivalent(u.rad / u.s)
    assert omega_lh.unit.is_equivalent(u.rad / u.s)
    left_hand_side = omega_lh ** -2
    right_hand_side = 1 / (omega_ci ** 2 + omega_pi ** 2) + omega_ci ** -1 * omega_ce ** -1
    assert np.isclose(left_hand_side.value, right_hand_side.value)

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=5e19 * u.m ** -3, ion='asdfasd')

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=-5e19 * u.m ** -3, ion='asdfasd')

    with pytest.raises(ValueError):
        lower_hybrid_frequency(
            np.nan * u.T, n_i=-5e19 * u.m ** -3, ion='asdfasd')

    with pytest.warns(u.UnitsWarning):
        assert lower_hybrid_frequency(1.3, 1e19) == lower_hybrid_frequency(1.3 * u.T,
                                                                           1e19 * u.m ** -3)
