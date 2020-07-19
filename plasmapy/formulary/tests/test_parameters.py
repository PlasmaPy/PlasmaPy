"""Tests for functions that calculate plasma parameters."""

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c, e, m_e, m_p, mu0
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.parameters import (
    Alfven_speed,
    Bohm_diffusion,
    betaH_,
    cs_,
    cwp_,
    DB_,
    Debye_length,
    Debye_number,
    gyrofrequency,
    gyroradius,
    Hall_parameter,
    inertial_length,
    ion_sound_speed,
    kappa_thermal_speed,
    lambdaD_,
    lower_hybrid_frequency,
    magnetic_energy_density,
    magnetic_pressure,
    mass_density,
    nD_,
    oc_,
    plasma_frequency,
    pmag_,
    pth_,
    rc_,
    rho_,
    rhoc_,
    thermal_pressure,
    thermal_speed,
    ub_,
    upper_hybrid_frequency,
    va_,
    vth_,
    vth_kappa_,
    wc_,
    wp_,
    wlh_,
    wuh_,
)
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.exceptions import (
    PhysicsError,
    PhysicsWarning,
    RelativityError,
    RelativityWarning,
)
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

B = 1.0 * u.T
Z = 1
ion = "p"
m_i = m_p
n_i = 5e19 * u.m ** -3
n_e = Z * 5e19 * u.m ** -3
rho = n_i * m_i + n_e * m_e
T_e = 1e6 * u.K
T_i = 1e6 * u.K
k_1 = 3e1 * u.m ** -1
k_2 = 3e7 * u.m ** -1

B_arr = np.array([0.001, 0.002]) * u.T
B_nanarr = np.array([0.001, np.nan]) * u.T
B_allnanarr = np.array([np.nan, np.nan]) * u.T

rho_arr = np.array([5e-10, 2e-10]) * u.kg / u.m ** 3
rho_infarr = np.array([np.inf, 5e19]) * u.m ** -3
rho_negarr = np.array([-5e19, 6e19]) * u.m ** -3

T_arr = np.array([1e6, 2e6]) * u.K
T_nanarr = np.array([1e6, np.nan]) * u.K
T_nanarr2 = np.array([np.nan, 2e6]) * u.K
T_allnanarr = np.array([np.nan, np.nan]) * u.K
T_negarr = np.array([1e6, -5151.0]) * u.K

V = 25.2 * u.m / u.s
V_arr = np.array([25, 50]) * u.m / u.s
V_nanarr = np.array([25, np.nan]) * u.m / u.s
V_allnanarr = np.array([np.nan, np.nan]) * u.m / u.s

mu = m_p.to(u.u).value


class Test_mass_density:
    r"""Test the mass_density function in parameters.py."""

    def test_particleless(self):
        with pytest.raises(ValueError):
            mass_density(1 * u.m ** -3)

    def test_wrong_units(self):
        with pytest.raises(u.UnitTypeError):
            mass_density(1 * u.J)

    def test_handle_nparrays(self):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(mass_density)


# Assertions below that are in CGS units with 2-3 significant digits
# are generally from the NRL Plasma Formulary.


def test_Alfven_speed():
    r"""Test the Alfven_speed function in parameters.py."""

    # TODO: break this test up until multiple tests

    assert np.isclose(
        Alfven_speed(1 * u.T, 1e-8 * u.kg * u.m ** -3).value,
        8920620.580763856,
        rtol=1e-6,
    )

    V_A = Alfven_speed(B, n_i)
    assert np.isclose(V_A.value, (B / np.sqrt(mu0 * n_i * (m_p + m_e))).si.value)

    assert Alfven_speed(B, rho) == Alfven_speed(B, n_i)

    assert Alfven_speed(B, rho).unit.is_equivalent(u.m / u.s)

    assert Alfven_speed(B, rho) == Alfven_speed(-B, rho)

    assert Alfven_speed(B, 4 * rho) == 0.5 * Alfven_speed(B, rho)

    assert Alfven_speed(2 * B, rho) == 2 * Alfven_speed(B, rho)

    # Case when Z=1 is assumed
    with pytest.warns(RelativityWarning):
        assert np.isclose(
            Alfven_speed(5 * u.T, 5e19 * u.m ** -3, ion="H+"),
            Alfven_speed(5 * u.T, 5e19 * u.m ** -3, ion="p"),
            atol=0 * u.m / u.s,
            rtol=1e-3,
        )

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
        Alfven_speed(np.array([5, 6, 7]) * u.T, np.array([5, 6]) * u.m ** -3)

    assert np.isnan(Alfven_speed(B_nanarr, rho_arr)[1])

    with pytest.raises(ValueError):
        Alfven_speed(B_arr, rho_negarr)

    with pytest.raises(u.UnitTypeError):
        Alfven_speed(5 * u.A, n_i, ion="p")

    with pytest.raises(TypeError):
        Alfven_speed(B, 5, ion="p")

    with pytest.raises(u.UnitsError):
        Alfven_speed(B, 5 * u.m ** -2, ion="p")

    with pytest.raises(InvalidParticleError):
        Alfven_speed(B, n_i, ion="spacecats")

    with pytest.warns(RelativityWarning):  # relativistic
        Alfven_speed(5e1 * u.T, 5e19 * u.m ** -3, ion="p")

    with pytest.raises(RelativityError):  # super-relativistic
        Alfven_speed(5e8 * u.T, 5e19 * u.m ** -3, ion="p")

    with pytest.raises(ValueError):
        Alfven_speed(0.001 * u.T, -5e19 * u.m ** -3, ion="p")

    assert np.isnan(Alfven_speed(np.nan * u.T, 1 * u.m ** -3, ion="p"))

    assert np.isnan(Alfven_speed(1 * u.T, np.nan * u.m ** -3, ion="p"))

    with pytest.raises(RelativityError):
        assert Alfven_speed(np.inf * u.T, 1 * u.m ** -3, ion="p") == np.inf * u.m / u.s

    with pytest.raises(RelativityError):
        assert Alfven_speed(-np.inf * u.T, 1 * u.m ** -3, ion="p") == np.inf * u.m / u.s

    with pytest.warns(u.UnitsWarning):
        assert Alfven_speed(1.0, n_i) == Alfven_speed(1.0 * u.T, n_i)

    Alfven_speed(1 * u.T, 5e19 * u.m ** -3, ion="p")
    # testing for user input z_mean
    testMeth1 = Alfven_speed(1 * u.T, 5e19 * u.m ** -3, ion="p", z_mean=0.8).si.value
    testTrue1 = 3084015.75214846
    errStr = f"Alfven_speed() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1, testTrue1, atol=0.0, rtol=1e-6), errStr

    assert_can_handle_nparray(Alfven_speed)


def test_ion_sound_speed():
    r"""Test the ion_sound_speed function in parameters.py."""

    assert np.isclose(
        ion_sound_speed(
            T_i=1.3232 * u.MK, T_e=1.831 * u.MK, ion="p", gamma_e=1, gamma_i=3
        ).value,
        218816.06086407552,
    )

    assert np.isclose(
        ion_sound_speed(
            T_i=1.3232 * u.MK,
            T_e=1.831 * u.MK,
            n_e=n_e,
            k=k_1,
            ion="p",
            gamma_e=1,
            gamma_i=3,
        ).value,
        218816.06086407552,
    )

    assert np.isclose(
        ion_sound_speed(
            T_i=1.3232 * u.MK,
            T_e=1.831 * u.MK,
            n_e=n_e,
            k=k_2,
            ion="p",
            gamma_e=1,
            gamma_i=3,
        ).value,
        552.3212936293337,
    )

    assert np.isclose(
        ion_sound_speed(
            T_i=0.88 * u.MK,
            T_e=1.28 * u.MK,
            n_e=n_e,
            k=0 * u.m ** -1,
            ion="p",
            gamma_e=1.2,
            gamma_i=3.4,
        ).value,
        193328.52857788358,
    )

    # case when Z=1 is assumed
    # assert ion_sound_speed(T_i=T_i, T_e=T_e, ion='p+') == ion_sound_speed(T_i=T_i, T_e=T_e,
    # ion='H-1')

    assert ion_sound_speed(
        T_i=T_i, T_e=0 * u.K, n_e=n_e, k=k_1, ion="p+"
    ).unit.is_equivalent(u.m / u.s)

    with pytest.raises(RelativityError):
        ion_sound_speed(T_i=T_i, T_e=T_e, n_e=n_e, k=k_1, gamma_i=np.inf)

    with pytest.warns(PhysicsWarning):
        ion_sound_speed(T_i=T_i, T_e=T_e, n_e=n_e)

    with pytest.warns(PhysicsWarning):
        ion_sound_speed(T_i=T_i, T_e=T_e, k=k_1)

    with pytest.raises(u.UnitTypeError):
        ion_sound_speed(
            T_i=np.array([5, 6, 5]) * u.K,
            T_e=np.array([3, 4]) * u.K,
            n_e=np.array([5, 6, 5]) * u.m ** -3,
            k=np.array([3, 4]) * u.m ** -3,
        )

    with pytest.raises(TypeError):  # Is this test right??????
        ion_sound_speed(5 * u.T)

    with pytest.raises(TypeError):
        ion_sound_speed("p")

    with pytest.raises(PhysicsError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, gamma_i=0.9999)

    with pytest.raises(PhysicsError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, gamma_e=0.9999)

    with pytest.raises(TypeError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, gamma_e="sdjklsf")

    with pytest.raises(TypeError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, gamma_i="fsdfas")

    with pytest.raises(InvalidParticleError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, ion="cupcakes")

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=-np.abs(T_i), T_e=0 * u.K)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, n_e=-np.abs(n_e), k=k_1)

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_i, T_e=0 * u.K, n_e=n_e, k=-np.abs(k_1))

    with pytest.warns(RelativityWarning):
        ion_sound_speed(T_i=5e11 * u.K, T_e=0 * u.K)

    with pytest.raises(RelativityError):
        ion_sound_speed(T_i=5e19 * u.K, T_e=0 * u.K)

    with pytest.raises(u.UnitTypeError):
        ion_sound_speed(T_i=5 * u.A, T_e=0 * u.K, n_e=n_e, k=k_1)

    assert np.isnan(ion_sound_speed(T_i=T_nanarr, T_e=0 * u.K, n_e=n_e, k=k_1)[1])

    assert np.isnan(ion_sound_speed(T_e=T_nanarr, T_i=0 * u.K, n_e=n_e, k=k_1)[1])

    with pytest.raises(ValueError):
        ion_sound_speed(T_i=T_negarr, T_e=0 * u.K, n_e=n_e, k=k_1)

    with pytest.raises(ValueError):
        ion_sound_speed(T_e=T_negarr, T_i=0 * u.K, n_e=n_e, k=k_1)

    with pytest.warns(u.UnitsWarning):
        assert ion_sound_speed(
            T_e=1.2e6, T_i=0 * u.K, n_e=n_e, k=k_1
        ) == ion_sound_speed(T_e=1.2e6 * u.K, T_i=0 * u.K, n_e=n_e, k=k_1)

    with pytest.warns(u.UnitsWarning):
        assert ion_sound_speed(
            T_i=1.3e6, T_e=0 * u.K, n_e=n_e, k=k_1
        ) == ion_sound_speed(T_i=1.3e6 * u.K, T_e=0 * u.K, n_e=n_e, k=k_1)

    ion_sound_speed(T_e=1.2e6 * u.K, T_i=0 * u.K, n_e=n_e, k=k_1)
    # testing for user input z_mean
    testMeth1 = ion_sound_speed(
        T_e=1.2e6 * u.K, T_i=0 * u.K, n_e=n_e, k=0 * u.m ** -1, z_mean=0.8
    ).si.value
    testTrue1 = 89018.09
    errStr = f"ion_sound_speed() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1, testTrue1, atol=0.0, rtol=1e-6), errStr

    assert_can_handle_nparray(ion_sound_speed)


def test_thermal_speed():
    r"""Test the thermal_speed function in parameters.py"""
    assert thermal_speed(T_e).unit.is_equivalent(u.m / u.s)

    assert thermal_speed(T_e) > thermal_speed(T_e, "p")

    # The NRL Plasma Formulary uses a definition of the electron
    # thermal speed that differs by a factor of sqrt(2).
    assert np.isclose(thermal_speed(1 * u.MK).value, 5505694.743141063)

    with pytest.raises(u.UnitTypeError):
        thermal_speed(5 * u.m)

    with pytest.raises(ValueError):
        thermal_speed(-T_e)

    with pytest.warns(RelativityWarning):
        thermal_speed(1e9 * u.K)

    with pytest.raises(RelativityError):
        thermal_speed(5e19 * u.K)

    with pytest.warns(u.UnitsWarning):
        assert thermal_speed(1e5) == thermal_speed(1e5 * u.K)

    assert thermal_speed(T_i, particle="p").unit.is_equivalent(u.m / u.s)

    # The NRL Plasma Formulary uses a definition of the particle thermal
    # speed that differs by a factor of sqrt(2).
    assert np.isclose(
        thermal_speed(1 * u.MK, particle="p").si.value, 128486.56960876315
    )

    # Explicitly check all three modes and dimensionalities
    # ndim = 1
    assert np.isclose(thermal_speed(T_e, method="most_probable", ndim=1).si.value, 0.0)

    # Regression tests start here!
    assert np.isclose(
        thermal_speed(T_e, method="rms", ndim=1).si.value, 3893114.2008620175
    )

    assert np.isclose(
        thermal_speed(T_e, method="mean_magnitude", ndim=1).si.value, 3106255.714310189
    )

    # ndim = 2
    assert np.isclose(
        thermal_speed(T_e, method="most_probable", ndim=2).si.value, 3893114.2008620175
    )

    assert np.isclose(
        thermal_speed(T_e, method="rms", ndim=2).si.value, 5505694.902726359
    )

    assert np.isclose(
        thermal_speed(T_e, method="mean_magnitude", ndim=2).si.value, 4879295.066124102
    )

    # ndim = 3
    assert np.isclose(
        thermal_speed(T_e, method="most_probable", ndim=3).si.value, 5505694.902726359
    )

    assert np.isclose(
        thermal_speed(T_e, method="rms", ndim=3).si.value, 6743071.595560921
    )

    assert np.isclose(
        thermal_speed(T_e, method="mean_magnitude", ndim=3).si.value, 6212511.428620378
    )

    # Case when Z=1 is assumed
    assert thermal_speed(T_i, particle="p") == thermal_speed(T_i, particle="H-1+")

    assert thermal_speed(1 * u.MK, particle="e+") == thermal_speed(1 * u.MK)

    with pytest.raises(u.UnitTypeError):
        thermal_speed(5 * u.m, particle="p")

    with pytest.raises(ValueError):
        thermal_speed(-T_e, particle="p")

    with pytest.warns(RelativityWarning):
        thermal_speed(1e11 * u.K, particle="p")

    with pytest.raises(RelativityError):
        thermal_speed(1e14 * u.K, particle="p")

    with pytest.raises(InvalidParticleError):
        thermal_speed(T_i, particle="asdfasd")

    with pytest.warns(u.UnitsWarning):
        assert thermal_speed(1e6, particle="p") == thermal_speed(
            1e6 * u.K, particle="p"
        )

    assert np.isclose(
        thermal_speed(1e6 * u.K, method="mean_magnitude").si.value, 6212510.3969422
    )

    assert np.isclose(
        thermal_speed(1e6 * u.K, method="rms").si.value, 6743070.475775486
    )

    # Test invalid method
    with pytest.raises(ValueError):
        thermal_speed(T_i, method="sadks")

    # Test invalid ndim
    with pytest.raises(ValueError):
        thermal_speed(T_i, ndim=4)

    assert_can_handle_nparray(thermal_speed)


def test_thermal_pressure():
    assert thermal_pressure(T_e, n_i).unit.is_equivalent(u.Pa)

    # TODO: may be array issues with arg "mass"
    assert_can_handle_nparray(thermal_pressure)


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
            kappa_thermal_speed(self.T_e, self.kappaInvalid, particle=self.particle)
        return

    def test_invalid_method(self):
        """
        Checks if function raises error when invalid method is passed as an
        argument.
        """
        with pytest.raises(ValueError):
            kappa_thermal_speed(
                self.T_e, self.kappa, particle=self.particle, method="invalid"
            )
        return

    def test_probable1(self):
        """
        Tests if expected value is returned for a set of regular inputs.
        """
        known1 = kappa_thermal_speed(
            self.T_e, self.kappa, particle=self.particle, method="most_probable"
        )
        errStr = (
            f"Kappa thermal velocity should be {self.probable1True} "
            f"and not {known1.si.value}."
        )
        assert np.isclose(known1.value, self.probable1True, rtol=1e-8, atol=0.0), errStr
        return

    def test_rms1(self):
        """
        Tests if expected value is returned for a set of regular inputs.
        """
        known1 = kappa_thermal_speed(
            self.T_e, self.kappa, particle=self.particle, method="rms"
        )
        errStr = (
            f"Kappa thermal velocity should be {self.rms1True} "
            f"and not {known1.si.value}."
        )
        assert np.isclose(known1.value, self.rms1True, rtol=1e-8, atol=0.0), errStr
        return

    def test_mean1(self):
        """
        Tests if expected value is returned for a set of regular inputs.
        """
        known1 = kappa_thermal_speed(
            self.T_e, self.kappa, particle=self.particle, method="mean_magnitude"
        )
        errStr = (
            f"Kappa thermal velocity should be {self.mean1True} "
            f"and not {known1.si.value}."
        )
        assert np.isclose(known1.value, self.mean1True, rtol=1e-8, atol=0.0), errStr
        return

    def test_handle_nparrays(self, kwargs=None):
        """Test for ability to handle numpy array quantities"""
        if kwargs is None:
            kwargs = {"kappa": 2}
        assert_can_handle_nparray(kappa_thermal_speed, kwargs=kwargs)


def test_gyrofrequency():
    r"""Test the gyrofrequency function in parameters.py."""

    assert gyrofrequency(B).unit.is_equivalent(u.rad / u.s)

    assert gyrofrequency(B, to_hz=True).unit.is_equivalent(u.Hz)

    assert np.isclose(gyrofrequency(1 * u.T).value, 175882008784.72018)

    assert np.isclose(gyrofrequency(2.4 * u.T).value, 422116821083.3284)

    assert np.isclose(gyrofrequency(1 * u.T, to_hz=True).value, 27992490076.528206)

    assert np.isclose(gyrofrequency(2.4 * u.T, signed=True).value, -422116821083.3284)

    assert np.isclose(gyrofrequency(1 * u.G).cgs.value, 1.76e7, rtol=1e-3)

    with pytest.raises(TypeError):
        gyrofrequency(u.m)

    with pytest.raises(u.UnitTypeError):
        gyrofrequency(u.m * 1)

    assert np.isnan(gyrofrequency(B_nanarr)[-1])

    # The following is a test to check that equivalencies from astropy
    # are working.
    omega_ce = gyrofrequency(2.2 * u.T)
    f_ce = (omega_ce / (2 * np.pi)) / u.rad
    f_ce_use_equiv = omega_ce.to(u.Hz, equivalencies=[(u.cy / u.s, u.Hz)])
    assert np.isclose(f_ce.value, f_ce_use_equiv.value)

    with pytest.warns(u.UnitsWarning):
        assert gyrofrequency(5.0) == gyrofrequency(5.0 * u.T)

    assert gyrofrequency(B, particle=ion).unit.is_equivalent(u.rad / u.s)

    assert np.isclose(gyrofrequency(1 * u.T, particle="p").value, 95788335.834874)

    assert np.isclose(gyrofrequency(2.4 * u.T, particle="p").value, 229892006.00369796)

    assert np.isclose(gyrofrequency(1 * u.G, particle="p").cgs.value, 9.58e3, rtol=2e-3)

    assert gyrofrequency(-5 * u.T, "p") == gyrofrequency(5 * u.T, "p")

    # Case when Z=1 is assumed
    # assert gyrofrequency(B, particle='p+') == gyrofrequency(B, particle='H-1')

    assert gyrofrequency(B, particle="e+") == gyrofrequency(B)

    with pytest.warns(u.UnitsWarning):
        gyrofrequency(8, "p")

    with pytest.raises(u.UnitTypeError):
        gyrofrequency(5 * u.m, "p")

    with pytest.raises(InvalidParticleError):
        gyrofrequency(8 * u.T, particle="asdfasd")

    with pytest.warns(u.UnitsWarning):
        # TODO this should be WARNS, not RAISES. and it's probably still raised
        assert gyrofrequency(5.0, "p") == gyrofrequency(5.0 * u.T, "p")

    gyrofrequency(1 * u.T, particle="p")
    # testing for user input Z
    testMeth1 = gyrofrequency(1 * u.T, particle="p", Z=0.8).si.value
    testTrue1 = 76630665.79318453
    errStr = f"gyrofrequency() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1, testTrue1, atol=0.0, rtol=1e-5), errStr

    assert_can_handle_nparray(gyrofrequency, kwargs={"signed": True})

    assert_can_handle_nparray(gyrofrequency, kwargs={"signed": False})


def test_gyroradius():
    r"""Test the gyroradius function in parameters.py."""

    assert gyroradius(B, T_i=T_e).unit.is_equivalent(u.m)

    assert gyroradius(B, Vperp=25 * u.m / u.s).unit.is_equivalent(u.m)

    Vperp = 1e6 * u.m / u.s
    Bmag = 1 * u.T
    omega_ce = gyrofrequency(Bmag)
    analytical_result = (Vperp / omega_ce).to(
        u.m, equivalencies=u.dimensionless_angles()
    )
    assert gyroradius(Bmag, Vperp=Vperp) == analytical_result

    with pytest.raises(TypeError):
        gyroradius(u.T)

    with pytest.raises(u.UnitTypeError):
        gyroradius(5 * u.A, Vperp=8 * u.m / u.s)

    with pytest.raises(u.UnitTypeError):
        gyroradius(5 * u.T, Vperp=8 * u.m)

    with pytest.raises(ValueError):
        gyroradius(np.array([5, 6]) * u.T, Vperp=np.array([5, 6, 7]) * u.m / u.s)

    assert np.isnan(gyroradius(np.nan * u.T, Vperp=1 * u.m / u.s))

    with pytest.raises(ValueError):
        gyroradius(3.14159 * u.T, T_i=-1 * u.K)

    with pytest.warns(u.UnitsWarning):
        assert gyroradius(1.0, Vperp=1.0) == gyroradius(
            1.0 * u.T, Vperp=1.0 * u.m / u.s
        )

    with pytest.warns(u.UnitsWarning):
        assert gyroradius(1.1, T_i=1.2) == gyroradius(1.1 * u.T, T_i=1.2 * u.K)

    with pytest.raises(ValueError):
        gyroradius(1.1 * u.T, Vperp=1 * u.m / u.s, T_i=1.2 * u.K)

    with pytest.raises(u.UnitTypeError):
        gyroradius(1.1 * u.T, Vperp=1.1 * u.m, T_i=1.2 * u.K)

    assert gyroradius(B, particle="p", T_i=T_i).unit.is_equivalent(u.m)

    assert gyroradius(B, particle="p", Vperp=25 * u.m / u.s).unit.is_equivalent(u.m)

    # Case when Z=1 is assumed
    assert np.isclose(
        gyroradius(B, particle="p", T_i=T_i),
        gyroradius(B, particle="H+", T_i=T_i),
        atol=1e-6 * u.m,
    )

    gyroPos = gyroradius(B, particle="p", Vperp=V)
    gyroNeg = gyroradius(B, particle="p", Vperp=-V)
    assert gyroPos == gyroNeg

    Vperp = 1e6 * u.m / u.s
    Bmag = 1 * u.T
    omega_ci = gyrofrequency(Bmag, particle="p")
    analytical_result = (Vperp / omega_ci).to(
        u.m, equivalencies=u.dimensionless_angles()
    )
    assert gyroradius(Bmag, particle="p", Vperp=Vperp) == analytical_result

    T2 = 1.2 * u.MK
    B2 = 123 * u.G
    particle2 = "alpha"
    Vperp2 = thermal_speed(T2, particle=particle2)
    gyro_by_vperp = gyroradius(B2, particle="alpha", Vperp=Vperp2)
    assert gyro_by_vperp == gyroradius(B2, particle="alpha", T_i=T2)

    explicit_positron_gyro = gyroradius(1 * u.T, particle="positron", T_i=1 * u.MK)
    assert explicit_positron_gyro == gyroradius(1 * u.T, T_i=1 * u.MK)

    with pytest.raises(TypeError):
        gyroradius(u.T, particle="p", Vperp=8 * u.m / u.s)

    with pytest.raises(ValueError):
        gyroradius(B, particle="p", T_i=-1 * u.K)

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

    with pytest.raises(u.UnitTypeError):
        gyroradius(1.1 * u.T, particle="p", Vperp=1.1 * u.m, T_i=1.2 * u.K)

    with pytest.raises(u.UnitTypeError):
        gyroradius(1.1 * u.T, particle="p", Vperp=1.2 * u.m, T_i=1.1 * u.K)


class Test_gyroradius:

    # some custom numpy array tests here, because of the T_i / Vperp situation
    def test_handle_numpy_array(self):
        # Tests to verify that can handle Quantities with numpy array as the value:
        assert gyroradius(B_arr, Vperp=V_arr)[0] == gyroradius(B_arr[0], Vperp=V_arr[0])
        assert gyroradius(B_arr, T_i=T_arr)[0] == gyroradius(B_arr[0], T_i=T_arr[0])

    def test_handle_mixed_Qarrays(self):
        # If both Vperp or Ti are input as Qarrays, but only one of the two is valid
        # at each element, then that's fine, the function should work:
        assert gyroradius(B_arr, Vperp=V_nanarr, T_i=T_nanarr2)[0] == gyroradius(
            B_arr[0], Vperp=V_nanarr[0], T_i=T_nanarr2[0]
        )

    def test_raise_two_valid_inputs(self):
        # If both Vperp or Ti are nan-less, Qarrays or not, should raise ValueError:
        with pytest.raises(ValueError):
            gyroradius(B_arr, Vperp=V, T_i=T_arr)
        with pytest.raises(ValueError):
            gyroradius(B_arr, Vperp=V_arr, T_i=T_i)

    def test_all_valid_and_one_valid(self):
        # If one of (Vperp, Ti) is a valid and one is Qarray with at least one valid, ValueError:
        with pytest.raises(ValueError):
            gyroradius(B_arr, Vperp=V, T_i=T_nanarr)
        with pytest.raises(ValueError):
            gyroradius(B_arr, Vperp=V_nanarr, T_i=T_i)

    def test_scalar_and_nan_qarray(self):
        # If either Vperp or Ti is a valid scalar and the other is a Qarray of all nans,
        # should do something valid and not raise a ValueError
        assert np.all(np.isfinite(gyroradius(B_arr, Vperp=V, T_i=T_allnanarr)))
        assert np.all(np.isfinite(gyroradius(B_arr, Vperp=V_allnanarr, T_i=T_i)))

    def test_keeps_arguments_unchanged(self):
        Vperp1 = u.Quantity([np.nan, 1], unit=u.m / u.s)
        Vperp2 = u.Quantity([np.nan, 1], unit=u.m / u.s)  # an exact copy
        T_i = u.Quantity([1, np.nan], unit=u.K)

        gyroradius(B_arr, Vperp=Vperp1, T_i=T_i)
        assert_quantity_allclose(Vperp1, Vperp2)


def test_plasma_frequency():
    r"""Test the plasma_frequency function in parameters.py."""

    assert plasma_frequency(n_e).unit.is_equivalent(u.rad / u.s)

    assert plasma_frequency(n_e, to_hz=True).unit.is_equivalent(u.Hz)

    assert np.isclose(plasma_frequency(1 * u.cm ** -3).value, 5.64e4, rtol=1e-2)

    assert np.isclose(
        plasma_frequency(1 * u.cm ** -3, particle="N").value, 3.53e2, rtol=1e-1
    )

    assert np.isclose(
        plasma_frequency(1 * u.cm ** -3, particle="N", to_hz=True).value,
        56.19000195094519,
    )

    with pytest.raises(TypeError):
        plasma_frequency(u.m ** -3)

    with pytest.raises(u.UnitTypeError):
        plasma_frequency(5 * u.m ** -2)

    assert np.isnan(plasma_frequency(np.nan * u.m ** -3))

    with pytest.warns(u.UnitsWarning):
        assert plasma_frequency(1e19) == plasma_frequency(1e19 * u.m ** -3)

        assert plasma_frequency(n_i, particle="p").unit.is_equivalent(u.rad / u.s)

    # Case where Z=1 is assumed
    assert plasma_frequency(n_i, particle="H-1+") == plasma_frequency(n_i, particle="p")

    assert np.isclose(
        plasma_frequency(mu * u.cm ** -3, particle="p").value, 1.32e3, rtol=1e-2
    )

    with pytest.raises(ValueError):
        plasma_frequency(n=5 * u.m ** -3, particle="sdfas")

    with pytest.warns(u.UnitsWarning):
        plasma_freq_no_units = plasma_frequency(1e19, particle="p")
        assert plasma_freq_no_units == plasma_frequency(1e19 * u.m ** -3, particle="p")

    plasma_frequency(1e17 * u.cm ** -3, particle="p")
    # testing for user input z_mean
    testMeth1 = plasma_frequency(1e17 * u.cm ** -3, particle="p", z_mean=0.8).si.value
    testTrue1 = 333063562455.4028
    errStr = f"plasma_frequency() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1, testTrue1, atol=0.0, rtol=1e-6), errStr

    assert_can_handle_nparray(plasma_frequency)


def test_Debye_length():
    r"""Test the Debye_length function in parameters.py."""

    assert Debye_length(T_e, n_e).unit.is_equivalent(u.m)

    assert np.isclose(Debye_length(1 * u.eV, 1 * u.cm ** -3).value, 7.43, atol=0.005)

    with pytest.warns(u.UnitsWarning):
        Debye_length(5, 5 * u.m ** -3)

    with pytest.raises(u.UnitTypeError):
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

    assert_can_handle_nparray(Debye_length)


def test_Debye_number():
    r"""Test the Debye_number function in parameters.py."""

    assert Debye_number(T_e, n_e).unit.is_equivalent(u.dimensionless_unscaled)

    T_e_eV = T_e.to(u.eV, equivalencies=u.temperature_energy())
    assert np.isclose(Debye_number(T_e, n_e).value, Debye_number(T_e_eV, n_e).value)

    assert np.isclose(Debye_number(1 * u.eV, 1 * u.cm ** -3).value, 1720862385.43342)

    with pytest.warns(u.UnitsWarning):
        Debye_number(T_e, 4)

    with pytest.raises(ValueError):
        Debye_number(None, n_e)

    with pytest.raises(u.UnitTypeError):
        Debye_number(5 * u.m, 5 * u.m ** -3)

    with pytest.raises(u.UnitTypeError):
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

    assert_can_handle_nparray(Debye_number)


def test_inertial_length():
    r"""Test the inertial_length function in parameters.py."""

    assert inertial_length(n_i, particle="p").unit.is_equivalent(u.m)

    assert np.isclose(
        inertial_length(mu * u.cm ** -3, particle="p").cgs.value, 2.28e7, rtol=0.01
    )

    inertial_length_electron_plus = inertial_length(5.351 * u.m ** -3, particle="e+")
    assert inertial_length_electron_plus == inertial_length(
        5.351 * u.m ** -3, particle="e"
    )

    assert inertial_length(n_i, particle="p") == inertial_length(n_i, particle="p")

    with pytest.warns(u.UnitsWarning):
        inertial_length(4, particle="p")

    with pytest.raises(u.UnitTypeError):
        inertial_length(4 * u.m ** -2, particle="p")

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m ** -3, particle="p")

    with pytest.raises(InvalidParticleError):
        inertial_length(n_i, particle=-135)

    with pytest.warns(u.UnitsWarning):
        inertial_length_no_units = inertial_length(1e19, particle="p")
        assert inertial_length_no_units == inertial_length(
            1e19 * u.m ** -3, particle="p"
        )

    assert inertial_length(n_e, "e-").unit.is_equivalent(u.m)

    assert np.isclose(
        inertial_length(1 * u.cm ** -3, "e-").cgs.value, 5.31e5, rtol=1e-3
    )

    with pytest.warns(u.UnitsWarning):
        inertial_length(5, "e-")

    with pytest.raises(u.UnitTypeError):
        inertial_length(5 * u.m, "e-")

    with pytest.raises(ValueError):
        inertial_length(-5 * u.m ** -3, "e-")

    with pytest.warns(u.UnitsWarning):
        assert inertial_length(1e19, "e-") == inertial_length(1e19 * u.m ** -3, "e-")

    assert_can_handle_nparray(inertial_length)


def test_magnetic_pressure():
    r"""Test the magnetic_pressure function in parameters.py."""

    assert magnetic_pressure(B_arr).unit.is_equivalent(u.Pa)

    assert magnetic_pressure(B).unit.is_equivalent(u.Pa)

    assert magnetic_pressure(B).unit.name == "Pa"

    assert magnetic_pressure(B).value == magnetic_energy_density(B).value

    assert magnetic_pressure(B) == magnetic_energy_density(B.to(u.G))

    assert np.isclose(magnetic_pressure(B).value, 397887.35772973835)

    with pytest.warns(u.UnitsWarning):
        magnetic_pressure(5)

    with pytest.raises(u.UnitTypeError):
        magnetic_pressure(5 * u.m)

    assert np.isnan(magnetic_pressure(np.nan * u.T))

    with pytest.raises(ValueError):
        magnetic_pressure(5j * u.T)

    assert np.isnan(magnetic_pressure(B_nanarr)[-1])

    with pytest.warns(u.UnitsWarning):
        assert magnetic_pressure(22.2) == magnetic_pressure(22.2 * u.T)

    assert_can_handle_nparray(magnetic_pressure)


def test_magnetic_energy_density():
    r"""Test the magnetic_energy_density function in parameters.py."""

    assert magnetic_energy_density(B_arr).unit.is_equivalent(u.J / u.m ** 3)

    assert magnetic_energy_density(B).unit.is_equivalent("J / m3")

    assert magnetic_energy_density(B).value == magnetic_pressure(B).value

    assert_quantity_allclose(
        magnetic_energy_density(2 * B), 4 * magnetic_energy_density(B)
    )

    assert_quantity_allclose(magnetic_energy_density(B).value, 397887.35772973835)

    assert_quantity_allclose(
        magnetic_energy_density(B), magnetic_energy_density(B.to(u.G))
    )

    assert isinstance(magnetic_energy_density(B_arr), u.Quantity)

    with pytest.warns(u.UnitsWarning):
        magnetic_energy_density(5)

    with pytest.raises(u.UnitTypeError):
        magnetic_energy_density(5 * u.m)

    assert np.isnan(magnetic_energy_density(np.nan * u.T))

    with pytest.raises(ValueError):
        magnetic_energy_density(5j * u.T)

    assert np.isnan(magnetic_energy_density(B_nanarr)[-1])

    with pytest.warns(u.UnitsWarning):
        assert magnetic_energy_density(22.2) == magnetic_energy_density(22.2 * u.T)

    assert_can_handle_nparray(magnetic_energy_density)


def test_upper_hybrid_frequency():
    r"""Test the upper_hybrid_frequency function in parameters.py."""

    omega_uh = upper_hybrid_frequency(B, n_e=n_e)
    omega_uh_hz = upper_hybrid_frequency(B, n_e=n_e, to_hz=True)
    omega_ce = gyrofrequency(B)
    omega_pe = plasma_frequency(n=n_e)
    assert omega_ce.unit.is_equivalent(u.rad / u.s)
    assert omega_pe.unit.is_equivalent(u.rad / u.s)
    assert omega_uh.unit.is_equivalent(u.rad / u.s)
    assert omega_uh_hz.unit.is_equivalent(u.Hz)
    left_hand_side = omega_uh ** 2
    right_hand_side = omega_ce ** 2 + omega_pe ** 2
    assert np.isclose(left_hand_side.value, right_hand_side.value)

    assert np.isclose(omega_uh_hz.value, 69385868857.90918)

    with pytest.raises(ValueError):
        upper_hybrid_frequency(5 * u.T, n_e=-1 * u.m ** -3)

    with pytest.warns(u.UnitsWarning):
        assert upper_hybrid_frequency(1.2, 1.3) == upper_hybrid_frequency(
            1.2 * u.T, 1.3 * u.m ** -3
        )

    with pytest.warns(u.UnitsWarning):
        assert upper_hybrid_frequency(1.4 * u.T, 1.3) == upper_hybrid_frequency(
            1.4, 1.3 * u.m ** -3
        )

    assert_can_handle_nparray(upper_hybrid_frequency)


def test_lower_hybrid_frequency():
    r"""Test the lower_hybrid_frequency function in parameters.py."""

    ion = "He-4 1+"
    omega_ci = gyrofrequency(B, particle=ion)
    omega_pi = plasma_frequency(n=n_i, particle=ion)
    omega_ce = gyrofrequency(B)
    omega_lh = lower_hybrid_frequency(B, n_i=n_i, ion=ion)
    omega_lh_hz = lower_hybrid_frequency(B, n_i=n_i, ion=ion, to_hz=True)
    assert omega_ci.unit.is_equivalent(u.rad / u.s)
    assert omega_pi.unit.is_equivalent(u.rad / u.s)
    assert omega_ce.unit.is_equivalent(u.rad / u.s)
    assert omega_lh.unit.is_equivalent(u.rad / u.s)
    left_hand_side = omega_lh ** -2
    right_hand_side = (
        1 / (omega_ci ** 2 + omega_pi ** 2) + omega_ci ** -1 * omega_ce ** -1
    )
    assert np.isclose(left_hand_side.value, right_hand_side.value)

    assert np.isclose(omega_lh_hz.value, 299878691.3223296)

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=5e19 * u.m ** -3, ion="asdfasd")

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=-5e19 * u.m ** -3, ion="asdfasd")

    with pytest.raises(ValueError):
        lower_hybrid_frequency(np.nan * u.T, n_i=-5e19 * u.m ** -3, ion="asdfasd")

    with pytest.warns(u.UnitsWarning):
        assert lower_hybrid_frequency(1.3, 1e19) == lower_hybrid_frequency(
            1.3 * u.T, 1e19 * u.m ** -3
        )
    assert_can_handle_nparray(lower_hybrid_frequency)


def test_Bohm_diffusion():
    r"""Test Mag_Reynolds in dimensionless.py"""

    T_e = 5000 * u.K
    B = 10 * u.T

    assert (Bohm_diffusion(T_e, B)).unit == u.m ** 2 / u.s

    with pytest.warns(u.UnitsWarning):
        Bohm_diffusion(5000, B)

    with pytest.raises(u.UnitTypeError):
        Bohm_diffusion(2.2 * u.kg, B)


def test_parameters_aliases():
    r"""Test all aliases defined in parameters.py"""

    assert rho_ is mass_density
    assert va_ is Alfven_speed
    assert cs_ is ion_sound_speed
    assert vth_ is thermal_speed
    assert pth_ is thermal_pressure
    assert vth_kappa_ is kappa_thermal_speed
    assert betaH_ is Hall_parameter
    assert oc_ is gyrofrequency
    assert wc_ is gyrofrequency
    assert rc_ is gyroradius
    assert rhoc_ is gyroradius
    assert wp_ is plasma_frequency
    assert lambdaD_ is Debye_length
    assert nD_ is Debye_number
    assert cwp_ is inertial_length
    assert pmag_ is magnetic_pressure
    assert ub_ is magnetic_energy_density
    assert wuh_ is upper_hybrid_frequency
    assert wlh_ is lower_hybrid_frequency
    assert DB_ is Bohm_diffusion
