"""Tests for functions that uses Distribution functions."""

import numpy as np
import pytest
import scipy.integrate as spint

from astropy import units as u
from astropy.constants import k_B, m_e

from plasmapy.formulary.distribution import (
    kappa_velocity_1D,
    kappa_velocity_3D,
    Maxwellian_1D,
    Maxwellian_speed_1D,
    Maxwellian_speed_2D,
    Maxwellian_speed_3D,
    Maxwellian_velocity_2D,
    Maxwellian_velocity_3D,
)
from plasmapy.formulary.speeds import kappa_thermal_speed, thermal_speed


class Test_Maxwellian_1D:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T_e = 30000 * u.K
        cls.v = 1e5 * u.m / u.s
        cls.v_drift = 1000000 * u.m / u.s
        cls.v_drift2 = 0 * u.m / u.s
        cls.v_drift3 = 1e5 * u.m / u.s
        cls.start = -5000
        cls.stop = -cls.start
        cls.dv = 10000 * u.m / u.s
        cls.v_vect = np.arange(cls.start, cls.stop, dtype="float64") * cls.dv
        cls.particle = "e"
        cls.vTh = thermal_speed(cls.T_e, particle=cls.particle, method="most_probable")
        cls.distFuncTrue = 5.851627151617136e-07

    def test_max_noDrift(self):
        """
        Checks maximum value of distribution function is in expected place,
        when there is no drift applied.
        """
        max_index = Maxwellian_1D(
            self.v_vect, T=self.T_e, particle=self.particle, v_drift=0 * u.m / u.s
        ).argmax()
        assert np.isclose(self.v_vect[max_index].value, 0.0)

    def test_max_drift(self):
        """
        Checks maximum value of distribution function is in expected place,
        when there is drift applied.
        """
        max_index = Maxwellian_1D(
            self.v_vect, T=self.T_e, particle=self.particle, v_drift=self.v_drift
        ).argmax()
        assert np.isclose(self.v_vect[max_index].value, self.v_drift.value)

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # converting vTh to unitless
        vTh = self.vTh.si.value
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * vTh
        # integrating, this should be close to 1
        integ = spint.quad(
            Maxwellian_1D,
            -infApprox,
            infApprox,
            args=(self.T_e, self.particle, 0, vTh, "unitless"),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from quad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = "Integral of distribution function should be 1 and not {integVal}."
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_std(self):
        """
        Tests standard deviation of function?
        """
        std = (
            Maxwellian_1D(self.v_vect, T=self.T_e, particle=self.particle)
            * self.v_vect**2
            * self.dv
        ).sum()
        std = np.sqrt(std)
        T_distri = (std**2 / k_B * m_e).to(u.K)
        assert np.isclose(T_distri.value, self.T_e.value)

    def test_units_no_vTh(self):
        """
        Tests distribution function with units, but not passing vTh.
        """
        distFunc = Maxwellian_1D(
            v=self.v, T=self.T_e, particle=self.particle, units="units"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_units_vTh(self):
        """
        Tests distribution function with units and passing vTh.
        """
        distFunc = Maxwellian_1D(
            v=self.v, T=self.T_e, vTh=self.vTh, particle=self.particle, units="units"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_unitless_no_vTh(self):
        """
        Tests distribution function without units, and not passing vTh.
        """
        # converting T to SI then stripping units
        T_e = self.T_e.to(u.K, equivalencies=u.temperature_energy())
        T_e = T_e.si.value
        distFunc = Maxwellian_1D(
            v=self.v.si.value, T=T_e, particle=self.particle, units="unitless"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_unitless_vTh(self):
        """
        Tests distribution function without units, and with passing vTh.
        """
        # converting T to SI then stripping units
        T_e = self.T_e.to(u.K, equivalencies=u.temperature_energy())
        T_e = T_e.si.value
        distFunc = Maxwellian_1D(
            v=self.v.si.value,
            T=T_e,
            vTh=self.vTh.si.value,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0 with units. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = Maxwellian_1D(
            v=self.v,
            T=self.T_e,
            particle=self.particle,
            v_drift=self.v_drift2,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift_units(self):
        """
        Testing vdrifts with values
        """
        testVal = ((self.vTh**2 * np.pi) ** (-1 / 2)).si.value
        distFunc = Maxwellian_1D(
            v=self.v,
            T=self.T_e,
            particle=self.particle,
            v_drift=self.v_drift3,
            units="units",
        )
        errStr = f"Distribution function should be {testVal} and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr


class Test_Maxwellian_speed_1D:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 1.0 * u.eV
        cls.particle = "H+"
        # get thermal velocity and thermal velocity squared
        cls.vTh = thermal_speed(cls.T, particle=cls.particle, method="most_probable")
        cls.v = 1e5 * u.m / u.s
        cls.v_drift = 0 * u.m / u.s
        cls.v_drift2 = 1e5 * u.m / u.s
        cls.distFuncTrue = 1.72940389716217e-27
        cls.distFuncDrift = 2 * (cls.vTh**2 * np.pi) ** (-1 / 2)

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # setting up integration from 0 to 10*vTh
        xData1D = np.arange(0, 10.01, 0.01) * self.vTh
        yData1D = Maxwellian_speed_1D(v=xData1D, T=self.T, particle=self.particle)
        # integrating, this should be close to 1
        integ = spint.trapz(y=yData1D, x=xData1D)
        exceptStr = "Integral of distribution function should be 1."
        assert np.isclose(integ.value, 1), exceptStr

    def test_units_no_vTh(self):
        """
        Tests distribution function with units, but not passing vTh.
        """
        distFunc = Maxwellian_speed_1D(
            v=self.v, T=self.T, particle=self.particle, units="units"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_units_vTh(self):
        """
        Tests distribution function with units and passing vTh.
        """
        distFunc = Maxwellian_speed_1D(
            v=self.v, T=self.T, vTh=self.vTh, particle=self.particle, units="units"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_unitless_no_vTh(self):
        """
        Tests distribution function without units, and not passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_speed_1D(
            v=self.v.si.value, T=T, particle=self.particle, units="unitless"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_unitless_vTh(self):
        """
        Tests distribution function without units, and with passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_speed_1D(
            v=self.v.si.value,
            T=T,
            vTh=self.vTh.si.value,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0 with units. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = Maxwellian_speed_1D(
            v=self.v,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift_units(self):
        """
        Testing vdrifts with values
        """
        distFunc = Maxwellian_speed_1D(
            v=self.v,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift2,
            units="units",
        )
        errStr = f"Distribution function should be 0.0 and not {distFunc}."
        assert np.isclose(
            distFunc.value, self.distFuncDrift.value, rtol=1e-5, atol=0.0
        ), errStr


class Test_Maxwellian_velocity_2D:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 1.0 * u.eV
        cls.particle = "H+"
        # get thermal velocity and thermal velocity squared
        cls.vTh = thermal_speed(cls.T, particle=cls.particle, method="most_probable")
        cls.vx = 1e5 * u.m / u.s
        cls.vy = 1e5 * u.m / u.s
        cls.vx_drift = 0 * u.m / u.s
        cls.vy_drift = 0 * u.m / u.s
        cls.vx_drift2 = 1e5 * u.m / u.s
        cls.vy_drift2 = 1e5 * u.m / u.s
        cls.distFuncTrue = 7.477094598799251e-55

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # converting vTh to unitless
        vTh = self.vTh.si.value
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * vTh
        # integrating, this should be close to 1
        integ = spint.dblquad(
            Maxwellian_velocity_2D,
            -infApprox,
            infApprox,
            lambda y: -infApprox,
            lambda y: infApprox,
            args=(self.T, self.particle, 0, 0, vTh, "unitless"),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from dblquad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = "Integral of distribution function should be 1 and not {integVal}."
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_units_no_vTh(self):
        """
        Tests distribution function with units, but not passing vTh.
        """
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx, vy=self.vy, T=self.T, particle=self.particle, units="units"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_units_vTh(self):
        """
        Tests distribution function with units and passing vTh.
        """
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx,
            vy=self.vy,
            T=self.T,
            vTh=self.vTh,
            particle=self.particle,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_unitless_no_vTh(self):
        """
        Tests distribution function without units, and not passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx.si.value,
            vy=self.vy.si.value,
            T=T,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_unitless_vTh(self):
        """
        Tests distribution function without units, and with passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx.si.value,
            vy=self.vy.si.value,
            T=T,
            vTh=self.vTh.si.value,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0 with units. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx,
            vy=self.vy,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift,
            vy_drift=self.vy_drift,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift_units(self):
        """
        Testing vdrifts with values
        """
        testVal = ((self.vTh**2 * np.pi) ** (-1)).si.value
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx,
            vy=self.vy,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
            units="units",
        )
        errStr = f"Distribution function should be {testVal} and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr


@pytest.mark.slow
class Test_Maxwellian_speed_2D:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 1.0 * u.eV
        cls.particle = "H+"
        # get thermal velocity and thermal velocity squared
        cls.vTh = thermal_speed(cls.T, particle=cls.particle, method="most_probable")
        cls.v = 1e5 * u.m / u.s
        cls.v_drift = 0 * u.m / u.s
        cls.v_drift2 = 1e5 * u.m / u.s
        cls.distFuncTrue = 2.2148166449365907e-26

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # setting up integration from 0 to 10*vTh
        xData1D = np.arange(0, 10.001, 0.001) * self.vTh
        yData1D = Maxwellian_speed_2D(v=xData1D, T=self.T, particle=self.particle)
        # integrating, this should be close to 1
        integ = spint.trapz(y=yData1D, x=xData1D)
        exceptStr = "Integral of distribution function should be 1."
        assert np.isclose(integ.value, 1), exceptStr

    def test_units_no_vTh(self):
        """
        Tests distribution function with units, but not passing vTh.
        """
        distFunc = Maxwellian_speed_2D(
            v=self.v, T=self.T, particle=self.particle, units="units"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_units_vTh(self):
        """
        Tests distribution function with units and passing vTh.
        """
        distFunc = Maxwellian_speed_2D(
            v=self.v, T=self.T, vTh=self.vTh, particle=self.particle, units="units"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_unitless_no_vTh(self):
        """
        Tests distribution function without units, and not passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_speed_2D(
            v=self.v.si.value, T=T, particle=self.particle, units="unitless"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_unitless_vTh(self):
        """
        Tests distribution function without units, and with passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_speed_2D(
            v=self.v.si.value,
            T=T,
            vTh=self.vTh.si.value,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0 with units. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = Maxwellian_speed_2D(
            v=self.v,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift_units(self):
        """
        Testing vdrifts with values
        """
        with pytest.raises(NotImplementedError):
            distFunc = Maxwellian_speed_2D(
                v=self.v,
                T=self.T,
                particle=self.particle,
                v_drift=self.v_drift2,
                units="units",
            )


@pytest.mark.slow
class Test_Maxwellian_velocity_3D:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 1.0 * u.eV
        cls.particle = "H+"
        # get thermal velocity and thermal velocity squared
        cls.vTh = thermal_speed(cls.T, particle=cls.particle, method="most_probable")
        cls.vx = 1e5 * u.m / u.s
        cls.vy = 1e5 * u.m / u.s
        cls.vz = 1e5 * u.m / u.s
        cls.vx_drift = 0 * u.m / u.s
        cls.vy_drift = 0 * u.m / u.s
        cls.vz_drift = 0 * u.m / u.s
        cls.vx_drift2 = 1e5 * u.m / u.s
        cls.vy_drift2 = 1e5 * u.m / u.s
        cls.vz_drift2 = 1e5 * u.m / u.s
        cls.distFuncTrue = 6.465458269306909e-82

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # converting vTh to unitless
        vTh = self.vTh.si.value
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * vTh
        # integrating, this should be close to 1
        integ = spint.tplquad(
            Maxwellian_velocity_3D,
            -infApprox,
            infApprox,
            lambda z: -infApprox,
            lambda z: infApprox,
            lambda z, y: -infApprox,
            lambda z, y: infApprox,
            args=(self.T, self.particle, 0, 0, 0, vTh, "unitless"),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from tplquad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = "Integral of distribution function should be 1 and not {integVal}."
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_units_no_vTh(self):
        """
        Tests distribution function with units, but not passing vTh.
        """
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            particle=self.particle,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_units_vTh(self):
        """
        Tests distribution function with units and passing vTh.
        """
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            vTh=self.vTh,
            particle=self.particle,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_unitless_no_vTh(self):
        """
        Tests distribution function without units, and not passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx.si.value,
            vy=self.vy.si.value,
            vz=self.vz.si.value,
            T=T,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_unitless_vTh(self):
        """
        Tests distribution function without units, and with passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx.si.value,
            vy=self.vy.si.value,
            vz=self.vz.si.value,
            T=T,
            vTh=self.vTh.si.value,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0 with units. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift,
            vy_drift=self.vy_drift,
            vz_drift=self.vz_drift,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift_units(self):
        """
        Testing vdrifts with values
        """
        testVal = ((self.vTh**2 * np.pi) ** (-3 / 2)).si.value
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
            vz_drift=self.vz_drift2,
            units="units",
        )
        errStr = f"Distribution function should be {testVal} and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr


class Test_Maxwellian_speed_3D:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 1.0 * u.eV
        cls.particle = "H+"
        # get thermal velocity and thermal velocity squared
        cls.vTh = thermal_speed(cls.T, particle=cls.particle, method="most_probable")
        cls.v = 1e5 * u.m / u.s
        cls.v_drift = 0 * u.m / u.s
        cls.v_drift2 = 1e5 * u.m / u.s
        cls.distFuncTrue = 1.8057567503860518e-25

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # setting up integration from 0 to 10*vTh
        xData1D = np.arange(0, 10.01, 0.01) * self.vTh
        yData1D = Maxwellian_speed_3D(v=xData1D, T=self.T, particle=self.particle)
        # integrating, this should be close to 1
        integ = spint.trapz(y=yData1D, x=xData1D)
        exceptStr = "Integral of distribution function should be 1."
        assert np.isclose(integ.value, 1), exceptStr

    def test_units_no_vTh(self):
        """
        Tests distribution function with units, but not passing vTh.
        """
        distFunc = Maxwellian_speed_3D(
            v=self.v, T=self.T, particle=self.particle, units="units"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_units_vTh(self):
        """
        Tests distribution function with units and passing vTh.
        """
        distFunc = Maxwellian_speed_3D(
            v=self.v, T=self.T, vTh=self.vTh, particle=self.particle, units="units"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_unitless_no_vTh(self):
        """
        Tests distribution function without units, and not passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_speed_3D(
            v=self.v.si.value, T=T, particle=self.particle, units="unitless"
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_unitless_vTh(self):
        """
        Tests distribution function without units, and with passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = Maxwellian_speed_3D(
            v=self.v.si.value,
            T=T,
            vTh=self.vTh.si.value,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0 with units. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = Maxwellian_speed_3D(
            v=self.v,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift_units(self):
        """
        Testing vdrifts with values
        """
        with pytest.raises(NotImplementedError):
            distFunc = Maxwellian_speed_3D(
                v=self.v,
                T=self.T,
                particle=self.particle,
                v_drift=self.v_drift2,
                units="units",
            )


class Test_kappa_velocity_1D:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T_e = 30000 * u.K
        cls.kappa = 4
        cls.kappaInvalid = 3 / 2
        cls.v = 1e5 * u.m / u.s
        cls.v_drift = 1000000 * u.m / u.s
        cls.v_drift2 = 0 * u.m / u.s
        cls.v_drift3 = 1e5 * u.m / u.s
        cls.start = -5000
        cls.stop = -cls.start
        cls.dv = 10000 * u.m / u.s
        cls.v_vect = np.arange(cls.start, cls.stop, dtype="float64") * cls.dv
        cls.particle = "e"
        cls.vTh = kappa_thermal_speed(cls.T_e, kappa=cls.kappa, particle=cls.particle)
        cls.distFuncTrue = 6.637935187755855e-07

    def test_invalid_kappa(self):
        """
        Checks if function raises error when kappa <= 3/2 is passed as an
        argument.
        """
        with pytest.raises(ValueError):
            kappa_velocity_1D(
                v=self.v,
                T=self.T_e,
                kappa=self.kappaInvalid,
                particle=self.particle,
                units="units",
            )

    def test_max_noDrift(self):
        """
        Checks maximum value of distribution function is in expected place,
        when there is no drift applied.
        """
        max_index = kappa_velocity_1D(
            self.v_vect,
            T=self.T_e,
            kappa=self.kappa,
            particle=self.particle,
            v_drift=0 * u.m / u.s,
        ).argmax()
        assert np.isclose(self.v_vect[max_index].value, 0.0)

    def test_max_drift(self):
        """
        Checks maximum value of distribution function is in expected place,
        when there is drift applied.
        """
        max_index = kappa_velocity_1D(
            self.v_vect,
            T=self.T_e,
            kappa=self.kappa,
            particle=self.particle,
            v_drift=self.v_drift,
        ).argmax()
        assert np.isclose(self.v_vect[max_index].value, self.v_drift.value)

    def test_maxwellian_limit(self):
        """
        Tests the limit of large kappa to see if kappa distribution function
        converges to Maxwellian.
        """
        return

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # converting vTh to unitless
        vTh = self.vTh.si.value
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * vTh
        # integrating, this should be close to 1
        integ = spint.quad(
            kappa_velocity_1D,
            -infApprox,
            infApprox,
            args=(self.T_e, self.kappa, self.particle, 0, vTh, "unitless"),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from quad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = "Integral of distribution function should be 1 and not {integVal}."
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_std(self):
        """
        Tests standard deviation of function?
        """
        std = (
            kappa_velocity_1D(
                self.v_vect, T=self.T_e, kappa=self.kappa, particle=self.particle
            )
            * self.v_vect**2
            * self.dv
        ).sum()
        std = np.sqrt(std)
        T_distri = (std**2 / k_B * m_e).to(u.K)
        assert np.isclose(T_distri.value, self.T_e.value)

    def test_units_no_vTh(self):
        """
        Tests distribution function with units, but not passing vTh.
        """
        distFunc = kappa_velocity_1D(
            v=self.v,
            T=self.T_e,
            kappa=self.kappa,
            particle=self.particle,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_units_vTh(self):
        """
        Tests distribution function with units and passing vTh.
        """
        distFunc = kappa_velocity_1D(
            v=self.v,
            T=self.T_e,
            kappa=self.kappa,
            vTh=self.vTh,
            particle=self.particle,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_unitless_no_vTh(self):
        """
        Tests distribution function without units, and not passing vTh.
        """
        # converting T to SI then stripping units
        T_e = self.T_e.to(u.K, equivalencies=u.temperature_energy())
        T_e = T_e.si.value
        distFunc = kappa_velocity_1D(
            v=self.v.si.value,
            T=T_e,
            kappa=self.kappa,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_unitless_vTh(self):
        """
        Tests distribution function without units, and with passing vTh.
        """
        # converting T to SI then stripping units
        T_e = self.T_e.to(u.K, equivalencies=u.temperature_energy())
        T_e = T_e.si.value
        distFunc = kappa_velocity_1D(
            v=self.v.si.value,
            T=T_e,
            kappa=self.kappa,
            vTh=self.vTh.si.value,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0 with units. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = kappa_velocity_1D(
            v=self.v,
            T=self.T_e,
            kappa=self.kappa,
            particle=self.particle,
            v_drift=self.v_drift2,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift_units(self):
        """
        Testing vdrifts with values
        """
        testVal = 6.755498543630533e-07
        distFunc = kappa_velocity_1D(
            v=self.v,
            T=self.T_e,
            kappa=self.kappa,
            particle=self.particle,
            v_drift=self.v_drift3,
            units="units",
        )
        errStr = f"Distribution function should be {testVal} and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr


@pytest.mark.slow
class Test_kappa_velocity_3D:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 1.0 * u.eV
        cls.kappa = 4
        cls.kappaInvalid = 3 / 2
        cls.particle = "H+"
        # get thermal velocity and thermal velocity squared
        cls.vTh = kappa_thermal_speed(cls.T, kappa=cls.kappa, particle=cls.particle)
        cls.vx = 1e5 * u.m / u.s
        cls.vy = 1e5 * u.m / u.s
        cls.vz = 1e5 * u.m / u.s
        cls.vx_drift = 0 * u.m / u.s
        cls.vy_drift = 0 * u.m / u.s
        cls.vz_drift = 0 * u.m / u.s
        cls.vx_drift2 = 1e5 * u.m / u.s
        cls.vy_drift2 = 1e5 * u.m / u.s
        cls.vz_drift2 = 1e5 * u.m / u.s
        cls.distFuncTrue = 1.1847914288918793e-22

    def test_invalid_kappa(self):
        """
        Checks if function raises error when kappa <= 3/2 is passed as an
        argument.
        """
        with pytest.raises(ValueError):
            kappa_velocity_3D(
                vx=self.vx,
                vy=self.vy,
                vz=self.vz,
                T=self.T,
                kappa=self.kappaInvalid,
                particle=self.particle,
                units="units",
            )

    #    def test_maxwellian_limit(self):
    #        """
    #        Tests the limit of large kappa to see if kappa distribution function
    #        converges to Maxwellian.
    #        """
    #        kappaLarge = 100
    #        kappaDistFunc = kappa_velocity_3D(vx=self.vx,
    #                                          vy=self.vy,
    #                                          vz=self.vz,
    #                                          T=self.T,
    #                                          kappa=kappaLarge,
    #                                          particle=self.particle,
    #                                          vx_drift=self.vx_drift2,
    #                                          vy_drift=self.vy_drift2,
    #                                          vz_drift=self.vz_drift2,
    #                                          units="units")
    #        Teff =  self.T
    #        maxwellDistFunc = Maxwellian_velocity_3D(vx=self.vx,
    #                                                 vy=self.vy,
    #                                                 vz=self.vz,
    #                                                 T=Teff,
    #                                                 particle=self.particle,
    #                                                 vx_drift=self.vx_drift2,
    #                                                 vy_drift=self.vy_drift2,
    #                                                 vz_drift=self.vz_drift2,
    #                                                 units="units")
    #        errStr = (f"Distribution function should be {maxwellDistFunc} "
    #                  f"and not {kappaDistFunc}.")
    #        assert np.isclose(kappaDistFunc.value,
    #                          maxwellDistFunc.value,
    #                          rtol=1e-5,
    #                          atol=0.0), errStr
    #
    #        return

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # converting vTh to unitless
        vTh = self.vTh.si.value
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * vTh
        # integrating, this should be close to 1
        integ = spint.tplquad(
            kappa_velocity_3D,
            -infApprox,
            infApprox,
            lambda z: -infApprox,
            lambda z: infApprox,
            lambda z, y: -infApprox,
            lambda z, y: infApprox,
            args=(self.T, self.kappa, self.particle, 0, 0, 0, vTh, "unitless"),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from tplquad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = "Integral of distribution function should be 1 and not {integVal}."
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_units_no_vTh(self):
        """
        Tests distribution function with units, but not passing vTh.
        """
        distFunc = kappa_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            kappa=self.kappa,
            particle=self.particle,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_units_vTh(self):
        """
        Tests distribution function with units and passing vTh.
        """
        distFunc = kappa_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            kappa=self.kappa,
            vTh=self.vTh,
            particle=self.particle,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_unitless_no_vTh(self):
        """
        Tests distribution function without units, and not passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = kappa_velocity_3D(
            vx=self.vx.si.value,
            vy=self.vy.si.value,
            vz=self.vz.si.value,
            T=T,
            kappa=self.kappa,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_unitless_vTh(self):
        """
        Tests distribution function without units, and with passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        T = T.si.value
        distFunc = kappa_velocity_3D(
            vx=self.vx.si.value,
            vy=self.vy.si.value,
            vz=self.vz.si.value,
            T=T,
            kappa=self.kappa,
            vTh=self.vTh.si.value,
            particle=self.particle,
            units="unitless",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0 with units. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = kappa_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            kappa=self.kappa,
            particle=self.particle,
            vx_drift=self.vx_drift,
            vy_drift=self.vy_drift,
            vz_drift=self.vz_drift,
            units="units",
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift_units(self):
        """
        Testing vdrifts with values
        """
        testVal = 1.2376545373917465e-13
        distFunc = kappa_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            kappa=self.kappa,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
            vz_drift=self.vz_drift2,
            units="units",
        )
        errStr = f"Distribution function should be {testVal} and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr
