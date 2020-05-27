"""Tests for functions that uses Distribution functions."""

import numpy as np
import pytest
from astropy import units as u
import scipy.integrate as spint
import time

from astropy.constants import m_p, m_e, c, mu0, k_B, e, eps0
from astropy.units.core import UnitTypeError
from ..distribution import (
    Maxwellian_1D,
    Maxwellian_speed_1D,
    Maxwellian_velocity_2D,
    Maxwellian_speed_2D,
    Maxwellian_velocity_3D,
    Maxwellian_speed_3D,
    kappa_velocity_1D,
    kappa_velocity_3D,
)
from ..parameters import thermal_speed, kappa_thermal_speed

# test class for Maxwellian_1D (velocity) function:


class Test_Maxwellian_1D(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T_e = 30000 * u.K
        self.v = 1e5 * u.m / u.s
        self.v_drift = 1000000 * u.m / u.s
        self.v_drift2 = 0 * u.m / u.s
        self.v_drift3 = 1e5 * u.m / u.s
        self.start = -5000
        self.stop = -self.start
        self.dv = 10000 * u.m / u.s
        self.v_vect = np.arange(self.start, self.stop, dtype="float64") * self.dv
        self.particle = "e"
        self.vTh = thermal_speed(
            self.T_e, particle=self.particle, method="most_probable"
        )
        self.distFuncTrue = 5.851627151617136e-07

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
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * self.vTh.si.value

        def Maxwellian_1D_wrap(v, T, particle="e", v_drift=0, vTh=np.nan):
            v = v * (u.m / u.s)
            return Maxwellian_1D(v, T, particle, v_drift, vTh).si.value

        # integrating, this should be close to 1
        integ = spint.quad(
            Maxwellian_1D_wrap,
            -infApprox,
            infApprox,
            args=(self.T_e, self.particle, 0, self.vTh),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from quad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = (
            "Integral of distribution function should be 1 " f"and not {integVal}."
        )
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_std(self):
        """
        Tests standard deviation of function?
        """
        std = (
            Maxwellian_1D(self.v_vect, T=self.T_e, particle=self.particle)
            * self.v_vect ** 2
            * self.dv
        ).sum()
        std = np.sqrt(std)
        T_distri = (std ** 2 / k_B * m_e).to(u.K)
        assert np.isclose(T_distri.value, self.T_e.value)

    def test_vTh_is_nan(self):
        """
        Tests distribution function with units, but not passing vTh.
        """
        distFunc = Maxwellian_1D(
            v=self.v, T=self.T_e, particle=self.particle
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_set_vTh(self):
        """
        Tests distribution function with passing vTh.
        """
        distFunc = Maxwellian_1D(
            v=self.v, T=self.T_e, vTh=self.vTh, particle=self.particle
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

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
        testVal = ((self.vTh ** 2 * np.pi) ** (-1 / 2)).si.value
        distFunc = Maxwellian_1D(
            v=self.v,
            T=self.T_e,
            particle=self.particle,
            v_drift=self.v_drift3,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_velocity_wrong_unit(self):
        """
        Testing velocity passed with wrong unit
        """

        with pytest.raises(UnitTypeError):
            Maxwellian_1D(
                v=self.v * u.m,
                T=self.T_e,
                particle=self.particle,
                v_drift=self.v_drift3,
            )

    def test_temperature_wrong_unit(self):
        """
        Testing temperature passed with wrong unit
        """
        with pytest.raises(UnitTypeError):
            Maxwellian_1D(
                v=self.v,
                T=0 * u.C,
                particle=self.particle,
                v_drift=self.v_drift3,
            )

    def test_convert_v_unit(self):
        testVal = ((self.vTh ** 2 * np.pi) ** (-1 / 2)).si.value
        distFunc = Maxwellian_1D(
            v=self.v * 100 / u.m * u.cm,
            T=self.T_e,
            particle=self.particle,
            v_drift=self.v_drift3,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr


    def test_convert_drift_units(self):
        testVal = ((self.vTh ** 2 * np.pi) ** (-1 / 2)).si.value
        distFunc = Maxwellian_1D(
            v=self.v,
            T=self.T_e,
            particle=self.particle,
            v_drift=self.v_drift3 * 100 / u.m * u.cm,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_temperature_units(self):
        testVal = ((self.vTh ** 2 * np.pi) ** (-1 / 2)).si.value
        distFunc = Maxwellian_1D(
            v=self.v,
            T=self.T_e.to(u.eV, equivalencies=u.temperature_energy()),
            particle=self.particle,
            v_drift=self.v_drift3,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

# test class for Maxwellian_speed_1D function
class Test_Maxwellian_speed_1D(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T = 1.0 * u.eV
        self.particle = "H+"
        # get thermal velocity and thermal velocity squared
        self.vTh = thermal_speed(self.T, particle=self.particle, method="most_probable")
        self.v = 1e5 * u.m / u.s
        self.v_drift = 0 * u.m / u.s
        self.v_drift2 = 1e5 * u.m / u.s
        self.distFuncTrue = 1.72940389716217e-27
        self.distFuncDrift = 2 * (self.vTh ** 2 * np.pi) ** (-1 / 2)

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # setting up integration from 0 to 10*vTh
        xData1D = np.arange(0, 10.01, 0.01) * self.vTh.si.value

        def Maxwellian_speed_1D_wrap(v, T, particle="e", v_drift=0, vTh=np.nan):
            v = v * (u.m / u.s)
            return Maxwellian_speed_1D(v, T, particle, v_drift, vTh).si.value

        yData1D = Maxwellian_speed_1D_wrap(v=xData1D, T=self.T, particle=self.particle)
        # integrating, this should be close to 1
        integ = spint.trapz(y=yData1D, x=xData1D)
        exceptStr = "Integral of distribution function should be 1."
        assert np.isclose(integ, 1), exceptStr

    def test_units_no_vTh(self):
        """
        Tests distribution function with not passing vTh.
        """
        distFunc = Maxwellian_speed_1D(
            v=self.v, T=self.T, particle=self.particle
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
        Tests distribution function with passing vTh.
        """
        distFunc = Maxwellian_speed_1D(
            v=self.v, T=self.T, vTh=self.vTh, particle=self.particle
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = Maxwellian_speed_1D(
            v=self.v,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift,
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
        )
        errStr = f"Distribution function should be 0.0 " f"and not {distFunc}."
        assert np.isclose(
            distFunc.value, self.distFuncDrift.value, rtol=1e-5, atol=0.0
        ), errStr

    def test_convert_v_unit(self):
        distFunc = Maxwellian_speed_1D(
            v=self.v * 100 / u.m * u.cm,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift2,
        )
        errStr = f"Distribution function should be 0.0 " f"and not {distFunc}."
        assert np.isclose(
            distFunc.value, self.distFuncDrift.value, rtol=1e-5, atol=0.0
        ), errStr

    def test_convert_drift_units(self):
        distFunc = Maxwellian_speed_1D(
            v=self.v,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift2 * 100 / u.m * u.cm,
        )
        errStr = f"Distribution function should be 0.0 " f"and not {distFunc}."
        assert np.isclose(
            distFunc.value, self.distFuncDrift.value, rtol=1e-5, atol=0.0
        ), errStr

    def test_convert_temperature_units(self):
        distFunc = Maxwellian_speed_1D(
            v=self.v,
            T=self.T.to(u.eV, equivalencies=u.temperature_energy()),
            particle=self.particle,
            v_drift=self.v_drift2,
        )
        errStr = f"Distribution function should be 0.0 " f"and not {distFunc}."
        assert np.isclose(
            distFunc.value, self.distFuncDrift.value, rtol=1e-5, atol=0.0
        ), errStr


# test class for Maxwellian_velocity_2D function
class Test_Maxwellian_velocity_2D(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T = 1.0 * u.eV
        self.particle = "H+"
        # get thermal velocity and thermal velocity squared
        self.vTh = thermal_speed(self.T, particle=self.particle, method="most_probable")
        self.vx = 1e5 * u.m / u.s
        self.vy = 1e5 * u.m / u.s
        self.vx_drift = 0 * u.m / u.s
        self.vy_drift = 0 * u.m / u.s
        self.vx_drift2 = 1e5 * u.m / u.s
        self.vy_drift2 = 1e5 * u.m / u.s
        self.distFuncTrue = 7.477094598799251e-55

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * self.vTh.si.value

        def Maxwellian_velocity_2D_wrap(
            vx, vy, T, particle, vx_drift, vy_drift, vTh
        ):
            return Maxwellian_velocity_2D(
                vx, vy, T, particle, vx_drift, vy_drift, vTh
            ).si.value

        # integrating, this should be close to 1
        integ = spint.dblquad(
            Maxwellian_velocity_2D,
            -infApprox,
            infApprox,
            lambda y: -infApprox,
            lambda y: infApprox,
            args=(self.T, self.particle, 0, 0, self.vTh),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from dblquad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = (
            "Integral of distribution function should be 1 " f"and not {integVal}."
        )
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_no_vTh(self):
        """
        Tests distribution function with not passing vTh.
        """
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx, vy=self.vy, T=self.T, particle=self.particle
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_vTh(self):
        """
        Tests distribution function with passing vTh.
        """
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx,
            vy=self.vy,
            T=self.T,
            vTh=self.vTh,
            particle=self.particle,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_zero_drift(self):
        """
        Testing inputting drift equal to 0. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx,
            vy=self.vy,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift,
            vy_drift=self.vy_drift,
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
        testVal = ((self.vTh ** 2 * np.pi) ** (-1)).si.value
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx,
            vy=self.vy,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_v_unit(self):
        testVal = ((self.vTh ** 2 * np.pi) ** (-1)).si.value
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx * 100 / u.m * u.cm,
            vy=self.vy,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_drift_units(self):
        testVal = ((self.vTh ** 2 * np.pi) ** (-1)).si.value
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx,
            vy=self.vy,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift2 * 100  / u.m * u.cm,
            vy_drift=self.vy_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_temperature_units(self):
        testVal = ((self.vTh ** 2 * np.pi) ** (-1)).si.value
        distFunc = Maxwellian_velocity_2D(
            vx=self.vx,
            vy=self.vy,
            T=self.T.to(u.eV, equivalencies=u.temperature_energy()),
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_wrong_drift_unit(self):
        """
        Testing vdrifts with wrong units
        """
        with pytest.raises(UnitTypeError):
            Maxwellian_velocity_2D(
                vx=self.vx,
                vy=self.vy,
                T=self.T,
                particle=self.particle,
                vx_drift=0 * u.m / u.s / u.s,
                vy_drift=self.vy_drift2,
            )


# test class for Maxwellian_speed_2D function
class Test_Maxwellian_speed_2D(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T = 1.0 * u.eV
        self.particle = "H+"
        # get thermal velocity and thermal velocity squared
        self.vTh = thermal_speed(self.T, particle=self.particle, method="most_probable")
        self.v = 1e5 * u.m / u.s
        self.v_drift = 0 * u.m / u.s
        self.v_drift2 = 1e5 * u.m / u.s
        self.distFuncTrue = 2.2148166449365907e-26

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # setting up integration from 0 to 10*vTh
        xData1D = np.arange(0, 10.001, 0.001) * self.vTh.si.value

        def Maxwellian_speed_2D_wrap(v, T, particle="e", v_drift=0, vTh=np.nan):
            v = v * (u.m / u.s)
            return Maxwellian_speed_2D(v, T, particle, v_drift, vTh).si.value

        yData1D = Maxwellian_speed_2D_wrap(v=xData1D, T=self.T, particle=self.particle)
        # integrating, this should be close to 1
        integ = spint.trapz(y=yData1D, x=xData1D)
        exceptStr = "Integral of distribution function should be 1."
        assert np.isclose(integ, 1), exceptStr

    def test_no_vTh(self):
        """
        Tests distribution function with not passing vTh.
        """
        distFunc = Maxwellian_speed_2D(
            v=self.v, T=self.T, particle=self.particle
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_set_vTh(self):
        """
        Tests distribution function with passing vTh.
        """
        distFunc = Maxwellian_speed_2D(
            v=self.v, T=self.T, vTh=self.vTh, particle=self.particle
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_no_vTh(self):
        """
        Tests distribution function without passing vTh.
        """
        # converting T to SI then stripping units
        T = self.T.to(u.K, equivalencies=u.temperature_energy())
        distFunc = Maxwellian_speed_2D(
            v=self.v, T=T, particle=self.particle
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0), errStr

    def test_zero_drift(self):
        """
        Testing inputting drift equal to 0. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = Maxwellian_speed_2D(
            v=self.v,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift,
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
            )

    def test_convert_v_unit(self):
        distFunc = Maxwellian_speed_2D(
            v=self.v * 100 / u.m * u.cm,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_convert_drift_units(self):
        distFunc = Maxwellian_speed_2D(
            v=self.v,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift * 100 / u.m * u.cm,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_convert_temperature_units(self):
        distFunc = Maxwellian_speed_2D(
            v=self.v,
            T=self.T.to(u.eV, equivalencies=u.temperature_energy()),
            particle=self.particle,
            v_drift=self.v_drift * 100 / u.m * u.cm,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr


# # test class for Maxwellian_velocity_3D function
class Test_Maxwellian_velocity_3D(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T = 1.0 * u.eV
        self.particle = "H+"
        # get thermal velocity and thermal velocity squared
        self.vTh = thermal_speed(self.T, particle=self.particle, method="most_probable")
        self.vx = 1e5 * u.m / u.s
        self.vy = 1e5 * u.m / u.s
        self.vz = 1e5 * u.m / u.s
        self.vx_drift = 0 * u.m / u.s
        self.vy_drift = 0 * u.m / u.s
        self.vz_drift = 0 * u.m / u.s
        self.vx_drift2 = 1e5 * u.m / u.s
        self.vy_drift2 = 1e5 * u.m / u.s
        self.vz_drift2 = 1e5 * u.m / u.s
        self.distFuncTrue = 6.465458269306909e-82

    @pytest.mark.slow
    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # converting vTh to unitless
        vTh = self.vTh.si.value
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * self.vTh.si.value

        def Maxwellian_velocity_3D_wrap(
            vx, vy, vz, T, particle="e", vx_drift=0, vy_drift=0, vz_drift=0, vTh=np.nan
        ):
            vx = vx * (u.m / u.s)
            vy = vy * (u.m / u.s)
            vz = vz * (u.m / u.s)
            return Maxwellian_velocity_3D(
                vx, vy, vz, T, particle, vx_drift, vy_drift, vz_drift, vTh
            ).si.value

        integ = spint.tplquad(
            Maxwellian_velocity_3D_wrap,
            -infApprox,
            infApprox,
            lambda z: -infApprox,
            lambda z: infApprox,
            lambda z, y: -infApprox,
            lambda z, y: infApprox,
            args=(self.T, self.particle, 0, 0, 0, vTh),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from tplquad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = (
            "Integral of distribution function should be 1 " f"and not {integVal}."
        )
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_no_vTh(self):
        """
        Tests distribution function without passing vTh.
        """
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            particle=self.particle,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_set_vTh(self):
        """
        Tests distribution function with passing vTh.
        """
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            vTh=self.vTh,
            particle=self.particle,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_zero_drift_units(self):
        """
        Testing inputting drift equal to 0. These should just
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
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift(self):
        """
        Testing vdrifts with values
        """
        testVal = ((self.vTh ** 2 * np.pi) ** (-3 / 2)).si.value
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
            vz_drift=self.vz_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_wrong_drift_unit(self):
        """
        Testing vdrifts with wrong units
        """
        with pytest.raises(UnitTypeError):
            Maxwellian_velocity_3D(
                vx=self.vx,
                vy=self.vy,
                vz=self.vz,
                T=self.T,
                particle=self.particle,
                vx_drift=0 * u.m / u.s / u.s,
                vy_drift=self.vy_drift2,
                vz_drift=self.vy_drift2,
            )

    def test_convert_v_unit(self):
        testVal = ((self.vTh ** 2 * np.pi) ** (-3 / 2)).si.value
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx * 100 / u.m * u.cm,
            vy=self.vy,
            vz=self.vz * 100 / u.m * u.cm,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
            vz_drift=self.vz_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_drift_units(self):
        testVal = ((self.vTh ** 2 * np.pi) ** (-3 / 2)).si.value
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            particle=self.particle,
            vx_drift=self.vx_drift2 * 100 / u.m * u.cm,
            vy_drift=self.vy_drift2 * 100 / u.m * u.cm,
            vz_drift=self.vz_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_temperature_units(self):
        testVal = ((self.vTh ** 2 * np.pi) ** (-3 / 2)).si.value
        distFunc = Maxwellian_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T.to(u.eV, equivalencies=u.temperature_energy()),
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
            vz_drift=self.vz_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr


# test class for Maxwellian_speed_3D function
class Test_Maxwellian_speed_3D(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T = 1.0 * u.eV
        self.particle = "H+"
        # get thermal velocity and thermal velocity squared
        self.vTh = thermal_speed(self.T, particle=self.particle, method="most_probable")
        self.v = 1e5 * u.m / u.s
        self.v_drift = 0 * u.m / u.s
        self.v_drift2 = 1e5 * u.m / u.s
        self.distFuncTrue = 1.8057567503860518e-25

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # setting up integration from 0 to 10*vTh
        xData1D = np.arange(0, 10.01, 0.01) * self.vTh.si.value

        def Maxwellian_speed_3D_wrap(v, T, particle="e", v_drift=0, vTh=np.nan):
            v = v * (u.m / u.s)
            return Maxwellian_speed_3D(v, T, particle, v_drift, vTh).si.value

        yData1D = Maxwellian_speed_3D_wrap(v=xData1D, T=self.T, particle=self.particle)
        # integrating, this should be close to 1
        integ = spint.trapz(y=yData1D, x=xData1D)
        exceptStr = "Integral of distribution function should be 1."
        assert np.isclose(integ, 1), exceptStr

    def test_no_vTh(self):
        """
        Tests distribution function without units passing vTh.
        """
        distFunc = Maxwellian_speed_3D(
            v=self.v, T=self.T, particle=self.particle
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_set_vTh(self):
        """
        Tests distribution function with units and passing vTh.
        """
        distFunc = Maxwellian_speed_3D(
            v=self.v, T=self.T, vTh=self.vTh, particle=self.particle
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

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
            )

    def test_convert_v_unit(self):
        distFunc = Maxwellian_speed_3D(
            v=self.v * 100 / u.m * u.cm,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_convert_drift_units(self):
        distFunc = Maxwellian_speed_3D(
            v=self.v,
            T=self.T,
            particle=self.particle,
            v_drift=self.v_drift * 100 / u.m * u.cm,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_convert_temperature_units(self):
        distFunc = Maxwellian_speed_3D(
            v=self.v,
            T=self.T.to(u.eV, equivalencies=u.temperature_energy()),
            particle=self.particle,
            v_drift=self.v_drift * 100 / u.m * u.cm,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

# kappa
# test class for kappa_velocity_1D function:
class Test_kappa_velocity_1D(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T_e = 30000 * u.K
        self.kappa = 4
        self.kappaInvalid = 3 / 2
        self.v = 1e5 * u.m / u.s
        self.v_drift = 1000000 * u.m / u.s
        self.v_drift2 = 0 * u.m / u.s
        self.v_drift3 = 1e5 * u.m / u.s
        self.start = -5000
        self.stop = -self.start
        self.dv = 10000 * u.m / u.s
        self.v_vect = np.arange(self.start, self.stop, dtype="float64") * self.dv
        self.particle = "e"
        self.vTh = kappa_thermal_speed(
            self.T_e, kappa=self.kappa, particle=self.particle
        )
        self.distFuncTrue = 6.637935187755855e-07

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

    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # converting vTh to unitless
        vTh = self.vTh.si.value
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * vTh

        def kappa_velocity_1D_wrap(v, T, kappa, particle="e", v_drift=0, vTh=np.nan):
            v = v * (u.m / u.s)
            return kappa_velocity_1D(v, T, kappa, particle, v_drift, vTh).si.value

        # integrating, this should be close to 1
        integ = spint.quad(
            kappa_velocity_1D_wrap,
            -infApprox,
            infApprox,
            args=(self.T_e, self.kappa, self.particle, 0, vTh),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from quad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = (
            "Integral of distribution function should be 1 " f"and not {integVal}."
        )
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_std(self):
        """
        Tests standard deviation of function?
        """
        std = (
            kappa_velocity_1D(
                self.v_vect, T=self.T_e, kappa=self.kappa, particle=self.particle
            )
            * self.v_vect ** 2
            * self.dv
        ).sum()
        std = np.sqrt(std)
        T_distri = (std ** 2 / k_B * m_e).to(u.K)
        assert np.isclose(T_distri.value, self.T_e.value)

    def test_no_vTh(self):
        """
        Tests distribution function without passing vTh.
        """
        distFunc = kappa_velocity_1D(
            v=self.v,
            T=self.T_e,
            kappa=self.kappa,
            particle=self.particle,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_vTh(self):
        """
        Tests distribution function with passing vTh.
        """
        distFunc = kappa_velocity_1D(
            v=self.v,
            T=self.T_e,
            kappa=self.kappa,
            vTh=self.vTh,
            particle=self.particle,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_zero_drift(self):
        """
        Testing inputting drift equal to 0. These should just
        get passed and not have extra units applied to them.
        """
        distFunc = kappa_velocity_1D(
            v=self.v,
            T=self.T_e,
            kappa=self.kappa,
            particle=self.particle,
            v_drift=self.v_drift2,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_value_drift(self):
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
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_v_unit(self):
        testVal = 6.755498543630533e-07
        distFunc = kappa_velocity_1D(
            v=self.v * 100 / u.m * u.cm,
            T=self.T_e,
            kappa=self.kappa,
            particle=self.particle,
            v_drift=self.v_drift3,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_drift_units(self):
        testVal = 6.755498543630533e-07
        distFunc = kappa_velocity_1D(
            v=self.v,
            T=self.T_e,
            kappa=self.kappa,
            particle=self.particle,
            v_drift=self.v_drift3 * 100 / u.m * u.cm,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_temperature_units(self):
        testVal = 6.755498543630533e-07
        distFunc = kappa_velocity_1D(
            v=self.v,
            T=self.T_e.to(u.eV, equivalencies=u.temperature_energy()),
            kappa=self.kappa,
            particle=self.particle,
            v_drift=self.v_drift3,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr


# test class for kappa_velocity_3D function
class Test_kappa_velocity_3D(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests """
        self.T = 1.0 * u.eV
        self.kappa = 4
        self.kappaInvalid = 3 / 2
        self.particle = "H+"
        # get thermal velocity and thermal velocity squared
        self.vTh = kappa_thermal_speed(self.T, kappa=self.kappa, particle=self.particle)
        self.vx = 1e5 * u.m / u.s
        self.vy = 1e5 * u.m / u.s
        self.vz = 1e5 * u.m / u.s
        self.vx_drift = 0 * u.m / u.s
        self.vy_drift = 0 * u.m / u.s
        self.vz_drift = 0 * u.m / u.s
        self.vx_drift2 = 1e5 * u.m / u.s
        self.vy_drift2 = 1e5 * u.m / u.s
        self.vz_drift2 = 1e5 * u.m / u.s
        self.distFuncTrue = 1.1847914288918793e-22

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
    #                                          vz_drift=self.vz_drift2)
    #        Teff =  self.T
    #        maxwellDistFunc = Maxwellian_velocity_3D(vx=self.vx,
    #                                                 vy=self.vy,
    #                                                 vz=self.vz,
    #                                                 T=Teff,
    #                                                 particle=self.particle,
    #                                                 vx_drift=self.vx_drift2,
    #                                                 vy_drift=self.vy_drift2,
    #                                                 vz_drift=self.vz_drift2)
    #        errStr = (f"Distribution function should be {maxwellDistFunc} "
    #                  f"and not {kappaDistFunc}.")
    #        assert np.isclose(kappaDistFunc.value,
    #                          maxwellDistFunc.value,
    #                          rtol=1e-5,
    #                          atol=0.0), errStr
    #
    #        return

    @pytest.mark.slow
    def test_norm(self):
        """
        Tests whether distribution function is normalized, and integrates to 1.
        """
        # converting vTh to unitless
        vTh = self.vTh.si.value
        # setting up integration from -10*vTh to 10*vTh, which is close to Inf
        infApprox = 10 * vTh

        def kappa_velocity_3D_wrap(
            vx, vy, vz, T, kappa, particle="e", vx_drift=0, vy_drift=0, vz_drift=0, vTh=np.nan
        ):
            vx = vx * (u.m / u.s)
            vy = vy * (u.m / u.s)
            vz = vz * (u.m / u.s)
            return kappa_velocity_3D(
                vx, vy, vz, T, kappa, particle, vx_drift, vy_drift, vz_drift, vTh
            ).si.value

        # integrating, this should be close to 1
        integ = spint.tplquad(
            kappa_velocity_3D_wrap,
            -infApprox,
            infApprox,
            lambda z: -infApprox,
            lambda z: infApprox,
            lambda z, y: -infApprox,
            lambda z, y: infApprox,
            args=(self.T, self.kappa, self.particle, 0, 0, 0, vTh),
            epsabs=1e0,
            epsrel=1e0,
        )
        # value returned from tplquad is (integral, error), we just need
        # the 1st
        integVal = integ[0]
        exceptStr = (
            "Integral of distribution function should be 1 " f"and not {integVal}."
        )
        assert np.isclose(integVal, 1, rtol=1e-3, atol=0.0), exceptStr

    def test_no_vTh(self):
        """
        Tests distribution function without passing vTh.
        """
        distFunc = kappa_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            kappa=self.kappa,
            particle=self.particle,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

    def test_set_vTh(self):
        """
        Tests distribution function with passing vTh.
        """
        distFunc = kappa_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            kappa=self.kappa,
            vTh=self.vTh,
            particle=self.particle,
        )
        errStr = (
            f"Distribution function should be {self.distFuncTrue} "
            f"and not {distFunc}."
        )
        assert np.isclose(
            distFunc.value, self.distFuncTrue, rtol=1e-5, atol=0.0
        ), errStr

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
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_v_unit(self):
        testVal = 1.2376545373917465e-13
        distFunc = kappa_velocity_3D(
            vx=self.vx * 100 / u.m * u.cm,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            kappa=self.kappa,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2,
            vz_drift=self.vz_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_drift_units(self):
        testVal = 1.2376545373917465e-13
        distFunc = kappa_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T,
            kappa=self.kappa,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2 * 100 / u.m * u.cm,
            vz_drift=self.vz_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr

    def test_convert_temperature_units(self):
        testVal = 1.2376545373917465e-13
        distFunc = kappa_velocity_3D(
            vx=self.vx,
            vy=self.vy,
            vz=self.vz,
            T=self.T.to(u.eV, equivalencies=u.temperature_energy()),
            kappa=self.kappa,
            particle=self.particle,
            vx_drift=self.vx_drift2,
            vy_drift=self.vy_drift2 * 100 / u.m * u.cm,
            vz_drift=self.vz_drift2,
        )
        errStr = f"Distribution function should be {testVal} " f"and not {distFunc}."
        assert np.isclose(distFunc.value, testVal, rtol=1e-5, atol=0.0), errStr
