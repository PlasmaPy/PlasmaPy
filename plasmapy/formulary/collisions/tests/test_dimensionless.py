import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.collisions.dimensionless import (
    coupling_parameter,
    Knudsen_number,
)
from plasmapy.utils import exceptions
from plasmapy.utils.exceptions import CouplingWarning
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray


class Test_coupling_parameter:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 11604 * u.K
        cls.n_e = 1e21 * u.cm**-3
        cls.particles = ("e", "p")
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.V = 1e4 * u.km / u.s
        cls.True1 = 2.3213156755481195
        cls.True_zmean = 10.689750083758698
        cls.True_quantum = 0.3334662805238162

    def test_symmetry(self):
        result = coupling_parameter(self.T, self.n_e, self.particles)
        resultRev = coupling_parameter(self.T, self.n_e, self.particles[::-1])
        assert result == resultRev

    def test_known1(self):
        """
        Test for known value.
        """
        methodVal = coupling_parameter(
            self.T,
            self.n_e,
            self.particles,
            z_mean=np.nan * u.dimensionless_unscaled,
            V=np.nan * u.m / u.s,
            method="classical",
        )
        testTrue = np.isclose(self.True1, methodVal, rtol=1e-1, atol=0.0)
        errStr = f"Coupling parameter should be {self.True1} and not {methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 * (1 + 1e-15)
        methodVal = coupling_parameter(
            self.T,
            self.n_e,
            self.particles,
            z_mean=np.nan * u.dimensionless_unscaled,
            V=np.nan * u.m / u.s,
            method="classical",
        )
        testTrue = not np.isclose(methodVal, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"Coupling parameter value test gives {methodVal} and "
            f"should not be equal to {fail1}."
        )
        assert testTrue, errStr

    def test_zmean(self):
        """
        Test value obtained when arbitrary z_mean is passed
        """
        methodVal = coupling_parameter(
            self.T,
            self.n_e,
            self.particles,
            z_mean=self.z_mean,
            V=np.nan * u.m / u.s,
            method="classical",
        )
        testTrue = np.isclose(self.True_zmean, methodVal, rtol=1e-1, atol=0.0)
        errStr = f"Coupling parameter should be {self.True_zmean} and not {methodVal}."
        assert testTrue, errStr

    # TODO vector z_mean
    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"method": "classical"},
            {"method": "quantum"},
        ],
    )
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans, kwargs):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            coupling_parameter, insert_some_nans, insert_all_nans, kwargs
        )

    def test_quantum(self):
        """
        Testing quantum method for coupling parameter.
        """
        methodVal = coupling_parameter(
            self.T, self.n_e, self.particles, method="quantum"
        )
        testTrue = np.isclose(self.True_quantum, methodVal, rtol=1e-1, atol=0.0)
        errStr = (
            f"Coupling parameter should be {self.True_quantum} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_kwarg_method_error(self):
        """Testing kwarg `method` fails is not 'classical' or 'quantum'"""
        with pytest.raises(ValueError):
            coupling_parameter(self.T, self.n_e, self.particles, method="not a method")


class Test_Knudsen_number:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.length = 1 * u.nm
        cls.T = 11604 * u.K
        cls.n_e = 1e17 * u.cm**-3
        cls.particles = ("e", "p")
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.V = 1e4 * u.km / u.s
        cls.True1 = 440.4757187793204

    def test_symmetry(self):
        with pytest.warns(CouplingWarning):
            result = Knudsen_number(self.length, self.T, self.n_e, self.particles)
            resultRev = Knudsen_number(
                self.length, self.T, self.n_e, self.particles[::-1]
            )
        assert result == resultRev

    def test_known1(self):
        """
        Test for known value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Knudsen_number(
                self.length,
                self.T,
                self.n_e,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = np.isclose(self.True1, methodVal, rtol=1e-1, atol=0.0)
        errStr = f"Knudsen number should be {self.True1} and not {methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 * (1 + 1e-15)
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Knudsen_number(
                self.length,
                self.T,
                self.n_e,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = not np.isclose(methodVal, fail1, rtol=0.0, atol=1e-16)
        errStr = (
            f"Knudsen number value test gives {methodVal} and "
            f"should not be equal to {fail1}."
        )
        assert testTrue, errStr

    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(Knudsen_number, insert_some_nans, insert_all_nans, {})
