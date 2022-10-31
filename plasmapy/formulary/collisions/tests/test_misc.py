import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.collisions.misc import mobility, Spitzer_resistivity
from plasmapy.utils import exceptions
from plasmapy.utils.exceptions import CouplingWarning
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray


class Test_Spitzer_resistivity:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 11604 * u.K
        cls.n = 1e12 * u.cm**-3
        cls.particles = ("e", "p")
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.V = 1e4 * u.km / u.s
        cls.True1 = 1.2665402649805445e-3
        cls.True_zmean = 0.00020264644239688712

    def test_symmetry(self):
        result = Spitzer_resistivity(self.T, self.n, self.particles)
        resultRev = Spitzer_resistivity(self.T, self.n, self.particles[::-1])
        assert result == resultRev

    def test_known1(self):
        """
        Test for known value.
        """
        methodVal = Spitzer_resistivity(
            self.T,
            self.n,
            self.particles,
            z_mean=np.nan * u.dimensionless_unscaled,
            V=np.nan * u.m / u.s,
            method="classical",
        )
        testTrue = np.isclose(self.True1, methodVal.si.value, rtol=1e-1, atol=0.0)
        errStr = f"Spitzer resistivity should be {self.True1} and not {methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 * (1 + 1e-15)
        methodVal = Spitzer_resistivity(
            self.T,
            self.n,
            self.particles,
            z_mean=np.nan * u.dimensionless_unscaled,
            V=np.nan * u.m / u.s,
            method="classical",
        )
        testTrue = not np.isclose(methodVal.si.value, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"Spitzer resistivity value test gives {methodVal} and "
            f"should not be equal to {fail1}."
        )
        assert testTrue, errStr

    def test_zmean(self):
        """Testing Spitzer when z_mean is passed."""
        methodVal = Spitzer_resistivity(
            self.T,
            self.n,
            self.particles,
            z_mean=self.z_mean,
            V=np.nan * u.m / u.s,
            method="classical",
        )
        testTrue = np.isclose(self.True_zmean, methodVal.si.value, rtol=1e-1, atol=0.0)
        errStr = f"Spitzer resistivity should be {self.True_zmean} and not {methodVal}."
        assert testTrue, errStr

    # TODO vector z_mean
    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            Spitzer_resistivity, insert_some_nans, insert_all_nans, {}
        )


class Test_mobility:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 11604 * u.K
        cls.n_e = 1e17 * u.cm**-3
        cls.particles = ("e", "p")
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.V = 1e4 * u.km / u.s
        cls.True1 = 0.13066090887074902
        cls.True_zmean = 0.32665227217687254

    def test_symmetry(self):
        with pytest.warns(CouplingWarning):
            result = mobility(self.T, self.n_e, self.particles)
            resultRev = mobility(self.T, self.n_e, self.particles[::-1])
        assert result == resultRev

    def test_known1(self):
        """
        Test for known value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = mobility(
                self.T,
                self.n_e,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = np.isclose(self.True1, methodVal.si.value, rtol=1e-1, atol=0.0)
        errStr = f"Mobility should be {self.True1} and not {methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 * (1 + 1e-15)
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = mobility(
                self.T,
                self.n_e,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = not np.isclose(methodVal.si.value, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"Mobility value test gives {methodVal} and "
            f"should not be equal to {fail1}."
        )
        assert testTrue, errStr

    def test_zmean(self):
        """Testing mobility when z_mean is passed."""
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = mobility(
                self.T,
                self.n_e,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = np.isclose(self.True_zmean, methodVal.si.value, rtol=1e-1, atol=0.0)
        errStr = f"Mobility should be {self.True_zmean} and not {methodVal}."
        assert testTrue, errStr

    # TODO vector z_mean
    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(mobility, insert_some_nans, insert_all_nans, {})
