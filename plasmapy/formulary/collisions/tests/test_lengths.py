import numpy as np
import pytest

from astropy import units as u

from plasmapy.formulary.collisions import coulomb, lengths
from plasmapy.utils import exceptions
from plasmapy.utils.exceptions import CouplingWarning
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray


class Test_impact_parameter_perp:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 11604 * u.K
        cls.particles = ("e", "p")
        cls.V = 1e4 * u.km / u.s
        cls.True1 = 7.200146594293746e-10

    def test_symmetry(self):
        result = lengths.impact_parameter_perp(self.T, self.particles)
        resultRev = lengths.impact_parameter_perp(self.T, self.particles[::-1])
        assert result == resultRev

    def test_known1(self):
        """
        Test for known value.
        """
        methodVal = lengths.impact_parameter_perp(
            self.T, self.particles, V=np.nan * u.m / u.s
        )
        testTrue = np.isclose(self.True1, methodVal.si.value, rtol=1e-1, atol=0.0)
        errStr = (
            "Distance of closest approach for 90 degree Coulomb "
            f"collision, impact_parameter_perp, should be {self.True1} and "
            f"not {methodVal}."
        )
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 + 1e-15
        methodVal = lengths.impact_parameter_perp(
            self.T, self.particles, V=np.nan * u.m / u.s
        )
        testTrue = not np.isclose(methodVal.si.value, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"impact_parameter_perp value test gives {methodVal} and "
            f"should not be equal to {fail1}."
        )
        assert testTrue, errStr

    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            lengths.impact_parameter_perp, insert_some_nans, insert_all_nans, {}
        )

    assert np.isclose(
        coulomb.Coulomb_logarithm(1 * u.eV, 5 * u.m**-3, ("e", "e")),
        coulomb.Coulomb_logarithm(11604.5220 * u.K, 5 * u.m**-3, ("e", "e")),
    )


class Test_impact_parameter:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 11604 * u.K
        cls.T_arr = np.array([1, 2]) * u.eV
        cls.n_e = 1e17 * u.cm**-3
        cls.n_e_arr = np.array([1e17, 2e17]) * u.cm**-3
        cls.particles = ("e", "p")
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.V = 1e4 * u.km / u.s
        cls.True1 = np.array([7.200146594293746e-10, 2.3507660003984624e-08])

    def test_symmetry(self):
        result = lengths.impact_parameter(self.T, self.n_e, self.particles)
        resultRev = lengths.impact_parameter(self.T, self.n_e, self.particles[::-1])
        assert result == resultRev

    def test_known1(self):
        """
        Test for known value.
        """
        methodVal = lengths.impact_parameter(
            self.T,
            self.n_e,
            self.particles,
            z_mean=np.nan * u.dimensionless_unscaled,
            V=np.nan * u.m / u.s,
            method="classical",
        )
        bmin, bmax = methodVal
        methodVal = bmin.si.value, bmax.si.value
        testTrue = np.allclose(self.True1, methodVal, rtol=1e-1, atol=0.0)
        errStr = f"Impact parameters should be {self.True1} and not {methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 * (1 + 1e-15)
        methodVal = lengths.impact_parameter(
            self.T,
            self.n_e,
            self.particles,
            z_mean=np.nan * u.dimensionless_unscaled,
            V=np.nan * u.m / u.s,
            method="classical",
        )
        bmin, bmax = methodVal
        methodVal = bmin.si.value, bmax.si.value
        testTrue = not np.allclose(methodVal, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"Impact parameter value test gives {methodVal} and "
            f"should not be equal to {fail1}."
        )
        assert testTrue, errStr

    def test_bad_method(self):
        """Testing failure when invalid method is passed."""
        with pytest.raises(ValueError):
            lengths.impact_parameter(
                self.T,
                self.n_e,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="meow",
            )

    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"method": "classical"},
            {"method": "GMS-1"},
            {"method": "GMS-2", "z_mean": 1.0 * u.dimensionless_unscaled},
            {"method": "GMS-3"},
            {"method": "GMS-4"},
            {"method": "GMS-5", "z_mean": 1.0 * u.dimensionless_unscaled},
            {"method": "GMS-6", "z_mean": 1.0 * u.dimensionless_unscaled},
        ],
    )
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans, kwargs):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            lengths.impact_parameter, insert_some_nans, insert_all_nans, kwargs
        )

    @pytest.mark.parametrize(
        "n_e_shape,T_shape",
        # Scalar T
        [
            ((2, 3, 5), (1,)),
            # Scalar n
            ((1,), (2, 3, 5)),
            # Both arrays of equal size
            ((2, 3, 5), (2, 3, 5)),
            # Higher dimensional test
            ((2, 3, 5, 4, 2), (2, 3, 5, 4, 2)),
        ],
    )
    def test_extend_output_for_array_input(self, n_e_shape, T_shape):
        """
        Test to verify that if either/or T and n_e are arrays, the
        resulting bmin and bmax have the correct shapes.

        This is necessary in addition to test_handle_nparrays to ensure
        that the output arrays are extended correctly.
        """

        output_shape = T_shape if len(T_shape) >= len(n_e_shape) else n_e_shape

        n_e = self.n_e * np.ones(n_e_shape)
        T = self.T * np.ones(T_shape)

        bmin, bmax = lengths.impact_parameter(T, n_e, self.particles)

        msg = f"wrong shape for {n_e.shape = } and {T.shape = }"

        assert bmin.shape == output_shape, f"Bmin {msg}"
        assert bmax.shape == output_shape, f"Bmax {msg}"


class Test_mean_free_path:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.T = 11604 * u.K
        cls.n_e = 1e17 * u.cm**-3
        cls.particles = ("e", "p")
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.V = 1e4 * u.km / u.s
        cls.True1 = 4.4047571877932046e-07

    def test_symmetry(self):
        with pytest.warns(CouplingWarning):
            result = lengths.mean_free_path(self.T, self.n_e, self.particles)
            resultRev = lengths.mean_free_path(self.T, self.n_e, self.particles[::-1])
        assert result == resultRev

    def test_known1(self):
        """
        Test for known value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = lengths.mean_free_path(
                self.T,
                self.n_e,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = np.isclose(self.True1, methodVal.si.value, rtol=1e-1, atol=0.0)
        errStr = f"Mean free path should be {self.True1} and not {methodVal}."
        assert testTrue, errStr

    def test_fail1(self):
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 * (1 + 1e-15)
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = lengths.mean_free_path(
                self.T,
                self.n_e,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="classical",
            )
        testTrue = not np.isclose(methodVal.si.value, fail1, rtol=1e-16, atol=0.0)
        errStr = (
            f"Mean free path value test gives {methodVal} and "
            f"should not be equal to {fail1}."
        )
        assert testTrue, errStr

    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            lengths.mean_free_path, insert_some_nans, insert_all_nans, {}
        )
