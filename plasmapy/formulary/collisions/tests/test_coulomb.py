import numpy as np
import pytest

from astropy import units as u
from astropy.constants import c
from astropy.tests.helper import assert_quantity_allclose

import plasmapy.particles.exceptions

from plasmapy.formulary.collisions.coulomb import Coulomb_logarithm
from plasmapy.utils import exceptions
from plasmapy.utils.exceptions import CouplingWarning
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray


class Test_Coulomb_logarithm:
    @classmethod
    def setup_class(cls):
        """initializing parameters for tests"""
        cls.temperature1 = 10 * 11604 * u.K
        cls.T_arr = np.array([1, 2]) * u.eV
        cls.density1 = 1e20 * u.cm**-3
        cls.n_arr = np.array([1e20, 2e20]) * u.cm**-3
        cls.temperature2 = 1 * 11604 * u.K
        cls.density2 = 1e23 * u.cm**-3
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.particles = ("e", "p")
        cls.ls_min_interp = 3.4014290066940966
        cls.gms1 = 3.4014290066940966
        cls.ls_min_interp_negative = -3.4310536971592493
        cls.gms1_negative = -3.4310536971592493
        cls.ls_full_interp = 3.6349941014645157
        cls.gms2 = 3.6349941014645157
        cls.ls_full_interp_negative = -1.379394033464292
        cls.gms2_negative = -1.379394033464292
        cls.ls_clamp_mininterp = 3.4014290066940966
        cls.gms3 = 3.4014290066940966
        cls.ls_clamp_mininterp_negative = 2
        cls.gms3_negative = 2
        cls.ls_clamp_mininterp_non_scalar = (2, 2)
        cls.gms3_non_scalar = (2, 2)
        cls.hls_min_interp = 3.401983996820073
        cls.gms4 = 3.401983996820073
        cls.hls_min_interp_negative = 0.0005230791851781715
        cls.gms4_negative = 0.0005230791851781715
        cls.hls_max_interp = 3.7196690506837693
        cls.gms5 = 3.7196690506837693
        cls.hls_max_interp_negative = 0.03126832674323108
        cls.gms5_negative = 0.03126832674323108
        cls.hls_full_interp = 3.635342040477818
        cls.gms6 = 3.635342040477818
        cls.hls_full_interp_negative = 0.030720859361047514
        cls.gms6_negative = 0.030720859361047514

    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"method": "classical"},
            {"method": "ls"},
            {"method": "ls_min_interp"},
            {"method": "GMS-1"},
            {"method": "ls_full_interp", "z_mean": 1.0 * u.dimensionless_unscaled},
            {"method": "GMS-2", "z_mean": 1.0 * u.dimensionless_unscaled},
            {"method": "ls_clamp_mininterp"},
            {"method": "GMS-3"},
            {"method": "hls_min_interp"},
            {"method": "GMS-4"},
            {"method": "hls_max_interp", "z_mean": 1.0 * u.dimensionless_unscaled},
            {"method": "GMS-5", "z_mean": 1.0 * u.dimensionless_unscaled},
            {"method": "hls_full_interp", "z_mean": 1.0 * u.dimensionless_unscaled},
            {"method": "GMS-6", "z_mean": 1.0 * u.dimensionless_unscaled},
        ],
    )
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans, kwargs):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            Coulomb_logarithm, insert_some_nans, insert_all_nans, kwargs
        )

    def test_unknown_method(self):
        """Test that function will raise ValueError on non-existent method"""
        with pytest.raises(ValueError):
            Coulomb_logarithm(
                self.T_arr[0],
                self.n_arr[0],
                self.particles,
                method="welcome our new microsoft overlords",
            )

    def test_handle_invalid_V(self):
        """Test that V default, V = None, and V = np.nan all give the same result"""
        with pytest.warns(CouplingWarning):
            methodVal_0 = Coulomb_logarithm(
                self.T_arr[0],
                self.n_arr[0],
                self.particles,
                z_mean=1 * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
            )
            methodVal_1 = Coulomb_logarithm(
                self.T_arr[0],
                self.n_arr[0],
                self.particles,
                z_mean=1 * u.dimensionless_unscaled,
                V=None,
            )
            methodVal_2 = Coulomb_logarithm(
                self.T_arr[0],
                self.n_arr[0],
                self.particles,
                z_mean=1 * u.dimensionless_unscaled,
            )
        assert_quantity_allclose(methodVal_0, methodVal_1)
        assert_quantity_allclose(methodVal_0, methodVal_2)

    def test_handle_zero_V(self):
        """Test that V == 0 returns a PhysicsError"""
        with pytest.raises(exceptions.PhysicsError):
            Coulomb_logarithm(
                self.T_arr[0],
                self.n_arr[0],
                self.particles,
                z_mean=1 * u.dimensionless_unscaled,
                V=0 * u.m / u.s,
            )

    def test_handle_V_arraysizes(self):
        """Test that different sized V input array gets handled by _boilerplate"""
        with pytest.warns(CouplingWarning):
            methodVal_0 = Coulomb_logarithm(
                self.T_arr[0],
                self.n_arr[0],
                self.particles,
                z_mean=1 * u.dimensionless_unscaled,
                V=np.array([np.nan, 3e7]) * u.m / u.s,
            )
            methodVal_1 = Coulomb_logarithm(
                self.T_arr[1],
                self.n_arr[0],
                self.particles,
                z_mean=1 * u.dimensionless_unscaled,
                V=np.array([1e7, np.nan]) * u.m / u.s,
            )
            methodVal_2 = Coulomb_logarithm(
                self.T_arr,
                self.n_arr[0],
                self.particles,
                z_mean=1 * u.dimensionless_unscaled,
                V=np.array([np.nan, np.nan]) * u.m / u.s,
            )
        assert_quantity_allclose(methodVal_0[0], methodVal_2[0])
        assert_quantity_allclose(methodVal_1[1], methodVal_2[1])

    def test_symmetry(self):
        with pytest.warns(CouplingWarning):
            lnLambda = Coulomb_logarithm(
                self.temperature1, self.density2, self.particles
            )
            lnLambdaRev = Coulomb_logarithm(
                self.temperature1, self.density2, self.particles[::-1]
            )
        assert lnLambda == lnLambdaRev

    def test_Chen_Q_machine(self):
        """
        Tests whether Coulomb logarithm gives value consistent with
        Chen's Introduction to Plasma Physics and Controlled Fusion
        section 5.6.2 Q-machine example.
        """
        T = 0.2 * u.eV
        T = T.to(u.K, equivalencies=u.temperature_energy())
        n = 1e15 * u.m**-3
        # factor of np.log(2) corrects for different definitions of thermal
        # velocity. Chen uses v**2 = k * T / m  whereas we use
        # v ** 2 = 2 * k * T / m
        lnLambdaChen = 9.1 + np.log(2)
        lnLambda = Coulomb_logarithm(T, n, ("e", "p"))
        testTrue = np.isclose(lnLambda, lnLambdaChen, rtol=1e-1, atol=0.0)
        errStr = (
            "Q-machine value of Coulomb logarithm should be "
            f"{lnLambdaChen} and not {lnLambda}."
        )
        assert testTrue, errStr

    def test_Chen_lab(self):
        """
        Tests whether Coulomb logarithm gives value consistent with
        Chen's Introduction to Plasma Physics and Controlled Fusion
        section 5.6.2 lab plasma example.
        """
        T = 2 * u.eV
        T = T.to(u.K, equivalencies=u.temperature_energy())
        n = 1e17 * u.m**-3
        # factor of np.log(2) corrects for different definitions of thermal
        # velocity. Chen uses v**2 = k * T / m  whereas we use
        # v ** 2 = 2 * k * T / m
        lnLambdaChen = 10.2 + np.log(2)
        lnLambda = Coulomb_logarithm(T, n, ("e", "p"))
        testTrue = np.isclose(lnLambda, lnLambdaChen, rtol=1e-1, atol=0.0)
        errStr = (
            "Lab plasma value of Coulomb logarithm should be "
            f"{lnLambdaChen} and not {lnLambda}."
        )
        assert testTrue, errStr

    def test_Chen_torus(self):
        """
        Tests whether Coulomb logarithm gives value consistent with
        Chen's Introduction to Plasma Physics and Controlled Fusion
        section 5.6.2 torus example.
        """
        T = 100 * u.eV
        T = T.to(u.K, equivalencies=u.temperature_energy())
        n = 1e19 * u.m**-3
        # factor of np.log(2) corrects for different definitions of thermal
        # velocity. Chen uses v**2 = k * T / m  whereas we use
        # v ** 2 = 2 * k * T / m
        lnLambdaChen = 13.7 + np.log(2)
        lnLambda = Coulomb_logarithm(T, n, ("e", "p"))
        testTrue = np.isclose(lnLambda, lnLambdaChen, rtol=1e-1, atol=0.0)
        errStr = (
            "Torus value of Coulomb logarithm should be "
            f"{lnLambdaChen} and not {lnLambda}."
        )
        assert testTrue, errStr

    def test_Chen_fusion(self):
        """
        Tests whether Coulomb logarithm gives value consistent with
        Chen's Introduction to Plasma Physics and Controlled Fusion
        section 5.6.2 fusion reactor example.
        """
        T = 1e4 * u.eV
        T = T.to(u.K, equivalencies=u.temperature_energy())
        n = 1e21 * u.m**-3
        # factor of np.log(2) corrects for different definitions of thermal
        # velocity. Chen uses v**2 = k * T / m  whereas we use
        # v ** 2 = 2 * k * T / m
        lnLambdaChen = 16 + np.log(2)
        with pytest.warns(exceptions.RelativityWarning):
            lnLambda = Coulomb_logarithm(T, n, ("e", "p"))
        testTrue = np.isclose(lnLambda, lnLambdaChen, rtol=1e-1, atol=0.0)
        errStr = (
            "Fusion reactor value of Coulomb logarithm should be "
            f"{lnLambdaChen} and not {lnLambda}."
        )
        assert testTrue, errStr

    def test_Chen_laser(self):
        """
        Tests whether Coulomb logarithm gives value consistent with
        Chen's Introduction to Plasma Physics and Controlled Fusion
        section 5.6.2 laser plasma example.
        """
        T = 1e3 * u.eV
        T = T.to(u.K, equivalencies=u.temperature_energy())
        n = 1e27 * u.m**-3
        # factor of np.log(2) corrects for different definitions of thermal
        # velocity. Chen uses v**2 = k * T / m  whereas we use
        # v ** 2 = 2 * k * T / m
        lnLambdaChen = 6.8 + np.log(2)
        with pytest.warns(exceptions.RelativityWarning):
            lnLambda = Coulomb_logarithm(T, n, ("e", "p"))
        testTrue = np.isclose(lnLambda, lnLambdaChen, rtol=1e-1, atol=0.0)
        errStr = (
            "Laser plasma value of Coulomb logarithm should be "
            f"{lnLambdaChen} and not {lnLambda}."
        )
        assert testTrue, errStr

    def test_ls_min_interp(self):
        """
        Test for first version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="ls_min_interp",
            )
        testTrue = np.isclose(methodVal, self.ls_min_interp, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for ls_min_interp should be "
            f"{self.ls_min_interp} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS1(self):
        """
        Test for first version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="GMS-1",
            )
        testTrue = np.isclose(methodVal, self.gms1, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-1 should be "
            f"{self.gms1} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_ls_min_interp_negative(self):
        """
        Test for first version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks for when
        a negative (invalid) Coulomb logarithm is returned.
        """
        with pytest.warns(exceptions.CouplingWarning, match="depends on weak coupling"):
            Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="ls_min_interp",
            )

    def test_GMS1_negative(self):
        """
        Test for first version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks for when
        a negative (invalid) Coulomb logarithm is returned.
        """
        with pytest.warns(exceptions.CouplingWarning, match="depends on weak coupling"):
            Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=np.nan * u.dimensionless_unscaled,
                V=np.nan * u.m / u.s,
                method="GMS-1",
            )

    def test_ls_full_interp(self):
        """
        Test for second version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="ls_full_interp",
            )
        testTrue = np.isclose(methodVal, self.ls_full_interp, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-2 should be "
            f"{self.ls_full_interp} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS2(self):
        """
        Test for second version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-2",
            )
        testTrue = np.isclose(methodVal, self.gms2, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-2 should be "
            f"{self.gms2} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_ls_full_interp_negative(self):
        """
        Test for second version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks for when
        a negative (invalid) Coulomb logarithm is returned.
        """
        with pytest.warns(exceptions.CouplingWarning, match="depends on weak coupling"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="ls_full_interp",
            )

    def test_GMS2_negative(self):
        """
        Test for second version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks for when
        a negative (invalid) Coulomb logarithm is returned.
        """
        with pytest.warns(exceptions.CouplingWarning, match="depends on weak coupling"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-2",
            )

    def test_ls_clamp_mininterp(self):
        """
        Test for third version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="ls_clamp_mininterp",
            )
        testTrue = np.isclose(methodVal, self.ls_clamp_mininterp, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for ls_clamp_minterp should be "
            f"{self.ls_clamp_mininterp} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS3(self):
        """
        Test for third version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-3",
            )
        testTrue = np.isclose(methodVal, self.gms3, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-3 should be "
            f"{self.gms3} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_ls_clamp_mininterp_negative(self):
        """
        Test for third version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        a positive value is returned whereas the classical Coulomb
        logarithm would return a negative value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="ls_clamp_mininterp",
            )
        testTrue = np.isclose(
            methodVal, self.ls_clamp_mininterp_negative, rtol=1e-6, atol=0.0
        )
        errStr = (
            f"Coulomb logarithm for GMS-3 should be "
            f"{self.ls_clamp_mininterp_negative} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS3_negative(self):
        """
        Test for third version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        a positive value is returned whereas the classical Coulomb
        logarithm would return a negative value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-3",
            )
        testTrue = np.isclose(methodVal, self.gms3_negative, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-3 should be "
            f"{self.gms3_negative} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_ls_clamp_mininterp_non_scalar_density(self):
        """
        Test for third version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        passing in a collection of density values returns a
        collection of Coulomb logarithm values.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                10 * 1160 * u.K,
                (1e23 * u.cm**-3, 1e20 * u.cm**-3),
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="ls_clamp_mininterp",
            )
        testTrue = np.isclose(
            methodVal, self.ls_clamp_mininterp_non_scalar, rtol=1e-6, atol=0.0
        )
        errStr = (
            f"Coulomb logarithm for GMS-3 should be "
            f"{self.ls_clamp_mininterp_non_scalar} and not {methodVal}."
        )
        assert testTrue.all(), errStr

    def test_GMS3_non_scalar_density(self):
        """
        Test for third version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        passing in a collection of density values returns a
        collection of Coulomb logarithm values.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                10 * 1160 * u.K,
                (1e23 * u.cm**-3, 1e20 * u.cm**-3),
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-3",
            )
        testTrue = np.isclose(methodVal, self.gms3_non_scalar, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-3 should be "
            f"{self.gms3_non_scalar} and not {methodVal}."
        )
        assert testTrue.all(), errStr

    def test_hls_min_interp(self):
        """
        Test for fourth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="hls_min_interp",
            )
        testTrue = np.isclose(methodVal, self.hls_min_interp, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-4 should be "
            f"{self.hls_min_interp} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS4(self):
        """
        Test for fourth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-4",
            )
        testTrue = np.isclose(methodVal, self.gms4, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-4 should be "
            f"{self.gms4} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_hls_min_interp_negative(self):
        """
        Test for fourth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        a positive value is returned whereas the classical Coulomb
        logarithm would return a negative value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="hls_min_interp",
            )
        testTrue = np.isclose(
            methodVal, self.hls_min_interp_negative, rtol=1e-6, atol=0.0
        )
        errStr = (
            f"Coulomb logarithm for GMS-4 should be "
            f"{self.hls_min_interp_negative} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS4_negative(self):
        """
        Test for fourth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        a positive value is returned whereas the classical Coulomb
        logarithm would return a negative value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-4",
            )
        testTrue = np.isclose(methodVal, self.gms4_negative, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-4 should be "
            f"{self.gms4_negative} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_hls_max_interp(self):
        """
        Test for fifth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="hls_max_interp",
            )
        testTrue = np.isclose(methodVal, self.hls_max_interp, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-5 should be "
            f"{self.hls_max_interp} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS5(self):
        """
        Test for fifth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-5",
            )
        testTrue = np.isclose(methodVal, self.gms5, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-5 should be "
            f"{self.gms5} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_hls_max_interp_negative(self):
        """
        Test for fifth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        a positive value is returned whereas the classical Coulomb
        logarithm would return a negative value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="hls_max_interp",
            )
        testTrue = np.isclose(
            methodVal, self.hls_max_interp_negative, rtol=1e-6, atol=0.0
        )
        errStr = (
            f"Coulomb logarithm for GMS-5 should be "
            f"{self.hls_max_interp_negative} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS5_negative(self):
        """
        Test for fifth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        a positive value is returned whereas the classical Coulomb
        logarithm would return a negative value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-5",
            )
        testTrue = np.isclose(methodVal, self.gms5_negative, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-5 should be "
            f"{self.gms5_negative} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_hls_full_interp(self):
        """
        Test for sixth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="hls_full_interp",
            )
        testTrue = np.isclose(methodVal, self.hls_full_interp, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-6 should be "
            f"{self.hls_full_interp} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS6(self):
        """
        Test for sixth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002).
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature1,
                self.density1,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-6",
            )
        testTrue = np.isclose(methodVal, self.gms6, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-6 should be "
            f"{self.gms6} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_hls_full_interp_negative(self):
        """
        Test for sixth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        a positive value is returned whereas the classical Coulomb
        logarithm would return a negative value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="hls_full_interp",
            )
        testTrue = np.isclose(
            methodVal, self.hls_full_interp_negative, rtol=1e-6, atol=0.0
        )
        errStr = (
            f"Coulomb logarithm for GMS-6 should be "
            f"{self.hls_full_interp_negative} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_GMS6_negative(self):
        """
        Test for sixth version of Coulomb logarithm from Gericke,
        Murillo, and Schlanges PRE (2002). This checks whether
        a positive value is returned whereas the classical Coulomb
        logarithm would return a negative value.
        """
        with pytest.warns(exceptions.PhysicsWarning, match="strong coupling effects"):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                z_mean=self.z_mean,
                V=np.nan * u.m / u.s,
                method="GMS-6",
            )
        testTrue = np.isclose(methodVal, self.gms6_negative, rtol=1e-6, atol=0.0)
        errStr = (
            f"Coulomb logarithm for GMS-6 should be "
            f"{self.gms6_negative} and not {methodVal}."
        )
        assert testTrue, errStr

    def test_ls_full_interp_zmean_error(self):
        """
        Tests whether ls_full_interp raises z_mean error when a z_mean is not
        provided.
        """
        with pytest.raises(ValueError):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                method="ls_full_interp",
            )

    def test_GMS2_zmean_error(self):
        """
        Tests whether GMS-2 raises z_mean error when a z_mean is not
        provided.
        """
        with pytest.raises(ValueError):
            methodVal = Coulomb_logarithm(
                self.temperature2, self.density2, self.particles, method="GMS-2"
            )

    def test_hls_max_interp_zmean_error(self):
        """
        Tests whether hls_max_interp raises z_mean error when a z_mean is not
        provided.
        """
        with pytest.raises(ValueError):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                method="hls_max_interp",
            )

    def test_GMS5_zmean_error(self):
        """
        Tests whether GMS-5 raises z_mean error when a z_mean is not
        provided.
        """
        with pytest.raises(ValueError):
            methodVal = Coulomb_logarithm(
                self.temperature2, self.density2, self.particles, method="GMS-5"
            )

    def test_hls_full_interp_zmean_error(self):
        """
        Tests whether hls_full_interp raises z_mean error when a z_mean is not
        provided.
        """
        with pytest.raises(ValueError):
            methodVal = Coulomb_logarithm(
                self.temperature2,
                self.density2,
                self.particles,
                method="hls_full_interp",
            )

    def test_GMS6_zmean_error(self):
        """
        Tests whether GMS-6 raises z_mean error when a z_mean is not
        provided.
        """
        with pytest.raises(ValueError):
            methodVal = Coulomb_logarithm(
                self.temperature2, self.density2, self.particles, method="GMS-6"
            )

    def test_relativity_warn(self):
        """Tests whether relativity warning is raised at high velocity."""
        with pytest.warns(exceptions.RelativityWarning):
            Coulomb_logarithm(1e5 * u.K, 1 * u.m**-3, ("e", "p"), V=0.9 * c)

    def test_relativity_error(self):
        """Tests whether relativity error is raised at light speed."""
        with pytest.raises(exceptions.RelativityError):
            Coulomb_logarithm(1e5 * u.K, 1 * u.m**-3, ("e", "p"), V=1.1 * c)

    def test_unit_conversion_error(self):
        """
        Tests whether unit conversion error is raised when arguments
        are given with incorrect units.
        """
        with pytest.raises(u.UnitTypeError):
            Coulomb_logarithm(
                1e5 * u.g, 1 * u.m**-3, ("e", "p"), V=29979245 * u.m / u.s
            )

    def test_single_particle_error(self):
        """
        Tests whether an error is raised if only a single particle is given.
        """
        with pytest.raises(ValueError):
            Coulomb_logarithm(1 * u.K, 5 * u.m**-3, "e")

    def test_invalid_particle_error(self):
        """
        Tests whether an error is raised when an invalid particle name
        is given.
        """
        with pytest.raises(plasmapy.particles.exceptions.InvalidParticleError):
            Coulomb_logarithm(1 * u.K, 5 * u.m**-3, ("e", "g"))

    n_e = np.array([1e9, 1e9, 1e24]) * u.cm**-3
    T = np.array([1e2, 1e7, 1e8]) * u.K
    Lambda = np.array([5.97, 21.66, 6.69])
    particles = ("e", "p")


# class Test_Coulomb_cross_section:
