import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.collisions.misc import (
    Bethe_stopping,
    Spitzer_resistivity,
    mobility,
)
from plasmapy.formulary.relativity import RelativisticBody
from plasmapy.particles.atomic import stopping_power
from plasmapy.particles.particle_class import Particle
from plasmapy.utils import exceptions
from plasmapy.utils._pytest_helpers import assert_can_handle_nparray
from plasmapy.utils.data.downloader import _API_CONNECTION_ESTABLISHED
from plasmapy.utils.exceptions import CouplingWarning, PhysicsWarning

check_database_connection = pytest.mark.skipif(
    not _API_CONNECTION_ESTABLISHED, reason="failed to connect to data repository"
)


@check_database_connection
@pytest.mark.parametrize(
    (
        "material",
        "rho",
        "I",
        "n_e",
        "T2",
    ),
    [
        # TODO: Add more materials
        (
            "ALUMINUM",
            2.7 * u.g / u.cm**3,
            166 * u.eV,
            7.8 * 10**23 * u.cm**-3,
            1.0 * u.MeV,
        )
    ],
)
def test_Bethe_stopping(material, rho, I, n_e, T2) -> None:  # noqa: E741
    """
    The NIST PSTAR and ASTAR databases should have stopping powers similar
    to those calculated by the Bethe formula beyond a material-dependent
    energy threshold.
    """

    energy_space = np.logspace(start=np.log10(T2.to(u.MeV).value), stop=2) * u.MeV
    particle = RelativisticBody(Particle("p+"), kinetic_energy=energy_space)
    Bethe_stopping_space = Bethe_stopping(I, n_e, particle.velocity, 1)
    _, NIST_stopping_space_density = stopping_power(
        Particle("p+"), material, energy_space, component="electronic"
    )
    NIST_stopping_space = NIST_stopping_space_density * rho

    assert np.isclose(Bethe_stopping_space, NIST_stopping_space).all()


class Test_Spitzer_resistivity:
    @classmethod
    def setup_class(cls) -> None:
        """Initializing parameters for tests"""
        cls.T = 11604 * u.K
        cls.n = 1e12 * u.cm**-3
        cls.particles = ("e", "p")
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.V = 1e4 * u.km / u.s
        cls.True1 = 1.2665402649805445e-3
        cls.True_zmean = 0.00020264644239688712

    def test_symmetry(self) -> None:
        result = Spitzer_resistivity(self.T, self.n, self.particles)
        resultRev = Spitzer_resistivity(self.T, self.n, self.particles[::-1])
        assert result == resultRev

    def test_known1(self) -> None:
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

    def test_fail1(self) -> None:
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

    def test_zmean(self) -> None:
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

    # TODO: vector z_mean
    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans) -> None:
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(
            Spitzer_resistivity, insert_some_nans, insert_all_nans, {}
        )


class Test_mobility:
    @classmethod
    def setup_class(cls) -> None:
        """Initializing parameters for tests"""
        cls.T = 11604 * u.K
        cls.n_e = 1e17 * u.cm**-3
        cls.particles = ("e", "p")
        cls.z_mean = 2.5 * u.dimensionless_unscaled
        cls.V = 1e4 * u.km / u.s
        cls.True1 = 0.13066090887074902
        cls.True_zmean = 0.32665227217687254

    def test_symmetry(self) -> None:
        with pytest.warns(CouplingWarning):  # noqa: PT031
            result = mobility(self.T, self.n_e, self.particles)
            resultRev = mobility(self.T, self.n_e, self.particles[::-1])
        assert result == resultRev

    def test_known1(self) -> None:
        """
        Test for known value.
        """
        with pytest.warns(PhysicsWarning, match="strong coupling effects"):
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

    def test_fail1(self) -> None:
        """
        Tests if test_known1() would fail if we slightly adjusted the
        value comparison by some quantity close to numerical error.
        """
        fail1 = self.True1 * (1 + 1e-15)
        with pytest.warns(PhysicsWarning, match="strong coupling effects"):
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
            f"Mobility value test gives {methodVal} and should not be equal to {fail1}."
        )
        assert testTrue, errStr

    def test_zmean(self) -> None:
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

    # TODO: vector z_mean
    @pytest.mark.parametrize("insert_some_nans", [[], ["V"]])
    @pytest.mark.parametrize("insert_all_nans", [[], ["V"]])
    def test_handle_nparrays(self, insert_some_nans, insert_all_nans) -> None:
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(mobility, insert_some_nans, insert_all_nans, {})
