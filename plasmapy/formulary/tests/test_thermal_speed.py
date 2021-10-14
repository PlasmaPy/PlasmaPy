"""
Module for testing functionality associated with calculating the
thermal speed.


- `~plasmapy.formulary.parameters.thermal_speed`
- `~plasmapy.formulary.parameters.thermal_speed_lite`
- `~plasmapy.formulary.parameters.thermal_speed_coefficients`
- `~plasmapy.formulary.parameters.vth_`
"""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.parameters import (
    kappa_thermal_speed,
    thermal_speed,
    thermal_speed_coefficients,
    thermal_speed_lite,
    vth_,
    vth_kappa_,
)
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.exceptions import RelativityError, RelativityWarning
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

T_e = 1e6 * u.K
T_i = 1e6 * u.K


@pytest.mark.parametrize(
    "alias, parent",
    [(vth_, thermal_speed), (vth_kappa_, kappa_thermal_speed)],
)
def test_aliases(alias, parent):
    assert alias is parent


class TestThermalSpeedCoefficients:
    ...


class TestThermalSpeed:
    """
    Test class for functionality around calculating the thermal speed.  This
    covers the functionality

        - `~plasmapy.formulary.parameters.thermal_speed`
        - `~plasmapy.formulary.parameters.thermal_speed_lite`
        - `~plasmapy.formulary.parameters.thermal_speed_coefficients`
    """
    _bound_attrs = [
        ("lite", thermal_speed_lite),
        ("coefficients", thermal_speed_coefficients),
    ]

    @pytest.mark.parametrize("alias", [vth_])
    def test_aliases(self, alias):
        """Test alias is properly defined."""
        assert alias is thermal_speed

    @pytest.mark.parametrize(
        "bound_name, bound_attr",
        _bound_attrs,
    )
    def test_lite_function_binding(self, bound_name, bound_attr):
        """Test expected attributes are bound correctly."""
        assert hasattr(thermal_speed, bound_name)
        attr = getattr(thermal_speed, bound_name)
        assert attr is bound_attr

    def test_lite_function_marking(self):
        """
        Test thermal_speed is marked as having a Lite-Function and that its
        __has_lite_func__ attribute is defined correctly.
        """
        assert hasattr(thermal_speed, "__has_lite_func__")
        assert len(thermal_speed.__has_lite_func__) == len(self._bound_attrs)

        bound_names, bound_qualnames = zip(*thermal_speed.__has_lite_func__)
        for name, attr in self._bound_attrs:
            assert name in bound_names

            index = bound_names.index(name)
            qualname = f"{attr.__module__}.{attr.__name__}"
            assert qualname == bound_qualnames[index]

    @pytest.mark.parametrize(
        "method, ndim, expected",
        [
            ("most_probable", 1, 0),
            ("most_probable", 2, 1),
            ("most_probable", 3, np.sqrt(2)),
            ("rms", 1, 1),
            ("rms", 2, np.sqrt(2)),
            ("rms", 3, np.sqrt(3)),
            ("mean_magnitude", 1, np.sqrt(2 / np.pi)),
            ("mean_magnitude", 2, np.sqrt(np.pi / 2)),
            ("mean_magnitude", 3, np.sqrt(8 / np.pi)),
            ("nrl", 1, 1),
            ("nrl", 2, 1),
            ("nrl", 3, 1),
        ],
    )
    def test_thermal_speed_coefficient_values(self, method, ndim, expected):
        """Test values returned by thermal_speed_coefficients."""
        assert np.isclose(
            thermal_speed_coefficients(method=method, ndim=ndim), expected
        )

    @pytest.mark.parametrize(
        "method, ndim, _raises",
        [
            ("most_probably", -1, ValueError),
            ("most_probably", 4, ValueError),
            ("most_probably", "not an int", ValueError),
            ("wrong method", 3, ValueError),
            (5, 1, ValueError),
            ({"wrong": 1}, "wrong", TypeError),
        ],
    )
    def test_thermal_speed_coefficients_raises(self, method, ndim, _raises):
        """Test raise Exception cases for thermal_speed_coefficients."""
        with pytest.raises(_raises):
            thermal_speed_coefficients(method=method, ndim=ndim)

    def test_thermal_speed_lite(self):
        ...

    @pytest.mark.parametrize(
        "inputs",
        [
            dict(T=5 * u.eV, particle=Particle("p"), method= "most_probable", ndim=3),
            dict(T=3000 * u.K, particle=Particle("e"), method="nrl", ndim=2),
            dict(
                T=5000 * u.K, particle=Particle("He+"), method="mean_magnitude", ndim=1
            ),
            dict(T=1 * u.eV, particle=Particle("Ar+"), method="rms", ndim=3),
        ],
    )
    def test_normal_vs_lite_values(self, inputs):
        """
        Test that thermal_speed and thermal_speed_lite calulate the same values
        for the same inputs.
        """
        T_unitless = inputs["T"].to(u.K, equivalencies=u.temperature_energy()).value
        m_unitless = inputs["particle"].mass.value

        normal = thermal_speed(**inputs)
        coeff = thermal_speed_coefficients(method=inputs["method"], ndim=inputs["ndim"])
        lite = thermal_speed_lite(T=T_unitless, mass=m_unitless, coeff=coeff)
        assert np.isclose(normal.value, lite)


def test_thermal_speed():
    r"""Test the thermal_speed function in parameters.py"""
    assert thermal_speed(T_e, "e-").unit.is_equivalent(u.m / u.s)

    assert thermal_speed(T_e, "e-") > thermal_speed(T_e, "p")

    # The NRL Plasma Formulary uses a definition of the electron
    # thermal speed that differs by a factor of sqrt(2).
    assert np.isclose(thermal_speed(1 * u.MK, "e-").value, 5505694.743141063)

    with pytest.raises(u.UnitTypeError):
        thermal_speed(5 * u.m, "e-")

    with pytest.raises(ValueError):
        thermal_speed(-T_e, "e-")

    with pytest.warns(RelativityWarning):
        thermal_speed(1e9 * u.K, "e-")

    with pytest.raises(RelativityError):
        thermal_speed(5e19 * u.K, "e-")

    with pytest.warns(u.UnitsWarning):
        assert thermal_speed(1e5, "e-") == thermal_speed(1e5 * u.K, "e-")

    assert thermal_speed(T_i, particle="p").unit.is_equivalent(u.m / u.s)

    # The NRL Plasma Formulary uses a definition of the particle thermal
    # speed that differs by a factor of sqrt(2).
    assert np.isclose(
        thermal_speed(1 * u.MK, particle="p").si.value, 128486.56960876315
    )

    # Explicitly check all three modes and dimensionalities
    # ndim = 1
    assert np.isclose(
        thermal_speed(T_e, "e-", method="most_probable", ndim=1).si.value, 0.0
    )

    # Regression tests start here!
    assert np.isclose(
        thermal_speed(T_e, "e-", method="rms", ndim=1).si.value, 3893114.2008620175
    )

    assert np.isclose(
        thermal_speed(T_e, "e-", method="mean_magnitude", ndim=1).si.value,
        3106255.714310189,
    )

    # ndim = 2
    assert np.isclose(
        thermal_speed(T_e, "e-", method="most_probable", ndim=2).si.value,
        3893114.2008620175,
    )

    assert np.isclose(
        thermal_speed(T_e, "e-", method="rms", ndim=2).si.value, 5505694.902726359
    )

    assert np.isclose(
        thermal_speed(T_e, "e-", method="mean_magnitude", ndim=2).si.value,
        4879295.066124102,
    )

    # ndim = 3
    assert np.isclose(
        thermal_speed(T_e, "e-", method="most_probable", ndim=3).si.value,
        5505694.902726359,
    )

    assert np.isclose(
        thermal_speed(T_e, "e-", method="rms", ndim=3).si.value, 6743071.595560921
    )

    assert np.isclose(
        thermal_speed(T_e, "e-", method="mean_magnitude", ndim=3).si.value,
        6212511.428620378,
    )

    # Case when Z=1 is assumed
    assert thermal_speed(T_i, particle="p") == thermal_speed(T_i, particle="H-1+")

    assert thermal_speed(1 * u.MK, particle="e+") == thermal_speed(1 * u.MK, "e-")

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
        thermal_speed(1e6 * u.K, "e-", method="mean_magnitude").si.value,
        6212510.3969422,
    )

    assert np.isclose(
        thermal_speed(1e6 * u.K, "e-", method="rms").si.value, 6743070.475775486
    )

    # Test invalid method
    with pytest.raises(ValueError):
        thermal_speed(T_i, "e-", method="sadks")

    # Test invalid ndim
    with pytest.raises(ValueError):
        thermal_speed(T_i, "e-", ndim=4)

    assert_can_handle_nparray(thermal_speed)


# test class for kappa_thermal_speed() function:
class Test_kappa_thermal_speed(object):
    @classmethod
    def setup_class(self):
        """initializing parameters for tests"""
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
