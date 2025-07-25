"""
Module for testing functionality associated with calculating the
plasma frequency.

- `~plasmapy.formulary.frequencies.plasma_frequency`
- `~plasmapy.formulary.frequencies.plasma_frequency_lite`
- `~plasmapy.formulary.frequencies.wp_`
"""

import astropy.units as u
import numpy as np
import pytest
from astropy.constants.si import m_p

from plasmapy.formulary.frequencies import plasma_frequency, plasma_frequency_lite, wp_
from plasmapy.particles._factory import _physical_particle_factory
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils._pytest_helpers import assert_can_handle_nparray


@pytest.mark.parametrize(
    ("alias", "parent"),
    [(wp_, plasma_frequency)],
)
def test_aliases(alias, parent) -> None:
    assert alias is parent


class TestPlasmaFrequency:
    """
    Test class for `plasmapy.formulary.frequencies.plasma_frequency`.

    Note: Testing of `plasma_frequency_lite` is done in a separate test
    class.
    """

    @pytest.mark.parametrize(
        ("bound_name", "bound_attr"),
        [("lite", plasma_frequency_lite)],
    )
    def test_lite_function_binding(self, bound_name: str, bound_attr) -> None:
        """Test expected attributes are bound correctly."""
        assert hasattr(plasma_frequency, bound_name)
        assert getattr(plasma_frequency, bound_name) is bound_attr

    def test_lite_function_marking(self) -> None:
        """
        Test plasma_frequency is marked as having a Lite-Function.
        """
        assert hasattr(plasma_frequency, "__bound_lite_func__")
        assert isinstance(plasma_frequency.__bound_lite_func__, dict)

        for bound_name, bound_origin in plasma_frequency.__bound_lite_func__.items():
            assert hasattr(plasma_frequency, bound_name)

            attr = getattr(plasma_frequency, bound_name)
            origin = f"{attr.__module__}.{attr.__name__}"
            assert origin == bound_origin

    @pytest.mark.parametrize(
        ("args", "kwargs", "_error"),
        [
            (("not a density", "e-"), {}, TypeError),
            ((5 * u.s, "e-"), {}, u.UnitTypeError),
            ((5 * u.m**-2, "e-"), {}, u.UnitTypeError),
            (
                (),
                {"n": 5 * u.m**-3, "particle": "not a particle"},
                InvalidParticleError,
            ),
        ],
    )
    def test_raises(self, args, kwargs, _error) -> None:
        """
        Test scenarios that cause plasma_frequency to raise an
        Exception.
        """
        with pytest.raises(_error):
            plasma_frequency(*args, **kwargs)

    @pytest.mark.parametrize(
        ("args", "kwargs", "_warning", "expected"),
        [
            (
                (1e19, "e-"),
                {},
                u.UnitsWarning,
                plasma_frequency(1e19 * u.m**-3, "e-"),
            ),
            ((1e19, "p"), {}, u.UnitsWarning, plasma_frequency(1e19 * u.m**-3, "p")),
        ],
    )
    def test_warns(self, args, kwargs, _warning, expected) -> None:
        """
        Test scenarios the cause plasma_frequency to issue a warning.
        """
        with pytest.warns(_warning):
            wp = plasma_frequency(*args, **kwargs)
        assert isinstance(wp, u.Quantity)
        assert wp.unit == u.rad / u.s

        if expected is not None:
            assert np.allclose(wp, expected)

    @pytest.mark.parametrize(
        ("args", "kwargs", "expected", "rtol"),
        [
            ((1 * u.cm**-3, "e-"), {}, 5.64e4, 1e-2),
            ((1 * u.cm**-3, "N+"), {}, 3.53e2, 1e-1),
            ((1e17 * u.cm**-3, "H-1"), {"Z": 0.8}, 333045427357.53955, 1e-6),
            (
                (5e19 * u.m**-3, "p"),
                {},
                plasma_frequency(5e19 * u.m**-3, particle="H-1+").value,
                1e-5,
            ),
            ((m_p.to(u.u).value * u.cm**-3,), {"particle": "p"}, 1.32e3, 1e-2),
        ],
    )
    def test_values(self, args, kwargs, expected, rtol) -> None:
        """Test various expected values."""
        wp = plasma_frequency(*args, **kwargs)

        assert isinstance(wp, u.Quantity)
        assert wp.unit == u.rad / u.s
        assert np.allclose(wp.value, expected, rtol=rtol)

    @pytest.mark.parametrize(
        ("args", "kwargs"),
        [((1 * u.cm**-3, "N+"), {}), ((1e12 * u.cm**-3,), {"particle": "p"})],
    )
    def test_to_hz(self, args, kwargs) -> None:
        """Test behavior of the ``to_hz`` keyword."""
        wp = plasma_frequency(*args, **kwargs)
        fp = plasma_frequency(*args, to_hz=True, **kwargs)

        assert isinstance(fp, u.Quantity)
        assert fp.unit == u.Hz
        assert fp.value == wp.value / (2.0 * np.pi)

    def test_nans(self) -> None:
        assert np.isnan(plasma_frequency(np.nan * u.m**-3, "e-"))

    def test_can_handle_numpy_arrays(self) -> None:
        assert_can_handle_nparray(plasma_frequency)


class TestPlasmaFrequencyLite:
    """Test class for `plasma_frequency_lite`."""

    @pytest.mark.parametrize(
        "inputs",
        [
            {"n": 1e12 * u.cm**-3, "particle": "e-"},
            {"n": 1e12 * u.cm**-3, "particle": "e-", "to_hz": True},
            {"n": 1e11 * u.cm**-3, "particle": "He", "Z": 0.8},
        ],
    )
    def test_normal_vs_lite_values(self, inputs) -> None:
        """
        Test that plasma_frequency and plasma_frequency_lite calculate
        the same values.
        """
        Z = inputs.get("Z", None)
        particle = _physical_particle_factory(inputs["particle"], Z=Z)

        inputs_unitless = {
            "n": inputs["n"].to(u.m**-3).value,
            "mass": particle.mass.value,
            "Z": np.abs(particle.charge_number),
        }

        if "to_hz" in inputs:
            inputs_unitless["to_hz"] = inputs["to_hz"]

        lite = plasma_frequency_lite(**inputs_unitless)

        normal = plasma_frequency(**inputs)
        assert np.allclose(normal.value, lite)
