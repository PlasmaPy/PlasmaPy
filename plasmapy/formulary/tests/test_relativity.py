"""Tests for functions in relativity.py."""

import astropy.units as u
import numpy as np
import pytest
from astropy.constants import c
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.relativity import (
    Lorentz_factor,
    RelativisticBody,
    relativistic_energy,
)
from plasmapy.particles import CustomParticle, electron, proton
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.exceptions import RelativityError


@pytest.mark.parametrize(
    ("speed", "expected"),
    [
        (0 * u.m / u.s, 1),
        (np.nan * u.m / u.s, np.nan),
        (c, np.inf),
        (-c, np.inf),
        (123456789 * u.m / u.s, 1.0973686),
        (-123456789 * u.m / u.s, 1.0973686),
        (555555555 * u.km / u.hr, 1.1664056),
        (0.99 * c, 7.088812),
        ([np.nan, 0, c.si.value] * u.m / u.s, [-np.nan, 1, np.inf]),
    ],
)
def test_Lorentz_factor(speed, expected) -> None:
    actual = Lorentz_factor(V=speed)
    assert u.allclose(actual, expected, equal_nan=True, rtol=1e-7)


@pytest.mark.parametrize(
    ("speed", "exception"),
    [
        (np.inf * u.m / u.s, RelativityError),
        (-np.inf * u.m / u.s, RelativityError),
        (1.00000001 * c, RelativityError),
        (-1.00000001 * c, RelativityError),
        (299792458 * u.kg / u.s, u.UnitTypeError),
    ],
)
def test_Lorentz_factor_exceptions(speed, exception) -> None:
    with pytest.raises(exception):
        Lorentz_factor(speed)


@pytest.mark.parametrize(
    ("speed", "warning"), [(2.2, u.UnitsWarning), (np.nan, u.UnitsWarning)]
)
def test_Lorentz_factor_warnings(speed, warning) -> None:
    with pytest.warns(warning):
        Lorentz_factor(speed)


@pytest.mark.parametrize(
    ("velocity", "mass", "expected"),
    [
        (123456789 * u.m / u.s, 1 * u.kg, 9.86265694e16 * u.J),
        (-123456789 * u.m / u.s, 1 * u.kg, 9.86265694e16 * u.J),
        (5e6 * u.m / u.s, 0 * u.kg, 0 * u.J),
        (0 * u.m / u.s, 1 * u.kg, 1 * u.kg * c**2),
        (np.nan * u.m / u.s, 1 * u.kg, np.nan * u.J),
        (100 * u.m / u.s, np.nan * u.kg, np.nan * u.J),
        (5e6 * u.m / u.s, np.inf * u.kg, np.inf * u.J),
        ([123456789, np.nan] * u.m / u.s, 1 * u.kg, [9.86265694e16, np.nan] * u.J),
        (123456789 * u.m / u.s, [1, 2] * u.kg, [9.86265694e16, 1.972531388e17] * u.J),
        ([1, 2] * u.Mm / u.s, [1e-15, 1e-16] * u.kg, [89.87601788, 8.98775179] * u.J),
    ],
)
def test_relativistic_energy(velocity, mass, expected) -> None:
    actual = relativistic_energy(particle=mass, V=velocity)
    assert u.allclose(actual, expected, rtol=1e-6, atol=1e-6 * u.J, equal_nan=True)
    assert expected.unit == u.J


@pytest.mark.parametrize(
    ("velocity", "mass", "exception"),
    [
        (1.00000001 * c, 1 * u.kg, RelativityError),
        (-1.00000001 * c, 1 * u.kg, RelativityError),
        (0 * c, -1 * u.kg, InvalidParticleError),
    ],
)
def test_relativistic_energy_exceptions(velocity, mass, exception) -> None:
    with pytest.raises(exception):
        relativistic_energy(V=velocity, particle=mass)


@pytest.mark.parametrize(
    ("velocity", "mass", "warning"),
    [
        (2.2, 5 * u.kg, u.UnitsWarning),
    ],
)
def test_relativistic_energy_warnings(velocity, mass, warning) -> None:
    with pytest.warns(warning):
        relativistic_energy(V=velocity, particle=mass)


proton_at_half_c_inputs = [
    ("v_over_c", 0.5),
    ("velocity", 0.5 * c),
    ("lorentz_factor", 1.1547005383792517),
    ("total_energy", 1.7358354725115025e-10 * u.J),
    ("kinetic_energy", 2.3255785652637692e-11 * u.J),
    ("momentum", 2.8950619440057805e-19 * u.kg * u.m / u.s),
]


@pytest.mark.parametrize(("attribute", "expected"), proton_at_half_c_inputs)
@pytest.mark.parametrize(("parameter", "argument"), proton_at_half_c_inputs)
def test_relativistic_body(parameter, argument, attribute, expected) -> None:
    """
    Test that when we create `RelativisticBody` instances for each of
    the different velocity-like arguments, that each of the resulting
    `RelativisticBody` attributes end up providing the correct value.
    """

    if parameter == "velocity":
        parameter = "V"

    kwargs = {"particle": proton, parameter: argument}

    relativistic_body = RelativisticBody(**kwargs)
    actual = getattr(relativistic_body, attribute)

    assert_quantity_allclose(actual, expected, rtol=1e-9)


@pytest.mark.parametrize(("attr_to_set", "set_value"), proton_at_half_c_inputs)
@pytest.mark.parametrize(("attr_to_test", "expected"), proton_at_half_c_inputs)
def test_relativistic_body_setters(
    attr_to_set, set_value, attr_to_test, expected
) -> None:
    """Test setting RelativisticBody attributes."""
    relativistic_body = RelativisticBody(proton, v_over_c=0.1)
    setattr(relativistic_body, attr_to_set, set_value)

    actual = getattr(relativistic_body, attr_to_test)

    assert u.isclose(expected, actual, rtol=1e-8), (
        f"When setting {attr_to_set} to {set_value!r} in a "
        f"RelativisticBody instance, the value of {attr_to_test} was "
        f"expected to be {expected!r}, but was instead {actual!r}."
    )


@pytest.mark.parametrize("particle", [electron, proton])
def test_relativistic_body_mass_energy(particle) -> None:
    """Test `RelativisticBody.mass_energy`."""
    relativistic_body = RelativisticBody(particle, v_over_c=0)
    expected = particle.mass * c**2

    actual = relativistic_body.mass_energy

    assert u.isclose(actual, expected, rtol=1e-9)


@pytest.mark.parametrize(
    ("kwargs", "exception"),
    [
        ({"V": 299792459 * (u.m / u.s)}, RelativityError),
        ({"v_over_c": 1.00001}, RelativityError),
        ({"total_energy": -1 * u.J}, ValueError),
        ({"kinetic_energy": -1 * u.J}, ValueError),
        ({"lorentz_factor": 0.99}, ValueError),
        ({"lorentz_factor": 0}, ValueError),
        ({"lorentz_factor": -1}, ValueError),
        ({"lorentz_factor": -1, "v_over_c": 0.5}, ValueError),
        ({"total_energy": 1 * u.J, "momentum": 1 * u.kg * u.m / u.s}, ValueError),
        ({}, ValueError),
        ({"lorentz_factor": "wrong type"}, TypeError),
        ({"lorentz_factor": 3 * u.m / u.s}, u.UnitConversionError),
    ],
)
def test_relativistic_body_init_exceptions(kwargs, exception) -> None:
    """
    Test that `RelativisticBody` raises the appropriate exceptions
    during instantiation.
    """
    with pytest.raises(exception):
        RelativisticBody(proton, **kwargs)


def test_relativistic_body_equality() -> None:
    """Test that a `RelativisticBody` instance equals itself."""
    relativistic_body = RelativisticBody(particle=proton, v_over_c=0.34)
    assert relativistic_body == relativistic_body  # noqa: PLR0124


@pytest.mark.parametrize(
    ("this", "that"),
    [
        (RelativisticBody("p+", v_over_c=0.2), RelativisticBody("p+", v_over_c=0.5)),
        (RelativisticBody("p+", v_over_c=0.2), RelativisticBody("e-", v_over_c=0.2)),
        (RelativisticBody("p+", V=c / 2), "different type"),
    ],
)
def test_relativistic_body_inequalities(this, that) -> None:
    """Test the inequality properties of `RelativisticBody`."""
    assert this != that


def test_relativistic_body_inequality_with_different_velocities() -> None:
    """
    Test that `RelativisticBody` instances are not equal when the
    velocity provided to them is different.
    """
    slower_body = RelativisticBody(particle=proton, v_over_c=0.23)
    faster_body = RelativisticBody(particle=proton, v_over_c=0.24)
    assert slower_body != faster_body


def test_relativistic_body_inequality_with_different_particles() -> None:
    """
    Test that `RelativisticBody` instances are not equal when the
    particles are different.
    """
    relativistic_proton = RelativisticBody(particle=proton, v_over_c=0.23)
    relativistic_electron = RelativisticBody(particle=electron, v_over_c=0.23)
    assert relativistic_proton != relativistic_electron


@pytest.mark.parametrize(
    "attr",
    [
        "lorentz_factor",
        "v_over_c",
        "velocity",
        "momentum",
        "total_energy",
        "kinetic_energy",
    ],
)
def test_relativistic_body_defined_using_mass(attr) -> None:
    """Test that a RelativisticBody can be provided a mass as the particle."""
    V = c / 2

    relativistic_proton = RelativisticBody(proton, V=V)
    expected = getattr(relativistic_proton, attr)

    relativistic_custom_particle = RelativisticBody(proton.mass, V=V)
    actual = getattr(relativistic_custom_particle, attr)

    assert u.isclose(actual, expected, rtol=1e-9)


@pytest.mark.xfail(reason="RelativisticBody does not yet accept nan velocities")
def test_relativistic_body_nan_velocity() -> None:
    """
    Test that RelativisticBody can be created with no velocity defined,
    and then have the velocity be nan.
    """
    relativistic_body = RelativisticBody("p+")
    assert np.isnan(relativistic_body.velocity)


def test_relativistic_body_for_custom_particle() -> None:
    """Test that `RelativisticBody` can be created using a `CustomParticle`."""
    mass = 1e-27 * u.kg
    custom_particle = CustomParticle(mass=mass)
    velocity = 0 * u.m / u.s
    relativistic_custom_particle = RelativisticBody(
        particle=custom_particle, V=velocity
    )
    assert u.isclose(relativistic_custom_particle.lorentz_factor, 1, rtol=1e-9)


def test_relativistic_body_with_particle_list() -> None:
    """
    Test that `RelativisticBody` can be instantiated with multiple
    particles.
    """
    particles = ["p+", "e-"]
    relativistic_particles = RelativisticBody(particle=particles, V=0 * u.m / u.s)
    np.testing.assert_allclose(relativistic_particles.lorentz_factor, 1)


def test_relativistic_body_with_multiple_velocities() -> None:
    """
    Test that `RelativisticBody` can be instantiated with an array of
    velocities.
    """
    velocities = np.array([0, 0.5]) * c
    relativistic_particles = RelativisticBody("p+", V=velocities)
    np.testing.assert_allclose(
        relativistic_particles.lorentz_factor, [1, 1.1547005383792517]
    )


def test_relativistic_body_with_multiple_particles_and_velocities() -> None:
    """
    Test that `RelativisticBody` can be instantiated with multiple
    particles and multiple velocities.
    """
    velocities = np.array([0.5, 0.7]) * c
    particles = ["p+", "e-"]
    relativistic_particles = RelativisticBody(particle=particles, V=velocities)
    np.testing.assert_allclose(relativistic_particles.velocity, velocities)


@pytest.mark.parametrize("function", [repr, str])
def test_relativistic_body_into_string(function) -> None:
    """Test that `repr` and `str` work on RelativisticBody."""
    relativistic_body = RelativisticBody("p+", V=5.0 * u.m / u.s)
    expected = "RelativisticBody(p+, 5.0 m / s)"
    actual = function(relativistic_body)
    assert actual == expected
