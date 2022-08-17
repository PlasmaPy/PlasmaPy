"""Tests for functions in relativity.py."""

import numpy as np
import pytest

from astropy import units as u
from astropy.constants import c

from plasmapy.formulary.relativity import (
    Lorentz_factor,
    relativistic_energy,
    RelativisticBody,
)
from plasmapy.particles import CustomParticle, electron, proton
from plasmapy.utils.code_repr import call_string
from plasmapy.utils.exceptions import RelativityError


def test_Lorentz_factor():
    r"""Test Lorentz_factor in relativity.py"""

    V = 123456789 * u.m / u.s
    assert np.isclose(Lorentz_factor(V), (1 / np.sqrt(1 - V**2 / c**2)).value)
    assert Lorentz_factor(-V) == Lorentz_factor(V)

    assert np.isclose(Lorentz_factor(0 * u.m / u.s), 1.0)
    assert Lorentz_factor(c) == np.inf

    V_arr = np.array([987532.0, 299792458]) * u.m / u.s
    gamma_arr = Lorentz_factor(V_arr)
    assert np.isclose(gamma_arr[0], (1 / np.sqrt(1 - V_arr[0] ** 2 / c**2)).value)
    assert gamma_arr[1] == np.inf

    assert (
        Lorentz_factor(3 * u.m / u.s) * u.dimensionless_unscaled
    ).unit == u.dimensionless_unscaled

    with pytest.raises(RelativityError):
        Lorentz_factor(1.0000000001 * c)

    with pytest.raises(ValueError), pytest.warns(u.UnitsWarning):
        Lorentz_factor(299792459)

    with pytest.warns(u.UnitsWarning):
        Lorentz_factor(2.2)

    with pytest.raises(u.UnitTypeError):
        Lorentz_factor(4 * u.kg)


def test_relativistic_energy():
    r"""Test relativistic_energy in relativity.py"""

    v = 123456789 * u.m / u.s
    m = 1 * u.kg
    assert np.isclose(
        relativistic_energy(m, v).value,
        ((1 / np.sqrt(1 - v**2 / c**2)) * m * c**2).value,
    )
    assert relativistic_energy(m, -v) == relativistic_energy(m, v)

    assert np.isclose(relativistic_energy(m, 0 * u.m / u.s).value, (m * c**2).value)
    assert relativistic_energy(m, c) == np.inf

    V_arr = np.array([987532.0, 299792458]) * u.m / u.s
    Energy_arr = relativistic_energy(m, V_arr)
    assert np.isclose(
        Energy_arr[0].value,
        ((1 / np.sqrt(1 - V_arr[0] ** 2 / c**2)) * m * c**2).value,
    )
    assert Energy_arr[1] == np.inf

    assert relativistic_energy(2 * u.kg, 3 * u.m / u.s).unit == u.J

    with pytest.raises(RelativityError):
        relativistic_energy(m, 1.0000000001 * c)

    with pytest.raises(RelativityError), pytest.warns(u.UnitsWarning):
        relativistic_energy(1, 299792459)

    with pytest.warns(u.UnitsWarning):
        relativistic_energy(m, 2.2)

    with pytest.raises(u.UnitTypeError):
        relativistic_energy(m, 4 * u.kg)

    with pytest.raises(ValueError):
        relativistic_energy(-m, v)


@pytest.fixture
def ultra_relativistic_proton():
    """Representing an ultra high energy cosmic ray (UHECR)."""
    return RelativisticBody(particle="p+", kinetic_energy=np.longdouble(1e21) * u.eV)


def test_ultra_relativistic_proton_lorentz_factor(ultra_relativistic_proton):
    assert ultra_relativistic_proton.lorentz_factor >= 1e6


def test_uhech_v_over_c(ultra_relativistic_proton):
    assert u.isclose(ultra_relativistic_proton.v_over_c, 1, atol=1e-10)


def test_ultra_relativistic_proton_speed(ultra_relativistic_proton):
    assert u.isclose(ultra_relativistic_proton.velocity, c, rtol=1e-14)


proton_at_half_c_inputs = [
    ("v_over_c", 0.5),
    ("velocity", 0.5 * c),
    ("lorentz_factor", 1.1547005383792517),
    ("total_energy", 1.7358354725115025e-10 * u.J),
    ("kinetic_energy", 2.3255785652637692e-11 * u.J),
    ("momentum", 2.8950619440057805e-19 * u.kg * u.m / u.s),
]


@pytest.mark.parametrize("attr, expected", proton_at_half_c_inputs)
@pytest.mark.parametrize("parameter, argument", proton_at_half_c_inputs)
def test_relativistic_body(parameter, argument, attr, expected):
    """Test attributes of RelativisticBody."""

    if parameter == "velocity":
        parameter = "V"

    kwargs = {"particle": proton, parameter: argument}
    call_str = call_string(RelativisticBody, kwargs=kwargs)

    relativistic_body = RelativisticBody(**kwargs)
    actual = getattr(relativistic_body, attr)

    actual_unit = getattr(actual, "unit", None)
    expected_unit = getattr(expected, "unit", None)

    assert actual_unit == expected_unit, (
        f"{call_str}.{attr}.unit was expected to give "
        f"{expected_unit}, but instead gave {actual_unit}."
    )

    assert u.isclose(actual, expected, rtol=1e-6), (
        f"{call_str}.{attr} was expected to be {expected}, but instead "
        f"gave {actual}."
    )


@pytest.mark.parametrize("attr_to_set, set_value", proton_at_half_c_inputs)
@pytest.mark.parametrize("attr_to_test, expected", proton_at_half_c_inputs)
def test_relativistic_body_setters(attr_to_set, set_value, attr_to_test, expected):
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
def test_relativistic_body_mass_energy(particle):
    """Test `RelativisticBody.mass_energy`."""
    relativistic_body = RelativisticBody(particle, v_over_c=0)
    expected = particle.mass * c**2

    actual = relativistic_body.mass_energy

    assert u.isclose(actual, expected, rtol=1e-9)


@pytest.mark.parametrize(
    "kwargs, exception",
    [
        ({"V": 299792459 * (u.m / u.s)}, RelativityError),
        ({"v_over_c": 1.00001}, RelativityError),
        ({"total_energy": -1 * u.J}, ValueError),
        ({"kinetic_energy": -1 * u.J}, ValueError),
        ({"lorentz_factor": 0.99}, ValueError),
        ({"lorentz_factor": 0}, ValueError),
        ({"lorentz_factor": -1}, ValueError),
        ({"lorentz_factor": -1}, ValueError),
        ({"lorentz_factor": -1, "v_over_c": 0.5}, ValueError),
        ({"total_energy": 1 * u.J, "momentum": 1 * u.kg * u.m / u.s}, ValueError),
        ({}, ValueError),
    ],
)
def test_relativistic_body_init_exceptions(kwargs, exception):
    """
    Test that `RelativisticBody` raises the appropriate exceptions
    during instantiation.
    """
    with pytest.raises(exception):
        RelativisticBody(proton, **kwargs)


def test_relativistic_body_equality():
    """Test that a `RelativisticBody` instance equals itself."""
    relativistic_body = RelativisticBody(particle=proton, v_over_c=0.34)
    assert relativistic_body == relativistic_body


def test_relativistic_body_inequality_with_different_velocities():
    """
    Test that `RelativisticBody` instances are not equal when the
    velocity provided to them is different.
    """
    slower_body = RelativisticBody(particle=proton, v_over_c=0.23)
    faster_body = RelativisticBody(particle=proton, v_over_c=0.24)
    assert slower_body != faster_body


def test_relativistic_body_inequality_with_different_particles():
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
def test_relativistic_body_defined_using_mass(attr):
    """Test that a RelativisticBody can be provided a mass as the particle."""
    V = c / 2

    relativistic_proton = RelativisticBody(proton, V=V)
    expected = getattr(relativistic_proton, attr)

    relativistic_custom_particle = RelativisticBody(proton.mass, V=V)
    actual = getattr(relativistic_custom_particle, attr)

    assert u.isclose(actual, expected, rtol=1e-9)


@pytest.mark.xfail(reason="RelativisticBody does not accept nan velocities yet")
def test_relativistic_body_nan_velocity():
    """
    Test that RelativisticBody can be created with no velocity defined,
    and then have the velocity be nan
    """
    relativistic_body = RelativisticBody("p+")
    assert np.isnan(relativistic_body.velocity)
