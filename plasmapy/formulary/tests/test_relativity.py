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
from plasmapy.particles import proton
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
    return RelativisticBody(particle="p+", kinetic_energy=np.float128(1e21) * u.eV)


def test_ultra_relativistic_proton_lorentz_factor(ultra_relativistic_proton):
    assert ultra_relativistic_proton.lorentz_factor >= 1e6


def test_uhech_v_over_c(ultra_relativistic_proton):
    assert u.isclose(ultra_relativistic_proton.v_over_c, 1, atol=1e-10)


def test_ultra_relativistic_proton_speed(ultra_relativistic_proton):
    assert u.isclose(ultra_relativistic_proton.velocity, c, rtol=1e-14)


proton_at_half_c_inputs = [
    ("v_over_c", 0.5),
    ("V", 0.5 * c),
    ("lorentz_factor", 1.1547005383792517),
    ("total_energy", 1.7358354725115025e-10 * u.J),
    ("kinetic_energy", 2.3255785652637692e-11 * u.J),
    ("momentum", 2.8950619440057805e-19 * u.kg * u.m / u.s),
]

# Need to add a test for mass_energy, which can't be used as an input


@pytest.mark.parametrize("attr, expected", proton_at_half_c_inputs)
@pytest.mark.parametrize("parameter, argument", proton_at_half_c_inputs)
def test_relativistic_body(parameter, argument, attr, expected):
    """Test attributes of RelativisticBody."""
    kwargs = {"particle": proton, parameter: argument}
    call_str = call_string(RelativisticBody, kwargs=kwargs)

    if attr == "V":
        attr = "velocity"

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
