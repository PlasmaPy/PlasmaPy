"""Tests for functionality contained in `plasmapy.formulary.misc`."""

import astropy.units as u
import numpy as np
import pytest

from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.misc import (
    Bohm_diffusion,
    DB_,
    magnetic_energy_density,
    magnetic_pressure,
    mass_density,
    pmag_,
    pth_,
    rho_,
    thermal_pressure,
    ub_,
)
from plasmapy.particles import Particle
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

B = 1.0 * u.T
B_arr = np.array([0.001, 0.002]) * u.T
B_nanarr = np.array([0.001, np.nan]) * u.T

n_i = 5e19 * u.m**-3

T_e = 1e6 * u.K


@pytest.mark.parametrize(
    "alias, parent",
    [
        (DB_, Bohm_diffusion),
        (ub_, magnetic_energy_density),
        (pmag_, magnetic_pressure),
        (rho_, mass_density),
        (pth_, thermal_pressure),
    ],
)
def test_aliases(alias, parent):
    """Test all aliases defined in misc.py"""
    assert alias is parent


class Test_mass_density:
    r"""Test the mass_density function in misc.py."""

    @pytest.mark.parametrize(
        "args, kwargs, conditional",
        [
            ((-1 * u.kg * u.m**-3, "He"), {}, pytest.raises(ValueError)),
            ((-1 * u.m**-3, "He"), {}, pytest.raises(ValueError)),
            (("not a Quantity", "He"), {}, pytest.raises(TypeError)),
            ((1 * u.m**-3,), {}, pytest.raises(TypeError)),
            ((1 * u.J, "He"), {}, pytest.raises(u.UnitTypeError)),
            ((1 * u.m**-3, None), {}, pytest.raises(TypeError)),
            (
                (1 * u.m**-3, "He"),
                {"z_ratio": "not a ratio"},
                pytest.raises(TypeError),
            ),
        ],
    )
    def test_raises(self, args, kwargs, conditional):
        with conditional:
            mass_density(*args, **kwargs)

    @pytest.mark.parametrize(
        "args, kwargs, expected",
        [
            ((1.0 * u.g * u.m**-3, ""), {}, 1.0e-3 * u.kg * u.m**-3),
            ((5.0e12 * u.cm**-3, "He"), {}, 3.32323849e-8 * u.kg * u.m**-3),
            (
                (5.0e12 * u.cm**-3, Particle("He")),
                {},
                3.32323849e-8 * u.kg * u.m**-3,
            ),
            (
                (5.0e12 * u.cm**-3, "He"),
                {"z_ratio": 0.5},
                1.66161925e-08 * u.kg * u.m**-3,
            ),
            (
                (5.0e12 * u.cm**-3, "He"),
                {"z_ratio": -0.5},
                1.66161925e-08 * u.kg * u.m**-3,
            ),
        ],
    )
    def test_values(self, args, kwargs, expected):
        assert np.isclose(mass_density(*args, **kwargs), expected)

    def test_handle_nparrays(self):
        """Test for ability to handle numpy array quantities"""
        assert_can_handle_nparray(mass_density)


def test_thermal_pressure():
    assert thermal_pressure(T_e, n_i).unit.is_equivalent(u.Pa)

    # TODO: may be array issues with arg "mass"
    assert_can_handle_nparray(thermal_pressure)


def test_magnetic_pressure():
    r"""Test the magnetic_pressure function in misc.py."""

    assert magnetic_pressure(B_arr).unit.is_equivalent(u.Pa)

    assert magnetic_pressure(B).unit.is_equivalent(u.Pa)

    assert magnetic_pressure(B).unit.name == "Pa"

    assert magnetic_pressure(B).value == magnetic_energy_density(B).value

    assert magnetic_pressure(B) == magnetic_energy_density(B.to(u.G))

    assert np.isclose(magnetic_pressure(B).value, 397887.35772973835)

    with pytest.warns(u.UnitsWarning):
        magnetic_pressure(5)

    with pytest.raises(u.UnitTypeError):
        magnetic_pressure(5 * u.m)

    assert np.isnan(magnetic_pressure(np.nan * u.T))

    with pytest.raises(ValueError):
        magnetic_pressure(5j * u.T)

    assert np.isnan(magnetic_pressure(B_nanarr)[-1])

    with pytest.warns(u.UnitsWarning):
        assert magnetic_pressure(22.2) == magnetic_pressure(22.2 * u.T)

    assert_can_handle_nparray(magnetic_pressure)


def test_magnetic_energy_density():
    r"""Test the magnetic_energy_density function in misc.py."""

    assert magnetic_energy_density(B_arr).unit.is_equivalent(u.J / u.m**3)

    assert magnetic_energy_density(B).unit.is_equivalent("J / m3")

    assert magnetic_energy_density(B).value == magnetic_pressure(B).value

    assert_quantity_allclose(
        magnetic_energy_density(2 * B), 4 * magnetic_energy_density(B)
    )

    assert_quantity_allclose(magnetic_energy_density(B).value, 397887.35772973835)

    assert_quantity_allclose(
        magnetic_energy_density(B), magnetic_energy_density(B.to(u.G))
    )

    assert isinstance(magnetic_energy_density(B_arr), u.Quantity)

    with pytest.warns(u.UnitsWarning):
        magnetic_energy_density(5)

    with pytest.raises(u.UnitTypeError):
        magnetic_energy_density(5 * u.m)

    assert np.isnan(magnetic_energy_density(np.nan * u.T))

    with pytest.raises(ValueError):
        magnetic_energy_density(5j * u.T)

    assert np.isnan(magnetic_energy_density(B_nanarr)[-1])

    with pytest.warns(u.UnitsWarning):
        assert magnetic_energy_density(22.2) == magnetic_energy_density(22.2 * u.T)

    assert_can_handle_nparray(magnetic_energy_density)


def test_Bohm_diffusion():
    r"""Test Mag_Reynolds in dimensionless.py"""

    T_e = 5000 * u.K
    B = 10 * u.T

    assert (Bohm_diffusion(T_e, B)).unit == u.m**2 / u.s

    with pytest.warns(u.UnitsWarning):
        Bohm_diffusion(5000, B)

    with pytest.raises(u.UnitTypeError):
        Bohm_diffusion(2.2 * u.kg, B)
