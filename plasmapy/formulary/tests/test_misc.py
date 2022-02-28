"""Tests for functionality contained in `plasmapy.formulary.misc`."""

import pytest

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
def test_parameters_aliases(alias, parent):
    """Test all aliases defined in parameters.py"""
    assert alias is parent
