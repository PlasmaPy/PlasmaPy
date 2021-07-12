"""Adds custom functionality to pytest"""
# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    import os

    if "PLASMAPY_PLOT_TESTS" not in os.environ:
        matplotlib.use("Agg")

import datetime


# coverage : ignore
def pytest_configure(config):
    """Adds @pytest.mark.slow annotation for marking slow tests for optional skipping"""
    config.addinivalue_line(
        "markers",
        (
            "slow: mark test as slow to run. Used to mark tests that execute in more than 1000x the "
            "time of the median test. Tests marked slow may be skipped with 'pytest -m 'not slow'' "
            "or exclusively executed with 'pytest -m slow'."
        ),
    )


import pytest
import astropy.units as u

from plasmapy.plasma.symbolicequilibrium import SymbolicEquilibrium
from plasmapy.plasma.simplefluxsurface import SimpleFluxSurface


@pytest.fixture(scope="module")
def equilibrium():
    import plasmaboundaries

    params = plasmaboundaries.ITER.copy()
    equilibrium = SymbolicEquilibrium(**params, B0=5.2, config="single-null")
    return equilibrium


@pytest.fixture(scope="module", params=[-0.01, -0.02, "simple"])
def flux_surface(request, equilibrium):
    psi_value = request.param
    if psi_value == "simple":
        return SimpleFluxSurface(
            equilibrium.aspect_ratio,
            safety_factor=1,
            major_radius=6.2 * u.m,
            minor_radius=2 * u.m,
            axial_elongation=equilibrium.elongation,
            axial_toroidal_field=5.3 * u.T,
            q0 = 3.0,
            axial_safety_factor = 3.0,
        )

    return equilibrium.get_flux_surface(psi_value)
