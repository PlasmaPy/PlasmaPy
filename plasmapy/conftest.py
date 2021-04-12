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

from plasmapy.plasma.symbolicequilibrium import SymbolicEquilibrium


@pytest.fixture(scope="module")
def equilibrium():
    import plasmaboundaries

    params = plasmaboundaries.ITER.copy()
    equilibrium = SymbolicEquilibrium(**params, B0=5.2, config="single-null")
    return equilibrium


@pytest.fixture(scope="module")
def flux_surface(equilibrium, psi_value=-0.01):
    return equilibrium.get_flux_surface(psi_value)


import os

from hypothesis import settings, Verbosity

settings.register_profile("ci", max_examples=1000)
settings.register_profile("dev", max_examples=10)
settings.register_profile("debug", max_examples=10, verbosity=Verbosity.verbose)
settings.load_profile(os.getenv(u"HYPOTHESIS_PROFILE", "default"))
