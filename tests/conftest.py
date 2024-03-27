"""Adds custom functionality to pytest."""

import os

from hypothesis import Verbosity, settings


# Force MPL to use non-gui backends for testing.
try:
    import matplotlib as mpl
except ImportError:
    pass
else:
    import os

    if "PLASMAPY_PLOT_TESTS" not in os.environ:
        mpl.use("Agg")


def pytest_configure(config) -> None:  # coverage: ignore
    """Adds @pytest.mark.slow annotation for marking slow tests for optional skipping."""
    config.addinivalue_line(
        "markers",
        (
            "slow: mark test as slow to run. Used to mark tests that execute in more than 1000x the "
            "time of the median test. Tests marked slow may be skipped with 'pytest -m 'not slow'' "
            "or exclusively executed with 'pytest -m slow'."
        ),
    )


settings.register_profile("ci", max_examples=1000)
settings.register_profile("dev", max_examples=10)
settings.register_profile("debug", max_examples=10, verbosity=Verbosity.verbose)
settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "default"))
