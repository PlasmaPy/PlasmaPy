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


def pytest_configure(config):  # coverage: ignore
    """Adds @pytest.mark.slow annotation for marking slow tests for optional skipping."""
    config.addinivalue_line(
        "markers",
        (
            "slow: mark test as slow to run. Used to mark tests that execute in more than 1000x the "
            "time of the median test. Tests marked slow may be skipped with 'pytest -m 'not slow'' "
            "or exclusively executed with 'pytest -m slow'."
        ),
    )
