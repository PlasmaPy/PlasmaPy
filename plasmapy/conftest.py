"""Adds custom functionality to pytest."""

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib as mpl
except ImportError:
    pass
else:
    import os

    if "PLASMAPY_PLOT_TESTS" not in os.environ:
        mpl.use("Agg")
