"""Adds custom functionality to pytest."""

import os

from hypothesis import Verbosity, settings

settings.register_profile("ci", max_examples=1000)
settings.register_profile("dev", max_examples=10)
settings.register_profile("debug", max_examples=10, verbosity=Verbosity.verbose)
settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "default"))

collect_ignore = ["noxfile.py"]

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib as mpl
except ImportError:
    pass
else:
    if "PLASMAPY_PLOT_TESTS" not in os.environ:
        mpl.use("Agg")
