# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._base_init import *
# ----------------------------------------------------------------------------

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

__name__ = "plasmapy"

__doc__ = ("A community-developed and community-driven open source "
           "core Python package for plasma physics.")

class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple((int(val) for val in "3.6".split('.'))):
    raise UnsupportedPythonError("plasmapy does not support Python < {}".format(3.6))


if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.
    from . import atomic
    from . import classes
    from . import constants
    from . import diagnostics
    from . import mathematics
    from . import physics
    from . import utils

