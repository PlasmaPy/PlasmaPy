"""
PlasmaPy is a community-developed and community-driven open source core
Python package for plasma physics.
"""

# Check that PlasmaPy is being imported from a recent enough version of
# Python so that otherwise an appropriate ImportError with a useful
# error message will be raised.  This check should be done before
# anything else to avoid error messages that obfuscate the source of the
# problem.  This file must be Python 2 compliant enough to avoid raising
# any SyntaxError exceptions, and thus cannot have f-strings.

import sys

if sys.version_info < (3, 6):
    raise ImportError(
        "PlasmaPy requires Python version 3.6 or higher, but is being "
        "called from Python version {}.".format(sys.version.split()[0]))

# Place all imports that require Python 3.6+ **after** the Python
# version check.

from . import import_helpers  # noqa

# TODO: Create _minimum_versions from requirements/requirements.txt
# All of this could be put into import_helpers._check_versions.

_minimum_versions = {
    'numpy': '1.13',
    'astropy': '2.0',
    'scipy': '0.19',
    'matplotlib': '2.0',
    'mpmath': '1.0',
    'lmfit': '0.9.7',
    'roman': '1.4',
    'colorama': '0.3',
    'cython': '0.26',
}

import_helpers._check_versions(_minimum_versions)

__name__ = "plasmapy"

__doc__ = ("A community-developed and community-driven open source "
           "core Python package for plasma physics.")

# The file version.py is created by installing PlasmaPy with setup.py
# using functionality from astropy_helpers.  If this has not been run,
# then we can set the default values to None.

try:
    from .version import version as __version__
    from .version import githash as _githash
except (ImportError, ModuleNotFoundError):
    __version__ = None
    _githash = None

try:
    from . import atomic
    from . import classes
    from . import constants
    from . import diagnostics
    from . import mathematics
    from . import physics
    from . import utils
except ImportError:
    raise ImportError("Unable to load PlasmaPy subpackages.")

# Allow astropy.units to be imported as plasmapy.units. This is the only
# place in the code where units should not be abbreviated as u.

try:
    from astropy import units  # do not abbreviate
except ImportError:
    raise ImportError("Unable to import astropy.units as a PlasmaPy submodule")

# Clean up the top-level namespace

del _minimum_versions, sys, import_helpers
