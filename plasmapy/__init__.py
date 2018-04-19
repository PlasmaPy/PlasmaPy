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

# All imports that require Python 3.6+ should be placed after the Python
# version check.

from . import utils

__name__ = "plasmapy"

__doc__ = ("A community-developed and community-driven open source core "
           "Python package for plasma physics.")

_minimum_versions = {
    'numpy': '1.13',
    'astropy': '2.0',
    'scipy': '0.19',
    }

utils.check_versions(_minimum_versions)

# The file version.py is created by installing PlasmaPy with setup.py
# using functionality from astropy_helpers.  If this has not been run,
# then we will not create the __version__ attribute.

try:
    from .version import version as __version__
    from .version import githash as _githash
    from .version import cython_version as _cython_version
except ImportError:
    pass

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

# Allow astropy.units to be imported from PlasmaPy. This is the
# only place in the code where units should not be abbreviated as u.

try:
    from astropy import units
except ImportError:
    raise ImportError("Unable to import astropy.units as a PlasmaPy submodule")

# Clean up PlasmaPy's top-level namespace

del _minimum_versions, sys
