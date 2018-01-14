# This file must be Python 2 compliant enough so that the appropriate
# ImportError is raised when PlasmaPy is being imported from Python 2.
# This precludes the usage of some features such as f-strings.  Helper
# functions should be put into utils/import_helpers.py except when
# Python 2 compliance is needed to raise the correct ImportError.

import sys


def _split_version(version):
    """Separate a string including digits separated by periods into a
    tuple of integers."""
    return tuple(int(ver) for ver in version.split('.'))


__name__ = "plasmapy"

__doc__ = ("A community-developed and community-driven open source core "
           "Python package for plasma physics.")

__minimum_python_version__ = '3.6'

__minimum_versions__ = {
    'numpy': '1.13',
    'astropy': '2.0',
    'scipy': '0.19',
    }

if sys.version_info < _split_version(__minimum_python_version__):
    raise ImportError(
        "PlasmaPy requires Python version {} or higher, but is being called "
        "from Python version {}."
        .format(__minimum_python_version__, sys.version.split()[0]))

from . import utils

utils.check_versions(__minimum_versions__)

# The file version.py is created by installing PlasmaPy with setup.py
# using functionality from astropy_helpers.  If this has not been run,
# then we will not create the __version__ attribute.

try:
    from .version import version as __version__
except ImportError:
    pass

try:
    from .classes import Plasma
    from . import classes
    from . import constants
    from . import atomic
    from . import mathematics
    from . import physics
    from . import diagnostics
    from . import utils
except ImportError:
   raise ImportError("Unable to load PlasmaPy subpackages.")

try:
    from astropy import units
except ImportError:
    raise ImportError("Unable to import astropy.units as a PlasmaPy submodule")

# A more extensive and thoughtful method for cleaning up our top-level
# namespace is in Astropy's __init__.py (see also pull request #210).

del (__minimum_python_version__, __minimum_versions__)
