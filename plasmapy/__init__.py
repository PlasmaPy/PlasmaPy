from ._metadata import (
    name as __name__,
    version as __version__,
    description as __doc__,
    author as __author__,
)

import sys
import warnings

__minimum_python_version__ = '3.6'
__minimum_numpy_version__ = '1.13.0'
__minimum_astropy_version__ = '2.0.0'


def _split_version(version):
    return tuple(int(ver) for ver in version.split('.'))


def _min_required_version(required, current):  # coveralls: ignore
    r""" Return `True` if the current version meets the required minimum
        version and `False` if not/ if not installed.

        Right now `required` and `current` are just '.' separated strings
        but it would be good to make this more general and accept modules.
    """
    return _split_version(current) >= _split_version(required)


def _check_numpy_version():  # coveralls: ignore
    r""" Make sure numpy in installed and meets the minimum version requirements
    """
    required_version = False
    np_ver = None

    try:
        from numpy import __version__ as np_ver
        required_version = _min_required_version(__minimum_numpy_version__,
                                                 np_ver)
    except ImportError:
        pass

    if not required_version:
        ver_error = ("Numpy {} or above is required for PlasmaPy. The "
                     "currently installed version is {}"
                     ).format(__minimum_numpy_version__, np_ver)
        raise ImportError(ver_error)


def _check_astropy_version():  # coveralls: ignore
    r""" Make sure astropy in installed and meets the minimum version requirements
    """
    required_version = False
    ap_ver = None

    try:
        from astropy import __version__ as ap_ver
        required_version = _min_required_version(__minimum_astropy_version__,
                                                 ap_ver)
    except ImportError:
        pass

    if not required_version:
        ver_error = ("Astropy {} or above is required for PlasmaPy. The "
                     "currently installed version is {}"
                     ).format(__minimum_astropy_version__, ap_ver)
        raise ImportError(ver_error)


if (sys.version_info < _split_version(__minimum_python_version__)):  # coveralls: ignore
    warnings.warn("PlasmaPy does not support Python 3.5 and below")

_check_numpy_version()
_check_astropy_version()

try:
    from .classes import Plasma
    from . import classes
    from . import constants
    from . import atomic
    from . import math
    from . import physics
    from . import utils
except ImportError:  # coveralls: ignore
    raise ImportError("Unable to load PlasmaPy subpackages.")
