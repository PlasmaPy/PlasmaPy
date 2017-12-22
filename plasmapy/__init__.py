import sys
import warnings
import importlib

__minimum_python_version__ = '3.6'

__minimum_versions__ = {
    'numpy': '1.13',
    'astropy': '2.0',
    }


def _split_version(version):
    """Separate a string including digits separated by periods into a
    tuple of integers."""
    return tuple(int(ver) for ver in version.split('.'))


def _check_python(minimum_python_version):
    """Raises an ImportError if the version of Python is not at least
    the version given by a string representing the minimum version
    number."""
    required_python = _split_version(minimum_python_version)
    current_python = sys.version_info

    if current_python < required_python:
        raise ImportError(f"PlasmaPy requires Python {minimum_python_version} "
                          "or newer.") from None


def _check_versions(minimum_versions):
    """Raises an ImportError if a dependent package is not installed
    and at the required version number, or provides a warning if the
    version of the dependent package cannot be found."""

    for module_name in minimum_versions.keys():
        minimum_version = minimum_versions[module_name]

        try:
            module = importlib.import_module(module_name)
            module_version = module.__version__
        except ImportError:
            raise ImportError(f"Unable to import {module_name} while "
                              "importing PlasmaPy.") from None
        except AttributeError:
            warnings.warn(f"{module_name}.__version__ was not found while "
                          "importing PlasmaPy", UserWarning)
        else:
            if minimum_version > module_version:
                raise ImportError(
                    f"{module_name} {minimum_version} or newer is required "
                    "for PlasmaPy. The currently installed version is "
                    f"{module_version}.") from None


_check_python(__minimum_python_version__)
_check_versions(__minimum_versions__)

try:
    from .classes import Plasma
    from . import classes
    from . import constants
    from . import atomic
    from . import mathematics
    from . import physics
    from . import utils
except ImportError as exc:
    raise ImportError("Unable to load PlasmaPy subpackages.") from exc
