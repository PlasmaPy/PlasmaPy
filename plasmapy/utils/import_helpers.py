import importlib
import warnings
from distutils.version import LooseVersion


def check_versions(minimum_versions):
    """Raises an ImportError if a dependent package is not installed
    and at the required version number, or provides a warning if the
    version of the dependent package cannot be found."""

    for module_name in minimum_versions.keys():
        minimum_version = LooseVersion(minimum_versions[module_name])

        try:
            module = importlib.import_module(module_name)
            module_version = LooseVersion(module.__version__)
        except ImportError:
            raise ImportError(f"Unable to import {module_name} while "
                              "importing PlasmaPy.") from None
        except AttributeError:  # coveralls: ignore
            warnings.warn(f"{module_name} version {minimum_version.vstring} "
                          "is required for PlasmaPy.  However, the version of "
                          f"{module_name} could not be determined to check if "
                          "this requirement is met.")
        else:
            if minimum_version > module_version:
                raise ImportError(
                    f"{module_name} {minimum_version} or newer is required "
                    "for PlasmaPy. The currently installed version is "
                    f"{module_version}.") from None
