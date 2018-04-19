import importlib
import warnings
import distutils.version as dv


def _split_version(version):
    """Separate a string including digits separated by periods into a
    tuple of integers."""
    return tuple(int(ver) for ver in version.split('.'))


def _check_versions(minimum_versions):
    """
    Raise an `ImportError` if a dependent package is not installed and
    at the required version number, or issue a
    `~plasmapy.utils.PlasmaPyWarning` if the version of the dependent
    package cannot be found.
    """
    for module_name in minimum_versions.keys():
        minimum_version = dv.LooseVersion(minimum_versions[module_name])
        try:
            module = importlib.import_module(module_name)
            module_version = dv.LooseVersion(module.__version__)
        except ImportError:
            raise ImportError(
                f"Unable to import {module_name} while importing PlasmaPy.") from None
        except ModuleNotFoundError:
            raise ImportError(
                f"Unable to find {module_name} while importing PlasmaPy.") from None
        except AttributeError:  # coveralls: ignore
            warnings.warn(
                f"{module_name} version {minimum_version.vstring} "
                "is required for PlasmaPy.  However, the version of "
                f"{module_name} could not be determined to check if "
                "this requirement is met.", UserWarning)
        else:
            if minimum_version > module_version:
                raise ImportError(
                    f"{module_name} {minimum_version} or newer is required "
                    "for PlasmaPy. The currently installed version is "
                    f"{module_version}.") from None
