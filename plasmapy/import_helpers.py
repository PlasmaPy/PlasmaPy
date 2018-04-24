import importlib
import warnings
import distutils.version as dv
from os.path import dirname

requirements_dir = dirname(dirname(__file__)) + '/requirements/'


def split_version(version):
    """
    Separate a string including digits separated by periods into a
    tuple of integers.
    """
    return tuple(int(ver) for ver in version.split('.'))


def get_minimum_versions(file='requirements.txt'):
    """
    Get the minimum versions from `requirements.txt` in the
    `requirements` directory of `plasmapy`.
    """
    requirements = open(requirements_dir + file)
    minimum_versions = {}

    for line in requirements.readlines():
        if ("(>=") in line:
            package = line.split('(>=')[0].strip()
            minimum_version = line.split('(>=')[1].split(')')[0].strip()
        elif line.isalnum():
            package = line.strip()
            minimum_version = '0.1'
        else:
            raise ImportError(f"Invalid package requirement in requirements.txt: {line}")
        minimum_versions[package] = minimum_version

    requirements.close()

    return minimum_versions


def check_versions(minimum_versions=None):
    """
    Raise an `ImportError` if a dependent package is not installed and
    at the required version number, or issue a `UserWarning` if the
    version of the dependent package cannot be found.
    """
    if not minimum_versions:
        minimum_versions = get_minimum_versions()

    for module_name in minimum_versions.keys():
        minimum_version = dv.LooseVersion(minimum_versions[module_name])
        try:
            module = importlib.import_module(module_name)
            module_version = dv.LooseVersion(module.__version__)
        except (ImportError, ModuleNotFoundError):
            raise ImportError(
                f"Unable to import PlasmaPy because the required "
                f"package {module_name} cannot be imported.") from None
        except AttributeError:  # coveralls: ignore
            warnings.warn(
                f"{module_name} version {minimum_version.vstring} "
                "is required for PlasmaPy.  However, the version of "
                f"{module_name} could not be determined to check if "
                "this requirement is met.", UserWarning)
        else:  # coveralls: ignore
            if minimum_version > module_version:
                raise ImportError(
                    f"PlasmaPy requires {module_name} {minimum_version}"
                    f" or newer, but {module_name} {module_version} "
                    f"is currently installed.") from None
