"""Tests for the consistency of requirements."""

import os
import setuptools
import toml

from typing import Dict, List

import plasmapy

base_directory = os.path.realpath(f"{plasmapy.__path__[0]}/..")
requirements_directory = f"{base_directory}/requirements"
requirements_prefixes = ("build", "docs", "extras", "install", "tests")


def read_requirements_txt_file(prefix: str, requirements_directory: str) -> List[str]:
    """
    Read in a text file containing requirements.

    Parameters
    ----------
    prefix : str
        The prefix to a filename, such as `"build"` for :file:`build.txt`.

    requirements_directory : str
        The path to the directory containing the requirements file.

    Returns
    -------
    list of str
        A `list` containing the lines of the requirements file,
        excluding lines that do not start with an alphabetic character
        (like comments).
    """
    filename = f"{requirements_directory}/{prefix}.txt"
    with open(filename) as file:
        lines_of_file = file.readlines()
    return [line.strip() for line in lines_of_file if line[0].isalpha()]


def get_requirements_from_txt() -> Dict[str, str]:
    """
    Get the requirements from the .txt files in the requirements
    directory.
    """
    return {
        prefix: read_requirements_txt_file(prefix, requirements_directory)
        for prefix in requirements_prefixes
    }


def get_requirements_from_setup_cfg() -> Dict[str, str]:
    """Get the requirements that are contained in setup.cfg."""
    configuration = setuptools.config.read_configuration(f"{base_directory}/setup.cfg")
    return {
        "docs": configuration["options"]["extras_require"]["docs"],
        "extras": configuration["options"]["extras_require"]["extras"],
        "install": configuration["options"]["install_requires"],
        "tests": configuration["options"]["extras_require"]["tests"],
    }


def get_requirements_from_pyproject_toml() -> Dict[str, str]:
    """Get the requirements that are contained in pyproject.toml."""
    pyproject_toml = toml.load(f"{base_directory}/pyproject.toml")
    return {"build": pyproject_toml["build-system"]["requires"]}
