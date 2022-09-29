"""Tests for the consistency of requirements."""

import pathlib
import pytest
import setuptools
import tomli

from typing import Dict, Set

import plasmapy

base_directory = pathlib.Path(__file__).parents[2]
requirements_directory = base_directory / "requirements"
requirements_prefixes = ("build", "docs", "extras", "install", "tests")


def read_requirements_txt_file(prefix: str, requirements_directory: str) -> Set[str]:
    """
    Read in a text file containing requirements.

    Parameters
    ----------
    prefix : str
        The prefix to a filename, such as ``"build"`` for ``build.txt``.

    requirements_directory : str
        The path to the directory containing the requirements file.

    Returns
    -------
    list of str
        A `list` containing the lines of the requirements file,
        excluding lines that do not start with an alphabetic character
        (like comments).
    """
    file = (requirements_directory / prefix).with_suffix(".txt")
    lines_of_file = file.read_text().splitlines()
    return {line.strip() for line in lines_of_file if line[0].isalpha()}


def get_requirements_from_txt() -> Dict[str, Set[str]]:
    """
    Get the requirements from the ``.txt`` files in the ``requirements``
    directory.
    """
    return {
        prefix: read_requirements_txt_file(prefix, requirements_directory)
        for prefix in requirements_prefixes
    }


def get_requirements_from_setup_cfg() -> Dict[str, Set[str]]:
    """Get the requirements that are contained in ``setup.cfg``."""
    configuration = setuptools.config.read_configuration(f"{base_directory}/setup.cfg")
    d = {
        key: set(configuration["options"]["extras_require"][key])
        for key in "docs extras tests".split()
    }
    d["install"] = set(configuration["options"]["install_requires"])
    return d


def get_requirements_from_pyproject_toml() -> Dict[str, Set[str]]:
    """Get the requirements that are contained in pyproject.toml."""
    with open(f"{base_directory}/pyproject.toml", "rb") as pyproject:
        pyproject_toml = tomli.load(pyproject)
    return {"build": set(pyproject_toml["build-system"]["requires"])}


requirements_from_txt = get_requirements_from_txt()
requirements_from_pyproject_toml = get_requirements_from_pyproject_toml()
requirements_from_setup_cfg = get_requirements_from_setup_cfg()

requirements_table = [
    (
        "build",
        requirements_from_txt["build"],
        "requirements/build.txt",
        requirements_from_pyproject_toml["build"],
        "pyproject.toml",
    ),
    (
        "install",
        requirements_from_txt["install"],
        "requirements/install.txt",
        requirements_from_setup_cfg["install"],
        "setup.cfg",
    ),
    (
        "extras",
        requirements_from_txt["extras"],
        "requirements/extras.txt",
        requirements_from_setup_cfg["extras"],
        "setup.cfg",
    ),
    (
        "docs",
        requirements_from_txt["docs"],
        "requirements/docs.txt",
        requirements_from_setup_cfg["docs"] - requirements_from_txt["extras"],
        "setup.cfg",
    ),
    (
        "tests",
        requirements_from_txt["tests"],
        "requirements/tests.txt",
        requirements_from_setup_cfg["tests"] - requirements_from_txt["extras"],
        "setup.cfg",
    ),
]


@pytest.mark.parametrize("prefix, req1, file1, req2, file2", requirements_table)
def test_consistency_of_requirements(
    prefix: str, req1: Set[str], file1: str, req2: Set[str], file2: str
):
    """
    Check that the ``prefix``-category of requirements in ``req1`` (from
    ``file1``) match the requirements in ``req2`` (from ``file2``).

    Note that in ``setup.cfg``, the ``docs`` and ``tests`` requirements
    also include ``extras``, so these should be subtracted out in the
    parametrization for this test.
    """
    error_messages = []

    if requirements_in_req1_but_not_req2 := req1 - req2:
        error_messages.append(
            f"The following {prefix} requirements are in {file1} but "
            f"not {file2}: {requirements_in_req1_but_not_req2}."
        )

    if requirements_in_rec2_but_not_rec1 := req2 - req1:
        error_messages.append(
            f"The following {prefix} requirements are in {file2} but "
            f"not {file1}: {requirements_in_rec2_but_not_rec1}."
        )

    if errmsg := "".join(error_messages).strip():
        pytest.fail(errmsg)
