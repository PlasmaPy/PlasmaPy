"""Tests of documentation infrastructure."""

import glob
import os

from pathlib import Path
from typing import Optional


def _change_directories_to_namespaces(path: str) -> list[str]:
    directories = [str(d[0]).removeprefix("../") for d in os.walk(path)]

    return [
        d.replace("/", ".").strip()
        for d in directories
        if not d.endswith(("tests")) and "." not in d and "/_" not in d
    ]


def _change_python_files_to_namespaces(path: str) -> list[str]:
    modules = list(glob.glob(f"{path}/**/[a-zA-Z]*.py"))

    return [
        m.removeprefix("../").removesuffix(".py").replace("/", ".")
        for m in modules
        if "tests/test_" not in m and "/_" not in m
    ]


def _get_modules(path: str) -> list[str]:

    subpackages = [
        *_change_directories_to_namespaces("../../plasmapy/"),
        *_change_directories_to_namespaces("plasmapy_sphinx/"),
    ]

    modules = [
        *_change_python_files_to_namespaces("../../plasmapy/"),
        *_change_python_files_to_namespaces("plasmapy_sphinx/")
    ]

    return sorted(subpackages + modules)


def _get_rst_files(
        docpath: str = ".",
        exclude: Optional[set[str]] = None,
) -> list[Path]:

    globs = glob.glob(f"{docpath}/**/[a-z][A-Z]*.rst")
    return



def _get_automodapi_directives(docpath: str = ".") -> list[str]:



def test_api_static():
    ...
