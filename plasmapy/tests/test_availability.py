"""
Test for common discrepancies and inconsistencies in code related to
availability and testing.
"""

import pytest
from collections import namedtuple
import re
from ..utils import RunTestError
from ..import_helpers import requirements_dir
from os.path import dirname

filenames = [
    'automated-code-tests.txt',
    'automated-documentation-tests.txt',
    'environment.txt',
    'environment.yml',
    'requirements.txt',
]

Requirements = namedtuple('Requirements', [
    'file',
    'readlines',
    'lines',
    'packages',
])

files = {}

for filename in filenames:
    file = open(requirements_dir + filename, 'rt')
    readlines = file.readlines()
    lines = {line.strip('\n').rstrip(' ') for line in readlines}

    if filename[-4:] == '.txt':

        packages = {
            line.split('(')[0].strip() if '(' in line else line.strip() for line in lines
        }

    elif filename[-4:] == '.yml':

        packages = set()
        for line in lines:
            if ':' not in line and 'python' not in line and 'plasmapy' not in line and line != "":
                packages |= {line.replace(' - ', '').strip()}
    else:
        packages = None

    files[filename] = Requirements(file, readlines, lines, packages)


@pytest.mark.parametrize('filename', filenames)
def test_for_no_duplicate_lines(filename):
    assert len(files[filename].lines) == len(set(files[filename].lines)), \
        f"Duplicate lines in {filename}."


def test_requirements_issubset_automated_code_tests():
    """
    Test that all of the lines in requirements.txt are also in
    automated-code-tests.txt.  This function strictly tests both
    packages and the minimum versions of packages.
    """
    missing_lines = files['requirements.txt'].lines - files['automated-code-tests.txt'].lines
    if missing_lines:
        raise Exception(
            f"The following lines from requirements.txt are not in "
            f"automated-code-tests.txt: {missing_lines}")


def test_requirements_issubset_automated_documentation_tests():
    """
    Test that all of the packages in requirements.txt are also in
    automated-code-tests.txt.
    """
    missing_packages = \
        files['requirements.txt'].packages - files['automated-documentation-tests.txt'].packages
    if missing_packages:
        raise Exception(
            f"The following packages from requirements.txt are not in "
            f"automated-documentation-tests.txt: {missing_packages}")


def test_requirements_issubset_environment_txt():
    """
    Test that all of the packages in requirements.txt are also in
    environment.txt.
    """
    missing_packages = files['requirements.txt'].packages - files['environment.txt'].packages
    if missing_packages:
        raise Exception(
            f"The following packages from requirements.txt are missing "
            f"from environment.txt: {missing_packages}")


def test_requirements_issubset_environment_yml():
    """
    Test that all of the packages in requirements.yml are also in
    environment.yml.
    """
    missing_packages = files['requirements.txt'].packages - files['environment.yml'].packages
    if missing_packages:
        raise Exception(
            f"The following packages from requirements.txt are missing "
            f" from environment.yml: {missing_packages}")


def test__environment_txt_and_yml_files():
    """
    Test that environment.txt and environment.yml have the same
    packages.
    """
    symmetric_difference = files['environment.txt'].packages ^ files['environment.yml'].packages
    if symmetric_difference:
        raise Exception(
            f"The following packages are in one but not both of the "
            f"environment.txt and environment.yml requirements files: "
            f"{symmetric_difference}")


def test_no_formatted_string_literals():
    """
    Test that there are no f-strings in top-level __init__.py.

    PlasmaPy's top-level __init__.py should not have any f-strings
    because then a very unhelpful SyntaxError will be raised instead of
    a helpful ImportError when attempting to import PlasmaPy in Python
    3.5 or below.

    Because we do not test that importing PlasmaPy from Python 2.7 or
    3.5 will raise the correct exception, an error like this could
    potentially remain hidden for a while unless we remember to check it
    manually.  This test checks each line in __init__.py and makes sure
    that there are no f-strings.

    """
    init_filename = dirname(dirname(__file__)) + "/__init__.py"
    init_file = open(init_filename)
    lines = init_file.readlines()
    init_file.close()

    fstring_regex_pattern = r'[ \(\[\+\{]f' + "['" + '"]'

    errors = []
    for line, line_number in zip(lines, range(1, len(lines) + 1)):
        if re.search(fstring_regex_pattern, line) is not None:
            errors.append((line, line_number))

    if errors:
        errmsg = (
            f"The file {init_filename} contains formatted string "
            f"literals (f-strings). However, f-strings may not be used "
            f"in PlasmaPy's main __init__.py so that the correct "
            f"exception is raised when attempting to import from "
            f"an unsupported version of Python. The lines containing "
            f"f-strings are:\n")
        for (line, line_number) in errors:
            errmsg += f"[Line {line_number}] {line}\n"
        raise SyntaxError(errmsg)
