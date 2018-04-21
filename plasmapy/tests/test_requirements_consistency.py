"""
Test for common discrepancies and inconsistencies in requirements files.
"""

import pytest
from collections import namedtuple

from ..utils import RunTestError

from ..import_helpers import requirements_dir

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
])

files = {}

for filename in filenames:
    file = open(requirements_dir + filename, 'rt')
    readlines = file.readlines()
    lines = {line.strip('\n') for line in readlines}
    files[filename] = Requirements(file, readlines, lines)


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
        raise RunTestError(
            f"The following lines from requirements.txt are not in "
            f"automated-code-tests.txt: {missing_lines}")
