"""Tests for import_helpers.py"""

import pytest
import sys
import numpy

from ..import_helpers import (check_python, check_versions)


def _split_version(version):
    """Separate a string including digits separated by periods into a
    tuple of integers."""
    return tuple(int(ver) for ver in version.split('.'))


def _create_different_versions(current_version, time):
    """Create a list containing version number strings that are either
    newer or older than/same as the inputted version string."""

    major, minor, *extra = _split_version(current_version)

    if len(extra) == 1:
        patch = extra[0]
    else:
        patch = 0

    if time == 'older':
        increment = -1
    elif time == 'newer':
        increment = +1

    different_versions = [
        str(major) + '.' + str(minor + increment),
        str(major + increment) + '.' + str(minor),
        ]

    if time == 'after' or patch > 0:
        different_versions.append(
            str(major) + '.' + str(minor) + '.' + str(patch + increment))

    if time == 'before':
        different_versions.append(current_version)
        different_versions.append(str(major) + '.' + str(minor) + '.0')

    return different_versions


_python_version = sys.version.split()[0]


@pytest.mark.parametrize(['older_version', 'newer_version'],
                         zip(_create_different_versions(_python_version, 'older'),
                             _create_different_versions(_python_version, 'newer')))
def test_check_python(older_version, newer_version):
    """Test that check_python will raise an ImportError when
    minimum_python_version is newer than the current version, and will
    not raise an ImportError when minimum_python_version is the same
    as or older than the current version."""

    with pytest.raises(ImportError):
        check_python(minimum_python_version=newer_version)
        raise Exception(
            "check_python is not raising an exception when it should be "
            f"raising one, with python_version = {_python_version} and "
            f"newer_version = {newer_version}")

    try:
        check_python(minimum_python_version=older_version)
    except ImportError as e:
        raise ImportError(
            "check_python is raising an exception when it should not be "
            f"raising one, with python_version = {_python_version} and "
            f"older_version = {older_version}") from e


@pytest.mark.parametrize(['older_version', 'newer_version'],
                         zip(_create_different_versions(numpy.__version__, 'older'),
                             _create_different_versions(numpy.__version__, 'newer')))
def test_check_versions(older_version, newer_version):
    """Test that check_versions will raise an ImportError when the
    minimum version for NumPy is newer than the current version, and
    will not raise an ImportError when the minimum version is the same
    as or older than the current version."""

    newer_minimum_versions = {'numpy': newer_version}
    with pytest.raises(ImportError):
        check_versions(newer_minimum_versions)
        raise Exception(
            "check_versions is not raising an exception when it should be "
            f"raising one, with numpy.__version__ = {numpy.__version__} "
            f"and newer_version = {newer_version}.")

    older_minimum_versions = {'numpy': older_version}
    try:
        check_versions(older_minimum_versions)
    except ImportError as e:
        raise ImportError(
            "check_versions is raising an exception when it should not be "
            f"raising one, with numpy.__version__ = {numpy.__version__} "
            f"and older_version = {older_version}")
