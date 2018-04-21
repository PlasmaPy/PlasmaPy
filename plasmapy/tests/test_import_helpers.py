"""Tests for import_helpers.py"""

import pytest
import numpy

from ..import_helpers import check_versions, split_version


def _create_distinct_versions(current_version, time):
    """Create a list containing version number strings that are either
    newer or older than/same as the inputted version string."""

    major, minor, *extra = split_version(current_version)

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


@pytest.mark.parametrize('newer_version',
                         _create_distinct_versions(numpy.__version__, 'newer'))
def test_check_versions_newer(newer_version):
    """Test that check_versions will raise an ImportError when the
    minimum version for NumPy is newer than the current version.
    """

    newer_minimum_versions = {'numpy': newer_version}
    with pytest.raises(ImportError):
        check_versions(newer_minimum_versions)
        raise Exception(
            "check_versions is not raising an exception when it should be "
            f"raising one, with numpy.__version__ = {numpy.__version__} "
            f"and newer_version = {newer_version}.")


@pytest.mark.parametrize('older_version',
                         _create_distinct_versions(numpy.__version__, 'older'))
def test_check_versions_older(older_version):
    """Test that check_versions will not raise an ImportError when the minimum
    version is the same as or older than the current version.
    """

    older_minimum_versions = {'numpy': older_version}
    try:
        check_versions(older_minimum_versions)
    except ImportError as e:
        raise ImportError(
            "check_versions is raising an exception when it should not be "
            f"raising one, with numpy.__version__ = {numpy.__version__} "
            f"and older_version = {older_version}") from e
