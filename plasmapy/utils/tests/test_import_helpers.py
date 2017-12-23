"""Tests for import_helpers.py"""

import pytest
import sys

from ..import_helpers import (check_python, check_versions)


def test_check_python():
    """Test that check_python will raise an ImportError when
    minimum_python_version is newer than the current version, and will
    not raise an ImportError when minimum_python_version is the same
    as or older than the current version."""

    python_version = sys.version.split()[0]

    try:
        check_python(minimum_python_version=python_version)
    except ImportError as e:
        raise Exception("check_python is raising an exception for the version "
                        f"of Python currently in use ({python_version})")

    major_str, minor_str, patch_str = python_version.split('.')

    major, minor, patch = int(major_str), int(minor_str), int(patch_str)

    newer_versions = [
        str(major+1) + '.0',
        str(major) + '.' + str(minor+1),
        str(major) + '.' + str(minor) + '.' + str(patch+1),
        ]

    older_versions = [
        str(major) + '.' + str(minor) + '.0',
        str(major-1) + '.' + str(minor),
        str(major) + '.' + str(minor-1),
        ]

    for newer_version in newer_versions:
        with pytest.raises(ImportError):
            check_python(minimum_python_version=newer_version)
            raise Exception(
                "check_python is not raising an exception when it should be "
                f"raising one, with python_version = {python_version} and "
                f"newer_version = {newer_version}")

    for older_version in older_versions:
        try:
            check_python(minimum_python_version=older_version)
        except ImportError as e:
            raise ImportError("check_python is raising an exception when it "
                              "should not be raising one, with python_version"
                              f" = {python_version} and older_version = "
                              f"{older_version}")
