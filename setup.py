#!/usr/bin/env python
# https://github.com/pypa/pip/issues/7953#issuecomment-645133255
import site
import sys

from setuptools import setup

site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

# Get configuration information from all of the various subpackages.
# See the docstring for setup_helpers.update_package_files for more
# details.
setup(use_scm_version=True)
