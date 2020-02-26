#!/usr/bin/env python
import os
import sys
from itertools import chain
import builtins

from setuptools import setup
from setuptools.config import read_configuration

# Append cwd for pip 19
sys.path.append(os.path.abspath("."))

################################################################################
# Programmatically generate some extras combos.
################################################################################
extras = read_configuration("setup.cfg")["options"]["extras_require"]

# Dev is everything
extras["dev"] = list(chain(*extras.values()))

# All is everything but tests and docs
exclude_keys = ("tests", "docs", "dev")
ex_extras = dict(filter(lambda i: i[0] not in exclude_keys, extras.items()))
# Concatenate all the values together for 'all'
extras["all"] = list(chain.from_iterable(ex_extras.values()))

VERSION_TEMPLATE = """
# Note that we need to fall back to the hard-coded version if either
# setuptools_scm can't be imported or setuptools_scm can't determine the
# version, so we catch the generic 'Exception'.
try:
    from setuptools_scm import get_version
    __version__ = get_version(root='..', relative_to=__file__)
except Exception:
    __version__ = '{version}'
""".lstrip()

# Get configuration information from all of the various subpackages.
# See the docstring for setup_helpers.update_package_files for more
# details.
setup(
    extras_require=extras,
    use_scm_version={
        "write_to": os.path.join("plasmapy", "version.py"),
        "write_to_template": VERSION_TEMPLATE,
    },
)
