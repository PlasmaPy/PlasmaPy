#!/usr/bin/env python
import importlib
import os
import pkgutil
import sys
import plasmapy.addons

from itertools import chain
from setuptools import setup
from setuptools.config import read_configuration

# Append cwd for pip 19
sys.path.append(os.path.abspath("."))

################################################################################
# Programmatically generate some extras combos.
################################################################################
extras = read_configuration("setup.cfg")['options']['extras_require']

# Dev is everything
extras['dev'] = list(chain(*extras.values()))

# All is everything but tests and docs
exclude_keys = ("tests", "docs", "dev")
ex_extras = dict(filter(lambda i: i[0] not in exclude_keys, extras.items()))
# Concatenate all the values together for 'all'
extras['all'] = list(chain.from_iterable(ex_extras.values()))


# Determine entry_points for namespace pacakge plasmapy.addons
# - this needs to be so that the namespace is actually populated
# - copied from https://packaging.python.org/guides/creating-and-discovering-plugins/
def iter_namespace(ns_pkg):
    # Specifying the second argument (prefix) to iter_modules makes the
    # returned name an absolute name instead of a relative one. This allows
    # import_module to work without having to do additional modification to
    # the name.
    return pkgutil.iter_modules(ns_pkg.__path__, ns_pkg.__name__ + ".")


discovered_plugins = {
    name: importlib.import_module(name)
    for finder, name, ispkg
    in iter_namespace(plasmapy.addons)
}
extras['entry_points'] = discovered_plugins

# Get configuration information from all of the various subpackages.
# See the docstring for setup_helpers.update_package_files for more
# details.
setup(extras_require=extras, use_scm_version=True)
