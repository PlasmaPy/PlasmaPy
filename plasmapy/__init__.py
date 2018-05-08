# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._base_init import *
# ----------------------------------------------------------------------------

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

__name__ = "plasmapy"

__doc__ = ("A community-developed and community-driven open source "
           "core Python package for plasma physics.")

class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple((int(val) for val in "3.6".split('.'))):
    raise UnsupportedPythonError("plasmapy does not support Python < {}".format(3.6))


if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.
    from . import atomic
    from . import classes
    from . import constants
    from . import diagnostics
    from . import mathematics
    from . import physics
    from . import utils


def online_help(query):
    """
    Search the online PlasmaPy documentation for the given query from plasmapy.org
    Opens the results in the default web browser.
    Requires an active Internet connection.
    Redirects to Astropy.units in case of query 'unit' or 'units'

    Parameters
    ----------
    query : str
        The search query.
    """
    from urllib.parse import urlencode
    import webbrowser

    url = 'http://docs.plasmapy.org/en/stable/search.html?\
    {0}&check_keywords=yes&area=default'.format(urlencode({'q': query}))

    if(query.lower() in ('unit','units')):
        url = 'http://docs.astropy.org/en/stable/units/'

    webbrowser.open(url)
