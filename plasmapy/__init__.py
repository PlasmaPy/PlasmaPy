"""
PlasmaPy: A plasma physics Python package
================================================

Documentation is available in the docstrings,
online at https://docs.plasmapy.org (accessible also using
the ``plasmapy.online_help`` function).

Contents
--------
PlasmaPy provides the following functionality:

Subpackages
-----------
Each of these subpackages (except for `formulary` and `atomic`) requires an
explicit import, for example, via ``import plasmapy.diagnostics``.

::

 atomic                            --- Database for atoms, isotopes, ions...
 classes                           --- (WIP) classes used in multiple places
 data                              --- Data used for testing and examples
 diagnostics                       --- Experimental research data analysis
 formulary                         --- Plasma theory analysis formulae
 utils                             --- Various utilities

Utility tools
-------------
::

 online_help       --- Search the online documentation
 __version__       --- PlasmaPy version string
 __citation__      --- PlasmaPy citation instructions

"""
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from .version import version as __version__
from .version import githash as __githash__

from . import formulary
from . import atomic
# ----------------------------------------------------------------------------

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

__citation__ = (
    "Instructions on how to cite and acknowledge PlasmaPy are provided in the "
    "online documentation at: http://docs.plasmapy.org/en/latest/about/citation.html"
)

if sys.version_info < tuple((int(val) for val in "3.6".split('.'))):
    raise Exception("PlasmaPy does not support Python < {}".format(3.6))


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

    url = ('http://docs.plasmapy.org/en/stable/search.html?'
           '{0}&check_keywords=yes&area=default').format(urlencode({'q': query}))

    if(query.lower() in ('unit', 'units')):
        url = 'http://docs.astropy.org/en/stable/units/'

    webbrowser.open(url)


del sys
del test
