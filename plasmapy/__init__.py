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
Each of these subpackages requires an explicit import, for example,
via ``import plasmapy.physics``.

::

 atomic                            --- Database for atoms, isotopes, ions...
 classes                           --- (WIP) classes used in multiple places
 constants                         --- (WIP?) wrapper for astropy.constants.si
 data                              --- Data used for testing and examples
 diagnostics                       --- Experimental research data analysis
 mathematics                       --- General formulae used elsewhere
 physics                           --- Plasma theory functionality
 transport                         --- Transport theory functionality
 utils                             --- Various utilities

Utility tools
-------------
::

 test              --- Run PlasmaPy unit tests
 online_help       --- Search the online documentation
 __version__       --- PlasmaPy version string
 __citation__      --- PlasmaPy citation template

"""
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


if sys.version_info < tuple((int(val) for val in "3.6".split("."))):
    raise RuntimeError("plasmapy does not support Python < {}".format(3.6))


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

    url = (
        "http://docs.plasmapy.org/en/stable/search.html?" "{0}&check_keywords=yes&area=default"
    ).format(urlencode({"q": query}))

    if query.lower() in ("unit", "units"):
        url = "http://docs.astropy.org/en/stable/units/"

    webbrowser.open(url)


__citation__ = """@misc{plasmapy_community_2018_1238132,
  author       = {{PlasmaPy Community} and
                  Murphy, Nicholas A. and
                  Leonard, Andrew J. and
                  Sta\'nczak, Dominik and
                  Kozlowski, Pawel M. and
                  Langendorf, Samuel J. and
                  Haggerty, Colby C. and
                  Beckers, Jasper P. and
                  Mumford, Stuart J. and
                  Parashar, Tulasi N. and
                  Huang, Yi-Min},
  title        = {{PlasmaPy: an open source community-developed 
                   Python package for plasma physics}},
  month        = apr,
  year         = 2018,
  doi          = {10.5281/zenodo.1238132},
  url          = {https://doi.org/10.5281/zenodo.1238132}
}"""
