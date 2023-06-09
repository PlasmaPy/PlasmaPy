"""
Welcome to the `plasmapy` package, an open source community-developed
Python package for the plasma community. Documentation is available in
the docstrings and `online <https://docs.plasmapy.org>`__ (accessible
also using the :func:`plasmapy.online_help` function).
"""
__all__ = [
    "online_help",
    "analysis",
    "diagnostics",
    "dispersion",
    "formulary",
    "particles",
    "plasma",
    "simulation",
    "utils",
    "__version__",
    "__citation__",
]

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

if sys.version_info < (3, 9):  # coverage: ignore # noqa: UP036
    raise ImportError(
        f"This version of PlasmaPy does not support Python {sys.version.split()[0]}."
        "Please upgrade to a newer version."
    )

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from importlib.metadata import PackageNotFoundError

from plasmapy import (
    analysis,
    diagnostics,
    dispersion,
    formulary,
    particles,
    plasma,
    simulation,
    utils,
)

# define version
try:
    try:
        from plasmapy._dev.scm_version import version as __version__
    except ImportError:
        from plasmapy._version import version as __version__
except Exception:  # coverage: ignore  # noqa: BLE001
    # package is not installed
    __version__ = "0.0.0"

    from warnings import warn

    warn(
        "plasmapy.__version__ not generated (set to '0.0.0'). It looks like "
        "the installation's broken. Ask on Element!",
    )

    del warn

# ----------------------------------------------------------------------------
#: PlasmaPy citation instructions
__citation__ = (
    "Instructions on how to cite and acknowledge PlasmaPy are provided in the "
    "online documentation at: http://docs.plasmapy.org/en/stable/about/citation.html"
)


def online_help(query: str):  # coverage: ignore
    """
    Open a webpage containing a search page in `PlasmaPy's documentation`_,
    or another page that contains relevant online help.

    This function requires an active internet connection, and will open
    the page in the default web browser.

    Parameters
    ----------
    query : str
        The search query.
    """
    import webbrowser

    from urllib.parse import urlencode

    url = (
        "http://docs.plasmapy.org/en/stable/search.html?"
        "{}&check_keywords=yes&area=default"
    ).format(urlencode({"q": query}))

    if query.lower() in {"unit", "units", "quantity", "quantities"}:
        url = "http://docs.astropy.org/en/stable/units/"

    webbrowser.open(url)


del sys
