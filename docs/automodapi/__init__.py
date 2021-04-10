# The code here was adapted from v0.14.0.dev0 of sphinx_automodapi
__all__ = ["package_dir", "templates_dir"]

import os

from sphinx.application import Sphinx


package_dir = os.path.abspath(os.path.dirname(__file__))
templates_dir = os.path.join(package_dir, "templates")


def setup(app: Sphinx):
    from .automodapi import setup as setup_automodapi

    rtn = setup_automodapi(app)

    return rtn
