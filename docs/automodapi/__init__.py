# The code here was adapted from v0.14.0.dev0 of sphinx_automodapi
__all__ = ["package_dir"]

import os

from sphinx.application import Sphinx


package_dir = os.path.abspath(os.path.dirname(__file__))
templates_dir = os.path.join(package_dir, "templates")


from .automodapi import setup_automodapi
from .automodsumm import setup_autosummary


def setup(app: Sphinx):
    setup_automodapi(app)
    # setup_autosummary(app)

    return {"parallel_read_safe": True, "parallel_write_safe": True}
