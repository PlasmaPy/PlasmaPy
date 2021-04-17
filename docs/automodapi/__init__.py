# The code here was adapted from v0.14.0.dev0 of sphinx_automodapi

from sphinx.application import Sphinx

from . import automodapi, automodsumm, utils


def setup(app: Sphinx):
    rtn = automodapi.setup(app)
    return rtn
