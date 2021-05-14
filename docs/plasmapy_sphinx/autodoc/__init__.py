"""My docstring."""
from sphinx.application import Sphinx

from . import automodapi


def setup(app: Sphinx):
    rtn = automodapi.setup(app)

    return rtn
