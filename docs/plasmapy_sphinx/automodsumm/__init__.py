from sphinx.application import Sphinx

from . import core, generate


def setup(app: Sphinx):
    rtn = core.setup(app)
    return rtn
