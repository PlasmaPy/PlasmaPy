from sphinx.application import Sphinx

from . import core, generate


def setup(app: Sphinx):
    """
    Sphinx ``setup()`` function for setting up the :rst:dir:`automodsumm`
    functionality.
    """
    rtn = core.setup(app)
    return rtn
