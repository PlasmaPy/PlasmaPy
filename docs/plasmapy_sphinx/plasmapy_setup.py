from sphinx.application import Sphinx


def setup(app: Sphinx):
    from . import setup as setup_automodapi
    from .extras import setup as setup_extras
    from .documenters import setup as setup_documenters

    setup_extras(app)
    setup_documenters(app)
    rtn = setup_automodapi(app)

    return rtn
