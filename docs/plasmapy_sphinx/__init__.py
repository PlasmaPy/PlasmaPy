"""
This `sphinx` documentation package was highly influenced by
`Astropy's <https://www.astropy.org/>`_
`sphinx_automodapi <https://sphinx-automodapi.readthedocs.io/en/latest/index.html>`_
package. In fact,
the `~automdoapi.automodapi` and `~automodapi.automodsumm` directive defined in
this package are a adaption/evolutions of `sphinx_automodapi.automodapi`
and `sphinx_automodapi.automodsumm` directives.
"""
# The code here was adapted from v0.14.0.dev0 of sphinx_automodapi

from sphinx.application import Sphinx

from . import automodapi, automodsumm, utils


def setup(app: Sphinx):
    """The `sphinx` ``setup()`` function for the `plasmapy_sphinx` extension."""

    rtn = automodapi.setup(app)
    return rtn
