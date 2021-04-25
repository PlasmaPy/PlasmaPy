"""
This `sphinx` extension package was highly influenced by
`Astropy's sphinx_automodapi
<https://sphinx-automodapi.readthedocs.io/en/latest/index.html>`_
package. In fact,
the :rst:dir:`automodapi` and :rst:dir:`automodsumm` directives defined in
this package are an adaption/evolution of the `sphinx_automodapi.automodapi`
and the `sphinx_automodapi.automodsumm` directives.

Defined Directives
------------------

A directive (`ref
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
#rst-directives>`_) is a generic block of explicit markup.  Along with roles, it
is one of the extension mechanisms of reST and, thus, Sphinx,

- :rst:dir:`automodapi`
- :rst:dir:`automodsumm`

Defined Configuration Values
----------------------------

Configuration values are variables that can be defined in the your `conf.py`
that controls the default behavior of the functionality defined by `plasmapy_sphinx`.

- :ref:`Values for automodapi <automodapi-confvals>`
- :ref:`Values for automodsumm <automodsumm-confvals>`

API
---
"""
# The code here was adapted from v0.14.0.dev0 of sphinx_automodapi

from sphinx.application import Sphinx

from . import automodapi, automodsumm, extras, utils


def setup(app: Sphinx):
    """The `sphinx` ``setup()`` function for the `plasmapy_sphinx` extension."""

    rtn = automodapi.setup(app)
    return rtn
