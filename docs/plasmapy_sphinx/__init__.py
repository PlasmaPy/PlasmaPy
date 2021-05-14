"""
This `Sphinx <https://www.sphinx-doc.org/en/master/>`_ extension package contains
custom documentation functionality used to document
`PlasmaPy's packages <https://github.com/PlasmaPy>`_ .

.. contents:: Content
   :local:

Installation
------------

This package is currently not released in any form and can only be obtained by
installing `plasmapy` directly from its
`GitHub repository <https://github.com/plasmapy/plasmapy>`_.  We do plan to breakout
`plasmapy_sphinx` into its own repository and release it to https://pypi.org at
a future date.

Setup
-----

To enable of of `plasmapy_sphinx`'s functionality it needs to be added to the
:confval:`extensions` configuration value in your ``conf.py`` file.

.. code-block::

    extensions = ["plasmapy_sphinx"]

If you don't want all of `plasmapy_sphinx`'s functionality, then any module
containing a Sphinx ``setup()`` function can be activated independently.  For
example, you can just implement the :rst:dir:`automodapi` functionality by
setting up the `plasmapy_sphinx.autodoc.automodapi` module like...

.. code-block::

    extensions = ["plasmapy_sphinx.autodoc.automodapi"]

The Rundown
-----------

`plasmapy_sphinx.autodoc`
~~~~~~~~~~~~~~~~~~~~~~~~~

   The `plasmapy_sphinx.autodoc` sub-package contains functionality that extends
   `sphinx.ext.autodoc`.  The defined directives are registered as
   `sphinx.ext.autodoc` directives and, thus, will behave similarly to
   :rst:dir:`autoclass`, :rst:dir:`autofunction`, etc. and work with
   `sphinx.ext.autosummary` and `plasmapy_sphinx.automodsumm`.

`plasmapy_sphinx.automodsumm`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   The `plasmapy_sphinx.automodsumm` sub-package is an evolution of
   `sphinx.ext.autosummary` and was adapted from `sphinx_automodapi.automodsumm`.
   Unlike the :rst:dir:`autosummary` directive where you have to manually list
   all the objects/members to appear in the summary table, the
   :rst:dir:`automodsumm` directive will inspect a given module and automatically
   populate the summary table with the detected objects/members.

`plasmapy_sphinx.directives`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   This sub-package contains custom reStructuredText directives and roles that
   do not fall under `plasmapy_sphinx.autodoc` or `plasmapy_sphinx.automodsumm`.

"""
# TODO: delete the below docstring when ready
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
is one of the extension mechanisms of reST and, thus, Sphinx.

- :rst:dir:`automodapi`
- :rst:dir:`automodsumm`

Defined Configuration Values
----------------------------

Configuration values are variables that can be defined in the `conf.py` file
to control the default behavior of the functionality defined by
`plasmapy_sphinx`.

- :ref:`Values for automodapi <automodapi-confvals>`
- :ref:`Values for automodsumm <automodsumm-confvals>`

API
---
"""
# The code here was adapted from v0.14.0.dev0 of sphinx_automodapi

from sphinx.application import Sphinx

from . import autodoc, automodsumm, extras, utils


def setup(app: Sphinx):
    """The `sphinx` ``setup()`` function for the `plasmapy_sphinx` extension."""

    rtn = autodoc.setup(app)
    return rtn
