"""
This `Sphinx <https://www.sphinx-doc.org/en/master/>`__ extension package
contains custom documentation functionality used to document
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

To enable `plasmapy_sphinx`'s functionality it needs to be added to the
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
# The code here was adapted from v0.14.0.dev0 of sphinx_automodapi

from sphinx.application import Sphinx

from . import autodoc, automodsumm, directives, utils


def setup(app: Sphinx):
    """The `sphinx` ``setup()`` function for the `plasmapy_sphinx` extension."""

    # Note: automodsum is setup by autodoc.setup since it is needed for
    # autodoc.automodapi

    directives.setup(app)
    rtn = autodoc.setup(app)
    return rtn
