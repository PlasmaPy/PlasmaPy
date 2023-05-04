"""
This sub-package contains functionality that extends `sphinx.ext.autodoc`.

*This functionality was highly influenced by and adapted from*
`sphinx.ext.autodoc` *and* `sphinx_automodapi.automodapi`.

.. contents:: Content
   :local:

Defined Directives
------------------

A directive (`ref
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
#rst-directives>`_) is a generic block of explicit markup.  Along with roles, it
is one of the extension mechanisms of reST and, thus, Sphinx.

+--------------------------------+-----------------------------------------------------+
| Directive                      | Description                                         |
+================================+=====================================================+
| | :rst:dir:`automodapi`        | An `~sphinx.ext.autodoc` directive that             |
| | ``.. automodapi:: modname``  | auto-generates documentation for a given module     |
|                                | ``modname`` (i.e. sub-package or ``.py`` file) by   |
|                                | inspecting and summarizing the object contained in  |
|                                | the module.                                         |
+--------------------------------+-----------------------------------------------------+

Defined Configuration Values
----------------------------

Configuration values are variables that can be defined in the ``conf.py`` file
to control the default behavior Sphinx and Sphinx extension packages like
`plasmapy_sphinx`.

+-------------------------------------------+------------------------------------------+
| Configuration Value                       | Description                              |
+===========================================+==========================================+
| :confval:`automodapi_default_toctree_dir` | Default directory for placing stub files |
|                                           | requested by the :rst:dir:`automodapi`   |
|                                           | directive.                               |
+-------------------------------------------+------------------------------------------+
| :confval:`automodapi_group_order`         | The order :rst:dir:`automodapi` displays |
|                                           | its group :rst:dir:`automodsumm` tables. |
+-------------------------------------------+------------------------------------------+
| |with_diagrams|                           | Define groups that should include        |
|                                           | inheritance diagrams.                    |
+-------------------------------------------+------------------------------------------+
| |include_diagram|                         | Control if :rst:dir:`automodapi` should  |
|                                           | display inheritance diagrams by default. |
+-------------------------------------------+------------------------------------------+

.. |with_diagrams| replace::
   :confval:`automodapi_groups_with_inheritance_diagrams`
.. |include_diagram| replace::
   :confval:`automodapi_include_inheritance_diagram`

"""
from sphinx.application import Sphinx

from . import automodapi


def setup(app: Sphinx):
    """
    Sphinx ``setup()`` function for setting up all of the `plasmapy_sphinx.autodoc`
    functionality, this includes `plasmapy_sphinx.automodsumm` functionality.
    """
    rtn = automodapi.setup(app)

    return rtn
