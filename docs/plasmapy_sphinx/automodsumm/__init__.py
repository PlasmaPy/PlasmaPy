"""
This sub-package contains functionality that defines the :rst:dir:`automodsumm`
directive and the
`stub file generation
<https://www.sphinx-doc.org/en/master/usage/extensions/autosummary.html
#sphinx-autogen-generate-autodoc-stub-pages>`_
for items listed in :rst:dir:`automodsumm` tables.

*This functionality was highly influenced by and adapted from*
`sphinx.ext.autosummary` *and* `sphinx_automodapi.automodsumm`.

.. contents:: Content
   :local:

Defined Directives
------------------

+--------------------------------+-----------------------------------------------------+
| Directive                      | Description                                         |
+================================+=====================================================+
| | :rst:dir:`automodsumm`       | A directive that generates and auto-populates       |
| | ``.. automodsumm:: modname`` | `~sphinx.ext.autosummary` tables for a given module |
|                                | ``modname`` (i.e. sub-package or ``.py`` file).     |
+--------------------------------+-----------------------------------------------------+

Defined Configuration Values
----------------------------

+--------------------------------------------------+-----------------------------------+
| Configuration Value                              | Description                       |
+==================================================+===================================+
| :confval:`automodapi_custom_groups`              | Used to define custom groups to   |
|                                                  | be displayed by the               |
|                                                  | :rst:dir:`automodsumm` and        |
|                                                  | :rst:dir:`automodapi` directives  |
+--------------------------------------------------+-----------------------------------+
| :confval:`automodapi_generate_module_stub_files` | Used to control is stub files are |
|                                                  | by default generated to modules   |
|                                                  | (i.e. sub-packages and ``.py``    |
|                                                  | files).                           |
+--------------------------------------------------+-----------------------------------+

Connected Sphinx Events
-----------------------

+------------------------------+-----------------------------------------------+
| Event                        | Connected                                     |
+==============================+===============================================+
| :event:`builder-inited`      | |gendoc|                                      |
+------------------------------+-----------------------------------------------+
| :event:`autodoc-skip-member` | |skip_mem|                                    |
+------------------------------+-----------------------------------------------+

.. |gendoc| replace:: `~plasmapy_sphinx.automodsumm.generate.GenDocsFromAutomodsumm`
.. |skip_mem| replace::
   `~plasmapy_sphinx.automodsumm.generate.GenDocsFromAutomodsumm.event_handler__autodoc_skip_member`

"""

from sphinx.application import Sphinx

from . import core, generate


def setup(app: Sphinx):
    """
    Sphinx ``setup()`` function for setting up the :rst:dir:`automodsumm`
    functionality.
    """
    rtn = core.setup(app)
    return rtn
