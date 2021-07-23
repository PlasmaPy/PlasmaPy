"""
Functionality for declaring and cross-referencing Sphinx configuration values.
Configuration values are variables that can be defined in the ``conf.py`` file
to control the default behavior for Sphinx and extension packages like
`plasmapy_sphinx`.  This code was taken and adapted from the ``conf.py`` files
of the `sphinx` and `sphinx_rtd_theme` packages.

.. rst:directive:: .. confval:: name

    | A directive used for documenting `Sphinx configuration values
      <https://www.sphinx-doc.org/en/master/usage/configuration.html>`_.  While
      this directive is not provided by Sphinx, it is used by Sphinx and, thus,
      cross-linking is provided through the :confval:`intersphinx_mapping`
      configuration value with the `sphinx.ext.intersphinx` extension.

    .. rubric:: Optional Fields

    .. rst:directive:option:: type
       :type: string

       An optional flag that specifies the configuration value's data type.

    .. rst:directive:option:: default
       :type: string

       An optional flag that specifies the default value for the configuration value.

    .. rubric:: Example

    The following example illustrates how to document the ``dummy_value``
    configuration value.

    .. code-block:: rst

        .. confval:: dummy_value

           :type: `bool`
           :default: `True`

           This is an example documentation for the configuration value
           ``dummy_value``.

    The code renders like...

    .. confval:: dummy_value

       :type: `bool`
       :default: `True`

       This is an example documentation for the configuration value
       ``dummy_value``.

.. rst:role:: confval

    This role is provided for easy cross-linking to a configuration value's
    definition.  For example, doing ``:confval:`dummy_value``` will cross-link
    to the ``dummy_value`` configuration value like :confval:`dummy_value`.  Or,
    a link to Sphinx's ``intersphinx_mapping`` configuration value goes like
    ``:confval:`intersphinx_mapping``` -> :confval:`intersphinx_mapping`.

    *Linking to external packages is made possible when using*
    `sphinx.ext.intersphinx`.

"""
from sphinx.application import Sphinx
from sphinx.domains.python import PyField
from sphinx.locale import _
from sphinx.util.docfields import Field


def setup(app: Sphinx) -> None:
    """
    A `sphinx` ``setup()`` function setting up the :rst:dir:`confval` directive
    and :rst:role:`confval` role.
    """
    # this was taken from the sphinx and sphinx_rtd_theme conf.py files and creates
    # the documenting directive `.. confval::` and role `:confval:` for documenting
    # sphinx configuration variables
    app.add_object_type(
        directivename="confval",
        rolename="confval",
        objname="configuration value",
        indextemplate="pair: %s; configuration value",
        doc_field_types=[
            PyField(
                name="type",
                names=("type",),
                label=_("Type"),
                has_arg=False,
                bodyrolename="class",
            ),
            Field(
                name="default",
                names=("default",),
                label=_("Default"),
                has_arg=False,
            ),
        ],
    )
