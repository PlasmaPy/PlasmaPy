"""
This module contains extra (optional) Sphinx setup parameters that are not needed
to run the `plasmapy_sphinx` extension but are useful for documenting.  See
`plasmapy_sphinx.extras.setup` for further details.
"""

import re

from sphinx import addnodes
from sphinx.application import Sphinx
from sphinx.domains.python import PyField
from sphinx.locale import _
from sphinx.util.docfields import Field, GroupedField


def parse_event(env, sig, signode):
    """
    Used to set up the ``event`` directive and role for documenting Sphinx events.
    Taken from the ``conf.py`` file of `Sphinx's documentation
    <https://github.com/sphinx-doc/sphinx/blob/
    8653ceca0021f6ac6ff0aac6c26e2a455c6d4b21/doc/conf.py#L123-L138>`_.

    Parameters
    ----------
    env : sphinx.environment.BuildEnvironment
        Instance of the Sphinx's build environment.

    sig : str
        The "signature" given the the event directive or role.  For example,

        .. code-block:: rst

            .. event:: foo(bar)

            :event:`foo`

        in the directive case ``foo(bar)`` would be the signature and in the role
        case ``foo`` would be the signature.

    signode : sphinx.addnodes.desc_signature
        A `docutils` Node for the object signatures.
    """
    event_sig_re = re.compile(r"([a-zA-Z-]+)\s*\((.*)\)")

    m = event_sig_re.match(sig)
    if not m:
        signode += addnodes.desc_name(sig, sig)
        return sig
    name, args = m.groups()
    signode += addnodes.desc_name(name, name)
    plist = addnodes.desc_parameterlist()
    for arg in args.split(","):
        arg = arg.strip()
        plist += addnodes.desc_parameter(arg, arg)
    signode += plist
    return name


def setup(app: Sphinx):
    """
    A `sphinx` ``setup()`` function for the extension package `plasmapy_sphinx`
    to set up convenient, but optional, functionality.  This adds:

    1. ``confval`` standard directive and role for documenting and cross-linking
       Sphinx configuration values.
    2. ``event`` standard directive and role for documenting and cross-linking
       Sphinx events.
    """
    # this was taken from the sphinx and sphinx_rtd_theme conf.py files and creates
    # the documenting directive `.. confval::` and role `:confval:` for documenting
    # sphinx configuration variables
    app.add_object_type(
        "confval",
        "confval",
        objname="configuration value",
        indextemplate="pair: %s; configuration value",
        doc_field_types=[
            PyField(
                "type",
                label=_("Type"),
                has_arg=False,
                names=("type",),
                bodyrolename="class",
            ),
            Field(
                "default",
                label=_("Default"),
                has_arg=False,
                names=("default",),
            ),
        ],
    )

    # this was taken from the sphinx conf.py file and creates the documenting
    # directive `.. confval::` and role `:confval:` for documenting sphinx events
    app.add_object_type(
        directivename="event",
        rolename="event",
        indextemplate="pair: %s; event",
        parse_node=parse_event,
        doc_field_types=[
            GroupedField(
                "parameter",
                label="Parameters",
                names=("param",),
                can_collapse=True,
            )
        ],
    )
