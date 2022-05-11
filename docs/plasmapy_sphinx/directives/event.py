"""
Functionality for declaring and cross-referencing
`Sphinx events
<https://www.sphinx-doc.org/en/master/extdev/appapi.html#sphinx-core-events>`_.
Sphinx events occur at specific points in the Sphinx build.  When an event
is reached a signal is emitted with `sphinx.application.Sphinx.emit` that causes
the build to "pause" and allow any connected functionality to do the desired
processing.  New functionality can be connected to an event with
`sphinx.application.Sphinx.connect` and new events can be created using
`sphinx.application.Sphinx.add_event`.

The code contained in this module was taken and adapted from the ``conf.py`` file
of `Sphinx's documentation
<https://github.com/sphinx-doc/sphinx/blob/
8653ceca0021f6ac6ff0aac6c26e2a455c6d4b21/doc/conf.py#L123-L138>`_.

.. rst:directive:: .. event:: name (signature)

    A directive used for documenting Sphinx events.  While this directive is not
    provided by Sphinx, it is usd by Sphinx and, thus, cross-linking is provided
    through the :confval:`intersphinx_mapping` configuration value with the
    `sphinx.ext.intersphinx` extension.

    The directive is used like...

    .. code-block:: rst

       .. event:: name (signature)

          Event documentation.

    | where ``name`` is the event name and ``signature`` represents what the
      connected function's signature should look like.

    .. rubric:: Example

    The following example illustrates how to document the ``dummy-event`` event.

    .. code-block:: rst

       .. event:: dummy-event (app, node)

          :param app: Instance of the Sphinx applications.
          :param node: The pending node to be resolved.

          This is just a dummy event for demonstration purposes.

    The code renders like...

    .. event:: dummy-event (app, node)

       :param app: Instance of the Sphinx applications.
       :param node: The pending node to be resolved.

       This is just a dummy event for demonstration purposes.


.. rst:role:: event

   This role is provided for easy cross-linking to an event's definition.  For
   example, doing ``:event:`dummy-event``` will cross-link to the
   ``dummy-event`` definition like :event:`dummy-event`.  Or, a link to
   Sphinx's ``builder-inited`` event goes like ``:event:`builder-inited``` ->
   :event:`builder-inited`.

   *Linking to external packages is made possible when using*
   `sphinx.ext.intersphinx`.

"""

import re

from sphinx import addnodes
from sphinx.application import Sphinx
from sphinx.util.docfields import GroupedField


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
    event_sig_re = re.compile(r"([a-zA-Z-_]+)\s*\((.*)\)")

    match = event_sig_re.match(sig)
    if not match:
        signode += addnodes.desc_name(sig, sig)
        return sig
    name, args = match.groups()
    signode += addnodes.desc_name(name, name)
    plist = addnodes.desc_parameterlist()
    for arg in args.split(","):
        arg = arg.strip()
        plist += addnodes.desc_parameter(arg, arg)
    signode += plist
    return name


def setup(app: Sphinx) -> None:
    """
    A `sphinx` ``setup()`` function setting up the :rst:dir:`event` directive
    and :rst:role:`event` role.
    """
    # this was taken from the sphinx conf.py file and creates the documenting
    # directive `.. event::` and role `:event:` for documenting sphinx events
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
