"""
This module contains functionality for auto-generating the stub files related to
the :rst:dir:`automodapi` and :rst:dir:`automodsumm` directives.
"""
__all__ = ["AutomodsummEntry", "AutomodsummRenderer"]

import os

from jinja2 import TemplateNotFound
from sphinx.ext.autosummary.generate import (
    AutosummaryEntry,
    AutosummaryRenderer,
)
from typing import Dict, Union

from .utils import templates_dir

if False:
    # noqa
    # for annotation, does not need real import
    from sphinx.application import Sphinx
    from sphinx.builders import Builder


class AutomodsummEntry(AutosummaryEntry):
    """
    A typed version of `~collections.namedtuple` representing an stub file
    entry for :rst:dir:`automodsumm`.
    """


class AutomodsummRenderer(AutosummaryRenderer):
    """
    A helper class for retrieving and rendering :rst:dir:`automodsumm` templates
    when writing stub files.

    Parameters
    ----------

    app : `sphinx.application.Sphinx`
        Instance of the `sphinx` application.

    template_dir : str
        Path to a specified template directory.
    """

    def __init__(
        self, app: Union["Builder", "Sphinx"], template_dir: str = None,
    ) -> None:

        asumm_path = templates_dir
        relpath = os.path.relpath(asumm_path, start=app.srcdir)
        app.config.templates_path.append(relpath)
        super().__init__(app, template_dir)

    def render(self, template_name: str, context: Dict) -> str:
        """
        Render a template file.  The render will first search for the template in
        the path specified by the sphinx configuration value :confval:`templates_path`,
        then the `~plasmapy_sphinx.templates_dir, and finally the
        :rst:dir:`autosummary` templates directory.  Upon finding the template,
        the values from the ``context`` dictionary will inserted into the
        template and returned.

        Parameters
        ----------
        template_name : str
            Name of the template file.

        context: dict
            Dictionary of values to be rendered (inserted) into the template.
        """
        if not template_name.endswith(".rst"):
            # if does not have '.rst' then objtype likely given for template_name
            template_name += ".rst"

        template = None
        for name in [template_name, "base.rst"]:
            for _path in ["", "automodapi/", "autosummary/"]:
                try:
                    template = self.env.get_template(_path + name)
                    return template.render(context)
                except TemplateNotFound:
                    pass

        if template is None:
            raise TemplateNotFound
