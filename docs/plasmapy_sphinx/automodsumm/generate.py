"""
This module contains functionality for auto-generating the stub files related to
the :rst:dir:`automodapi` and :rst:dir:`automodsumm` directives.
"""
__all__ = ["AutomodsummEntry", "AutomodsummRenderer", "GenDocsFromAutomodsumm"]

import os
import re

from jinja2 import TemplateNotFound
from sphinx.ext.autodoc.mock import mock
from sphinx.ext.autosummary import get_rst_suffix, import_by_name, import_ivar_by_name
from sphinx.ext.autosummary.generate import (
    AutosummaryEntry,
    AutosummaryRenderer,
    generate_autosummary_content,
)
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.osutil import ensuredir
from typing import Any, Dict, List, Union

from ..utils import templates_dir

if False:
    # noqa
    # for annotation, does not need real import
    from sphinx.application import Sphinx
    from sphinx.builders import Builder


logger = logging.getLogger(__name__)


class AutomodsummEntry(AutosummaryEntry):
    """
    A typed version of `~collections.namedtuple` representing an stub file
    entry for :rst:dir:`automodsumm`.

    Parameters
    ----------
    name : `str`
        The objects fully qualified name of the object for which the stub file
        will be generated.

    path : `str`
        Absolute file path to the toctree directory.  This is where the stub
        file will be placed.

    recursive : `bool`
        Specifies if stub file for modules and and sub-packages should be
        generated.

    template : `str`
        Name of the template file to be used in generating the stub file.
    """


class AutomodsummRenderer(AutosummaryRenderer):
    """
    A helper class for retrieving and rendering :rst:dir:`automodsumm` templates
    when writing stub files.

    Parameters
    ----------

    app : `sphinx.application.Sphinx`
        Instance of the `sphinx` application.
    """

    def __init__(self, app: "Sphinx") -> None:

        # add plasmapy_sphinx templates directory to the overall templates path
        asumm_path = templates_dir
        relpath = os.path.relpath(asumm_path, start=app.srcdir)
        app.config.templates_path.append(relpath)

        super().__init__(app)

    def render(self, template_name: str, context: Dict) -> str:
        """
        Render a template file.  The render will first search for the template in
        the path specified by the sphinx configuration value :confval:`templates_path`,
        then the `~plasmapy_sphinx.utils.templates_dir`, and finally the
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
            for _path in ["", "automodsumm/", "autosummary/"]:
                try:
                    template = self.env.get_template(_path + name)
                    return template.render(context)
                except TemplateNotFound:
                    pass

        if template is None:
            raise TemplateNotFound


class GenDocsFromAutomodsumm:
    """
    Class used for stub file generation from :rst:dir:`automodapi` and
    :rst:dir:`automodsumm`.  An instance of the class is connected to the Sphinx
    event :event:`builder-inited`, which is emitted when the builder object is
    created.
    """

    _re = {
        "automodsumm": re.compile(r"^\n?(\s*)\.\.\s+automodsumm::\s*(\S+)\s*(?:\n|$)"),
        "automodapi": re.compile(r"^\n?(\s*)\.\.\s+automodapi::\s*(\S+)\s*(?:\n|$)"),
        "option": re.compile(r"^\n?(\s+):(\S*):\s*(\S.*|)\s*(?:\n|$)"),
        "currentmodule": re.compile(
            r"^\s*\.\.\s+(|\S+:)(current)?module::\s*([a-zA-Z0-9_.]+)\s*$"
        ),
    }
    """
    Dictionary of regular expressions used for string matching a read document
    and identify key directives.
    """

    app = None  # type: "Sphinx"
    """Instance of the Sphinx application."""

    logger = logger
    """
    Instance of the `~sphinx.util.logging.SphinxLoggerAdapter` for report during
    builds.
    """

    def __call__(self, app: "Sphinx"):
        """
        Scan through source files, check for the :rst:dir:`automodsumm` and
        :rst:dir:`automodapi` directives, and auto generate any associated
        stub files.

        Parameters
        ----------
        app :  `~sphinx.application.Sphinx`
            Instance of the Sphinx application.


        .. note:: Adapted from :func:`sphinx.ext.autosummary.process_generate_options`.
        """
        self.app = app
        genfiles = app.config.autosummary_generate

        if genfiles is True:
            env = app.builder.env
            genfiles = [
                env.doc2path(x, base=None)
                for x in env.found_docs
                if os.path.isfile(env.doc2path(x))
            ]
        elif genfiles is False:
            pass
        else:
            ext = list(app.config.source_suffix)
            genfiles = [
                genfile + (ext[0] if not genfile.endswith(tuple(ext)) else "")
                for genfile in genfiles
            ]

            for entry in genfiles[:]:
                if not os.path.isfile(os.path.join(app.srcdir, entry)):
                    self.logger.warning(
                        __(f"automodsumm_generate: file not found: {entry}")
                    )
                    genfiles.remove(entry)

        if not genfiles:
            return

        suffix = get_rst_suffix(app)
        if suffix is None:
            self.logger.warning(
                __(
                    "automodsumm generates .rst files internally. "
                    "But your source_suffix does not contain .rst. Skipped."
                )
            )
            return

        imported_members = app.config.autosummary_imported_members
        with mock(app.config.autosummary_mock_imports):
            self.generate_docs(
                genfiles,
                suffix=suffix,
                base_path=app.srcdir,
                imported_members=imported_members,
                overwrite=app.config.autosummary_generate_overwrite,
                encoding=app.config.source_encoding,
            )

    def generate_docs(
        self,
        source_filenames: List[str],
        output_dir: str = None,
        suffix: str = ".rst",
        base_path: str = None,
        imported_members: bool = False,
        overwrite: bool = True,
        encoding: str = "utf-8",
    ) -> None:
        """
        Generate and write stub files for objects defined in the :rst:dir:`automodapi`
        and :rst:dir:`automodsumm` directives.

        Parameters
        ----------

        source_filenames : List[str]
            A list of all filenames for with the :rst:dir:`automodapi` and
            :rst:dir:`automodsumm` directives will be searched for.

        output_dir : `str`
            Directory for which the stub files will be written to.

        suffix : `str`
            (Default ``".rst"``) Suffix given to the written stub files.

        base_path : `str`
            The common base path for the filenames listed in ``source_filenames``.
            This is typically the source directory of the Sphinx application.

        imported_members : `bool`
            (Default `False`) Set `True` to include imported members in the
            stub file documentation for *module* object types.

        overwrite : `bool`
            (Default `True`)  Will cause existing stub files to be overwritten.

        encoding : `str`
            (Default: ``"utf-8"``) Encoding for the written stub files.


        .. note::  Adapted from
                   :func:`sphinx.ext.autosummary.generate.generate_autosummary_docs`.
        """
        app = self.app

        _info = self.logger.info
        _warn = self.logger.warning

        showed_sources = list(sorted(source_filenames))
        _info(
            __(f"[automodsumm] generating stub files for {len(showed_sources)} sources")
        )

        if output_dir:
            _info(__(f"[automodsumm] writing to {output_dir}"))

        if base_path is not None:
            source_filenames = [
                os.path.join(base_path, filename) for filename in source_filenames
            ]

        template = AutomodsummRenderer(app)

        # read
        items = self.find_in_files(source_filenames)

        # keep track of new files
        new_files = []

        if app:
            filename_map = app.config.autosummary_filename_map
        else:
            filename_map = {}

        # write
        for entry in sorted(set(items), key=str):
            if entry.path is None:
                # The corresponding automodsumm:: directive did not have
                # a :toctree: option
                continue

            path = output_dir or os.path.abspath(entry.path)
            ensuredir(path)

            try:
                name, obj, parent, modname = import_by_name(entry.name)
                qualname = name.replace(modname + ".", "")
            except ImportError as e:
                try:
                    # try to import as an instance attribute
                    name, obj, parent, modname = import_ivar_by_name(entry.name)
                    qualname = name.replace(modname + ".", "")
                except ImportError:
                    _warn(__(f"[automodsumm] failed to import {entry.name}: {e}"))
                    continue

            context = {}
            if app:
                context.update(app.config.autosummary_context)

            content = generate_autosummary_content(
                name,
                obj,
                parent,
                template,
                entry.template,
                imported_members,
                app,
                entry.recursive,
                context,
                modname,
                qualname,
            )

            filename = os.path.join(path, filename_map.get(name, name) + suffix)
            if os.path.isfile(filename):
                with open(filename, encoding=encoding) as f:
                    old_content = f.read()

                if content == old_content:
                    continue
                elif overwrite:  # content has changed
                    with open(filename, "w", encoding=encoding) as f:
                        f.write(content)
                    new_files.append(filename)
            else:
                with open(filename, "w", encoding=encoding) as f:
                    f.write(content)
                new_files.append(filename)

        # descend recursively to new files
        if new_files:
            self.generate_docs(
                new_files,
                output_dir=output_dir,
                suffix=suffix,
                base_path=base_path,
                imported_members=imported_members,
                overwrite=overwrite,
            )

    def find_in_files(self, filenames: List[str]) -> List[AutomodsummEntry]:
        """
        Search files for the :rst:dir:`automodapi` and :rst:dir:`automodsumm`
        directives and generate a list of
        `~plasmapy_sphinx.automodsumm.generate.AutomodsummEntry`'s indicating which stub
        files need to be generated.

        Parameters
        ----------
        filenames : List[str]
            List of filenames to be searched.


        .. note:: Adapted from
                  :func:`sphinx.ext.autosummary.generate.find_autosummary_in_files`.
        """
        documented = []  # type: List[AutomodsummEntry]
        for filename in filenames:
            with open(filename, encoding="utf-8", errors="ignore") as f:
                lines = f.read().splitlines()
                documented.extend(self.find_in_lines(lines, filename=filename))
        return documented

    def find_in_lines(
        self,
        lines: List[str],
        filename: str = None,
    ) -> List[AutomodsummEntry]:
        """
        Search a list of strings for the :rst:dir:`automodapi` and
        :rst:dir:`automodsumm` directives and generate a list of
        `~plasmapy_sphinx.automodsumm.generate.AutomodsummEntry`'s indicating which stub
        files need to be generated.

        Parameters
        ----------
        lines : List[str]
            List of strings to be searched.

        filename : str
            The file from which ``lines`` came from.


        .. note:: Adapted from
                  :func:`sphinx.ext.autosummary.generate.find_autosummary_in_lines`.
        """

        from ..autodoc.automodapi import AutomodapiOptions
        from .core import AutomodsummOptions

        documented = []  # type: List[AutomodsummEntry]

        current_module = None
        modname = ""

        options = {}  # type: Dict[str, Any]

        _option_cls = None

        in_automodapi_directive = False
        gather_objs = False

        last_line = False
        nlines = len(lines)

        for ii, line in enumerate(lines):
            if ii == nlines - 1:
                last_line = True

            # looking for option `   :option: option_args`
            if in_automodapi_directive:
                match = self._re["option"].search(line)
                if match is not None:
                    option_name = match.group(2)
                    option_args = match.group(3)
                    try:
                        option_args = _option_cls.option_spec[option_name](option_args)
                        options[option_name] = option_args
                    except (KeyError, TypeError):
                        pass
                else:
                    # done reading options
                    in_automodapi_directive = False
                    gather_objs = True

                if last_line:
                    # end of lines reached
                    in_automodapi_directive = False
                    gather_objs = True

                if in_automodapi_directive:
                    continue

            # looking for `.. automodsumm:: <modname>`
            match = self._re["automodsumm"].search(line)
            if match is not None:
                in_automodapi_directive = True
                # base_indent = match.group(1)
                modname = match.group(2)

                if current_module is None or modname == current_module:
                    pass
                elif not modname.startswith(f"{current_module}."):
                    modname = f"{current_module}.{modname}"
                _option_cls = AutomodsummOptions
                self.logger.info(f"[automodsumm] {modname}")

                if last_line:
                    # end of lines reached
                    in_automodapi_directive = False
                    gather_objs = True
                else:
                    continue

            # looking for `.. automodapi:: <modname>`
            match = self._re["automodapi"].search(line)
            if match is not None:
                in_automodapi_directive = True
                # base_indent = match.group(1)
                modname = match.group(2)

                if current_module is None or modname == current_module:
                    pass
                elif not modname.startswith(f"{current_module}."):
                    modname = f"{current_module}.{modname}"
                _option_cls = AutomodapiOptions

                if last_line:
                    # end of lines reached
                    in_automodapi_directive = False
                    gather_objs = True
                else:
                    continue

            # looking for `.. py:currentmodule:: <current_module>`
            match = self._re["currentmodule"].search(line)
            if match is not None:
                current_module = match.group(3)
                continue

            # gather objects and update documented list
            if gather_objs:
                process_options = _option_cls(
                    self.app,
                    modname,
                    options,
                    docname=filename,
                    _warn=self.logger.warning,
                )
                options = {
                    "toctree": process_options.toctree["abspath"],
                    "template": process_options.options.get("template", None),
                    "recursive": process_options.options.get("recursive", False),
                }

                exclude_modules = (
                    not self.app.config.automodapi_generate_module_stub_files
                )
                obj_list = process_options.generate_obj_list(
                    exclude_modules=exclude_modules
                )

                for name in obj_list:
                    documented.append(
                        AutomodsummEntry(
                            name=name,
                            path=options["toctree"],
                            template=options["template"],
                            recursive=options["recursive"],
                        )
                    )

                self.logger.info(
                    f"[automodsumm stub file gen] collected {len(obj_list):4d} "
                    f"object(s) in '{modname}'"
                )

                # reset for next search
                options = {}
                gather_objs = False
                _option_cls = None

        return documented

    @staticmethod
    def event_handler__autodoc_skip_member(
        app: "Sphinx", what: str, name: str, obj: Any, skip: bool, options: dict
    ):  # noqa
        """
        Event handler for the Sphinx event :event:`autodoc-skip-member`.  This
        handler ensures the ``__call__`` method is documented if defined by the
        associated class.
        """
        if what != "method":
            return

        if name == "__call__":
            return False
        return
