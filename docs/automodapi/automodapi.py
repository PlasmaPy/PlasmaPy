import re
import sys
import warnings

from collections import OrderedDict
from docutils.parsers.rst import directives
from docutils.statemachine import StringList
from sphinx.application import Sphinx

try:
    from sphinx.deprecation import RemovedInSphinx50Warning
except ImportError:

    class RemovedInSphinx50Warning(PendingDeprecationWarning):
        pass


from sphinx.ext.autodoc import bool_option, ModuleDocumenter
from sphinx.locale import __
from sphinx.util import logging
from typing import Any, Callable, Dict, Optional, Union

from .automodsumm import AutomodsummOptions, option_str_list
from .utils import default_grouping_info

if sys.version_info >= (3, 0):
    text_type = str
else:
    text_type = unicode  # noqa

logger = logging.getLogger(__name__)


_option_spec = option_spec = {
    **ModuleDocumenter.option_spec,
    "groups": option_str_list,
    "exclude-groups": option_str_list,
    "no-groups": bool_option,
    # "group-order": option_str_list,
    # "merge-groups": bool_option,
    "skip": option_str_list,
    "include-heading": bool_option,
    "members": lambda x: None,
    "heading-chars": directives.unchanged,
    "toctree": directives.unchanged,
    "no-toctree": bool_option,
    "no-main-docstr": bool_option,
    "inheritance-diagram": bool_option,
    "no-inheritance-diagram": bool_option,
}  # type: Dict[str, Callable]


class AutomodapiOptions(AutomodsummOptions):
    option_spec = _option_spec
    logger = logger

    def __init__(
        self,
        app: Sphinx,
        modname: str,
        options: Dict[str, Any],
        docname: str = None,
        _warn: Callable = None,
    ):
        super().__init__(app, modname, options, docname=docname, _warn=_warn)

    def condition_options(self):
        super().condition_options()
        self.condition_heading_chars_option()
        self.condition_include_heading_option()
        self.condition_inheritance_diagram_option()

    def condition_toctree_option(self):
        if "no-toctree" in self.options and self.options["no-toctree"]:
            if "toctree" in self.options:
                del self.options["toctree"]
        elif "toctree" not in self.options:
            self.options["toctree"] = self.app.config.automodapi_toctreedirnm

        super().condition_toctree_option()

    def condition_heading_chars_option(self):
        non_alphanumerics = re.compile("[^0-9a-zA-Z]]+")
        heading_chars = self.options.get("heading-chars", None)
        if (
            heading_chars is None
            or heading_chars == ""
            or len(heading_chars) < 2
            or non_alphanumerics.fullmatch(heading_chars) is None
        ):
            self.options["heading-chars"] = "-^"

    def condition_include_heading_option(self):
        if "include-heading" not in self.options:
            self.options["include-heading"] = False

    def condition_inheritance_diagram_option(self):
        if "no-inheritance-diagram" in self.options:
            self.options["inheritance-diagram"] = False
            del self.options["no-inheritance-diagram"]
        elif "inheritance-diagram" in self.options:
            self.options["inheritance-diagram"] = False
        else:
            self.options[
                "inheritance-diagram"
            ] = self.app.config.automodapi_inheritance_diagram

    def condition_group_options(self):
        if "no-groups" in self.options and self.options["no-groups"]:
            self.options["groups"] = []
            if "exclude-groups" in self.options:
                del self.options["exclude-groups"]

            return

        super().condition_group_options()

    @property
    def options_for_automodsumm(self):
        options = {}

        asumm_opts = list(AutomodsummOptions.option_spec)
        for name in asumm_opts:
            if name in self.options and name not in ("groups", "exclude-groups"):
                val = self.options[name]
                if isinstance(val, list):
                    val = ", ".join(val)

                if name == "toctree" and self.toctree["original"] is not None:
                    # automodsumm requires path relative to the confdir but
                    # Automodsumm.review_toctree_option generates the path relative
                    # to the file location
                    options[name] = self.toctree["original"]
                else:
                    options[name] = val

        return options


class ModAPIDocumenter(ModuleDocumenter):
    objtype = "modapi"
    directivetype = "module"
    titles_allowed = True
    content_indent = ""
    logger = logger
    option_spec = _option_spec

    _automod_option_spec_names = set(ModuleDocumenter.option_spec)

    _grouping_info = default_grouping_info

    _templates = {
        "mod-heading": "\n".join(["{modname} {pkg_or_mod}", "{underline}", ""]),
        "heading": "\n".join(["{title}", "{underline}"]),
        "automodsumm": "\n".join([".. automodsumm:: {modname}", "   :groups: {group}"]),
        "options": "   :{option}: {opt_args}",
        "inheritance-diagram": "\n".join(
            [
                ".. inheritance-diagram:: {cls_list}",
                "   :private-bases:",
                "   :parts: 1",
            ],
        ),
    }

    @property
    def grouping_info(self) -> Dict[str, Dict[str, Union[str, None]]]:

        group_order = tuple(self.env.app.config.automodapi_group_order)
        custom_group_info = self.env.app.config.automod_custom_groups

        group_names = set(self._grouping_info) | set(custom_group_info)
        remainder = list(group_names - set(group_order))
        if len(remainder) > 0:
            group_order += tuple(sorted(remainder))

        _grouping_info = OrderedDict()
        for name in group_order:
            if name in self._grouping_info:
                _info = self._grouping_info[name]
            else:
                _info = custom_group_info[name]

            _grouping_info.update(
                {
                    name: {
                        "title": _info.get("title", None),
                        "descr": _info.get("descr", None),
                        "dunder": _info.get("dunder", None),
                    },
                },
            )

        return _grouping_info

    def generate_more_content(self, modname):
        app = self.env.app
        inheritance_groups = app.config.automodapi_groups_with_inheritance_diagrams

        lines = []

        option_processor = AutomodapiOptions(
            app,
            modname,
            self.options,
            _warn=self.logger.warning,
            docname=self.env.docname,
        )
        asumm_options = option_processor.options_for_automodsumm
        mod_objs = option_processor.mod_objs_option_filtered
        heading_char = (
            option_processor.options["heading-chars"][1]
            if option_processor.options["include-heading"]
            else option_processor.options["heading-chars"][0]
        )
        include_inheritance_diagram = option_processor.options["inheritance-diagram"]

        # scan thru default groups first
        for group, info in self.grouping_info.items():
            if group not in mod_objs:
                continue

            title = info["title"]

            underline = heading_char * len(title)

            # sub-heading
            lines.extend(
                self._templates["heading"]
                .format(title=title, underline=underline)
                .splitlines()
            )
            lines.append("")

            # add automodsumm directive
            lines.extend(
                self._templates["automodsumm"]
                .format(modname=modname, group=group)
                .splitlines()
            )

            # add options for automodsumm directive
            for name, val in asumm_options.items():
                if name == "toctree" and group == "modules":
                    continue

                lines.extend(
                    self._templates["options"]
                    .format(option=name, opt_args=val)
                    .splitlines()
                )
            lines.append("")

            # add inheritance-diagram
            if group in inheritance_groups and include_inheritance_diagram:
                cls_list = " ".join(mod_objs[group]["qualnames"])
                lines.extend(
                    self._templates["inheritance-diagram"]
                    .format(cls_list=cls_list)
                    .splitlines()
                )
                lines.append("")

        return lines

    def generate_heading(self, modname):
        app = self.env.app
        sourcename = self.get_sourcename()

        option_processor = AutomodapiOptions(
            app,
            modname,
            self.options,
            _warn=self.logger.warning,
            docname=self.env.docname,
        )
        if not option_processor.options["include-heading"]:
            return

        modname = re.escape(modname)

        if option_processor.pkg_or_module == "pkg":
            pkg_or_mod = "Package"
        else:
            pkg_or_mod = "Module"

        heading_char = option_processor.options["heading-chars"][0]
        underline = heading_char * (len(modname) + 1 + len(pkg_or_mod))

        heading_lines = (
            self._templates["mod-heading"]
            .format(modname=modname, pkg_or_mod=pkg_or_mod, underline=underline)
            .splitlines()
        )

        for line in heading_lines:
            self.add_line(line, source=sourcename)

    def add_content(
        self, more_content: Optional[StringList], no_docstring: bool = False
    ) -> None:
        """Add content from docstrings, attribute documentation and user."""
        if no_docstring:
            warnings.warn(
                "The 'no_docstring' argument to %s.add_content() is deprecated."
                % self.__class__.__name__,
                RemovedInSphinx50Warning,
                stacklevel=2,
            )

        no_docstring = self.options.get("no-main-docstr", False)

        # set sourcename and add content from attribute documentation
        sourcename = self.get_sourcename()
        if self.analyzer:
            attr_docs = self.analyzer.find_attr_docs()
            if self.objpath:
                key = (".".join(self.objpath[:-1]), self.objpath[-1])
                if key in attr_docs:
                    no_docstring = True
                    # make a copy of docstring for attributes to avoid cache
                    # the change of autodoc-process-docstring event.
                    docstrings = [list(attr_docs[key])]

                    for i, line in enumerate(self.process_doc(docstrings)):
                        self.add_line(line, sourcename, i)

        # add content from docstrings
        if not no_docstring:
            docstrings = self.get_doc()
            if docstrings is None:
                # Do not call autodoc-process-docstring on get_doc() returns None.
                pass
            else:
                if not docstrings:
                    # append at least a dummy docstring, so that the event
                    # autodoc-process-docstring is fired and can add some
                    # content if desired
                    docstrings.append([])
                for i, line in enumerate(self.process_doc(docstrings)):
                    self.add_line(line, sourcename, i)

        # add additional content (e.g. from document), if present
        if more_content:
            for line, src in zip(more_content.data, more_content.items):
                self.add_line(line, src[0], src[1])

    def generate(
        self,
        more_content: Optional[StringList] = None,
        real_modname: str = None,
        check_module: bool = False,
        all_members: bool = False,
    ) -> None:
        """Generate reST for the object given by *self.name*, and possibly for
        its members.

        If *more_content* is given, include that content. If *real_modname* is
        given, use that module name to find attribute docs. If *check_module* is
        True, only generate if the object is defined in the module name it is
        imported from. If *all_members* is True, document all members.
        """

        docname = self.env.docname
        if not docname.endswith(".rst"):
            docname += ".rst"

        if not self.parse_name():
            # need a module to import
            logger.warning(
                __(
                    f"don't know which module to import for autodocumenting "
                    f"{self.name} (try placing a 'module' or 'currentmodule' "
                    f"directive in the document, or giving an explicit module name)"
                ),
                type="autodoc",
            )
            return

        # now, import the module and get object to document
        if not self.import_object():
            return

        # If there is no real module defined, figure out which to use.
        # The real module is used in the module analyzer to look up the module
        # where the attribute documentation would actually be found in.
        # This is used for situations where you have a module that collects the
        # functions and classes of internal submodules.
        real_modname = real_modname or self.get_real_modname()  # type: str

        # Generate heading
        self.generate_heading(modname=real_modname)

        # Generate the 'more_content' (automodsumm and inheritance diagrams)
        more_content = StringList()
        more_lines = self.generate_more_content(modname=real_modname)
        for line in more_lines:
            more_content.append(line, source=docname)

        # generate
        super().generate(
            more_content=more_content,
            real_modname=real_modname,
            check_module=check_module,
            all_members=all_members,
        )


def setup(app: Sphinx):
    from .automodsumm import setup as setup_automodsumm

    rtn = setup_automodsumm(app)

    app.setup_extension("sphinx.ext.inheritance_diagram")

    app.add_autodocumenter(ModAPIDocumenter)

    if not hasattr(app.config, "automodapi_inheritance_diagram"):
        app.add_config_value("automodapi_inheritance_diagram", True, True)
    if not hasattr(app.config, "automodapi_toctreedirnm"):
        app.add_config_value("automodapi_toctreedirnm", "api", True)
    # app.add_config_value("automodapi_writereprocessed", False, True)

    app.add_config_value(
        "automodapi_group_order",
        ("modules", "classes", "exceptions", "warnings", "functions", "variables"),
        True,
    )
    app.add_config_value(
        "automodapi_groups_with_inheritance_diagrams",
        ["classes", "exceptions", "warnings"],
        True,
    )

    return rtn
