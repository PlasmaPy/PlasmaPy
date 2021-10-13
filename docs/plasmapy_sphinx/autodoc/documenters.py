import inspect
import sys

from docutils.parsers.rst import directives
from docutils.statemachine import StringList
from plasmapy.utils.decorators.helpers import LiteFuncTupleEntry, _litefunc_registry
from sphinx.application import Sphinx
from sphinx.ext.autodoc import (
    AttributeDocumenter,
    DataDocumenter,
    FunctionDocumenter,
    ModuleLevelDocumenter,
)
from sphinx.locale import __
from sphinx.pycode import ModuleAnalyzer, PycodeError
from sphinx.util import logging
from typing import Any, List, Optional, Tuple

if sys.version_info >= (3, 0):
    text_type = str
else:
    text_type = unicode  # noqa

logger = logging.getLogger(__name__)

# TODO: drop the named tuple entries (LiteFuncTupleEntry) for __has_litefunc__


class LiteDataDocumenter(DataDocumenter):
    # TODO: remove origin option in favor of the registry
    # TODO: create and link configuration value for registry

    objtype = "data-lite"
    directivetype = "data"
    priority = AttributeDocumenter.priority + 10
    option_spec = {
        **ModuleLevelDocumenter.option_spec,
        "origin": directives.unchanged,
    }
    logger = logger

    @classmethod
    def can_document_member(
        cls,
        member: Any,
        membername: str,
        isattr: bool,
        parent: Any,
    ) -> bool:
        return isinstance(parent, LiteFuncDocumenter)

    def get_signature(self):
        if callable(self.object):
            sig = inspect.signature(self.object).__str__()
            self.options["annotation"] = sig
        else:
            self.options["annotation"] = ""

    def set_option_origin(self):
        name = self.name.replace("::", ".")
        try:
            orign = _litefunc_registry[name]["origin"]
        except KeyError:
            orign = None

        self.options["origin"] = orign

    def format_signature(self, **kwargs: Any) -> str:

        if self.args is None:
            if callable(self.object):
                sig = inspect.signature(self.object).__str__()
                self.args = sig.strip("()")

        self.logger.info(f" -----    args (after) = {self.args}")

        return super().format_signature(**kwargs)

    def generate_more_content(self):
        lines = []

        origin = self.options.get("origin", None)
        # qualname = f"{self.modname}.{self.format_name()}"

        if origin is not None:
            lines.append(
                f"**Bound copy of** `~{origin}`."
            )
            lines.append("")

        # return lines
        return []

    def get_module_comment(self, attrname: str) -> Optional[List[str]]:
        modname = self.modname

        self.set_option_origin()

        origin = self.options["origin"]

        if origin:
            osplit = origin.split(".")
            attrname = osplit[-1]
            modname = ".".join(osplit[:-1])

        try:
            analyzer = ModuleAnalyzer.for_module(modname)
            analyzer.analyze()
            key = ('', attrname)
            if key in analyzer.attr_docs:
                return list(analyzer.attr_docs[key])
        except PycodeError:
            pass

        return None

    def get_doc(self, ignore: int = None) -> List[List[str]]:
        docstrings = super().get_doc(ignore=ignore)

        origin = self.options.get("origin", None)
        if origin is not None:
            docstrings.insert(0, [f"**Bound copy of** `~{origin}`.", ])

        return docstrings

    def should_suppress_directive_header(self) -> bool:
        try:
            sig = self.format_signature()
            if sig != "":
                return True
        except Exception:
            pass

        return super().should_suppress_directive_header()

    def resolve_name(
        self,
        modname: str,
        parents: Any,
        path: str,
        base: Any,
    ) -> Tuple[str, List[str]]:

        # self.logger.info(
        #     f"~~~@@@ \n"
        #     f"   modname = {modname}\n"
        #     f"   parents = {parents}\n"
        #     f"      path = {path}\n"
        #     f"      base = {base}")
        if modname is None:
            if path:
                modname = path.rstrip('.')
            else:
                # if documenting a toplevel object without explicit module,
                # it can be contained in another auto directive ...
                modname = self.env.temp_data.get('autodoc:module')
                # ... or in the scope of a module directive
                if not modname:
                    modname = self.env.ref_context.get('py:module')
                # ... else, it stays None, which means invalid

            modname, sep, lfunc = modname.rpartition(".")
            parents = [lfunc]

        return modname, parents + [base]

    def generate(
        self, more_content: Optional[StringList] = None,
        real_modname: str = None,
        check_module: bool = False,
        all_members: bool = False,
    ) -> None:
        docname = self.env.docname
        if not docname.endswith(".rst"):
            docname += ".rst"

        if not self.parse_name():
            # need a module to import
            logger.warning(
                __('don\'t know which module to import for autodocumenting '
                   '%r (try placing a "module" or "currentmodule" directive '
                   'in the document, or giving an explicit module name)') %
                self.name, type='autodoc')
            return

        # now, import the module and get object to document
        if not self.import_object():
            return

        # self.get_signature()

        # Generate the 'more_content' (automodsumm and inheritance diagrams)
        more_content = StringList()
        more_lines = self.generate_more_content()
        self.logger.info("\n".join(more_lines))
        for line in more_lines:
            more_content.append(line, source=docname)

        super().generate(
            more_content=more_content,
            real_modname=real_modname,
            check_module=check_module,
            all_members=all_members,
        )

        # lines = '\n'.join(self.directive.result)
        # self.logger.info(f"--- {self.objpath} ---")
        # self.logger.info(f"{lines}")
        # self.logger.info(f"  MODNAME = {self.modname}")
        # self.logger.info(f"  OBJPATH = {self.objpath}")


class LiteFuncDocumenter(FunctionDocumenter):
    objtype = "function-lite"
    directivetype = "function"
    logger = logger

    priority = 10 + FunctionDocumenter.priority

    _templates = {
        "note": "\n".join(
            [
                ".. note:: To increase usability of `~{func_name}` several "
                "attributes/functions are manually bound to this function.",
                "",
                "",
            ]
        ),
        "autosummary_setup": [
            ".. rubric:: Summary of Bound Attributes",
            "",
            ".. autosummary::",
            "",
        ],
        "autosummary_item": "   ~{item}",
        "autodata_heading": [
            ".. rubric:: Bound Attributes",
            "",
        ],
        "autodata": ".. autodata-lite:: {item}",
    }

    @classmethod
    def can_document_member(
        cls, member: Any, membername: str, isattr: bool, parent: Any,
    ) -> bool:
        rtn = super().can_document_member(member, membername, isattr, parent)

        has_litefunc = hasattr(member, "__has_litefunc__")
        if not has_litefunc:
            return False

        return rtn and has_litefunc

    @property
    def has_litefunc(self):
        return hasattr(self.object, "__has_litefunc__")

    @property
    def is_alias(self):
        if hasattr(self.module, "__aliases__"):
            return self.name.split(".")[-1] in self.module.__aliases__
        else:
            return False

    @property
    def __has_litefunc__(self) -> LiteFuncTupleEntry:
        return self.object.__has_litefunc__

    def generate_more_content(self):
        if not self.has_litefunc or self.is_alias:
            return []

        if len(self.__has_litefunc__) == 0:
            return []

        bound_names = []
        origins = []
        for item in self.__has_litefunc__:
            bound_names.append(item.name)
            origins.append(item.origin)
        __has_litefunc__ = sorted(
            zip(bound_names, origins), key=lambda x: str.casefold(x[0])
        )
        bound_names = []
        origins = []
        for bname, origin in __has_litefunc__:
            if bname == "lite":
                bound_names.insert(0, bname)
                origins.insert(0, origin)
                continue

            bound_names.append(bname)
            origins.append(origin)

        lines = []

        shortnme = self.name.split('.')[-1]
        qualname = f"{self.modname}.{shortnme}"

        lines.extend(self._templates["note"].format(func_name=qualname).splitlines())
        lines.extend(self._templates["autosummary_setup"])

        for bname in bound_names:
            lines.append(
                self._templates["autosummary_item"].format(item=f"{qualname}.{bname}")
            )

        lines.append("")

        if "noindex" in self.options and self.options["noindex"]:
            return lines

        lines.extend(self._templates["autodata_heading"])

        for bname, origin in zip(bound_names, origins):
            self.logger.info(f"--- {bname}, {origin}")
            lines.append(
                self._templates["autodata"].format(item=f"{qualname}.{bname}")
            )
            lines.append(f"   :origin: {origin}")
            lines.append("")

        return lines

    def generate(
        self,
        more_content: Optional[StringList] = None,
        real_modname: str = None,
        check_module: bool = False,
        all_members: bool = False,
    ) -> None:
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

        self.logger.info(
            f"[autolitefunc] Documenting \n"
            f"      modename = {self.modname}\n"
            f"       objpath = {self.objpath}\n"
            f"       objtype = {self.objtype}\n"
            f"          name = {self.name}\n"
            f"  real_modname = {real_modname}\n"
            f"        object = {self.object}\n"
            f"    lite func? = {self.has_litefunc}\n"
            f"         alias = {self.is_alias}\n"
            f"        module = {self.module}."
        )

        # Generate the 'more_content' (automodsumm and inheritance diagrams)
        more_content = StringList()
        more_lines = self.generate_more_content()
        self.logger.info("\n".join(more_lines))
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
    app.add_autodocumenter(LiteDataDocumenter)
    app.add_autodocumenter(LiteFuncDocumenter)
