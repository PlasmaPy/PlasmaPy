import sys

from docutils.statemachine import StringList
from plasmapy.utils.decorators.helpers import LiteFuncTuple
from sphinx.application import Sphinx
from sphinx.ext.autodoc import DataDocumenter, FunctionDocumenter
from sphinx.locale import __
from sphinx.util import logging
from typing import Any, Optional

if sys.version_info >= (3, 0):
    text_type = str
else:
    text_type = unicode  # noqa

logger = logging.getLogger(__name__)


class LiteFuncDocumenter(FunctionDocumenter):
    objtype = "litefunc"
    directivetype = "function"
    logger = logger

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
            "   :nosignatures:",
            "",
        ],
        "autosummary_item": "   ~{item}",
        "autodata_heading": [
            ".. rubric:: Bound Attributes",
            "",
        ],
        "autodata": ".. autodata:: {item}",
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
    def __has_litefunc__(self) -> LiteFuncTuple:
        return self.object.__has_litefunc__

    def generate_more_content(self):
        if not self.has_litefunc or self.is_alias:
            return []

        bound_attrs = list(self.__has_litefunc__.bound_attrs)
        if len(bound_attrs) == 0:
            return []

        bound_attrs = sorted(bound_attrs, key=str.casefold)
        if "lite" in bound_attrs:
            # ensure lite is always first
            bound_attrs.remove("lite")
            bound_attrs.insert(0, "lite")

        lines = []

        shortnme = self.name.split('.')[-1]
        qualname = f"{self.modname}.{shortnme}"

        lines.extend(self._templates["note"].format(func_name=qualname).splitlines())
        lines.extend(self._templates["autosummary_setup"])

        for bname in bound_attrs:
            lines.append(
                self._templates["autosummary_item"].format(item=f"{qualname}.{bname}")
            )

        lines.append("")
        lines.extend(self._templates["autodata_heading"])

        for bname in bound_attrs:
            lines.append(
                self._templates["autodata"].format(item=f"{qualname}.{bname}")
            )
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
    app.add_autodocumenter(LiteFuncDocumenter)
