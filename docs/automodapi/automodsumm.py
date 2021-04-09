import os
import re

from sphinx.application import Sphinx
from sphinx.ext.autodoc.mock import mock
from sphinx.ext.autosummary import (
    Autosummary,
    get_rst_suffix,
    import_by_name,
    import_ivar_by_name,
)
from sphinx.ext.autosummary.generate import (
    AutosummaryRenderer,
    AutosummaryEntry,
    generate_autosummary_content,
)
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.osutil import ensuredir
from typing import Any, Callable, Dict, List

from .utils import find_mod_objs, automod_groupings


logger = logging.getLogger(__name__)


def option_str_list(argument):
    if argument is None:
        raise ValueError("argument required but none supplied")
    else:
        return [s.strip() for s in argument.split(",")]


class AutomodsummOptions:
    option_spec = {
        **Autosummary.option_spec,
        "groups": option_str_list,
        "exclude-groups": option_str_list,
        "skip": option_str_list,
    }
    logger = logger

    def __init__(self, app: Sphinx, modname: str, options: Dict[str, Any], _warn: Callable = None):
        self.app = app
        self.modname = modname
        self.options = options
        self.warn = _warn if _warn is not None else self.logger.warning

    @property
    def mod_objs(self):
        return find_mod_objs(self.modname, app=self.app)

    @property
    def groupings(self):
        dgroups, cgroups = automod_groupings(self.app)

        return dgroups | cgroups

    @property
    def custom_groups_defs(self):
        return self.app.config.automod_custom_groups

    def generate_obj_list(self) -> List[str]:
        try:
            mod_objs = self.mod_objs
        except ImportError:
            mod_objs = {}
            self.warn(f"Could not import module {self.modname}")

        # define groupings to include
        allowed_args = self.groupings | {"all"}
        do_groups = self.groupings
        if "groups" in self.options:
            opt_args = set(self.options["groups"])

            unknown_args = opt_args - allowed_args
            if len(unknown_args) > 0:
                self.warn(
                    f"Option 'groups' has unrecognized arguments "
                    f"{unknown_args}. Ignoring."
                )
                opt_args = opt_args - unknown_args

            if "all" not in opt_args:
                do_groups = opt_args

        # exclude groupings
        if "exclude-groups" not in self.options:
            if "groups" not in self.options:
                opt_args = {"modules"}
            else:
                opt_args = set()
        else:
            opt_args = set(self.options["exclude-groups"])

        unknown_args = opt_args - allowed_args
        if len(unknown_args) > 0:
            self.warn(
                f"Option 'exclude-groups' has unrecognized arguments "
                f"{unknown_args}. Ignoring."
            )
            opt_args = opt_args - unknown_args
        elif "all" in opt_args:
            self.warn(
                f"Arguments of 'groups' and 'exclude-groups' results in "
                f"no content."
            )
            content = []
            return content

        do_groups = do_groups - opt_args

        # objects to skip
        skip_names = set()
        if "skip" in self.options:
            skip_names = set(self.options["skip"])

        # gather content
        content = []
        for group in do_groups:
            try:
                for name, qualname in zip(
                        mod_objs[group]["names"], mod_objs[group]["qualnames"]
                ):
                    if not (name in skip_names or qualname in skip_names):
                        content.append(qualname)
            except KeyError:
                pass

        return sorted(content)


class Automodsumm(Autosummary):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    has_content = False
    option_spec = AutomodsummOptions.option_spec.copy()

    logger = logger

    def run(self):
        # env = self.state.document.settings.env
        app = self.env.app
        env = self.env
        modname = self.arguments[0]

        # for some reason, even though ``currentmodule`` is substituted in,
        # sphinx doesn't necessarily recognize this fact.  So we just force
        # it internally, and that seems to fix things
        env.temp_data['py:module'] = modname
        env.ref_context['py:module'] = modname

        nodelist = []

        self.logger.info(f"[automodsumm] Collecting content for '{modname}'.")
        option_processor = AutomodsummOptions(
            app, modname, self.options, _warn=self.warn,
        )
        content = option_processor.generate_obj_list()
        for ii, modname in enumerate(content):
            if not modname.startswith("~"):
                content[ii] = "~" + modname
        self.content = content

        nodelist.extend(Autosummary.run(self))
        return nodelist

    def get_items(self, names):
        try:
            self.bridge.genopt['imported-members'] = True
        except AttributeError:  # Sphinx < 4.0
            self.genopt['imported-members'] = True
        return Autosummary.get_items(self, names)


class GenDocsFromAutomodsumm:
    """Needed so stub file are automatically generated."""
    option_spec = AutomodsummOptions.option_spec.copy()

    _re = {
        "automodsumm": re.compile(
            r"^\n?(\s*)\.\.\s+automodsumm2::\s*(\S+)\s*(?:\n|$)"
        ),
        "automodapi": re.compile(r"^\n?(\s*)\.\.\s+automodapi::\s*(\S+)\s*(?:\n|$)"),
        "option": re.compile(r"^\n?(\s+):(\S*):\s*(\S.*|)\s*(?:\n|$)"),
        "currentmodule": re.compile(
            r"^\s*\.\.\s+(|\S+:)(current)?module::\s*([a-zA-Z0-9_.]+)\s*$"
        ),
    }

    app = None  # type: Sphinx
    logger = logger

    def __call__(self, app: Sphinx):
        """
        This routine is adapted from
        :func:`sphinx.ext.autosummary.process_generate_options` to scan through
        the source files, check for the `automodsumm` directive, and auto
        generate any associated stub files.
        """
        self.app = app
        genfiles = app.config.autosummary_generate

        if genfiles is True:
            env = app.builder.env
            genfiles = [env.doc2path(x, base=None) for x in env.found_docs
                        if os.path.isfile(env.doc2path(x))]
        elif genfiles is False:
            pass
        else:
            ext = list(app.config.source_suffix)
            genfiles = [genfile + (ext[0] if not genfile.endswith(tuple(ext)) else '')
                        for genfile in genfiles]

            for entry in genfiles[:]:
                if not os.path.isfile(os.path.join(app.srcdir, entry)):
                    self.logger.warning(
                        __('automodsumm_generate: file not found: %s'),
                        entry,
                    )
                    genfiles.remove(entry)

        if not genfiles:
            return

        suffix = get_rst_suffix(app)
        if suffix is None:
            self.logger.warning(
                __('automodsumm generates .rst files internally. '
                   'But your source_suffix does not contain .rst. Skipped.')
            )
            return

        # from sphinx.ext.autosummary.generate import generate_autosummary_docs

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

    def find_mod_objs(self, modname: str):
        return find_mod_objs(modname, app=self.app)

    def generate_docs(
            self,
            source_filenames: List[str],
            output_dir: str = None,
            suffix: str = '.rst',
            base_path: str = None,
            imported_members: bool = False,
            overwrite: bool = True,
            encoding: str = 'utf-8'
    ) -> None:
        """
        This code was adapted from
        :func:`sphinx.ext.autosummary.generate.generate_autosummary_docs`.
        """
        app = self.app

        _info = self.logger.info
        _warn = self.logger.warning

        showed_sources = list(sorted(source_filenames))
        if len(showed_sources) > 20:
            showed_sources = showed_sources[:10] + ['...'] + showed_sources[-10:]
        _info(__('[automodsumm] generating autosummary for: %s') %
              ', '.join(showed_sources))

        if output_dir:
            _info(__('[automodsumm] writing to %s') % output_dir)

        if base_path is not None:
            source_filenames = [
                os.path.join(base_path, filename) for filename in source_filenames
            ]

        template = AutosummaryRenderer(app)

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
                    _warn(__('[automodsumm] failed to import %r: %s') % (entry.name, e))
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
                    with open(filename, 'w', encoding=encoding) as f:
                        f.write(content)
                    new_files.append(filename)
            else:
                with open(filename, 'w', encoding=encoding) as f:
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

    def find_in_files(self, filenames: List[str]) -> List[AutosummaryEntry]:
        """
        Adapted from :func:`sphinx.ext.autosummary.generate.find_autosummary_in_files`.

        Find out what items are documented in source/*.rst.
        """
        documented = []  # type: List[AutosummaryEntry]
        for filename in filenames:
            with open(filename, encoding='utf-8', errors='ignore') as f:
                lines = f.read().splitlines()
                documented.extend(self.find_in_lines(lines, filename=filename))
        return documented

    def find_in_lines(
            self, lines: List[str], module: str = None, filename: str = None,
    ) -> List[AutosummaryEntry]:
        """
        Adapted from :func:`sphinx.ext.autosummary.generate.find_autosummary_in_lines`.

        Find out what items appear in automodsumm:: directives in the given lines.
        """
        documented = []  # type: List[AutosummaryEntry]

        current_module = None
        modname = ""

        base_indent = ""
        default_options = {
            "recursive": False,
            "toctree": None,
            "template": None,
        }  # type: Dict[str, Any]
        options = default_options.copy()

        in_automodsumm = False
        gather_objs = False

        for line in lines:
            # looking for option `   :option: option_args`
            if in_automodsumm:
                match = self._re["option"].search(line)
                if match is not None:
                    option_name = match.group(2)
                    option_args = match.group(3)
                    try:
                        option_args = \
                            AutomodsummOptions.option_spec[option_name](option_args)
                    except KeyError:
                        pass
                    options[option_name] = option_args
                else:
                    # done reading options
                    in_automodsumm = False
                    gather_objs = True

            # looking for `.. automodsumm:: <modname>`
            match = self._re["automodsumm"].search(line)
            if match is not None:
                in_automodsumm = True
                base_indent = match.group(1)
                modname = match.group(2)
                if current_module is not None \
                        and not modname.startswith(f"{current_module}."):
                    modname = f"{current_module}.{modname}"
                continue

            # looking for `.. py:currentmodule:: <current_module>`
            match = self._re["currentmodule"].search(line)
            if match is not None:
                current_module = match.group(3)
                continue

            # gather objects and update documented list
            if gather_objs:
                process_options = AutomodsummOptions(
                    self.app, modname, options, _warn=self.logger.warning,
                )
                self.logger.info(f"Processing {modname} with options {options}")
                obj_list = process_options.generate_obj_list()
                for name in obj_list:
                    documented.append(
                        AutosummaryEntry(
                            name,
                            options["toctree"],
                            options["template"],
                            options["recursive"],
                        )
                    )

                options = default_options.copy()
                gather_objs = False

        return documented


def setup_autosummary(app: Sphinx):
    import sphinx_automodapi

    app.setup_extension("sphinx.ext.autosummary")

    app.config.templates_path.append(
        os.path.join(sphinx_automodapi.__path__[0], "templates")
    )
    # app.setup_extension(autodoc_enhancements.__name__)
    # app.setup_extension("sphinx.ext.inheritance_diagram")

    # app.add_directive("automod-diagram", Automoddiagram)
    app.add_directive("automodsumm2", Automodsumm)

    gendocs_from_automodsumm = GenDocsFromAutomodsumm()
    app.connect("builder-inited", gendocs_from_automodsumm)

    # app.add_config_value("automodsumm_writereprocessed", False, True)
    # app.add_config_value("automodsumm_inherited_members", False, "env")

    automod_custom_groups = {
        "aliases": {
            "title": "Aliases",
            "descr": "",
            "dunder": "__aliases__",
            "option_name": "aliases",
        },
        "lite-functions": {
            "title": "Lite Functions",
            "descr": "",
            "dunder": "__flites__",
            "option_name": "lite-funcs",
        }
    }
    app.add_config_value("automod_custom_groups", automod_custom_groups, True)

    return {"parallel_read_safe": True, "parallel_write_safe": True}
