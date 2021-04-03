import abc
import inspect
import io
import os
import re
import sys

from docutils.parsers.rst.directives import flag
from sphinx.application import Sphinx
from sphinx.ext.autosummary import Autosummary
from sphinx.util import logging


class AutosumAliases(Autosummary):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    has_content = False
    option_spec = dict(Autosummary.option_spec)

    def run(self):
        env = self.state.document.settings.env
        modname = self.arguments[0]

        nodelist = []

        try:
            alias_names, alias_qualnames, alias_objs = self.find_alias_objs(modname)
        except ImportError:
            alias_names = []
            self.warn(f"Could not import module {modname}")

        try:
            self.content = alias_names

            # for some reason, even though ``currentmodule`` is substituted in,
            # sphinx doesn't necessarily recognize this fact.  So we just force
            # it internally, and that seems to fix things
            env.temp_data['py:module'] = modname
            env.ref_context['py:module'] = modname

            # can't use super because Sphinx/docutils has trouble return
            # super(Autosummary,self).run()
            nodelist.extend(Autosummary.run(self))

            return nodelist
        finally:  # has_content = False for the AutosumAliases
            self.content = []

    def get_items(self, names):
        try:
            self.bridge.genopt['imported-members'] = True
        except AttributeError:  # Sphinx < 4.0
            self.genopt['imported-members'] = True
        return Autosummary.get_items(self, names)

    @staticmethod
    def find_alias_objs(modname: str):
        __import__(modname)
        mod = sys.modules[modname]

        # collect all aliases
        if hasattr(mod, "__aliases__"):
            pkg_items = [
                (obj_name, getattr(mod, obj_name)) for obj_name in mod.__aliases__
            ]
        else:
            pkg_items = []

        # filter out modules
        alias_names = []
        alias_objs = []
        for name, obj in pkg_items:
            if inspect.ismodule(obj):
                continue
            alias_names.append(name)
            alias_objs.append(obj)

        # get fully qualified names
        alias_qualnames = []
        for name, obj in zip(alias_names, alias_objs):
            if hasattr(obj, "__module__") and hasattr(obj, "__name__"):
                alias_qualnames.append(f"{obj.__module__}.{obj.__name__}")
            else:
                alias_qualnames.append(f"{modname}.{name}")

        return alias_names, alias_qualnames, alias_objs


class ProcessAutosumAliases:
    _lineendrex = r'(?:\n|$)'
    _hdrex = r'^\n?(\s*)\.\. autosum_aliases::\s*(\S+)\s*' + _lineendrex
    _oprex1 = r'(?:\1(\s+)\S.*' + _lineendrex + ')'
    _oprex2 = r'(?:\1\4\S.*' + _lineendrex + ')'
    _automodsummrex = re.compile(_hdrex + '(' + _oprex1 + '?' + _oprex2 + '*)',
                                 re.MULTILINE)

    @staticmethod
    def find_alias_objs(modname: str):
        return AutosumAliases.find_alias_objs(modname)

    def convert_to_autosummary_lines(self, filename: str, app: Sphinx):
        logger = logging.getLogger(__name__)
        filepath = os.path.join(app.builder.env.srcdir, filename)

        with io.open(filepath, encoding='utf8') as fr:
            filestr = fr.read()

        spl = self._automodsummrex.split(filestr)
        # 0th entry is the stuff before the first automodsumm line
        indent1s = spl[1::5]
        mods = spl[2::5]
        opssecs = spl[3::5]
        indent2s = spl[4::5]
        remainders = spl[5::5]

        # only grab automodsumm sections and convert them to autosummary with the
        # entries for all the public objects
        newlines = []

        # loop over all autosum_aliases in this document
        for i, (i1, i2, modnm, ops, rem) in enumerate(zip(indent1s, indent2s, mods,
                                                          opssecs, remainders)):
            allindent = i1 + ('    ' if i2 is None else i2)

            # filter out functions-only, classes-only, and ariables-only
            # options if present.
            option_lines = ops.split('\n')
            # toskip = []
            # allowedpkgnms = []
            # funcsonly = clssonly = varsonly = False
            # for i, ln in reversed(list(enumerate(option_lines))):
            #     if ':functions-only:' in ln:
            #         funcsonly = True
            #         del option_lines[i]
            #     if ':classes-only:' in ln:
            #         clssonly = True
            #         del option_lines[i]
            #     if ':variables-only:' in ln:
            #         varsonly = True
            #         del option_lines[i]
            #     if ':skip:' in ln:
            #         toskip.extend(_str_list_converter(ln.replace(':skip:', '')))
            #         del option_lines[i]
            #     if ':allowed-package-names:' in ln:
            #         allowedpkgnms.extend(
            #             _str_list_converter(ln.replace(':allowed-package-names:', '')))
            #         del option_lines[i]
            # if [funcsonly, clssonly, varsonly].count(True) > 1:
            #     msg = ('Defined more than one of functions-only, classes-only, '
            #            'and variables-only.  Skipping this directive.')
            #     lnnum = sum([spl[j].count('\n') for j in range(i * 5 + 1)])
            #     logger.warning('[automodsumm] ' + msg, (fn, lnnum))
            #     continue

            # Use the currentmodule directive so we can just put the local names
            # in the autosummary table.  Note that this doesn't always seem to
            # actually "take" in Sphinx's eyes, so in `Automodsumm.run`, we have to
            # force it internally, as well.
            newlines.extend([i1 + '.. currentmodule:: ' + modnm,
                             '',
                             '.. autosummary::'])
            newlines.extend(option_lines)

            # ols = True if len(allowedpkgnms) == 0 else allowedpkgnms
            for name, qualname, obj in zip(*self.find_alias_objs(modnm)):
                # if nm in toskip:
                #     continue
                # if funcsonly and not inspect.isroutine(obj):
                #     continue
                # if clssonly and not inspect.isclass(obj):
                #     continue
                # if varsonly and (inspect.isclass(obj) or inspect.isroutine(obj)):
                #     continue
                newlines.append(allindent + name)

        # add one newline at the end of the autosummary block
        newlines.append('')

        return newlines

    def generate_docs(
        self,
        lines,
        srcfn,
        app=None,
        suffix=".rst",
        base_path=None,
        builder=None,
        template_dir=None,
    ):
        from sphinx_automodapi.utils import \
            find_autosummary_in_lines_for_automodsumm as find_autosummary_in_lines
        from sphinx_automodapi.utils import cleanup_whitespace
        from sphinx.ext.autosummary import import_by_name, get_documenter
        from sphinx.jinja2glue import BuiltinTemplateLoader
        from sphinx.util.inspect import safe_getattr
        from sphinx.util.osutil import ensuredir
        from jinja2 import FileSystemLoader, TemplateNotFound
        from jinja2.sandbox import SandboxedEnvironment

        logger = logging.getLogger(__name__)

        # Create our own templating environment - here we use Astropy's
        # templates rather than the default autosummary templates, in order to
        # allow docstrings to be shown for methods.
        template_dirs = [
            os.path.join(os.path.dirname(__file__), 'templates'),
            os.path.join(base_path, '_templates'),
        ]
        if builder is not None:
            # allow the user to override the templates
            template_loader = BuiltinTemplateLoader()
            template_loader.init(builder, dirs=template_dirs)
        else:
            if template_dir:
                template_dirs.insert(0, template_dir)
            template_loader = FileSystemLoader(template_dirs)
        template_env = SandboxedEnvironment(loader=template_loader)

        # read
        items = find_autosummary_in_lines(lines, filename=srcfn)
        if len(items) > 0:
            msg = '[autosum_aliases] {1}: found {0} autosum_aliases entries to generate'
            logger.info(msg.format(len(items), srcfn))

        # remove possible duplicates
        items = list(set(items))

        # keep track of new files
        new_files = []

        # write
        for name, path, template_name, inherited_mem in sorted(items):

            if path is None:
                # The corresponding autosummary:: directive did not have
                # a :toctree: option
                continue

            path = os.path.abspath(os.path.join(base_path, path))
            ensuredir(path)

            try:
                import_by_name_values = import_by_name(name)
            except ImportError as e:
                logger.warning('[autosum_aliases] failed to import %r: %s' % (name, e))
                continue

            # if block to accommodate Sphinx's v1.2.2 and v1.2.3 respectively
            if len(import_by_name_values) == 3:
                name, obj, parent = import_by_name_values
            elif len(import_by_name_values) == 4:
                name, obj, parent, module_name = import_by_name_values

            fn = os.path.join(path, name + suffix)

            # skip it if it exists
            if os.path.isfile(fn):
                continue

            new_files.append(fn)

            f = open(fn, 'w')

            try:

                doc = get_documenter(app, obj, parent)

                if template_name is not None:
                    template = template_env.get_template(template_name)
                else:
                    tmplstr = 'autosummary_core/%s.rst'
                    try:
                        template = template_env.get_template(tmplstr % doc.objtype)
                    except TemplateNotFound:
                        template = template_env.get_template(tmplstr % 'base')

                def get_members_mod(obj, typ, include_public=None):
                    """
                    typ = None -> all
                    """
                    if include_public is None:
                        include_public = []
                    items = []
                    for name in dir(obj):
                        try:
                            documenter = get_documenter(app, safe_getattr(obj, name),
                                                        obj)
                        except AttributeError:
                            continue
                        if typ is None or documenter.objtype == typ:
                            items.append(name)
                    public = [x for x in items
                              if x in include_public or not x.startswith('_')]
                    return public, items

                def get_members_class(
                        obj, typ, include_public=None, include_base=False
                ):
                    """
                    typ = None -> all
                    include_base -> include attrs that are from a base class
                    """
                    if include_public is None:
                        include_public = []
                    items = []

                    # using dir gets all of the attributes, including the elements
                    # from the base class, otherwise use __slots__ or __dict__
                    if include_base:
                        names = dir(obj)
                    else:
                        # Classes deriving from an ABC using the `abc` module will
                        # have an empty `__slots__` attribute in Python 3, unless
                        # other slots were declared along the inheritance chain. If
                        # the ABC-derived class has empty slots, we'll use the
                        # class `__dict__` instead.
                        declares_slots = (
                                hasattr(obj, '__slots__') and
                                not (type(obj) is abc.ABCMeta and
                                     len(obj.__slots__) == 0)
                        )

                        if declares_slots:
                            names = tuple(getattr(obj, '__slots__'))
                        else:
                            names = getattr(obj, '__dict__').keys()

                    for name in names:
                        try:
                            documenter = get_documenter(app, safe_getattr(obj, name),
                                                        obj)
                        except AttributeError:
                            continue
                        if typ is None or documenter.objtype == typ:
                            items.append(name)
                        elif typ == 'attribute' and documenter.objtype == 'property':
                            # In Sphinx 2.0 and above, properties have a separate
                            # objtype, but we treat them the same here.
                            items.append(name)
                    public = [x for x in items
                              if x in include_public or not x.startswith('_')]
                    return public, items

                ns = {}

                if doc.objtype == 'module':
                    ns['members'] = get_members_mod(obj, None)
                    ns['functions'], ns['all_functions'] = \
                        get_members_mod(obj, 'function')
                    ns['classes'], ns['all_classes'] = \
                        get_members_mod(obj, 'class')
                    ns['exceptions'], ns['all_exceptions'] = \
                        get_members_mod(obj, 'exception')
                elif doc.objtype == 'class':
                    if inherited_mem is not None:
                        # option set in this specifc directive
                        include_base = inherited_mem
                    else:
                        # use default value
                        # include_base = inherited_members
                        include_base = False

                    api_class_methods = ['__init__', '__call__']
                    ns['members'] = get_members_class(obj, None,
                                                      include_base=include_base)
                    ns['methods'], ns['all_methods'] = \
                        get_members_class(obj, 'method', api_class_methods,
                                          include_base=include_base)
                    ns['attributes'], ns['all_attributes'] = \
                        get_members_class(obj, 'attribute',
                                          include_base=include_base)
                    ns['methods'].sort()
                    ns['attributes'].sort()

                parts = name.split('.')
                if doc.objtype in ('method', 'attribute'):
                    mod_name = '.'.join(parts[:-2])
                    cls_name = parts[-2]
                    obj_name = '.'.join(parts[-2:])
                    ns['class'] = cls_name
                else:
                    mod_name, obj_name = '.'.join(parts[:-1]), parts[-1]

                ns['fullname'] = name
                ns['module'] = mod_name
                ns['objname'] = obj_name
                ns['name'] = parts[-1]

                ns['objtype'] = doc.objtype
                ns['underline'] = len(obj_name) * '='

                # We now check whether a file for reference footnotes exists for
                # the module being documented. We first check if the
                # current module is a file or a directory, as this will give a
                # different path for the reference file. For example, if
                # documenting astropy.wcs then the reference file is at
                # ../wcs/references.txt, while if we are documenting
                # astropy.config.logging_helper (which is at
                # astropy/config/logging_helper.py) then the reference file is set
                # to ../config/references.txt
                if '.' in mod_name:
                    mod_name_dir = mod_name.split('.', 1)[1].replace('.', os.sep)
                else:
                    mod_name_dir = mod_name

                if (not os.path.isdir(os.path.join(base_path, mod_name_dir))
                        and os.path.isdir(os.path.join(base_path,
                                                       mod_name_dir.rsplit(os.sep, 1)[
                                                           0]))):
                    mod_name_dir = mod_name_dir.rsplit(os.sep, 1)[0]

                # We then have to check whether it exists, and if so, we pass it
                # to the template.
                if os.path.exists(
                        os.path.join(base_path, mod_name_dir, 'references.txt')):
                    # An important subtlety here is that the path we pass in has
                    # to be relative to the file being generated, so we have to
                    # figure out the right number of '..'s
                    ndirsback = path.replace(base_path, '').count(os.sep)
                    ref_file_rel_segments = ['..'] * ndirsback
                    ref_file_rel_segments.append(mod_name_dir)
                    ref_file_rel_segments.append('references.txt')
                    ns['referencefile'] = os.path.join(*ref_file_rel_segments).replace(
                        os.sep, '/')

                rendered = template.render(**ns)
                f.write(cleanup_whitespace(rendered))
            finally:
                f.close()


def process_autosum_aliases_generation(app: Sphinx):

    env = app.builder.env
    processor = ProcessAutosumAliases()

    # gather all files that need to be searched for `.. autosum_aliases::`
    files_to_search = []
    for docname in env.found_docs:
        filename = env.doc2path(docname)
        if os.path.isfile(filename):
            files_to_search.append(docname + os.path.splitext(filename)[1])

    # replace autosm_aliases lines with autosummry
    liness = []
    for filename in files_to_search:
        lines = processor.convert_to_autosummary_lines(filename, app)
        liness.append(lines)

    # generate
    for filename, lines in zip(files_to_search, liness):
        if len(lines) > 0:
            processor.generate_docs(
                lines,
                filename,
                app=app,
                builder=app.builder,
                base_path=app.srcdir,
            )


def setup(app: Sphinx):
    app.add_directive('autosum_aliases', AutosumAliases)
    app.connect("builder-inited", process_autosum_aliases_generation)
