"""The configuration file for building PlasmaPy's documentation."""

#!/usr/bin/env python3

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

# isort: off
sys.path.insert(0, os.path.abspath(".."))  # noqa: PTH100
sys.path.insert(0, os.path.abspath("."))  # noqa: PTH100
# isort: on

import cff_to_rst

from datetime import datetime
from pkg_resources import parse_version
from sphinx.application import Sphinx

# Generate author list from CITATION.cff

cff_to_rst.main()

from plasmapy import __version__ as release  # noqa: E402

# -- General configuration ------------------------------------------------
autosummary_generate = True
automodapi_custom_groups = {
    "aliases": {
        "title": "Aliases",
        "description": (
            """
            PlasmaPy provides :term:`aliases` of the most common plasma
            functionality for user convenience. Aliases in PlasmaPy are
            denoted with a trailing underscore (e.g., ``alias_``). For
            further details, please refer to the :ref:`contributor
            guide's section on aliases <aliases>`.
            """
        ),
        "dunder": "__aliases__",
    },
    "lite-functions": {
        "title": "Lite-Functions",
        "description": (
            """
            :term:`Lite-functions` are optimized versions of existing
            `plasmapy` functions that are intended for applications where
            computational efficiency matters most. Lite-functions accept
            numbers and NumPy arrays that are implicitly assumed to be
            in SI units, and do not accept |Quantity| objects as inputs.
            For further details, please refer to the :ref:`contributor
            guide's section on lite-functions <lite-functions>`.

            .. caution::

               Lite-functions do not include the safeguards that are
               included in most `plasmapy.formulary` functions. When
               using lite-functions, it is vital to double-check your
               implementation!
            """
        ),
        "dunder": "__lite_funcs__",
    },
}
automodapi_group_order = (
    "modules",
    "classes",
    "exceptions",
    "warnings",
    "functions",
    "aliases",
    "lite-functions",
    "variables",
)

# If your documentation needs a minimal Sphinx version, state it here.

needs_sphinx = "6.1.3"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones. When extensions are removed or added, please update the section
# in docs/doc_guide.rst on Sphinx extensions.

extensions = [
    "hoverxref.extension",
    "IPython.sphinxext.ipython_console_highlighting",
    "nbsphinx",
    "notfound.extension",
    "plasmapy_sphinx",
    "sphinx.ext.autodoc",
    "sphinx.ext.extlinks",
    "sphinx.ext.graphviz",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx_changelog",
    "sphinx_codeautolink",
    "sphinx_copybutton",
    "sphinx_gallery.load_style",
    "sphinx_issues",
    "sphinx_reredirects",
    "sphinx_tabs.tabs",
    "sphinxcontrib.bibtex",
]

# Configure sphinxcontrib-bibtex

bibtex_bibfiles = ["bibliography.bib"]
bibtex_default_style = "plain"
bibtex_reference_style = "author_year"
bibtex_cite_id = "{key}"

# Intersphinx generates automatic links to the documentation of objects
# in other packages. When mappings are removed or added, please update
# the section in docs/doc_guide.rst on references to other packages.

intersphinx_mapping = {
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "lmfit": ("https://lmfit.github.io/lmfit-py/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "numba": ("https://numba.readthedocs.io/en/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "pytest": ("https://docs.pytest.org/en/stable/", None),
    "python": ("https://docs.python.org/3/", None),
    "readthedocs": ("https://docs.readthedocs.io/en/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
    "sphinx_automodapi": (
        "https://sphinx-automodapi.readthedocs.io/en/latest/",
        None,
    ),
}

hoverxref_intersphinx = [
    "astropy",
    "lmfit",
    "numba",
    "numpy",
    "pandas",
    "pytest",
    "python",
    "readthedocs",
    "scipy",
    "sphinx",
    "sphinx_automodapi",
]

autoclass_content = "both"
autodoc_typehints_format = "short"

# Configure sphinx-issues

issues_github_path = "PlasmaPy/PlasmaPy"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The root toctree document.
root_doc = "index"

# General information about the project.
project = "PlasmaPy"
author = "PlasmaPy Community"
copyright = f"2015–{datetime.utcnow().year}, {author}"  # noqa: A001
language = "en"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The full version, including alpha/beta/rc tags.
#  Note: If plasmapy.__version__ can not be defined then it is set to 'unknown'.
#        However, release needs to be a semantic style version number, so set
#        the 'unknown' case to ''.
release = "" if release == "unknown" else release
pv = parse_version(release)
release = pv.public
version = ".".join(release.split(".")[:2])  # short X.Y version
revision = pv.local[1:] if pv.local is not None else ""

# This is added to the end of RST files — a good place to put substitutions to
# be used globally.
rst_epilog = ""
with open("common_links.rst") as cl:  # noqa: PTH123
    rst_epilog += cl.read()

rst_prolog = """
.. role:: py(code)
   :language: python

.. role:: bash(code)
   :language: bash
"""

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "notebooks/langmuir_samples",
    "**.ipynb_checkpoints",
    "plasmapy_sphinx",
    "common_links.rst",
    "**Untitled*",
]

html_extra_path = ["robots.txt"]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

default_role = "py:obj"

# Customizations for make linkcheck using regular expressions

linkcheck_allowed_redirects = {
    r"https://doi\.org/.+": r"https://.+",  # DOI links are more persistent
    r"https://docs.+\.org": r"https://docs.+\.org/en/.+",
    r"https://docs.+\.io": r"https://docs.+\.io/en/.+",
    r"https://docs.+\.com": r"https://docs.+\.com/en/.+",
    r"https://docs.+\.dev": r"https://docs.+\.dev/en/.+",
    r"https://en.wikipedia.org/wiki.+": "https://en.wikipedia.org/wiki.+",
    r"https://.+\.readthedocs\.io": r"https://.+\.readthedocs\.io/en/.+",
    r"https://www\.sphinx-doc\.org": r"https://www\.sphinx-doc\.org/en/.+",
    r"https://.+/github\.io": r"https://.+/github\.io/en/.+",
    r"https://.+": r".+(google|github).+[lL]ogin.+",  # some links require logins
    r"https://jinja\.palletsprojects\.com": r"https://jinja\.palletsprojects\.com/.+",
    r"https://pip\.pypa\.io": r"https://pip\.pypa\.io/en/.+",
    r"https://www.python.org/dev/peps/pep.+": "https://peps.python.org/pep.+",
}

# Hyperlinks for `make linkcheck` to ignore, such as links that point to
# setting options in PlasmaPy's GitHub account that require a login.

linkcheck_ignore = ["https://github.com/PlasmaPy/PlasmaPy/settings/secrets/actions"]

linkcheck_anchors = True
linkcheck_anchors_ignore = [
    "/room",
    r".+openastronomy.+",
    "L[0-9].+",
    "!forum/plasmapy",
]

redirects = {
    "contributing/install_dev": "../contributing/getting_ready.html",
    "development": "../contributing/",
    "development/changelog_guide": "../contributing/changelog_guide.html",
    "development/code_guide": "../contributing/code_guide.html",
    "development/doc_guide": "../contributing/doc_guide.html",
    "development/index": "../contributing/index.html",
    "development/install_dev": "../contributing/getting_ready.html",
    "development/release_guide": "../contributing/release_guide.html",
    "development/testing_guide": "../contributing/testing_guide.html",
    "whatsnew": "../changelog/",
    "whatsnew/0.1.0": "../changelog/0.1.0.html",
    "whatsnew/0.1.1": "../changelog/0.1.0.html",
    "whatsnew/0.2.0": "../changelog/0.1.0.html",
    "whatsnew/0.3.1": "../changelog/0.1.0.html",
    "whatsnew/0.4.0": "../changelog/0.1.0.html",
    "whatsnew/0.5.0": "../changelog/0.1.0.html",
    "whatsnew/0.6.0": "../changelog/0.1.0.html",
    "whatsnew/0.7.0": "../changelog/0.1.0.html",
    "whatsnew/0.8.1": "../changelog/0.1.0.html",
    "whatsnew/0.9.0": "../changelog/0.1.0.html",
    "whatsnew/0.9.1": "../changelog/0.1.0.html",
    "whatsnew/2023.1.0": "../changelog/2023.1.0.html",
    "whatsnew/index": "../changelog/index.html",
}

# Use a code highlighting style that meets the WCAG AA contrast standard
pygments_style = "default"

# adapted from sphinx-hoverxref conf.py
if os.environ.get("READTHEDOCS"):
    # Building on Read the Docs
    hoverxref_api_host = "https://readthedocs.org"

    if os.environ.get("PROXIED_API_ENDPOINT"):
        # Use the proxied API endpoint
        # - A RTD thing to avoid a CSRF block when docs are using a
        #   custom domain
        hoverxref_api_host = "/_"

hoverxref_tooltip_maxwidth = 600  # RTD main window is 696px
hoverxref_auto_ref = True
hoverxref_mathjax = True
hoverxref_sphinxtabs = True

# hoverxref has to be applied to these
hoverxref_domains = ["py", "cite"]
hoverxref_roles = ["confval", "term"]

hoverxref_role_types = {
    # roles with cite domain
    "p": "tooltip",
    "t": "tooltip",
    #
    # roles with py domain
    "attr": "tooltip",
    "class": "tooltip",
    "const": "tooltip",
    "data": "tooltip",
    "exc": "tooltip",
    "func": "tooltip",
    "meth": "tooltip",
    "mod": "tooltip",
    "obj": "tooltip",
    #
    # roles with std domain
    "confval": "tooltip",
    "hoverxref": "tooltip",
    "ref": "tooltip",
    "term": "tooltip",
}

# Using sphinx.ext.extlinks lets us simplify the process of creating
# links to commonly used external sites. The key of the extlink becomes
# a new role, and the corresponding tuple contains the base url and the
# caption. For example, we can now do :orcid:`0000-0000-0000-0000` and
# have a link create to the corresponding ORCID page. New roles should
# be added to rst-roles in tox.ini to avoid being caught by
# flake8-rst-docstrings.

extlinks = {
    "orcid": ("https://orcid.org/%s", "%s"),
    "wikipedia": ("https://en.wikipedia.org/wiki/%s", "%s"),
}

# Specify patterns to ignore when doing a nitpicky documentation build.
# These may include common expressions like "real number" as well as
# workarounds for nested inline literals as defined in docs/common_links.py

python_role = "py:.*"

nitpick_ignore_regex = [
    # Before adding patterns for type specifications in docstrings, note
    # that information on the *meaning* of a parameter should be
    # included in the parameter description instead.
    (python_role, "and"),
    (python_role, "array .*"),
    (python_role, "array_like"),
    (python_role, "callable"),
    # for defaults that are words, numbers, particle symbols like "p+"
    (python_role, r"default: [-\+]?\w+[-\+]?\.?\d*"),
    # for defaults that are lists, tuples, sets, dictionaries, and strings
    (python_role, r"default: ((\[.*\])|(\(.*\))|(\{.*\})|(\".*\")|(\'.*\'))"),
    # for defaults that are calls like Particle("p+") or items from a dictionary
    (python_role, r"default: \w+[\.\w]*[\(\[].*[\)\]]"),
    (python_role, "dictionary.*"),
    (python_role, "function"),
    (python_role, ".*integer.*"),
    (python_role, "iterable"),
    (python_role, "key"),
    (python_role, "keyword-only"),
    (python_role, ".* object"),
    (python_role, "optional"),
    (python_role, "or"),
    (python_role, "Real"),
    (python_role, ".*real number.*"),
    (python_role, ".*representation.*"),
    (python_role, "shape.*"),
    (python_role, r"u\..*"),
    (python_role, ".*Unit.*"),
    # pytest helpers
    (python_role, "_pytest.*"),
    (python_role, "Failed"),
    # charged_particle_radiography
    (python_role, "1"),
    (python_role, "2 ints"),
    (python_role, "a single int"),
    (python_role, "same"),
    (python_role, "Tuple of 1"),
    # thomson
    (python_role, "Ne"),
    (python_role, "Ni"),
    # utils
    (python_role, "docstring of"),
    (python_role, "validation specifications"),
    # for reST workarounds defined in docs/common_links.rst
    (python_role, "git"),
    (python_role, "h5py"),
    (python_role, "IPython.sphinxext.ipython_console_highlighting"),
    (python_role, "lmfit"),
    (python_role, "mpmath"),
    (python_role, "nbsphinx"),
    (python_role, "numba"),
    (python_role, "xarray"),
    # plasmapy_sphinx
    (python_role, "automod.*"),
    (python_role, "Builder"),
    (python_role, "docutils.*"),
    (python_role, "level"),
    (python_role, ".*member.*"),
    (python_role, "OptionSpec"),
    (python_role, "py"),
    (python_role, "[Ss]phinx.*"),  # also for reST workarounds in docs/common_links.rst
    # The following patterns still need to be fixed.
    (python_role, "json.decoder.JSONDecoder"),
    (python_role, "plasmapy.analysis.swept_langmuir.find_floating_potential"),
    (python_role, "plasmapy.particles.particle_collections"),
    (python_role, "plasmapy.utils.decorators.lite_func"),
]

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_logo = "./_static/with-text-light-190px.png"
html_theme_options = {
    "logo_only": True,
    #
    # TOC options
    #   https://sphinx-rtd-theme.readthedocs.io/en/stable/configuring.html#theme-options
    "includehidden": False,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# A list of prefixes that are ignored for sorting the Python module
# index (e.g., if this is set to ['foo.'], then foo.bar is shown under
# B, not F).
modindex_common_prefix = ["plasmapy."]

# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "PlasmaPydoc"

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',
    #
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',
    #
    # Additional stuff for the LaTeX preamble.
    # 'preamble': '',
    #
    # Latex figure (float) alignment
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (
        root_doc,
        "PlasmaPy.tex",
        "PlasmaPy Documentation",
        "PlasmaPy Community",
        "manual",
    )
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(root_doc, "plasmapy", "PlasmaPy Documentation", [author], 1)]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        root_doc,
        "PlasmaPy",
        "PlasmaPy Documentation",
        author,
        "PlasmaPy",
        "Python package for plasma physics",
        "Miscellaneous",
    )
]

html_favicon = "./_static/icon.ico"

# -- NBSphinx options

nbsphinx_thumbnails = {
    "notebooks/*": "_static/graphic-circular.png",
    "notebooks/*/*": "_static/graphic-circular.png",
    "notebooks/diagnostics/langmuir_analysis": (
        "_static/notebook_images/langmuir_analysis.png"
    ),
    "notebooks/formulary/magnetosphere": (
        "_static/notebook_images/mms.png"
    ),  # public domain
    "notebooks/getting_started/units": (
        "_static/notebook_images/astropy_logo_notext.png"
    ),  # CC BY-SA
    "notebooks/formulary/solar_plasma_beta": "_static/notebook_images/coronal_loops.png",
    "notebooks/plasma/grids_cartesian": (
        "_static/notebook_images/uniform_grid_thumbnail.png"
    ),
    "notebooks/plasma/grids_nonuniform": (
        "_static/notebook_images/nonuniform_grid_thumbnail.png"
    ),
}

# adapted from
# https://github.com/spatialaudio/nbsphinx/blob/58b8034dd9d7349c1b4ac3e7a7d6baa87ab2a6a9/doc/conf.py

# This is processed by Jinja2 and inserted before each notebook
nbsphinx_prolog = r"""
{% set docname = 'docs/' + env.doc2path(env.docname, base=None) %}
{% set nb_base = 'tree' if env.config.revision else 'blob' %}
{% set nb_where = env.config.revision if env.config.revision else 'main' %}

.. raw:: html

    <div class="admonition note">
      <p style="margin-bottom:0px">
        This page was generated by
        <a href="https://nbsphinx.readthedocs.io/">nbsphinx</a> from
        <a class="reference external" href="https://github.com/PlasmaPy/PlasmaPy/{{ nb_base|e }}/{{ nb_where|e }}/{{ docname|e }}">{{ docname|e }}</a>.
        <br>
        Interactive online version:
        <a href="https://mybinder.org/v2/gh/PlasmaPy/PlasmaPy/{{ nb_where|e }}/?filepath={{ docname|e }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a>.
      </p>
    </div>

.. raw:: latex

    \nbsphinxstartnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{The following section was generated from
    \sphinxcode{\sphinxupquote{\strut {{ docname | escape_latex }}}} \dotfill}}
"""


def setup(app: Sphinx) -> None:
    app.add_config_value("revision", "", rebuild=True)
    app.add_css_file("css/admonition_color_contrast.css")
    app.add_css_file("css/plasmapy.css", priority=600)
