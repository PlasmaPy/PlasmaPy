************************
Documentation guidelines
************************

Documentation that is up-to-date and understandable is vital to the
health of a software project. This page describes the documentation
requirements and guidelines to be followed during the development of
PlasmaPy and affiliated packages.

`PlasmaPy's online documentation`_ is hosted by `Read the Docs`_ and is
available at these locations.

* The documentation corresponding to the most recent official release
  is labelled ``stable`` and is found at
  `https://docs.plasmapy.org/en/stable/
  <https://docs.plasmapy.org/en/stable/>`_.
  The link `https://docs.plasmapy.org/ <https://docs.plasmapy.org/>`_
  redirects to ``stable``.

* The documentation corresponding to the ``main`` branch on
  `PlasmaPy's GitHub repository`_ is labelled ``latest`` and can be
  found at `https://docs.plasmapy.org/en/latest/
  <https://docs.plasmapy.org/en/latest/>`_.

A preview of the documentation is generated every time a pull request
is created or updated.  You can access this preview by scrolling down
to the checks at the bottom of a pull request, and clicking on
``Details`` next to ``docs/readthedocs.org:plasmapy``.

Building documentation
======================

To install all dependencies required to develop PlasmaPy on your local
computer, enter the top-level directory of the cloned repository and
run

.. code-block:: bash

  pip install tox -r requirements.txt

You can use `tox`_ to build the documentation from within the main
PlasmaPy repository directory by running

.. code-block:: bash

  tox -e build_docs

You can access the documentation landing page by opening
``docs/_build/html/index.html`` with your browser of choice.

When writing documentation, please make sure to fix any warnings that
arise.  To enforce this, the ``build_docs`` environment is set to fail
on encountering any warnings via the ``-W`` flag to ``sphinx-build``.

You can shorten the documentation build by running

.. code-block:: bash

  tox -e build_docs_no_examples

in order to build the documentation without executing the
:ref:`example notebooks <example_notebooks>`.  This command will pass
even if there are warnings.

If you have `make <https://www.gnu.org/software/make/>`_ installed,
then the documentation may be build by entering the ``docs/`` directory
and running

.. code-block:: bash

  make html

Documentation tools
===================

ReStructuredText
----------------

PlasmaPy's documentation is written using the `reStructuredText (ReST)
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
markup language. ReST is human readable when viewed within a
source code file or when printed out using `help`. ReST also contains
markup that allows the text to be transformed into `PlasmaPy's online
documentation`_. ReST files end in ``.rst``. Documentation contained
within ``.py`` files is written in ReST.

ReStructuredText Examples
~~~~~~~~~~~~~~~~~~~~~~~~~

Here we show some ReST that are commonly used in PlasmaPy.

This is an example of including headings for the document title, sections,
subsections, and so on.

.. code-block:: rst

  ==============
  Document title
  ==============

  Heading 1
  =========

  Heading 2
  ---------

  Heading 3
  ~~~~~~~~~

We can link to code objects by enclosing them in back ticks.

.. code-block:: rst

  Here is a reference to `plasmapy.particles` that will write out the
  full namespace when Sphinx generates the documentation and generates
  the link. Only the word "Particle" will show up if we prepend a
  tilde like in `~plasmapy.particles.particle_class.Particle`.

This linking will work for `python` commands as well as commonly used
packages like `numpy`, `astropy`, `scipy`, and `pandas`.  The full list of
`intersphinx <https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html>`_
mappings are defined in the ``intersphinx_mapping`` variable in
` ``docs/conf.py`` <https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/conf.py>`_.

Sphinx can format code blocks for Python and the Python console.

  .. code-block:: rst
    .. code-block:: python

      def sample_function():
          return 42

    .. code-block:: pycon

      >>> print(6 * 9)
      54

Here are some examples for linking to websites.

.. code-block:: rst

  Here is a link to `PlasmaPy's website <https://www.plasmapy.org>`_.

  We can link to PlasmaPy's latest documentation_ or `Python's website`_.

  .. _documentation: https://docs.plasmapy.org/en/latest/
  .. _`Python's documentation`: https://www.python.org/

Math can be written using `LaTeX <https://www.latex-project.org/>`_ commands

.. code-block:: rst

  .. math::

    \alpha = \beta + \gamma

Math can be in-line, like `` :math:`x` ``. Using unicode characters
makes math like `` :math:`α + β + γ` `` easier to read in source code.

Markdown
--------

A few of PlasmaPy's files are written using `Markdown
<https://www.markdownguide.org/>`_, such as README files and licenses
from other packages. Markdown is simpler but more limited than ReST.
Markdown files end with `.md`. Posts on GitHub are written in
`GitHub Flavored Markdown <https://github.github.com/gfm/>`_.
The following code block contains a few common examples of Markdown
formatting.

.. code-block:: markdown

  # Header 1

  ## Header 2

  Here is a link to [PlasmaPy's documentation](https://docs.plasmapy.org).

  We can make text **bold** or *italic*.

  We can write in-line code like `x = 1` or create a Python code block:

    ```Python
    y = 2
    z = 3
    ```

Sphinx
------

`Sphinx <https://www.sphinx-doc.org/>`_ is the software used to generate
`PlasmaPy's online documentation`_ from ReST files and Python docstrings.

Configuration
~~~~~~~~~~~~~

The ``docs/conf.py`` file contains the configuration information needed
to customize Sphinx behavior.
`Sphinx's documentation <https://www.sphinx-doc.org/>`_ lists the
`configuration options
<https://www.sphinx-doc.org/en/master/usage/configuration.html>`_ that
can be set.

Sphinx extensions
~~~~~~~~~~~~~~~~~

PlasmaPy documentation is built with the following Sphinx extensions:

* `sphinx.ext.autodoc
  <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_
  for including documentation from docstrings
* `sphinx.ext.intersphinx
  <https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html>`_
  for linking to other projects' documentation
* `sphinx.ext.graphviz
  <https://www.sphinx-doc.org/en/master/usage/extensions/graphviz.html>`_
  to allow `Graphviz <https://graphviz.org/>`_ graphs to be included
* `sphinx.ext.mathjax
  <https://www.sphinx-doc.org/en/master/usage/extensions/math.html#module-sphinx.ext.mathjax>`_
  for math rendering with `MathJax <https://www.mathjax.org/>`_
* `sphinx.ext.napoleon
  <https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html>`_
  for allowing NumPy style docstrings
* `sphinx.ext.todo
  <https://www.sphinx-doc.org/en/master/usage/extensions/todo.html>`_ to support
  ``todo`` directives
* `nbsphinx <https://nbsphinx.readthedocs.io>`_ for including
  `Jupyter`_ notebooks
* `sphinx_copybutton <https://sphinx-copybutton.readthedocs.io>`_ to add
  a "copy" button for code blocks
* `sphinx_gallery.load_style
  <https://sphinx-gallery.github.io/stable/advanced.html?highlight=load_style#using-only-sphinx-gallery-styles>`_
  for using sphinx-gallery styles
* IPython.sphinxext.ipython_console_highlighting
* `sphinx_changelog <https://sphinx-changelog.readthedocs.io>`_
  for rendering `towncrier`_ changelogs
* `plasmapy_sphinx` for customizations created for use in PlasmaPy

References to other packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Intersphinx <https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html>`_
allows the automatic generation of links to the documentation of
objects in other projects.  The mappings are defined in the
``intersphinx_mapping`` dictionary in ``docs/conf.py``, and include
`python`, `numpy`, `scipy`, `astropy`, `pandas`, `sphinx`, and `sphinx_automodapi`.

When we include ``` `astropy.units.Quantity` ``` in ReST documentation,
it will show up as `astropy.units.Quantity` and link to the appropriate
`object` in Astropy's documentation.

Substitutions
~~~~~~~~~~~~~

Some functions and classes are referred to repeatedly throughout the
documentation. ReST allows us to `define substitutions
<https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#substitution-definitions>`_.

.. code-block:: rst

    .. |Particle| replace:: `~plasmapy.particles.particle_class.Particle`

PlasmaPy has certain common substitutions pre-defined so that they can
be used throughout the documentation. For example, we can write
``|Quantity|`` instead of ``~astropy.units.Quantity``, and
``|Particle|`` instead of ``~plasmapy.particles.particle_class.Particle``.
For an up-to-date list of substitutions, please refer to the
`docs/common_links.rst`_ file.

Because substitutions are performed when Sphinx builds the
documentation, they will not be performed before `help` accesses the
docstring of an `object`. For example, when ``|Particle|`` is used in
a docstring, `help` will show it as ``|Particle|`` rather than
``~plasmapy.particles.particle_class.Particle``. Consequently,
substitutions should not be used in docstrings when it is important
that users have quick access to the full path of the `object` (such as
in the ``See Also`` section).


Writing documentation
=====================

Docstrings
----------

A docstring is a comment at the beginning of a function or another
object that provides information on how to use that function.
Docstrings begin with ``r"""`` (required when including backslashes,
such as using LaTeX code in equations) or ``"""``, and end with
``"""``.


.. code-block:: python

  def subtract(a, b, *, switch_order=False):
      r"""
      Return the difference between two integers. ← state what function does in 1–2 lines

      Add ∼1–3 sentences here for an extended summary of what the function
      does.

      Add ∼1–3 sentences here to clarify what the function does, if
      necessary. This extended summary is a good place to briefly define
      the quantity that is being returned.

      .. math::

          f(a, b) = a - b

      Parameters
      ----------
      a : `int`
          The left multiplicand.

      b : `int`
          The right multiplicand.

      switch_order : `bool`, optional, keyword-only
          If `True`, return :math:`a - b`. If `False`, then return
          :math:`b - a`. Defaults to `True`.

      Returns
      -------
      float
          The product of ``a`` and ``b``.

      Raises
      ------
      `TypeError`
          If ``a`` or ``b`` is not a `float`.

      Notes
      -----
      This section is used to provide extra information that cannot fit in
      the extended summary near the beginning of the docstring. This
      section should include a discussion of the physics behind a
      particular concept that should be understandable to someone who is
      taking their first plasma physics class. This section can also
      include a derivation of the quantity being calculated or a
      description of a particular algorithm.

      The next section contains example references to a journal article
      [1]_, a book [2]_, and a software package. Using a link with the
      digital object identifier (DOI) is helpful because of its permanence.
      We can also link to a website [3]_, though this is discouraged because

      References
      ----------
      .. [1] J. E. Foster, `Plasma-based water purification: Challenges and
         prospects for the future <https://doi.org/10.1063/1.4977921>`_,
         Physics of Plasmas, 22, 05501 (2017).

      .. [2] E. Gamma, R. Helm, R. Johnson, J. Vlissides, `Design Patterns:
         Elements of Reusable Object-Oriented Software
         <https://www.oreilly.com/library/view/design-patterns-elements/0201633612/>`_

      .. [3]

      Examples
      --------
      Include a few example usages of the function here.

      >>> from package.subpackage.module import subtract
      >>> subtract(9, 6)
      3
      >>> subtract(9, 6, switch_order=True)
      -3

      Here is an example where one line is too short.

      >>>

      PlasmaPy's test suite will check that these commands return the
      output that
      """
      if not isinstance(a, float) or not isinstance(b, float):
          raise TypeError("The arguments to multiply should be floats.")

      return b - a if switch_order else a - b

Documentation guidelines
========================

* All public functions, classes, and other objects should have a
  docstring.

* Documentation should be intended for

* Private functions, classes, and objects should generally have a
  docstring. These



*

Many words and software packages have more than one common acronym
  or spelling.

  -

Previewing documentation
========================

When a pull request is submitted to

.. Add picture of CI

Docstrings
==========

* All public classes, methods, and functions should have docstrings.

* PlasmaPy uses the `numpydoc`_ standard for docstrings.

* Docstrings must be raw string `literals
  <https://docs.python.org/3/reference/lexical_analysis.html#literals>`_
  if they contain backslashes. A raw string literal is denoted by
  having an ``r`` immediately precede quotes or triple quotes:

.. code-block:: python

   r"""
   I did not like unstable eigenfunctions at first, but then they
   grew on me.
   """

* Simple private functions may need only a one-line docstring.

Narrative Documentation
=======================

* Each subpackage must have narrative documentation describing its
  use.

.. _`Read the Docs`: https://readthedocs.org/
