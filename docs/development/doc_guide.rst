************************
Documentation guidelines
************************

Documentation that is up-to-date and understandable is vital to the
health of a software project.  This page describes the documentation
requirements and guidelines to be followed during the development of
PlasmaPy and affiliated packages.

PlasmaPy's online documentation is hosted by `Read the Docs
<https://readthedocs.org/>`_ and can be found online at:

* Most recent stable release:
  `https://docs.plasmapy.org <https://docs.plasmapy.org>`_ or
  `https://docs.plasmapy.org/en/stable/ <https://docs.plasmapy.org/en/stable/>`_

* Latest version on GitHub:
  `https://docs.plasmapy.org/en/latest/ <https://docs.plasmapy.org/en/latest/>`_

Markup languages
================

ReStructuredText
----------------

The documentation for PlasmaPy is written using `reStructuredText (ReST)
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
as its markup language. ReST is human readable when viewed within a
source code file or when printed out using `help`, and also contains
markup that allows the text to be transformed into `PlasmaPy's online
documentation <https://www.plasampy.org>`_. ReST is the markup language
used in docstrings and in the narrative documentation.  ReST files
end in ``.rst``.

This code block contains some ReST examples.

.. code-block:: rst

  ==============
  Document title
  ==============

  Here is a link to `PlasmaPy's website <https://www.plasmapy.org>`_.

  Heading 1
  =========

  Heading 1.1
  -----------

  Heading 1.1.1
  ~~~~~~~~~~~~~

  Here is a reference to `plasmapy.particles` that will write out the
  full namespace when Sphinx generates the documentation and generates
  the link. Only the word "Particle" will show up if we include write
  `~plasmapy.particles.particle_class.Particle`.

  Sphinx can format Python code blocks.

  .. code-block:: python

      def sample_function():
          return 42

  Math can be written using LaTeX commands.

  .. math::

      \alpha = \beta + \gamma

  Math can be in-line, like :math:`x`. Using unicode characters
  makes math like :math:`α + β + γ` easier to read.

Substitutions
~~~~~~~~~~~~~

There are some functions and classes that show up numerous times in the
documentation. Certain `common links`_ are defined in
``docs/common_links.rst``.  Instead of having to write
``~astropy.units.Quantity`` or
``~plasmapy.particles.particle_class.Particle``, we can write
``|Quantity|`` and ``Particle``.

Substitutions
~~~~~~~~~~~~~

The `common_links.rst
<https://github.com/PlasmaPy/PlasmaPy/blob/master/docs/common_links.rst>`_
file defines reStructuredText substitutions that can be used in
PlasmaPy's documentation.

.. code-block:: RST

    Instead of writing `~plasmapy.particles.particle_class.Particle`,
    we can write |Particle|.

Substitutions should be defined when they are used in multiple files.


Markdown
--------

A few of PlasmaPy's files are written using `Markdown
<https://www.markdownguide.org/>`_, such as README files and licenses
from other packages.  Markdown files end with `.md`.  Posts on GitHub
are written in `GitHub Flavored Markdown
<https://github.github.com/gfm/>`_.

Building documentation
======================

Sphinx
------

`Sphinx <https://www.sphinx-doc.org>`_ is the software package that is used to



PlasmaPy documentation is built with the following Sphinx extensions:

* `sphinx.ext.autodoc`
* `sphinx.ext.intersphinx`
* `sphinx.ext.graphviz`
* `sphinx.ext.mathjax`
* `sphinx.ext.napoleon`
* `sphinx.ext.todo`
* `nbsphinx`
* `sphinx_copybutton`
* `sphinx_gallery.load_style`
* `IPython.sphinxext.ipython_console_highlighting`
* `sphinx_changelog`
* `plasmapy_sphinx`




Documentation is built from the main branch on every commit pushed
to it.

Sphinx, the documentation generator of PlasmaPy, uses reStructuredText (reST)
as its markup language. A primer on reST is available at this `webpage
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
of Sphinx's website.

Using sphinx within the project
-------------------------------
To build docs locally, either:

* use `Tox <https://tox.readthedocs.io/en/latest/>`_ with ``tox -e build_docs`` from within the main PlasmaPy repository directory, or
* enter the ``docs`` directory and run ``make html``.

Afterwards, open ``docs/_build/html/index.html`` with your browser of choice.

Do try to solve warnings in documentation when writing your code. To enforce this,
The ``build_docs`` environment is set to fail on encountering any warnings via
the ``-W`` flag to ``sphinx-build``

.. note::
   The ``tox -e build_docs_no_examples`` command will build the documentation without
   executing the :ref:`example notebooks <example_notebooks>`. It will also
   pass with warnings.


Read the Docs
-------------

PlasmaPy's documentation is hosted on `Read the Docs`_.


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
          If `True`, return :math:`a - b`.  If `False`, then return
          :math:`b - a`.  Defaults to `True`.

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
  docstring.  These



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
  if they contain backslashes.  A raw string literal is denoted by
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

.. replace:: _`Read the Docs`: https://readthedocs.org/
