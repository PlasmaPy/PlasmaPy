*******************
Documentation Guide
*******************

Documentation that is up-to-date and understandable is vital to the
health of a software project. This page describes the documentation
requirements and guidelines to be followed during the development of
PlasmaPy and affiliated packages.

`PlasmaPy's online documentation`_ is hosted by `Read the Docs`_ and is
available at these locations.

* The documentation corresponding to the most recent official release
  is labelled ``stable`` and is found at
  `https://docs.plasmapy.org/ <https://docs.plasmapy.org/>`_ or
  `https://docs.plasmapy.org/en/stable/ <https://docs.plasmapy.org/en/stable/>`_.

* The documentation corresponding to the ``main`` branch on
  `PlasmaPy's GitHub repository`_ is labelled ``latest`` and can be
  found at `https://docs.plasmapy.org/en/latest/
  <https://docs.plasmapy.org/en/latest/>`_.

A preview of the documentation is generated every time a pull request
is created or updated. You can access this preview by scrolling down
to the checks at the bottom of a pull request, and clicking on
``Details`` next to ``docs/readthedocs.org:plasmapy``.

Building documentation
======================

There are two methods for building the documentation: make_ and tox_.
Using make_ will build the documentation based off of what is in the
current directory structure.  Using make_ is quicker for local builds
but requires that setting up all dependencies ahead of time.  Using
tox_ does not require setting up all dependencies ahead of time, but
is more computationally intensive since it builds the package and
environment before building the documentation.  Consequently, PlasmaPy
uses tox_ for building the documentation on continuous integration
testing platforms.

Prerequisites
-------------

To install all dependencies required to develop PlasmaPy on your local
computer, enter the top-level directory of the cloned repository and
run

.. code-block:: bash

   pip install -r requirements.txt

Building documentation with make
--------------------------------

If you have make_ installed, then you can build the documentation by
entering the ``docs/`` directory and running

.. code-block:: bash

   make html

Including the ``-j n`` flag in the ``make`` command will enable a
parallel build, where ``n`` is replaced with the number of processes.

To remove all files from a documentation build, run

.. code-block:: bash

   make clean

This command is needed when you make a change to a file that does not
trigger Sphinx_ to rebuild the file that you altered.

Building documentation with tox
-------------------------------

You can use tox_ to build the documentation from within the main
PlasmaPy repository directory by running

.. code-block:: bash

   tox -e build_docs

You can access the documentation landing page by opening
``docs/_build/html/index.html`` with your browser of choice.

When writing documentation, please make sure to fix any warnings that
arise. To enforce this, the ``build_docs`` environment is set to fail
on encountering any warnings via the ``-W`` flag to ``sphinx-build``.

You can shorten the documentation build by running

.. code-block:: bash

   tox -e build_docs_no_examples

in order to build the documentation without executing the
:ref:`example notebooks <example_notebooks>`. This command will pass
even if there are warnings.  The tox_ command may also include the
``--parallel`` flag to allow a parallel build.

Documentation tools
===================

ReStructuredText
----------------

PlasmaPy's documentation is written using the `reStructuredText (reST)
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
markup language. reST_ is human readable when viewed within a
source code file or when printed out using `help`. reST_ also contains
markup that allows the text to be transformed into `PlasmaPy's online
documentation`_. reST_ files end in ``.rst``. Documentation contained
within ``.py`` files are in the form of `docstrings
<https://en.wikipedia.org/wiki/Docstring>`_, which are written in reST_.

ReStructuredText Examples
~~~~~~~~~~~~~~~~~~~~~~~~~

Here we show some examples of reST_ that are commonly used in PlasmaPy.
Please refer to the documentation for Sphinx_ and reST_ for a list
of available directives and syntax.

This is an example of including headings for the document title,
sections, subsections, and so on. The lines surrounding each heading
are the same length as that heading.

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

We can link to code objects by enclosing them in single back ticks.

.. code-block:: rst

   Here is a reference to `plasmapy.particles` that will write out the
   full namespace when Sphinx generates the documentation and generates
   the link. Only the word "Particle" will show up if we prepend a
   tilde like in `~plasmapy.particles.particle_class.Particle`.

This reST_ block will render as:

   Here is a reference to `plasmapy.particles` that will write out the
   full namespace when Sphinx generates the documentation and generates
   the link. Only the word "Particle" will show up if we prepend a
   tilde like in `~plasmapy.particles.particle_class.Particle`.

This linking will work for Python_ commands as well as commonly used
packages like `numpy`, Astropy_, `scipy`, and `pandas`, as described in
the section on :ref:`external-references`.

Sphinx_ can format code blocks for Python and the Python console.

   .. code-block:: rst

      .. code-block:: python

         def sample_function():
             return 42

      .. code-block:: pycon

         >>> print(6 * 9)
         54

This reST block will render as:

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
   .. _`Python's website`: https://www.python.org/

Math can typically be written using
`LaTeX <https://www.latex-project.org/>`_ commands.

.. code-block:: rst

   .. math::

      \alpha = \beta + \gamma

This reST_ block will render as:

   .. math::

      \alpha = \beta + \gamma

Math can be in-line.

.. code-block:: rst

   An example of in-line math is :math:`x`. Using Unicode characters
   like :math:`α + β + γ` makes math easier to read in the source code.

This reST_ block will render as:

   An example of in-line math is :math:`x`. Using Unicode characters
   like :math:`α + β + γ` makes math easier to read in the source code.

Markdown
--------

A few of PlasmaPy's files are written using Markdown_, such as README
files and licenses from other packages. Markdown_ is simpler but more
limited than reST_. Markdown_ files end with ``.md``. Posts on GitHub
are written in `GitHub Flavored Markdown <https://github.github.com/gfm/>`_.
The following code block contains a few common examples of Markdown_
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

Sphinx_ is the software used to generate `PlasmaPy's online documentation`_
from reST_ files and Python docstrings.  It was originally created to
write Python's documentation and has become the de facto software for
documenting Python packages.  Almost all Python open-source packages
utilize Sphinx_ to generate their documentation.

Configuration
~~~~~~~~~~~~~

The `docs/conf.py`_ file contains the configuration information needed
to customize Sphinx_ behavior.
The documentation for Sphinx_ lists the `configuration options
<https://www.sphinx-doc.org/en/master/usage/configuration.html>`_
that can be set.

Sphinx extensions
~~~~~~~~~~~~~~~~~

PlasmaPy documentation is built with the following Sphinx_ extensions:

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
* `plasmapy_sphinx` for customizations created for use in PlasmaPy.
  Note that `plasmapy_sphinx` is expected to be broken out into its own
  package in the future.

These extensions are specified in :confval:`extensions` configuration
value in `docs/conf.py`_.

.. _external-references:

References to other packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Intersphinx_ allows the automatic generation of links to the
documentation of objects in other projects. This cross-package linking
is made possible with the `sphinx.ext.intersphinx` extension.

When we include ``astropy.units.Quantity`` in the documentation, it
will show up as `astropy.units.Quantity` with a link to the appropriate
page in Astropy documentation.  Similarly, ``~astropy.units.Quantity``
will show up as `~astropy.units.Quantity`.

The mappings are defined in the :confval:`intersphinx_mapping`
dictionary in `docs/conf.py`_, and include Python_, `numpy`, `scipy`,
Astropy_, `pandas`, Sphinx_, and `sphinx_automodapi`.

If a cross-link is not working as expected this is usually due to one
of the following reasons:

* A typo;
* The package not being defined in :confval:`intersphinx_mapping`, or
* The referenced source package not properly or fully indexing their
  own code, which is common in Python_ packages.

New source packages may be added, but please verify that references to
a function or class in that package show up correctly in PlasmaPy's
documentation. The name of the package does not always link as
expected.

Substitutions
~~~~~~~~~~~~~

Some functions and classes are referred to repeatedly throughout the
documentation. reST_ allows us to `define substitutions
<https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html#substitution-definitions>`_.

.. code-block:: rst

   .. |Particle| replace:: `~plasmapy.particles.particle_class.Particle`

PlasmaPy has certain common substitutions pre-defined so that they can
be used elsewhere in the documentation. For example, we can write
``|Quantity|`` instead of ``~astropy.units.Quantity``, and
``|Particle|`` instead of ``~plasmapy.particles.particle_class.Particle``.
For an up-to-date list of substitutions, please refer to the
`docs/common_links.rst`_ file.

Since substitutions are performed by Sphinx_ when the documentation
is built, any substitution used in docstrings will not show up when
using Python's `help` function (or the like). For example, when
``|Particle|`` is used in a docstring, `help` will show it as
``|Particle|`` rather than ``~plasmapy.particles.particle_class.Particle``.
Consequently, substitutions should not be used in docstrings when it is
important that users have quick access to the full path of the `object`
(such as in the ``See Also`` section).

Templating
~~~~~~~~~~

Sphinx_ uses the Jinja_ templating engine to generate HTML code. Jinja_
may be used within the documentation when templating is necessary.

Writing documentation
=====================

Docstrings
----------

A docstring is a comment at the beginning of a function or another
object that provides information on how to use that function.
Docstrings begin with ``r"""`` (required when including backslashes,
such as using LaTeX code in equations) or ``"""``, and end with
``"""``.

In order to improve readability and maintain consistency, PlasmaPy uses
the numpydoc_ standard for docstrings. Docstring conventions for
Python are more generally described in `PEP 257`_.

Example docstring
~~~~~~~~~~~~~~~~~

Here is an example docstring in the numpydoc_ format.

.. code-block:: python
   :caption: Example docstring

   import numpy as np
   import warnings

   def subtract(a, b, *, switch_order=False):
       r"""
       Compute the difference between two integers.

       Add ∼1–3 sentences here for an extended summary of what the function
       does. This extended summary is a good place to briefly define
       the quantity that is being returned.

       .. math::

          f(a, b) = a - b

      Parameters
      ----------
      a : `float`
          The left multiplicand.

      b : `float`
          The right multiplicand.

      switch_order : `bool`, optional, keyword-only
          If `True`, return :math:`a - b`. If `False`, then return
          :math:`b - a`. Defaults to `True`.

      Returns
      -------
      difference : float
          The difference between ``a`` and ``b``.

      Raises
      ------
      `ValueError`
          If ``a`` or ``b`` is `~numpy.inf`.

      Warns
      -----
      `UserWarning`
          If ``a`` or ``b`` is `~numpy.nan`.

      See Also
      --------
      add : Add two numbers.

      Notes
      -----
      The "Notes" section provides extra information that cannot fit in
      the extended summary near the beginning of the docstring. This
      section should include a discussion of the physics behind a
      particular concept that should be understandable to someone who is
      taking their first plasma physics class. This section can
      include a derivation of the quantity being calculated or a
      description of a particular algorithm.

      The next section contains example references to a journal article
      [1]_ and a book [2]_.

      References
      ----------
      .. [1] J. E. Foster, `Plasma-based water purification: Challenges and
         prospects for the future <https://doi.org/10.1063/1.4977921>`_,
         Physics of Plasmas, 22, 05501 (2017).

      .. [2] E. Gamma, R. Helm, R. Johnson, J. Vlissides, `Design Patterns:
         Elements of Reusable Object-Oriented Software
         <https://www.oreilly.com/library/view/design-patterns-elements/0201633612/>`_

      Examples
      --------
      Include a few example usages of the function here. Start with simple
      examples and then increase complexity when necessary.

      >>> from package.subpackage.module import subtract
      >>> subtract(9, 6)
      3

      Here is an example of a multi-line function call.

      >>> subtract(
      ...     9, 6, switch_order=True,
      ... )
      -3

      PlasmaPy's test suite will check that these commands provide the
      output that follows each function call.
      """
      if np.isinf(a) or np.isinf(b):
          raise ValueError("Cannot perform substraction operations involving infinity.")

      warnings.warn("The subtract function encountered a nan value.", UserWarning)

      return b - a if switch_order else a - b

Template docstring
~~~~~~~~~~~~~~~~~~

This template docstring may be copied into new functions. Usually only
some of the sections will be necessary for a particular function, but
any sections that are included should be in the order provided.

.. code-block:: python
  :caption: Docstring template

  def sample_function():
      r"""
      Compute ...

      Parameters
      ----------

      Returns
      -------

      Raises
      ------

      Warns
      -----

      See Also
      --------

      Notes
      -----

      References
      ----------

      Examples
      --------

      """

Documentation guidelines
========================

This section contains guidelines and best practices for writing
PlasmaPy documentation.

* Write documentation to be understandable to students taking their
  first course or beginning their first research project in plasma
  science. Include highly technical information only when necessary.

* Use technical jargon sparingly. Define technical jargon when
  necessary.

* Use the `active voice <https://en.wikipedia.org/wiki/Active_voice>`_
  in the present tense.

* Keep the documentation style consistent within a file or module, and
  preferably across all of PlasmaPy's documentation.

* Update code and corresponding documentation at the same time.

* Write sentences that are simple, concise, and direct rather than
  complicated, vague, or ambiguous. Prefer sentences with ≲ 20
  words.

* Avoid idioms, metaphors, and references that are specific to a
  particular culture.

* Many words and software packages have more than one common spelling
  or acronym. Use the spelling that is used in the file you are
  modifying, which is preferably the spelling used throughout
  `PlasmaPy's online documentation`_.

  * More generally, it is preferable to use the spelling that is used
    in `Python's documentation`_ or the spelling that is used most
    commonly.

  * Represent names and acronyms for a software package or language as
    they are represented in the documentation for each project. Common
    examples include "Python", "Astropy", "NumPy", and "reST".

* When referencing PlasmaPy functionality, write the full namespace
  path to where the functionality is defined, not where it is
  conveniently accessed. For example, write
  ``~plasmapy.formulary.parameters.Alfven_speed`` rather than
  ``~plasmapy.formulary.Alfven_speed``.

  This does not necessarily need to be done when referencing external
  packages, since each package may have their own standard.  For
  example, Astropy's |Quantity| class is defined in
  ``astropy.units.quantity.Quantity`` but is also indexed at
  ``~astropy.units.Quantity`` so either option will link to the same
  documentation.

* For readability, limit documentation line lengths to ≲ 72 characters
  (excluding leading spaces in docstrings). Longer line lengths may be
  used when necessary (e.g., for hyperlinks).

* Use indentations of 3 spaces for reST_ blocks.

* Avoid linking to websites that might disappear due to
  `link rot <https://en.wikipedia.org/wiki/Link_rot>`_ such as documents
  hosted on personal websites.

  * When including references, use a link that includes a
    `persistent identifier <https://en.wikipedia.org/wiki/Persistent_identifier>`_
    such as a digital object identifier (DOI) when one is available (e.g.,
    `https://doi.org/10.5281/zenodo.4602818
    <https://doi.org/10.5281/zenodo.4602818>`_).

  * `Wikipedia <https://www.wikipedia.org/>`_ articles may be linked to
    when they contain a well-developed and accurate description of a
    concept.

* Include both the original references for a topic as well as accessible
  pedagogical references. Prefer references that are open access over
  references that require purchase of a subscription or are behind a
  `paywall <https://en.wikipedia.org/wiki/Paywall>`_.

.. note::

   Emphasize important points with `admonitions
   <https://docutils.sourceforge.io/docs/ref/rst/directives.html#admonitions>`_
   like this one.

Docstring guidelines
--------------------

* All functions, classes, and objects that are part of PlasmaPy's
  public Application Programming Interface (API) must have a docstring
  that follows the numpydoc_ standard. Refer to the numpydoc_ standard
  for how to write docstrings for classes, class attributes, and
  constants.

* The short summary statement at the beginning of a docstring should be
  one line long, but may be longer if necessary.

* The extended summary that immediately follows the short summary should
  be ≲ 4 sentences long.  Any additional information should included in
  the "Notes" section.

* The short summary should start on the line immediately following the
  triple quotes.  There should not be any blank lines immediately before
  the closing triple quotes.

* The first line of the docstring for a function or method should begin
  with a word like "Calculate" or "Compute" and end with a period.

* The first line of an object that is not callable (for example, an
  attribute of a class decorated with `property`) should not begin with
  a verb and should end with a period.

* Keep the docstring indented at the same level as the ``r"""`` or
  ``"""`` that begins the docstring, except for reST_ constructs like
  lists, math, and code blocks. Use an indentation of four spaces more
  than the declaration of the object.

  .. code-block:: python

     def f():
         """This is indented four spaces relative to the `def` statement."""

* The first sentence of a docstring of a function should include a
  concise definition of the quantity being calculated, as in the
  following example.

  .. code-block:: python

     def beta(T, n, B):
         """Compute the ratio of thermal pressure to magnetic pressure."""

  When the definition of the quantity being calculated is unable to fit
  on ∼1–2 lines, include the definition in the extended summary instead.

  .. code-block:: python

     def beta(T, n, B):
         """
         Compute plasma beta.

         Plasma beta is the ratio of thermal pressure to magnetic pressure.
         """

* Put any necessary highly technical information in the "Notes" section
  of a docstring.

* Private code objects (e.g., code objects that begin with a single
  underscore, like ``_private_object``) should have docstrings. A
  docstring for a private code object may be a single line, and
  otherwise should be in numpydoc_ format.

  * Docstrings for private code objects do not get rendered in the
    online documentation, and should be intended for contributors.

* Dunder methods (e.g., code objects like ``__add__`` that begin and
  end with two underscores) only need docstrings if it is necessary to
  describe non-standard or potentially unexpected behavior. Custom
  behavior associated with dunder methods should be described in the
  class-level documentation.

  * Docstrings for most dunder methods are not rendered in the online
    documentation and should therefore be intended for contributors.

  * Docstrings for ``__init__`` and ``__new__`` are included in the
    class-level docstring when rendering documentation.

* When an attribute in a class has both a getter (which is decorated
  with `property`) and a ``setter`` decoration, then the getter and
  ``setter`` functionality should be documented in the docstring of
  the attribute decorated with ``@property``.

  .. code-block:: python

     class Person:
         @property
         def age(self):
             """Document both getter and setter here."""
             return self._age

         @age.setter
         def age(self, n):
             self._age = n

Narrative documentation guidelines
----------------------------------

* Each top-level subpackage in PlasmaPy must have corresponding
  narrative documentation.

* Use narrative documentation to describe how different functionality
  works together.

* Narrative documentation should be used when the full scope of some
  functionality cannot be adequately described within only the
  docstrings of that functionality.

* Use title case for page titles (e.g., "Title Case") and sentence case
  for all other headings (e.g., "Sentence case").

* When the narrative documentation does not index a subpackage
  (a directory) or module (a ``.py`` file) with ``automodule``,
  ``automodapi``, or the like, then it is required to create a stub
  file for that particular subpackage or module in `docs/api_static`_ .
  Here are the sample contents for a stub file for
  `plasmapy.particles.atomic`.  This file would be located at
  ``docs/api_static/plasmapy.particles.atomic.rst``.

  .. code-block:: rst

     :orphan:

     `plasmapy.particles.atomic`
     ===========================

     .. currentmodule:: plasmapy.particles.atomic

     .. automodapi::  plasmapy.particles.atomic

.. _`docs/api_static`: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/api_static/
.. _`docs/conf.py`: https://github.com/PlasmaPy/PlasmaPy/blob/main/docs/conf.py
.. _Jinja: https://jinja.palletsprojects.com/
.. _make: https://www.gnu.org/software/make/
.. _`PEP 257`: https://www.python.org/dev/peps/pep-0257/
.. _Sphinx: https://www.sphinx-doc.org/
