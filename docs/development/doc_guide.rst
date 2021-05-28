************************
Documentation Guidelines
************************

This document describes the documentation requirements and guidelines
to be followed during the development of PlasmaPy and affiliated
packages.  PlasmaPy's online documentation is hosted by
`Read the Docs <https://readthedocs.org/>`_ and

* Most recent stable release:
  `https://docs.plasmapy.org <https://docs.plasmapy.org>`_ or
  `https://docs.plasmapy.org/en/stable/ <https://docs.plasmapy.org/en/stable/>`_

* Latest version on GitHub: ``https://docs.plasmapy.org/en/latest/``

Essentials
==========

ReStructuredText
----------------

The documentation for PlasmaPy is written using `reStructuredText (ReST)
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_
as its markup language. ReST is human readable when viewed within a
source code file or when printed out using `help`, and also contains
markup that allows the text to be transformed into `PlasmaPy's online
documentation <https://www.plasampy.org>`_. ReST is the markup language
used in docstrings and in the narrative documentation.  ReST files
end in `.rst`.  Documentation enclosed in triple quotes (`""" ... """`)
within `.py` files is written in ReST.

Markdown
--------

A few of PlasmaPy's files are written using `Markdown
<https://www.markdownguide.org/>`_, such as README files and licenses
from other packages.  These files end with `.md`.  Markdown is the
default text for posts on `GitHub <https://github.com>`_.
`GitHub Flavored Markdown <GitHub Flavored Markdown>`_ contains

Sphinx
------


.. add plasmapy-sphinx later

Writing documentation
=====================

Docstrings
----------

The documentation for



Style guidelines
----------------

* All public functions, classes, and other objects should have a
  docstring.

*

Many words and software packages have more than one common acronym
  or spelling.

  -

Previewing documentation
========================

When a pull request is submitted to

.. Add picture of CI


Building documentation
======================
Documentation is built from the master branch on every commit pushed
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

Docstrings
==========

* All public classes, methods, and functions should have docstrings.

* PlasmaPy uses the `numpydoc standard for docstrings
  <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_\
  .

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
