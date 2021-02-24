************************
Documentation Guidelines
************************

This document describes the documentation requirements and guidelines
to be followed during the development of PlasmaPy and affiliated
packages.

Building documentation
======================
Documentation is built from the master branch on every commit pushed
to it.

Sphinx, the documentation generator of PlasmaPy, uses reStructuredText (reST) as its markup language. A primer on reST is available at this `webpage
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
