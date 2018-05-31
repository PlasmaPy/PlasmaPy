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


Using sphinx within the project
-------------------------------
To build docs locally, run ``python setup.py build_docs`` from
within the main PlasmaPy repository directory, then open
``docs/_build/index.html`` with your browser of choice.

Do try to solve warnings in documentation when writing your code.

Docstrings
==========

* All public classes, methods, and functions should have docstrings.

* PlasmaPy uses the `numpydoc
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  standard for docstrings.

* Docstrings must be raw string `literals
  <https://docs.python.org/3/reference/lexical_analysis.html#literals>`_
  if they contain backslashes.  A raw string literal is denoted by
  having an `r` immediately precede quotes or triple quotes:
  
.. code-block:: python

   r""" I did not like unstable eigenfunctions at first, but then they
   grew on me.
   
   """
    
* Simple functions may need only a one-line docstring.

Narrative Documentation
=======================

* Each subpackage must have narrative documentation describing its
  use.
