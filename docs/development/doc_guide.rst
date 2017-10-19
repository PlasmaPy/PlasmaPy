************************
Documentation Guidelines
************************

This document describes the documentation requirements and guidelines
to be followed during the development of PlasmaPy and affiliated
packages.

Docstrings
==========

* All public classes, methods, and functions should have docstrings.

* PlasmaPy uses the `numpydoc
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  standard for docstrings.

* Docstrings should generally be raw string `literals
  <https://docs.python.org/3/reference/lexical_analysis.html#literals>`_:
  
.. code-block:: python

   r""" I did not like unstable eigenfunctions at first, but then
   they grew on me.
   
   """
    
* Simple functions may need only a one-line docstring.

Narrative Documentation
=======================

* Each subpackage must have narrative documentation describing its
  use.
