**************************
PlasmaPy Coding Guidelines
**************************

This document describes the coding requirements and guidelines to be
followed during the development of PlasmaPy and any affiliated
packages.

Dependencies
============

* Code written for PlasmaPy must be compatible with Python 3.6 and
  later.

* PlasmaPy requires Astropy 2.0 or later, NumPy 1.13 or later, and
  matplotlib 2.0 or later.

Coding Style
============

* PlasmaPy follows the `PEP8 Style Guide for Python Code
  <http://www.python.org/dev/peps/pep-0008/>`_.  This style choice
  helps ensure that the code will be consistent and readable.

* The PEP 8 Speaks integration on GitHub will comment when there are
  any departures from the PEP 8 style guide.  PEP 8 compliance may be
  checked locally using the pep8 package.

* Departures from PEP 8 compliance should be used sparingly and only
  if there is a good reason.  A physics formula might be most readable
  if the line exceeds the 79 character limit, for example, if there
  are inconveniently placed parentheses that complicated indenting.
  However, departures from PEP 8 compliance should be considered a
  last resort.

* Follow the existing coding style within a subpackage.  

* Use standard abbreviations for imported packages when possible, such
  as ``import numpy as np``, ``import matplotlib as mpl``, ``import
  matplotlib.pyplot as plt``, and ``import astropy.units as u``.

* ``__init__.py`` files for modules should not contain any significant
  implementation code, but it can contain a docstring describing the
  module and code related to importing the module.  Any substantial
  functionality should be put into a separate file.g
  
* There should be at most one pun per 1284 lines of code.

Documentation
=============

* All public classes, methods, and functions should have docstrings
  using the numpydoc format.

* These docstrings should include usage examples.

Testing
=======

* Unit tests should be provided for all methods when possible.

* Bugs should be turned into test cases.  
  
* The Travis CI integration on GitHub runs tests whenever pull
  requests to PlasmaPy are updated.  The pytest module may used on a
  local computer.
  
* Tests are run frequently during code development, and slow tests may
  interrupt the flow of a contributor.  Tests should be efficient
  except as needed.

Warnings and Exceptions
=======================

* Debugging can be intensely frustrating when problems arise and the
  associated error messages do not provide useful information on the
  source of the problem.  Warnings and error messages must be helpful
  enough for new users to quickly understand any problems that arise.

* "Errors should never pass silently."  Users should be notified when
  problems arise by either issuing a warning or raising an exception.

* The exceptions raised by a method should be described in the
  method's docstring.  Documenting exceptions makes it easier for
  future developers to plan exception handling.
