**************************
PlasmaPy Coding Guidelines
**************************

This document describes the coding requirements and guidelines to be
followed during the development of PlasmaPy and affiliated packages.

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
    any departures from the PEP 8 style guide.

  * PEP 8 compliance may be checked locally using the pep8 package.

  * Departures from PEP 8 compliance should be used sparingly and only
    if there is a good reason.  A physics formula might be most
    readable if the line exceeds the 79 character limit, for example,
    if there are inconveniently placed parentheses that complicated
    indenting.  However, departures from PEP 8 compliance should be
    considered a last resort.

* Follow the existing coding style within a subpackage.  

* Use standard abbreviations for imported packages when possible, such
  as ``import numpy as np``, ``import matplotlib as mpl``, ``import
  matplotlib.pyplot as plt``, and ``import astropy.units as u``.

* ``__init__.py`` files for modules should not contain any significant
  implementation code, but it can contain a docstring describing the
  module and code related to importing the module.  Any substantial
  functionality should be put into a separate file.g
  
* There should be at most one pun per 1284 lines of code.

Commit Messages
===============

From `How to Write a Git Commit Message
<https://chris.beams.io/posts/git-commit/>`_:

* Separate subject from body with a blank line

* Limit the subject line to 50 characters

* Capitalize the subject line

* Do not end the subject line with a period

* Use the imperative mood in the subject line

* Wrap the body at 72 characters

* Use the body to explain what and why vs. how
  
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

Units
=====

* Code within PlasmaPy must use SI units to minimize the chance of
  ambiguity, and for consistency with the recognized international
  standard.  Physical formulae and expressions should be in base SI
  units.

  * Functions should not accept floats when an Astropy Quantity is
    expected.  In particular, functions should not accept floats and
    make the assumption that the value will be in SI units.  

  * A common convention among plasma physicists is to use
    electron-volts (eV) as a unit of temperature.  Strictly speaking,
    this unit corresponds not to temperature but is rather a measure
    of the thermal energy per particle.  Code within PlasmaPy must use
    the kelvin (K) as the unit of temperature to avoid unnecessary
    ambiguity.

* PlasmaPy uses the astropy.units package to give physical units to
  values.  

  * All units packages available in Python presently have some
    limitations, including incompatibility with some NumPy and SciPy
    functions.  These limitations are due to issues within NumPy
    itself.  Many of these limitations are being resolved, but require
    upstream fixes.

* Dimensionless units may be used when appropriate, such as for
  certain numerical simulations.  The conventions and normalizations
  should be clearly described in docstrings.

Equations and Physical Formulae
===============================

* If a quantity has several names, then the function name should be
  the one that provides the most physical insight into what the
  quantity represents.  For example, ``gyrofrequency`` indicates
  gyration, whereas ``Larmor_frequency`` indicates that this frequency
  is somehow related to someone named Larmor.  Similarly, using
  ``omega_ce`` as a function name will make the code less readable to
  people who are unfamiliar with this particular notation.

* Physical formulae should be inputted without first evaluating all of
  the physical constants.  For example, the following line of code
  obscures information about the physics being represented:

>>> omega_ce = 1.76e7*(B/units.G)*units.rad/units.s

  In contrast, the following line of code shows the exact formula
  which makes the code much more readable.

>>> omega_ce = (e * B) / (m_e * c)

  The origins of numerical coefficients in formulae should be
  documented.

* Docstrings should describe the physics associated with these
  quantities in ways that are understandable to students who are
  taking their first course in plasma physics while still being useful
  to experienced plasma physicists.

* SI units that were named after a person should not be capitalized
  except at the beginning of a sentence.

Angular Frequencies
===================

Unit conversions involving angles must be treated with care.  Angles
are dimensionless but do have units.  Angular velocity is often given
in units of radians per second, though dimensionally this is
equivalent to inverse seconds.  Astropy will treat radians
dimensionlessly when using the ``dimensionless_angles`` equivalency,
but ``dimensionless_angles`` does not account for the multiplicative
factor of ``2*pi`` that is used when converting between frequency (1 /
s) and angular frequency (rad / s).  An explicit way to do this
conversion is to set up an equivalency between cycles/s and Hz:

>>> from astropy import units
>>> f_ce = omega_ce.to(units.Hz, equivalencies=[(units.cy/units.s, units.Hz)])

However, ``dimensionless_angles`` does work when dividing a velocity
by an angular frequency to get a length scale:

>>> d_i = (c/omega_pi).to(units.m, equivalencies=units.dimensionless_angles())




