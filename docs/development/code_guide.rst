.. _code-development-guidelines:

************
Coding Guide
************

This page describes a set of common conventions, guidelines, and
strategies for contributing code to PlasmaPy. The goal of this guide is
not so much to enforce a common coding style, but rather to provide a
framework that helps us work together as a community. Having a
consistent coding style makes it less difficult to read code written by
a range of contributors. The particulars of the coding style do not
matter as much as consistency.

These guidelines are not rigid. There will be times when it is better to
partially rather than completely follow these guidelines. For example,
PlasmaPy usually follows the :pep:`8` style guide for Python code. This
guide indicates that variable names should be lowercase, with words
separated by underscores. However, the magnetic field is generally
denoted in uppercase as ``B`` to follow the conventions used in plasma
science.

.. important::

   These guidelines can (and should!) be refined by the community as we
   collectively learn new practices and our shared coding style changes.
   To propose a minor change to these guidelines, please submit a pull
   request. To propose a significant change, please raise an issue
   and/or bring this up at a community meeting for a discussion first.

Automatic code formatters
=========================

PlasmaPy makes heavy use of automatic code formatters such as black_ and
isort_. Using these tools helps us maintain a common code style without
having to spend excessive time worrying about the formatting of a
particular chunk of code or the order of import statements.

The use of automatic code formatters can be considered analogous to the
effects of using the `Imperius Curse <>`__ described in the novel *Harry
Potter and the Something or Other*:

   It was the most wonderful feeling. Harry felt a floating sensation as
   every thought and worry in his head was wiped gently away, leaving
   nothing but a vague, untraceable happiness. He stood there feeling
   immensely relaxed, only dimly aware of everyone watching him.

Names
=====

* Except as described below, use :pep:`8` conventions for naming
  variables, functions, classes, and constants.

  - Use lowercase words separated by underscores for function and
    variable names (e.g., ``function_name`` and ``variable_name``).

  - Use capitalized words without separators when naming a :term:`class`
    or an exception (e.g., ``ClassName`` or ``ExceptionName``). However,
    keep acronyms capitalized (e.g., ``MHDEquations``).

  - Use capital letters words separated by underscores for constants
    (e.g., ``CONSTANT`` or ``CONSTANT_NAME``).

* Use a capital letter for a :term:`parameter` when it matches the
  standard usage in plasma science.  For example, use ``B`` for magnetic
  field strength and ``T`` for temperature.

* Append ``_e`` to the name of a :term:`parameter` to indicate that it
  refers to electrons and ``_i`` to indicate that it refers to ions
  (e.g., ``T_e`` and ``T_i``).

* Python allows alphanumeric Unicode characters to be used in variable
  names (e.g., ``πλάσμα`` or ``φυσική``). These characters may be used
  for internal code when doing so improves readability (i.e. to match a
  symbol used in a paper or a standard symbol). Because non-ASCII
  characters are often difficult to enter on a keyboard, they should be
  avoided in sections of code that are under active development by
  multiple contributors. However, do not include non-ASCII characters in
  code that is part of the public API.

* If a quantity has several names, then the function name should be
  the one that provides the most physical insight into what the
  quantity represents.  For example, ``gyrofrequency`` indicates
  gyration, whereas ``Larmor_frequency`` indicates that this frequency
  is somehow related to someone named Larmor.  Similarly, using
  ``omega_ce`` as a function name will make the code less readable to
  people who are unfamiliar with this particular notation.

* Use names that are pronounceable and searchable.

* Avoid potentially ambiguous names such as ``temp`` and ``t``.

* To mark that an object is not part of PlasmaPy's public API, begin its
  name with a leading underscore (e.g., ``_private_variable``. Private
  variables should not be included in ``__all__``.

* In most situations, avoid single character variable names. Single
  character variable names may be used for standard plasma physics
  symbols (i.e., ``B``) or as indices in `for` loops (though more
  descriptive names are preferred).

* Intermediate variable names can provide additional context and
  meaning. For example, suppose we have a conditional operating on a
  complicated expression:

  .. code-block:: python

     if u[0] < x < u[1] and v[0] < y < v[1] and w[0] < z < w[1]: ...

  Defining an intermediate variable allows us to communicate the meaning
  and intent of the expression.

  .. code-block:: python

     point_in_grid_cell = u[0] < x < u[1] and v[0] < y < v[1] and w[0] < z < w[1]

     if point_in_grid_cell: ...

* Avoid unnecessary abbreviations, as these can make code more difficult
  to read. Clarity is more important than brevity, except for code that
  is frequently used interactively.

.. tip::

   Measure the length of a variable not by the number of characters, but
   rather by the time needed to understand its meaning.

Imports
=======

* Use absolute imports (e.g., ``from plasmapy.particles import Particle``).
  Relative imports (e.g., ``from ..particles import Particle``) are
  not recommended because they

* Avoid using star imports

Do not use star imports (e.g., ``from package.subpackage import *``)
  because

* Use standard abbreviations for imported packages.

  .. code-block::

     import numpy as np
     import astropy.units as u
     import astropy.constants as const
     import matplotlib.pyplot as plt
     import numba as nb
     import xarray as xr
     import pandas as pd

Units
=====

* PlasmaPy uses |astropy.units| to give physical units to values in the
  form of a |Quantity|.

* Use SI units within PlasmaPy, except when there is a strong
  justification to do otherwise.

  * Example notebooks may use non-SI units infrequently.

* Avoid using electron-volts as a unit of temperature within PlasmaPy,
  but allow arguments provided to a function

* Do not capitalize the names of units except at the beginning of a
  sentence, including when they are named after a person (except for
  "degree Celsius").

* Use operations between |Quantity| objects except when needed for
  performance.

.. _performance tips: https://docs.astropy.org/en/stable/units/index.html#performance-tips

Equations and physical formulae
===============================

* Physical formulae should be inputted without first evaluating all of
  the physical constants.  For example, the following line of code
  obscures information about the physics being represented:

>>> omega_ce = 1.76e7*(B/u.G)*u.rad/u.s   # doctest: +SKIP

  In contrast, the following line of code shows the exact formula
  which makes the code much more readable.

>>> omega_ce = (e * B) / (m_e * c)       # doctest: +SKIP

  The origins of numerical coefficients in formulae should be
  documented.

Temperature/energy equivalency
------------------------------



Comments
========



* Remove commented out code before merging a pull request.

Error messages
==============




Coding style
============

* Do not include any significant implementation code in
  :file:`__init__.py` files. Put any substantial functionality into a
  separate file.

* Use the `property` :term:`decorator` instead of getters and setters.

* Use formatted string literals (f-strings) instead of legacy formatting
  for strings.

  >>> package_name = "PlasmaPy"
  >>> print(f"The name of the package is {package_name}.")
  The name of the package is PlasmaPy.
  >>> print(f"{package_name=}")  # Python 3.8+ debugging shortcut
  package_name='PlasmaPy'
  >>> print(f"{package_name!r}")  # shortcut for f"{repr(package_name)}"
  'PlasmaPy'

* Do not use :term:`mutable` objects as default values in the function
  or method declaration. This can lead to unexpected behavior.

  .. code:: pycon

     >>> def function(l=[]):
     ...    l.append("x")
     ...    print(l)
     >>> function()
     ['x']
     >>> function()
     ['x', 'x']

* Limit usage of `lambda` functions to one-liners. For anything longer
  than that, use a nested function.

* Some plasma parameters depend on more than one |Quantity| of the same
  physical type. For example, when reading the following line of code,
  we cannot tell which is the electron temperature and which is the ion
  temperature without going to the function itself.

  .. code-block:: pycon

     f(1e6 * u.K, 2e6 * u.K)

  Spell out the :term:`parameter` names to improve readability and
  reduce the likelihood of errors.

  .. code-block::

     f(T_i = 1e6 * u.K, T_e = 2e6 * u.K)

  Similarly, when a function has parameters named ``T_e`` and ``T_i``,
  these parameters should be make :term:`keyword-only`.

.. note::

   Add the license for the google style guide, maybe?


Dependencies
============

* Follow the
