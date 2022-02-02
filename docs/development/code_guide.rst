.. _code-development-guidelines:

************
Coding Guide
************


Consistency

The particular coding style does not matter

.. important::

   Consistency improves the readability and understandability of code.

The specific conventions do not matter as much as consistency.

Names
=====

* If a quantity has several names, then the function name should be
  the one that provides the most physical insight into what the
  quantity represents.  For example, ``gyrofrequency`` indicates
  gyration, whereas ``Larmor_frequency`` indicates that this frequency
  is somehow related to someone named Larmor.  Similarly, using
  ``omega_ce`` as a function name will make the code less readable to
  people who are unfamiliar with this particular notation.

* Except as described below, use :pep:`8` conventions for naming
  variables, functions, classes, and constants.

  - Use lower case words separated by underscores for function and
    variable names (e.g., ``function_name`` and ``variable_name``).

  - Use capitalized words without separators when naming a :term:`class`
    or an :term:`exception` (e.g., ``ClassName`` or ``ExceptionName``).
    However, keep acronyms capitalized (e.g., ``MHDEquations``).

  - Use capital letters words separated by underscores for constants
    (e.g., ``CONSTANT`` or ``CONSTANT_NAME``).

* Use a capital letter for a :term:`parameter` when it
  (e.g., ``B`` for magnetic field strength and ``T`` for temperature).

* Append ``_e`` to the name of a :term:`parameter` to indicate that it
  refers to electrons and ``_i`` to indicate that it refers to ions
  (e.g., ``T_e`` and ``T_i``).

* Avoid non-ASCII characters in code that is part of the public API.

* Avoid potentially ambiguous names such as ``temp`` and ``t``.

* Avoid unnecessary abbreviations

.. tip::

   Measure the length of a variable not by the number of characters, but
   rather by the time needed to understand what the variable means.

Imports
=======

* Use absolute imports, such as

* Do not use star imports (e.g., ``from package.subpackage import *``)
  because

* Use standard abbreviations for imported packages.

  .. code-block::

     import numpy as np
     import astropy.units as u
     import matplotlib as mpl
     import matplotlib.pyplot as plt

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
