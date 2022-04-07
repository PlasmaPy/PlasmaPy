.. _code-development-guidelines:

************
Coding Guide
************

This guide describes common conventions, guidelines, and strategies for
contributing code to PlasmaPy. Having a shared coding style makes it
easier to understand code written by a range of contributors. The
purpose of this guide is not to provide a set of rigid guidelines that
must be adhered to, but rather to provide a common framework that helps
us to develop the package together as a community. The particulars of
the coding style are not as important as readability, maintainability,
and consistency. There will be times when it is better to partially
rather than completely follow these guidelines.

This guide can (and should!) be regularly refined by the PlasmaPy
community as we collectively learn new practices and our shared coding
style changes. Revisions to this guide can be proposed by submitting a
pull request or bringing up an idea at a community meeting.

Integrated development environments
===================================

.. Move this section to a page on getting set up to contribute?

An `integrated development environment`_ (IDE) is an application that
includes multiple tools needed for software development, such as a
source code editor, test automation tools, and a debugger. Commonly used
IDEs for Python include PyCharm_, `Visual Studio Code`_, and Atom_.

IDEs usually take some time and effort to learn, but eventually make
code development much easier. If are learning how to make a contribution
to an open source project for the first time, you might find it easier
to use a plain text editor that you are familiar with (e.g., Notepad++,
Sublime Text, emacs, or vi/vim) fo the time being. Alternatively, you
can do a web search for documentation or videos on :samp:`how to
contribute to a GitHub project with {name of IDE}`. In the long run, we
highly recommend taking the time to learn how to use an IDE.

.. _Atom: https://atom.io
.. _integrated development environment: https://en.wikipedia.org/wiki/Integrated_development_environment
.. _PyCharm: https://www.jetbrains.com/pycharm
.. _Visual Studio Code: https://code.visualstudio.com

Automatic code formatters
=========================

.. Move this to a page on getting set up to contribute?

PlasmaPy uses the pre-commit_ framework to perform validations and
automatically apply a consistent style to code contributions. Using
pre-commit_ helps us find errors and shortens code reviews. PlasmaPy's
pre-commit suite includes hooks such as:

* ``check-ast`` to verify that the Python code is valid.
* ``trailing-whitespace`` to remove trailing whitespace.
* black_ to format code.
* isort_ to sort imports.
* nbqa_ to format notebooks.

Most of the changes required by pre-commit_ can be applied
automatically. To apply these changes in a pull request, add a comment
that says ``pre-commit.ci autofix``. After doing this, be sure to `pull
the changes`_ from GitHub to your computer with ``git pull``.

.. _pull the changes: https://docs.github.com/en/get-started/using-git/getting-changes-from-a-remote-repository#pulling-changes-from-a-remote-repository

To enable pre-commit_ locally, open a terminal, enter the directory of
the PlasmaPy repository, and run:

.. code-block:: bash

   pip install pre-commit
   pre-commit install

Now suppose we added some trailing whitespace to :file:`some_file.py`
and attempted to commit it. If pre-commit_ has been installed, then the
``trailing-whitespace`` hook will cause pre-commit_ to fail while
modifying :file:`some_file.py` to remove the trailing whitespace.

.. code-block:: console

   $ git add some_file.py
   $ git commit -m "Add trailing whitespace"
   Trim Trailing Whitespace.................................................Failed
   - hook id: trailing-whitespace
   - exit code: 1
   - files were modified by this hook

At this point it will be necessary to run these two commands again to
commit the changes. The changes made by pre-commit_ will be unstaged and
thus could be seen by running ``git diff``. Sometimes pre-commit_ will
not be able to automatically fix the files, such as when there are
syntax errors in Python code. In these cases, the files will need to be
changed manually before running the ``git add`` and ``git commit``
commands again. Alternatively, the pre-commit_ hooks can be skipped
using :samp:`git commit --no-verify"` instead.

The pre-commit_ configuration is given in |.pre-commit-config.yaml|_.
After adding or updating pre-commit_ hooks, run the following command to
apply the changes to all files.

.. code-block:: bash

   pre-commit run --all-files

Names
=====

Names are the most fundamental means of communicating the intent and
purpose of code.

* Use :pep:`8` conventions for naming variables, functions, classes, and
  constants (except as described later in this section).

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

* PlasmaPy uses isort_ to sort import statements via a pre-commit_ hook.

* Use absolute imports (e.g., ``from plasmapy.particles import
  Particle``) rather than relative imports (e.g., ``from ..particles
  import Particle``).

* Avoid using star imports (e.g., ``from package.subpackage import *``)
  except in special situations.

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

  .. code-block:: pycon

     >>> import astropy.units as u
     >>> 5 * u.m / u.s
     <Quantity 5. m / s>

* Use SI units within PlasmaPy, except if there is strong justification
  otherwise.

  * Example notebooks may occasionally use non-SI units.

* Use |Unit| annotations for arguments and return values.

  .. code-block::

     @validate_quantities(
         n={"can_be_negative": False},
         validations_on_return={"equivalencies": u.dimensionless_angles()},
     )
     @particles.particle_input(require="charged")
     def inertial_length(n: u.m ** -3, ...) -> u.m:
         ...

* Use the `~plasmapy.utils.decorators.validators.validate_quantities`
  decorator to validate |Quantity| arguments and return values.

  .. code-block:: python

     from plasmapy.utils.decorators.validators import validate_quantities

     @validate_quantities

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
  the physical constants. For example, the following line of code
  obscures information about the physics being represented:

>>> omega_ce = 1.76e7*(B/u.G)*u.rad/u.s   # doctest: +SKIP

  In contrast, the following line of code shows the exact formula
  which makes the code much more readable.

>>> omega_ce = B * e / (m_e * c)       # doctest: +SKIP

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
