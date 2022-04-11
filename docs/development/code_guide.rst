.. _code-development-guidelines:

************
Coding Guide
************

This guide describes common conventions, guidelines, and strategies for
contributing code to PlasmaPy. The purpose of this guide is not to
provide a set of rigid guidelines that must be adhered to, but rather to
provide a common framework that helps us develop PlasmaPy together as a
community.

Having a shared coding style makes it easier to understand code written
by multiple contributors. The particulars of the coding style are not as
important as readability, maintainability, and consistency.

This guide can (and should!) be regularly refined by the PlasmaPy
community as we collectively learn new practices and our shared coding
style changes. Please feel free to propose revisions to this guide by
submitting a pull request or bringing up an idea at a community meeting.

Python resources
================

.. This section could be moved to a separate page on resources.

* `Python's documentation`_ is the ``GOTO`` place to look up different
  aspects of the Python_ language.

* The `Python Tutorial`_ provides a thorough introduction to the
  language.

* The `Python Standard Library`_ contains numerous components. Some of
  the ones used in PlasmaPy are:

  - `collections` — Container datatypes
  - `collections.abc` — Abstract base classes
  - `contextlib` — Utilities for ``with`` statement contexts
  - `dataclasses` — Data classes
  - `functools` — Higher-order functions and operations on callable objects
  - `glob` — Unix style pathname pattern expansion
  - `inspect` — Inspect live objects
  - `itertools` — Functions creating iterators for efficient looping
  - `json` — JSON encoder and decoder
  - `numbers` — Numeric abstract base classes, often used for type hint annotations
  - `os` — Miscellaneous operating system interfaces
  - `pdb` — The Python debugger
  - `re` — Regular expression operations
  - `string` — Common string operations
  - `sys` — System-specific parameters and functions
  - `typing` — Support for type hints
  - `warnings` — Warning control

.. _Python Standard Library: https://docs.python.org/3/library/
.. _Python Tutorial: https://docs.python.org/3/tutorial/index.html

Integrated development environments
===================================

.. Move this section to page on getting set up to contribute

An `integrated development environment`_ (IDE) is an application that
includes multiple tools needed for software development, such as a
source code editor, test automation tools, and a debugger. Commonly used
IDEs for Python include PyCharm_, `Visual Studio Code`_, and Atom_.

IDEs usually take some time and effort to learn, but eventually make
code development much easier. If are learning how to make a contribution
to an open source project for the first time, you might find it easier
to use a plain text editor that you are familiar with (e.g., Notepad++,
Sublime Text, emacs, or vi/vim) for the moment. Alternatively, you can
find tutorials or videos online for how to contribute to a GitHub
project with most common IDEs. In the long run, taking the time to learn
how to use an IDE is well worth it.

Pre-commit
==========

.. Move this section to page on getting set up to contribute

PlasmaPy uses the |pre-commit|_ framework to perform validations and
automatically apply a consistent style to code contributions. Using
|pre-commit|_ helps us find errors and shortens code reviews. PlasmaPy's
pre-commit suite includes hooks such as:

* ``check-ast`` to verify that the Python code is valid.
* ``trailing-whitespace`` to remove trailing whitespace.
* black_ to format code.
* isort_ to sort imports.
* nbqa_ to format notebooks.

Most of the changes required by |pre-commit|_ can be applied
automatically. To apply these changes in a pull request, add a comment
that says ``pre-commit.ci autofix``. After doing this, be sure to `pull
the changes`_ from GitHub to your computer with ``git pull``.

To enable |pre-commit|_ locally, open a terminal, enter the directory of
the PlasmaPy repository, and run:

.. code-block:: bash

   pip install pre-commit
   pre-commit install

Now suppose we added some trailing whitespace to :file:`some_file.py`
and attempted to commit it. If |pre-commit|_ has been installed, then
the ``trailing-whitespace`` hook will cause |pre-commit|_ to fail while
modifying :file:`some_file.py` to remove the trailing whitespace.

.. code-block:: console

   $ git add some_file.py
   $ git commit -m "Add trailing whitespace"
   Trim Trailing Whitespace.................................................Failed
   - hook id: trailing-whitespace
   - exit code: 1
   - files were modified by this hook

At this point it will be necessary to run these two commands again to
commit the changes. The changes made by |pre-commit|_ will be unstaged and
thus could be seen by running ``git diff``. Sometimes |pre-commit|_ will
not be able to automatically fix the files, such as when there are
syntax errors in Python code. In these cases, the files will need to be
changed manually before running the ``git add`` and ``git commit``
commands again. Alternatively, the |pre-commit|_ hooks can be skipped
using ``git commit --no-verify`` instead.

The |pre-commit|_ configuration is given in |.pre-commit-config.yaml|_.

After adding or updating |pre-commit|_ hooks, run the following command to
apply the changes to all files.

.. code-block:: bash

   pre-commit run --all-files

Names
=====

Names are our most fundamental means of communicating the intent and
purpose of code. Judicious choices of names can greatly improve the
understandability of code, while inadequate naming can obfuscate what
the code is supposed to be doing.

* Most IDEs have a built-in tool for simultaneously renaming all
  instances of a variable throughout a project, which can save a lot of
  time. For example, a `rename refactoring in PyCharm`_ can be done with
  :kbd:`Shift+F6` on Windows or Linux and :kbd:`⇧F6` or :kbd:`⌥⌘R` on
  macOS.

* Use :pep:`8` conventions for naming variables, functions, classes, and
  constants (except as described later in this list).

  - Use lowercase words separated by underscores for function and
    variable names (e.g., ``function_name`` and ``variable_name``).

  - Use capitalized words without separators when naming a :term:`class`
    or an exception (e.g., ``ClassName`` or ``ExceptionName``). However,
    keep acronyms capitalized (e.g., ``MHDEquations``).

  - Use capital letters words separated by underscores for constants
    (e.g., ``CONSTANT`` or ``CONSTANT_NAME``).

* Use a capital letter for a :term:`parameter` when it matches the
  standard usage in plasma science. For example, use ``B`` for magnetic
  field strength and ``T`` for temperature.

* Functions based on plasma parameters that are named after people may
  be capitalized (e.g., ``Alfven_speed`` and ``Debye_length``).

* Append ``_e`` to the name of a :term:`parameter` to indicate that it
  refers to electrons, ``_i`` to indicate that it refers to ions, and
  ``_p`` to indicate that it refers to protons (e.g., ``T_e`` and
  ``T_i``, and ``T_p``).

* Only ASCII_ characters should be used in code that is part of the
  public API_.

* Python allows alphanumeric Unicode characters to be used in object
  names (e.g., ``πλάσμα`` or ``φυσική``). These characters may be used
  for internal code when doing so improves readability (i.e. to match a
  symbol used in a paper or a standard symbol). Because non-\ ASCII_
  characters are often difficult to enter on a keyboard, they should be
  avoided in sections of code that are under active development by
  multiple contributors.

* If a plasma parameter has multiple names, then the function name
  should be the one that provides the most physical insight into what
  the quantity represents. For example, ``gyrofrequency`` indicates
  gyration, whereas ``Larmor_frequency`` indicates that this frequency
  is somehow related to someone named Larmor.

* Avoid naming functions by spelling out the name of the Greek
  character, as in

  * Similarly, using ``omega_ce`` as a function name will make the code
    less readable to people who are unfamiliar with this particular
    notation.

* Choose names that are pronounceable so that they are easier to
  remember and more compatible with screen reader (text-to-speech)
  technology.

* Choose names that are searchable (i.e. doing a web search for a name
  results in helpful

* Avoid potentially ambiguous names such as ``temp`` and ``t``.

* To mark that an object is not part of PlasmaPy's public API_, begin its
  name with a leading underscore (e.g., ``_private_variable``. Private
  variables should not be included in ``__all__``.

* In most situations, avoid single character variable names. Single
  character variable names may be used for standard plasma physics
  symbols (i.e., ``B``) or as indices in ``for`` loops (though more
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

* In general, avoid encoding type information in a variable name.

* Avoid unnecessary abbreviations, as these can make code more difficult
  to read. Clarity is more important than brevity, except for code that
  is used frequently and interactively (e.g., :command:`ls` or
  :command:`cd` in the Unix shell).

.. tip::

   Measure the length of a variable not by the number of characters, but
   rather by the time needed to understand its meaning.

   By this measure, |cggglm|_ is significantly longer than
   ``solve_gauss_markov_linear_model``.

.. _cggglm: http://www.netlib.org/lapack/explore-html/d9/d98/group__complex_o_t_h_e_reigen_ga4be128ffc05552459683f0aade5a7937.html
.. |cggglm| replace:: ``cggglm``

Imports
=======

* PlasmaPy uses isort_ to sort import statements via a |pre-commit|_
  hook.

* Use absolute imports (e.g., ``from plasmapy.particles import Particle``)
  rather than relative imports (e.g., ``from ..particles import Particle``).

* Avoid using star imports (e.g., ``from package.subpackage import *``)
  except in special situations.

* Importing a package, subpackage, or module rather than an individual
  code object has the benefit that the namespace provides helpful
  contextual information that can make code more understandable. For
  example, using ``json.loads`` is more understandable than using only
  ``loads``.

  For frequently used objects (e.g., |Particle|), using the full
  namespace will increase the clutter of the code without providing
  commensurately more information. This is also true for objects used as
  type hint annotations. For example, ``Optional[Union[Real, Complex]``
  is more understandable than
  ``typing.Optional[typing.Union[numbers.Real, numbers.Complex]]``.

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

* PlasmaPy uses |astropy.units|_ to give physical units to values in the
  form of a |Quantity|.

  .. code-block:: pycon

     >>> import astropy.units as u
     >>> 5 * u.m / u.s
     <Quantity 5. m / s>

  Using |astropy.units| improves compatibility with Python packages in
  adjacent fields such as astronomy and heliophysics.

* Use SI units within PlasmaPy, unless there is a strong justification
  to do otherwise. Example notebooks may use other unit systems.

* Use |Unit| annotations with the |validate_quantities| decorator to
  validate |Quantity| arguments and return values.

  .. code-block:: python

     from plasmapy.utils.decorators.validators import validate_quantities

     @validate_quantities(
        n={"can_be_negative": False},
        validations_on_return={"equivalencies": u.dimensionless_angles()},
     )
     def inertial_length(n: u.m ** -3, ...) -> u.m:
         ...

* Avoid using electron-volts as a unit of temperature within PlasmaPy
  because it is defined as a unit of energy. However, functions in
  `plasmapy.formulary` and elsewhere should accept temperatures in units
  of electron-volts, which can be done using |validate_quantities|.

* Non-standard unit conversions can be made using equivalencies_ such
  as `~astropy.units.temperature_energy`.

  .. code-block:: pycon

     >>> (1 * u.eV).to(u.K, equivalencies=u.temperature_energy())
     11604.518...

* Do not capitalize the names of units except at the beginning of a
  sentence, including when they are named after a person (except for
  "degree Celsius").

* Use operations between |Quantity| objects except when needed for
  performance.

* To improve performance in |Quantity| operations,

Particles
=========

* The |Particle| class provides an object-oriented interface for
  accessing basic particle data. |Particle| accepts
  :term:`particle-like` inputs.

  .. code-block:: pycon

     >>> from plasmapy.particles import Particle
     >>> alpha = Particle("He-4 2+")
     >>> alpha.mass
     <Quantity 6.6446...e-27 kg>
     >>> alpha.charge
     <Quantity 3.20435...e-19 C>

* The |particle_input| decorator can automatically transform a
  :term:`particle-like` :term:`argument` into a |Particle| instance, if
  the corresponding :term:`parameter` is decorated with |Particle|.

  .. code-block::

     from plasmapy.particles import particle_input, Particle

     @particle_input
     def recombine(ion: Particle):
          # ion is now a Particle instance
          return ion.recombine()

  The documentation for |particle_input| describes ways to ensure that
  the particle meets certain categorization criteria.

Equations and physical formulae
===============================

* Physical formulae should be inputted without first evaluating all of
  the physical constants. For example, the following line of code
  obscures information about the physics being represented:

  .. code-block:: pycon

     >>> omega_ce = 1.76e7*(B/u.G)*u.rad/u.s   # doctest: +SKIP

  In contrast, the following line of code shows the exact formula
  which makes the code much more readable.

  .. code-block:: pycon

     >>> omega_ce = B * e / (m_e * c)  # doctest: +SKIP

* The origins of numerical coefficients in formulae should generally be
  described in comments or in the docstring.

* References for equations should be included in the |bibliography|, as
  described in the |documentation guide|.

Coding style
============

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
     ...     l.append("x")
     ...     print(l)
     >>> function()
     ['x']
     >>> function()
     ['x', 'x']

* Use the `property` :term:`decorator` instead of getters and setters.

* Only use `lambda` functions for one-liners that are only used near
  where they are defined (e.g., when defining the default factory for a
  `~collections.defaultdict`). For anything longer than one line, define
  a define a function with ``def`` instead.

* Some plasma parameters depend on more than one |Quantity| of the same
  physical type. For example, when reading the following line of code,
  we cannot tell which is the electron temperature and which is the ion
  temperature without going to the function itself.

  .. code-block:: python

     f(1e6 * u.K, 2e6 * u.K)

  Spell out the :term:`parameter` names to improve readability and
  reduce the likelihood of errors.

  .. code-block:: python

     f(T_i = 1e6 * u.K, T_e = 2e6 * u.K)

  Similarly, when a function has parameters named ``T_e`` and ``T_i``,
  these parameters should be make :term:`keyword-only` to avoid
  ambiguity and reduce the chance of errors.

  .. code-block::

     def f(*, T_i, T_e):
         ...

* The ``__eq__`` and ``__ne__`` methods of a class should not raise
  exceptions. If the comparison for equality is being made between
  objects of different types, these methods should return `False`
  instead. This behavior is for consistency with operations like
  ``1 == "1"`` which will return `False`.

* List and dictionary comprehensions should be used for simple ``for``
  loops, like:

  .. code-block:: pycon

     >>> [x ** 2 for x in range(17) if x % 2 == 0]
     [0, 4, 16, 36, 64, 100, 144, 196, 256]
     >>> {x: x ** 2 for x in range(17) if x % 2 == 0}
     {0: 0, 2: 4, 4: 16, 6: 36, 8: 64, 10: 100, 12: 144, 14: 196, 16: 256}

* Avoid using global variables.

* Avoid putting any significant implementation code in
  :file:`__init__.py` files. Implementation details should be contained
  in a different file, and then imported into :file:`__init__.py`.

Aliases
=======

An :term:`alias` is an abbreviated version of a commonly used function
that is intended for interactive use. For example,
`~plasmapy.formulary.speeds.va_` is an alias to
`~plasmapy.formulary.speeds.Alfven_speed`.

* Aliases should only be defined for the most commonly used functions.

* An alias should be defined immediately after the original function.

* The name of an alias should in some way indicate what the alias is
  for. For example, `~plasmapy.formulary.lengths.cwp_` is a shortcut for
  for :math:`c/ω_p`\ .

* The name of an alias should end with a trailing underscore.

* Each alias should have a one-line docstring that refers users to the
  original function.

* The name of the main function should be included in ``__all__`` near
  the top of each module, and the name of the alias should be included
  in ``__aliases__``, which will then get appended to ``__all__``.

Here is a minimal example of an alias ``f_`` that would be for
``plasmapy.subpackage.module.function``.

.. code-block:: python

   __all__ = ["function"]
   __aliases__ = ["f_"]

   __all__ += __aliases__

   def function():
       ...

   f_ = function
   """Alias to `~plasmapy.subpackage.module.function`."""

Lite Functions
==============

Most functions in `plasmapy.formulary` use |astropy.units|_ to attach
units to values in the form of a |Quantity|, and also perform checks to
make sure that each :term:`argument` that is provided to a function is
valid. The use of |Quantity| operations and validations do not have a
noticeable performance penalty during typical interactive use, but the
performance penalty can become substantial for numerically intensive
applications.

A :term:`lite-function` is a lite weight version of another `plasmapy`
function. Most lite-functions are defined in `plasmapy.formulary`.

* The name of each lite-function should be the name of the original
  function with ``_lite`` appended at the end. For example,
  `~plasmapy.formulary.speeds.thermal_speed_lite` is the lite-function
  associated with `~plasmapy.formulary.speeds.thermal_speed`.

* Lite-functions assume SI units for all of arguments that represent
  physical quantities.

* Lite-functions should be defined immediately before the normal version
  of the function.

* Lite-functions are bound to their normal version as the ``lite``
  attribute using the `~plasmapy.utils.decorators.bind_lite_func`
  decorator.

* Each lite-function should be decorated with
  `~plasmapy.utils.decorators.preserve_signature`.

* A lite-function should usually be decorated with `numba.njit` (or the
  like) as a just-in-time compiler. If a decorator from `numba` is not
  able to be used, then it might be possible to use Cython_.

The following is a minimal implementation of a lite-function.

.. code-block:: python

   __all__ = ["function"]
   __lite_funcs__ = ["function_lite"]

   from numba import njit
   from numbers import Real
   from plasmapy.utils.decorators import bind_lite_func, preserve_signature

   __all__ += __lite_funcs__

   @preserve_signature
   @njit
   def function_lite(v: Real) -> Real:
       """
       The lite-function which accepts and returns real numbers in
       assumed SI units.
       """
       ...

   @bind_lite_func(function_lite)
   def function(v):
       """A function that accepts and returns Quantity arguments."""
       ...

Comments
========

A well-placed and well-written comment can prevent future frustrations.
However, comments are not inherent good. As code evolves, an
unmaintained comment may become outdated or get separated from the
section of code that it was meant to describe. Cryptic comments may end
up confusing contributors. In the worst case, an unmaintained comment
may contain inaccurate or misleading information.

* Refactor code to make it more readable, rather than explaining how it
  works.

* When a comment is used to define a variable, try renaming the
  variable to encode its meaning and intent.

  .. code-block:: python

     # collision frequency
     nu = 1e6 * u.s ** -1

     collision_frequency = 1e6 * u.s ** -1

* Use comments to communicate information that you wish you knew before
  starting to work on a particular section of code, including
  information that took some time to learn.

* Use comments to communicate information that the code cannot,
  such as why an alternative approach was *not* taken.

* Include enough contextual information that the comment will ideally
  make sense even if it is displaced from the code it was originally
  describing.

* Use comments to include references to books or articles that describe
  the equation, algorithm, or software design pattern that is being
  implemented.

* Include enough contextual information in the comment for a new user
  to be able to understand it.

* Remove commented out code before merging a pull request.

* When updating code, be sure to update the comments too!

* When a comment is used as the header for a section of code, consider
  extracting that section of code into. For example, we
  might start out with a function that includes multiple lines of code
  for each step.

  .. code-block:: python

     def analyze_experiment(data):
         # Step 1: calibrate the data
         ...
         # Step 2: normalize the data
         ...

  We can apply the `extract function refactoring pattern`_ by creating a
  separate function for each of these steps. The name of each function
  can often be very close to the comment.

  .. code-block:: python

     def calibrate_data(data):
         ...
         return calibrated_data

     def normalized_data(data):
         ...
         return normalized_data

     def analyze_experiment(data):
         calibrated_data = calibrate_data(data)
         normalized_data = normalize_data(calibrated_data)

  This refactoring strategy is appropriate for long functions where the
  different steps can be cleanly separated from each other. This pattern
  leads to functions that are shorter, more focused, more reusable, and
  easier to test. The original function no longer includes
  implementation details, and thus gives a high level view of what the
  function is doing. This pattern might not be appropriate if the
  different sections of code are intertwined with each other, like if
  both sections use the same intermediate variables.

Error messages
==============

Error messages are a vital but underappreciated form of documentation.
A good error message can help someone pinpoint the source of a problem
in seconds, while a cryptic or missing error message can lead to hours
of frustration.

* Use error messages to indicate the source of the problem while
  providing enough information for the user to fix it. Make sure that it
  is clear what the user should do next.

* Include diagnostic information when appropriate. For example, if an
  operation is being performed on every point in a grid, include the
  coordinates where the error happened.

* Write error messages that are concise when possible, as users often
  skim long error messages.

* Avoid including information that is irrelevant to the source of the
  problem.

* Write error messages in language that is plain enough to be
  understandable to someone who is undertaking their first research
  project.

  - If necessary, technical information may be placed after a plain
    language summary statement.

  - Alternatively, an error message may reference a docstring or a page
    in the narrative documentation.

* Write error messages that are friendly, supportive, and helpful. Error
  message should never be condescending or blame the user.

Requirements
============

* Package requirements are specified in multiple locations that need to
  be updated simultaneously.

  - The |requirements|_ directory contains multiple text files that
    contain build, installation, testing, documentation, and extra
    requirements.

  - The ``build-system.requires`` section of |pyproject.toml|_ includes
    the requirements for building PlasmaPy. This section should mirror
    :file:`requirements/build.txt`.

  - |setup.cfg|_ includes sections for the install, docs, tests, and
    extra requirements that should mirror the corresponding files in
    the |requirements|_ directory.

  - :file:`requirements/environment.yml` contains a Conda_ environment
    for PlasmaPy.

  - The :file:`tox.ini` file contains a testing environment for the
    minimal dependencies.

* Each release of `plasmapy` should support all minor versions of
  Python that have been released in the prior 42 months, and all minor
  versions of `numpy` that have been released in the last 24 months.
  This schedule was proposed in `NumPy Enhancement Proposal 29`_ for
  the scientific Python ecosystem, and has been adopted by upstream
  packages such as `numpy`, `matplotlib`, and `astropy`.

  .. tip::

     Tools like pyupgrade_ help automatically upgrade the code base to
     the minimum supported version of Python for the next release.

* In general, it is preferable to support minor releases of dependencies
  from the last ≲ 24 months, unless there is a new feature in a
  dependency that would be greatly beneficial for `plasmapy` development.

* Avoid setting maximum requirements such as ``sphinx <= 2.4.4`` because
  this can lead to version conflicts when PlasmaPy is installed
  alongside other packages. Instead it is preferable to update
  `plasmapy` to be compatible with the newest versions of each of its
  dependencies.

* Minor versions of Python are generally released in October of each
  year. However, it may take a few months before packages like NumPy_
  and Numba_ become compatible with the newest minor version of Python_.

.. _ASCII: https://en.wikipedia.org/wiki/ASCII
.. _Atom: https://atom.io
.. _equivalencies: https://docs.astropy.org/en/stable/units/equivalencies.html
.. _extract function refactoring pattern: https://refactoring.guru/extract-method
.. _integrated development environment: https://en.wikipedia.org/wiki/Integrated_development_environment
.. _nbqa: https://nbqa.readthedocs.io
.. _NumPy Enhancement Proposal 29: https://numpy.org/neps/nep-0029-deprecation_policy.html
.. _performance tips: https://docs.astropy.org/en/stable/units/index.html#performance-tips
.. _pull the changes: https://docs.github.com/en/get-started/using-git/getting-changes-from-a-remote-repository#pulling-changes-from-a-remote-repository
.. _PyCharm: https://www.jetbrains.com/pycharm
.. _pyupgrade: https://github.com/asottile/pyupgrade
.. _rename refactoring in PyCharm: https://www.jetbrains.com/help/pycharm/rename-refactorings.html
.. _Visual Studio Code: https://code.visualstudio.com
