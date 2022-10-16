.. _coding guide:

************
Coding Guide
************

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Introduction
============

This guide describes common conventions, guidelines, and strategies for
contributing code to PlasmaPy. The purpose of this guide is not to
provide a set of rigid guidelines that must be adhered to, but rather to
provide a common framework that helps us develop PlasmaPy together as a
community.

Having a shared coding style makes it easier to understand code written
by multiple contributors. The particulars of the coding style are not as
important as consistency, readability, and maintainability.

This guide can (and should!) be regularly refined by the PlasmaPy
community as we collectively learn new practices and our shared coding
style changes. Please feel free to propose revisions to this guide by
:ref:`submitting a pull request <code-contribution>` or by bringing up
an idea at a community meeting.

PlasmaPy generally follows the :pep:`8` style guide for Python code,
using auto-formatters such as black_ and isort_ that are executed using
pre-commit_.

Coding guidelines
=================

* Write short functions that do exactly one thing with no side effects.

* Use NumPy_ array options instead of ``for`` loops to make code more
  compact, readable, and performant.

* Instead of defining variables like ``a0``, ``a1``, & ``a2``, define
  these values in a collection such as an |ndarray| or a `list`.

* Use the `property` :term:`decorator` instead of getters and setters.

* Some plasma parameters depend on more than one |Quantity| of the same
  physical type. For example, when reading the following line of code,
  we cannot immediately tell which is the electron temperature and which
  is the ion temperature.

  .. code-block:: python

     f(1e6 * u.K, 2e6 * u.K)

  Spell out the :term:`parameter` names to improve readability and
  reduce the likelihood of errors.

  .. code-block:: python

     f(T_i = 1e6 * u.K, T_e = 2e6 * u.K)

  Similarly, when a function has parameters named ``T_e`` and ``T_i``,
  these parameters should be made |keyword-only| to avoid ambiguity and
  reduce the chance of errors.

  .. code-block:: python

     def f(*, T_i, T_e):
         ...

* The ``__eq__`` and ``__ne__`` methods of a class should not raise
  exceptions. If the comparison for equality is being made between
  objects of different types, these methods should return `False`
  instead. This behavior is for consistency with operations like
  :py:`1 == "1"` which will return `False`.

* Limit usage of ``lambda`` functions to one-liners, such as when
  defining the default factory of a `~collections.defaultdict`). For
  anything longer than one line, use ``def`` instead.

* List and dictionary comprehensions can be used for simple ``for``
  loops, like:

  .. code-block:: pycon

     >>> [x ** 2 for x in range(17) if x % 2 == 0]
     [0, 4, 16, 36, 64, 100, 144, 196, 256]

  A comprehension might be more readable when spread out over multiple
  lines.

  .. code-block:: pycon

     >>> {
     ...     x: x ** 2
     ...     for x in range(17)
     ...     if x % 2 == 0
     ... }
     {0: 0, 2: 4, 4: 16, 6: 36, 8: 64, 10: 100, 12: 144, 14: 196, 16: 256}

* Avoid putting any significant implementation code in
  :file:`__init__.py` files. Implementation details should be contained
  in a different file, and then imported into :file:`__init__.py`.

* Avoid defining global variables when possible.

* Use ``assert`` statements only in tests.

* Use formatted string literals (f-strings) instead of legacy formatting
  for strings.

  >>> package_name = "PlasmaPy"
  >>> print(f"The name of the package is {package_name}.")
  The name of the package is PlasmaPy.
  >>> print(f"{package_name=}")
  package_name='PlasmaPy'
  >>> print(f"{package_name!r}")  # shortcut for f"{repr(package_name)}"
  'PlasmaPy'

* Functions that accept |array_like| or |Quantity| inputs should accept
  and return |nan| (`not a number`_) values. This guideline applies when
  |nan| is the input as well as when |nan| values are included in an
  array.

  .. tip::

     Normally, :py:`numpy.nan == numpy.nan` evaluates to `False`, which
     complicates testing |nan| behavior. The ``equal_nan`` keyword of
     functions like `numpy.allclose` and `numpy.testing.assert_allclose`
     makes it so that |nan| is considered equal to itself.

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

* Use `pathlib` when working with paths to data files.

Names
=====

Names are our most fundamental means of communicating the intent and
purpose of code. Wisely chosen names can greatly improve the
understandability of code, while inadequate names can obfuscate what
the code is supposed to be doing.

* PlasmaPy generally uses the :pep:`8` conventions for variable names.

  - Use lowercase words separated by underscores for function and
    variable names (e.g., ``function_name`` and ``variable_name``).

  - Use capitalized words without separators when naming a class (e.g.,
    ``ClassName``), but keep acronyms capitalized (e.g.,
    ``MHDEquations``).

  - Use capital letters words separated by underscores when naming
    constants (e.g., ``CONSTANT`` or ``CONSTANT_NAME``).

  There are some situations in PlasmaPy which justify a departure from
  the :pep:`8` conventions.

  - Functions based on plasma parameters that are named after people may
    be capitalized (e.g., ``Alfven_speed``).

  - Capital letters may be used for a variable when it matches the
    standard usage in plasma science (e.g., ``B`` for magnetic field and
    ``T`` for temperature).

* Choose names that are pronounceable to make them more memorable and
  compatible with text-to-speech technology.

* Choose names will produce more relevant results when searching the
  internet.

* Avoid unnecessary abbreviations, as these make code harder to read.
  Prefer clarity over brevity, except for code that is used frequently
  and interactively (e.g., :command:`cd` or :command:`ls`).

  .. tip::

     Measure the length of a variable not by the number of characters,
     but rather by the time needed to understand its meaning.

     By this measure, ``cggglm`` is significantly longer than
     ``solve_gauss_markov_linear_model``.

* Avoid ambiguity. Does ``temp`` mean "temperature", "temporary", or
  "template"?

* Append ``_e`` to a variable name to indicate that it refers to
  electrons, ``_i`` for ions, and ``_p`` for protons (e.g., ``T_e``,
  ``T_i``, and ``T_p``).

* Only ASCII_ characters should be used in code that is part of the
  public API_.

* Python allows alphanumeric Unicode characters to be used in object
  names (e.g., ``πλάσμα`` or ``φυσική``). These characters may be used
  for *internal* code when doing so improves readability (i.e., to match
  a commonly used symbol) and in Jupyter_ notebooks.

* If a plasma parameter has multiple names, then use the name that
  provides the most physical insight. For example, ``gyrofrequency``
  indicates gyration but ``Larmor_frequency`` does not.

* It is *usually* preferable to name a variable after its name rather
  than its symbol.  An object named ``Debye_length`` is more broadly
  understandable and searchable than ``lambda_D``. However, there are
  some exceptions to this guideline.

  * Symbols used widely across plasma science can be used with low risk
    of confusion, such as :math:`T` for temperature or :math:`β` for
    plasma `~plasmapy.formulary.dimensionless.beta`.

  * Symbols that are defined in docstrings can be used with decreased
    likelihood of confusion.

  * Sometimes code that represents an equation will be more readable if
    the Unicode characters for the symbols are used, especially for
    complex equations. For someone who is familiar with the symbols,
    ``λ = c / ν`` will be more readable than ``lambda = c / nu`` or
    ``wavelength = speed_of_light / frequency``.

  * If an implementation is based on a journal article, then variable
    names may be based on the symbols used in that article. The article
    should be :ref:`cited <citation-instructions>` in the appropriate
    docstring so that it appears in the |bibliography|.

* To mark that an object is not part of PlasmaPy's public API_, begin
  its name with a leading underscore (e.g., ``_private_variable``).
  Private variables should not be included in ``__all__``.

* Avoid single character variable names except for standard plasma
  physics symbols (e.g., ``B``) or as indices in ``for`` loops.

* Avoid encoding type information in a variable name.

* Intermediate variable names can provide additional context and
  meaning. For example, suppose we have a conditional operating on a
  complicated expression:

  .. code-block:: python

     if u[0] < x < u[1] and v[0] < y < v[1] and w[0] < z < w[1]: ...

  Defining an intermediate variable allows us to communicate the meaning
  and intent of the expression.

  .. code-block:: python

     point_is_in_grid_cell = u[0] < x < u[1] and v[0] < y < v[1] and w[0] < z < w[1]

     if point_is_in_grid_cell:
         ...

  In ``for`` loops, this may take the form of assignment expressions
  with the walrus operator (``:=``).

.. tip::

   Most `integrated development environments <IDE>`_ (IDEs) have a
   built-in tool for simultaneously renaming a variable throughout a
   project. For example, a `rename refactoring in PyCharm
   <https://www.jetbrains.com/help/pycharm/rename-refactorings.html>`__
   can be done with :kbd:`Shift+F6` on Windows or Linux, and :kbd:`⇧F6`
   or :kbd:`⌥⌘R` on macOS.

Imports
=======

* Use standard abbreviations for imported packages:

  .. code-block:: python

     import numpy as np
     import astropy.units as u
     import astropy.constants as const
     import matplotlib.pyplot as plt
     import pandas as pd

* PlasmaPy uses isort_ to sort import statements via a |pre-commit|_
  hook.

* For infrequently used objects, import the package, subpackage, or
  module rather than the individual code object. Including more of the
  namespace provides contextual information that can make code easier to
  read. For example, ``json.loads`` is more readable than using only
  ``loads``.

* For frequently used objects (e.g., |Particle|) and type hint
  annotations (e.g., `~typing.Optional` and `~numbers.Real`), import the
  object directly instead of importing the package, subpackage, or
  module. Including more of the namespace would increase clutter and
  decrease readability without providing commensurately more
  information.

* Use absolute imports (e.g., :py:`from plasmapy.particles import
  Particle`) rather than relative imports (e.g., :py:`from ..particles
  import Particle`).

* Do not use star imports (e.g., :py:`from package.subpackage import *`),
  except in very limited situations.

Requirements
============

* Package requirements are specified in multiple locations that need to
  be updated simultaneously.

  - The |requirements|_ directory contains multiple text files that
    contain build, installation, testing, documentation, and extra
    requirements.

  - The ``build-system.requires`` section of |pyproject.toml|_ includes
    the requirements for building PlasmaPy. This section must mirror
    |requirements/build.txt|_.

  - |setup.cfg|_ includes sections for the install, docs, tests, and
    extra requirements that must mirror the corresponding files in
    the |requirements|_ directory.

  - |requirements/environment.yml|_ contains a Conda_ environment
    for PlasmaPy.

  - |tox.ini|_ contains a testing environment for the minimal
    dependencies.

* Each release of PlasmaPy should support all minor versions of
  Python that have been released in the prior 42 months, and all minor
  versions of NumPy_ that have been released in the last 24 months.
  This schedule was proposed in `NumPy Enhancement Proposal 29`_ for
  the scientific Python ecosystem, and has been adopted by upstream
  packages such as NumPy_, matplotlib_, and Astropy_.

  .. tip::

     Tools like pyupgrade_ help automatically upgrade the code base to
     the minimum supported version of Python for the next release.

* In general, it is preferable to support minor releases of dependencies
  from the last ≲ 24 months, unless there is a new feature in a
  dependency that would be greatly beneficial for `plasmapy` development.

* Do not set maximum requirements (e.g., ``numpy <= 1.22.3``) unless
  absolutely necessary. Maximum requirements can lead to version
  conflicts when installed alongside other packages. Instead, update
  PlasmaPy to become compatible with the latest versions of its
  dependencies. Similarly, do not require exact versions of packages
  (e.g., ``scipy == 1.5.3``).

* Minor versions of Python are generally released in October of each
  year. However, it may take a few months before packages like NumPy_
  and Numba_ become compatible with the newest minor version of Python_.

.. _code-contribution:

Branches, commits, and pull requests
====================================

Before making any changes, it is prudent to update your local
repository with the most recent changes from the development
repository:

.. code-block:: bash

  git fetch upstream

Changes to PlasmaPy should be made using branches.  It is usually best
to avoid making changes on your main branch so that it can be kept
consistent with the upstream repository.  Instead we can create a new
branch for the specific feature that you would like to work on:

.. code-block:: bash

  git branch *your-new-feature*

Descriptive branch names such as ``grad-shafranov`` or
``adding-eigenfunction-poetry`` are helpful, while vague names like
``edits`` are considered harmful.  After creating your branch locally,
let your fork of PlasmaPy know about it by running:

.. code-block:: bash

  git push --set-upstream origin *your-new-feature*

It is also useful to configure git so that only the branch you are
working on gets pushed to GitHub:

.. code-block:: bash

  git config --global push.default simple

Once you have set up your fork and created a branch, you are ready to
make edits to PlasmaPy.  Switch to your new branch by running:

.. code-block:: bash

  git checkout *your-new-feature*

Go ahead and modify files with your favorite text editor.  Be sure to
include tests and documentation with any new functionality.  We
recommend reading about `best practices for scientific computing
<https://doi.org/10.1371/journal.pbio.1001745>`_.  PlasmaPy uses the
`PEP 8 style guide for Python code
<https://www.python.org/dev/peps/pep-0008/>`_ and the `numpydoc format
for docstrings
<https://github.com/numpy/numpy/blob/main/doc/HOWTO_DOCUMENT.rst.txt>`_
to maintain consistency and readability.  New contributors should not
worry too much about precisely matching these styles when first
submitting a pull request, GitHub Actions will check pull requests
for :pep:`8` compatibility, and further changes to the style can be
suggested during code review.

You may periodically commit changes to your branch by running

.. code-block:: bash

  git add filename.py
  git commit -m "*brief description of changes*"

Committed changes may be pushed to the corresponding branch on your
GitHub fork of PlasmaPy using

.. code-block:: bash

  git push origin *your-new-feature*

or, more simply,

.. code-block:: bash

  git push

Once you have completed your changes and pushed them to the branch on
GitHub, you are ready to make a pull request.  Go to your fork of
PlasmaPy in GitHub.  Select "Compare and pull request".  Add a
descriptive title and some details about your changes.  Then select
"Create pull request".  Other contributors will then have a chance to
review the code and offer constructive suggestions.  You can continue
to edit the pull request by changing the corresponding branch on your
PlasmaPy fork on GitHub.  After a pull request is merged into the
code, you may delete the branch you created for that pull request.


Comments
========

A well-placed and well-written comment can prevent future frustrations.
However, comments are not inherently good. As code evolves, an
unmaintained comment may become outdated, or get separated from the
section of code that it was meant to describe. Cryptic and obsolete
comments may end up confusing contributors. In the worst case, an
unmaintained comment may contain inaccurate or misleading information
(hence the saying that "a comment is a lie waiting to happen").

.. important::

   The code we write should read like a book. The full meaning of code's
   functionality should be attainable by reading the code. Comments
   should only be used when the code itself cannot communicate its full
   meaning.

* Refactor code to make it more readable, rather than explaining how it
  works :cite:p:`wilson:2014`.

* Instead of using a comment to define a variable, rename the variable
  to encode its meaning and intent.  For example, code like:

  .. code-block:: python

     # collision frequency
     nu = 1e6 * u.s ** -1

  could be achieved with no comment by doing:

  .. code-block:: python

     collision_frequency = 1e6 * u.s ** -1

* Use comments to communicate information that you wish you knew before
  starting to work on a particular section of code, including
  information that took some time to learn.

* Use comments to communicate information that the code cannot,
  such as why an alternative approach was *not* taken.

* Use comments to include references to books or articles that describe
  the equation, algorithm, or software design pattern that is being
  implemented. Even better, include these references in docstrings.

* Provide enough contextual information in the comment for a new user
  to be able to understand it.

* Remove commented out code before merging a pull request.

* When updating code, be sure to review and update, if necessary, associated comments too!

* When a comment is used as the header for a section of code, consider
  extracting that section of code into its own function. For example, we
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
  can often be extracted directly from the comment.

  .. code-block:: python

     def calibrate_data(data):
         ...
         return calibrated_data

     def normalize_data(data):
         ...
         return normalized_data

     def analyze_experiment(data):
         calibrated_data = calibrate_data(data)
         normalized_data = normalize_data(calibrated_data)

  This refactoring pattern is appropriate for long functions where the
  different steps can be cleanly separated from each other. This pattern
  leads to functions that are shorter, more reusable, and easier to
  test. The original function contains fewer low-level implementation
  details and thus gives a higher level view of what the function is
  doing. This pattern reduces `cognitive complexity`_.

  The `extract function refactoring pattern`_ should be used
  judiciously, as taking it to an extreme and applying it at too fine of
  a scale can reduce readability and maintainability by producing overly
  fragmented code.

  .. hint::

     The `extract function refactoring pattern`_ might not be
     appropriate if the different sections of code are intertwined with
     each other (e.g., if both sections require the same intermediate
     variables). An alternative in such cases would be to create a class
     instead.

Error messages
==============

Error messages are a vital but underappreciated form of documentation.
A good error message can help someone pinpoint the source of a problem
in seconds, while a cryptic or missing error message can lead to hours
of frustration.

* Use error messages to indicate the source of the problem while
  providing enough information for the user to troubleshoot it. When
  possible, make it clear what the user should do next.

* Include diagnostic information when appropriate.  For example, if an
  error occurred at a single index in an array operation, then including
  the index where the error happened can help the user better understand
  the cause of the error.

* Write error messages that are concise when possible, as users often
  skim or skip long error messages.

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

Units
=====

PlasmaPy uses |astropy.units|_ to assign physical units to values in the
form of a |Quantity|.

.. code-block:: pycon

   >>> import astropy.units as u
   >>> 5 * u.m / u.s
   <Quantity 5. m / s>

Using |astropy.units|_ improves compatibility with Python packages in
adjacent fields such as astronomy and heliophysics. To get started with
|astropy.units|_, check out this `example notebook on units`_.

  .. caution::

     Some `scipy` functions silently drop units when used on |Quantity|
     instances.

* Only SI units should be used within PlasmaPy, unless there is a strong
  justification to do otherwise. Example notebooks may occasionally use
  other unit systems to show the flexibility of |astropy.units|_.

* Use operations between |Quantity| instances except when needed for
  performance. To improve performance in |Quantity| operations, check
  out `performance tips
  <https://docs.astropy.org/en/stable/units/index.html#performance-tips>`__
  for |astropy.units|_.

* Use unit annotations with the |validate_quantities| decorator to
  validate |Quantity| arguments and return values.

  .. code-block:: python

     from plasmapy.utils.decorators.validators import validate_quantities

     @validate_quantities(
        n={"can_be_negative": False},
        validations_on_return={"equivalencies": u.dimensionless_angles()},
     )
     def inertial_length(n: u.m ** -3, ...) -> u.m:
         ...

  .. caution::

     Recent versions of Astropy_ allow unit-aware |Quantity|
     annotations such as :py:`u.Quantity[u.m]`. However, these
     annotations are not yet compatible with |validate_quantities|.

* Avoid using electron-volts as a unit of temperature within PlasmaPy
  because it is defined as a unit of energy. However, functions in
  `plasmapy.formulary` and elsewhere should accept temperatures in units
  of electron-volts, which can be done using |validate_quantities|.

* Non-standard unit conversions can be made using equivalencies_ such
  as `~astropy.units.temperature_energy`.

  .. code-block:: pycon

     >>> (1 * u.eV).to(u.K, equivalencies=u.temperature_energy())
     11604.518...

* The names of SI units should not be capitalized except at the
  beginning of a sentence, including when they are named after a person.
  The sole exception is "degree Celsius".

Particles
=========

The |Particle| class provides an object-oriented interface for accessing
basic particle data. |Particle| accepts :term:`particle-like` inputs.

.. code-block:: pycon

   >>> from plasmapy.particles import Particle
   >>> alpha = Particle("He-4 2+")
   >>> alpha.mass
   <Quantity 6.6446...e-27 kg>
   >>> alpha.charge
   <Quantity 3.20435...e-19 C>

To get started with `plasmapy.particles`, check out this `example
notebook on particles`_.

* Avoid using implicit default particle assumptions for function
  arguments (see issue :issue:`453`).

* The |particle_input| decorator can automatically transform a
  |particle-like| |argument| into a |Particle|, |CustomParticle|, or
  |ParticleList| instance when the corresponding |parameter| is
  decorated with |ParticleLike|.

  .. code-block:: python

     from plasmapy.particles import particle_input, ParticleLike

     @particle_input
     def get_particle(particle: ParticleLike):
          return particle

  If we use ``get_particle`` on something |particle-like|, it will
  return the corresponding particle object.

  .. code-block:: pycon

     >>> return_particle("p+")
     Particle("p+")

  The documentation for |particle_input| describes ways to ensure that
  the particle meets certain categorization criteria.

Equations and Physical Formulae
===============================

* Physical formulae should be inputted without first evaluating all of
  the physical constants. For example, the following line of code
  obscures information about the physics being represented:

  .. code-block:: pycon

     >>> omega_ce = 1.76e7*(B/u.G)*u.rad/u.s  # doctest: +SKIP

  In contrast, the following line of code shows the exact formula
  which makes the code much more readable.

  .. code-block:: pycon

     >>> omega_ce = (e * B) / (m_e * c)  # doctest: +SKIP

  The origins of numerical coefficients in formulae should be
  documented.

* Docstrings should describe the physics associated with these
  quantities in ways that are understandable to students who are
  taking their first course in plasma physics while still being useful
  to experienced plasma physicists.

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

>>> from astropy import units as u
>>> f_ce = omega_ce.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])  # doctest: +SKIP

However, ``dimensionless_angles`` does work when dividing a velocity
by an angular frequency to get a length scale:

>>> d_i = (c/omega_pi).to(u.m, equivalencies=u.dimensionless_angles())  # doctest: +SKIP

.. _aliases:

Aliases
=======

An :term:`alias` is an abbreviated version of a commonly used function.
For example, `~plasmapy.formulary.speeds.va_` is an alias to
`~plasmapy.formulary.speeds.Alfven_speed`.

:term:`Aliases` are intended to give users the option for shortening
their code while maintaining some readability and explicit meaning. As
such, :term:`aliases` are given to functionality that already has a
widely-used symbol in plasma literature.

Here is a minimal example of an alias ``f_`` to ``function`` as would be
defined in :file:`plasmapy/subpackage/module.py`.

.. code-block:: python

   __all__ = ["function"]
   __aliases__ = ["f_"]

   __all__ += __aliases__

   def function():
       ...

   f_ = function
   """Alias to `~plasmapy.subpackage.module.function`."""

* Aliases should only be defined for functionality that already has a
  symbol that is widely used in the community's literature.  This is to
  ensure that the abbreviated function name is still widely
  understandable. For example, `~plasmapy.formulary.lengths.cwp_` is a
  shortcut for :math:`c/ω_p`\ .

* The name of an alias should end with a trailing underscore.

* An alias should be defined immediately after the original function.

* Each alias should have a one-line docstring that refers users to the
  original function.

* The name of the original function should be included in ``__all__``
  near the top of each module, and the name of the alias should be
  included in ``__aliases__``, which will then get appended to
  ``__all__``. This is done so both the :term:`alias` and the original
  function get properly documented.

* Aliases are intended for end users, and should not be used in PlasmaPy
  or other collaborative software development efforts because of
  reduced readability and searchability for someone new to plasma
  science.

.. _lite-functions:

Lite Functions
==============

Most functions in `plasmapy.formulary` accept |Quantity| instances as
arguments and use |validate_quantities| to verify that |Quantity|
arguments are valid. The use of |Quantity| operations and validations do
not noticeably impact performance during typical interactive use, but
the performance penalty can become significant for numerically intensive
applications.

A :term:`lite-function` is an optimized version of another `plasmapy`
function that accepts numbers and NumPy_ arrays in assumed SI units.
:term:`Lite-functions` skip all validations and instead prioritize
performance. Most :term:`lite-functions` are defined in
`plasmapy.formulary`.

.. caution::

   Unlike most `~plasmapy.formulary` functions, no validations are
   performed on the arguments provided to a :term:`lite-function` for
   the sake of computational efficiency. When using
   :term:`lite-functions`, it is vital to double-check your
   implementation!

Here is a minimal example of a :term:`lite-function` ``function_lite``
that corresponds to ``function`` as would be defined in
:file:`plasmapy/subpackage/module.py`.

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

* The name of each :term:`lite-function` should be the name of the
  original function with ``_lite`` appended at the end. For example,
  `~plasmapy.formulary.speeds.thermal_speed_lite` is the
  :term:`lite-function` associated with
  `~plasmapy.formulary.speeds.thermal_speed`.

* :term:`Lite-functions` assume SI units for all arguments that
  represent physical quantities.

* :term:`Lite-functions` should be defined immediately before the normal
  version of the function.

* :term:`Lite-functions` should be used by their associate non-lite
  counterpart, except for well reasoned exceptions. This is done to
  reduce code duplication.

* :term:`Lite-functions` are bound to their normal version as the
  ``lite`` attribute using the
  `~plasmapy.utils.decorators.lite_func.bind_lite_func` decorator. This
  allows the :term:`lite-function` to also be accessed like
  :py:`thermal_speed.lite()`.

* If a :term:`lite-function` is decorated with something like
  :py:`@njit`, then it should also be decorated with
  `~plasmapy.utils.decorators.helpers.preserve_signature`.  This
  preserves the function signature so interpreters can still
  give hints about function arguments.

* When possible, a :term:`lite-function` should incorporate `numba's
  just-in-time compilation
  <https://numba.pydata.org/numba-doc/latest/reference/jit-compilation.html>`__
  or utilize Cython_.  At a minimum any "extra" code beyond the raw
  calculation should be removed.

* The name of the original function should be included in ``__all__``
  near the top of each module, and the name of the :term:`lite-function`
  should be included in ``__lite_funcs__``, which will then get
  appended to ``__all__``. This is done so both the :term:`lite-function`
  and the original function get properly documented.

.. _example_notebooks:

Examples
========

.. _docs/notebooks: https://github.com/PlasmaPy/PlasmaPy/tree/main/docs/notebooks

Examples in PlasmaPy are written as Jupyter notebooks, taking advantage
of their mature ecosystems. They are located in `docs/notebooks`_. |nbsphinx|_
takes care of executing them at documentation build time and including them
in the documentation.

Please note that it is necessary to store notebooks with their outputs stripped
(use the "Edit -> Clear all" option in JupyterLab and the "Cell -> All Output -> Clear" option in the "classic" Jupyter Notebook). This accomplishes two goals:

1. helps with versioning the notebooks, as binary image data is not stored in
   the notebook
2. signals |nbsphinx|_ that it should execute the notebook.

.. note::

  In the future, verifying and running this step may be automated via a GitHub bot.
  Currently, reviewers should ensure that submitted notebooks have outputs stripped.

If you have an example notebook that includes packages unavailable in the
documentation building environment (e.g., ``bokeh``) or runs some heavy
computation that should not be executed on every commit, *keep the outputs in
the notebook* but store it in the repository with a ``preexecuted_`` prefix, e.g.,
:file:`preexecuted_full_3d_mhd_chaotic_turbulence_simulation.ipynb`.

Benchmarks
==========

.. _benchmarks: https://www.plasmapy.org/plasmapy-benchmarks
.. _benchmarks-repo: https://github.com/PlasmaPy/plasmapy-benchmarks
.. _asv: https://github.com/airspeed-velocity/asv
.. _asv-docs: https://asv.readthedocs.io/en/stable/

PlasmaPy has a set of `asv`_ benchmarks that monitor performance of its
functionalities.  This is meant to protect the package from performance
regressions. The benchmarks can be viewed at `benchmarks`_. They're
generated from results located in `benchmarks-repo`_. Detailed
instructions on writing such benchmarks can be found at `asv-docs`_.
Up-to-date instructions on running the benchmark suite will be located in
the README file of `benchmarks-repo`_.

Compatibility with Prior Versions of Python, NumPy, and Astropy
===============================================================

PlasmaPy releases will generally abide by the following standards,
which are adapted from `NEP 29`_ for the support of old versions of
Python_, NumPy_, and Astropy_.

* PlasmaPy should support at least the minor versions of Python
  initially released 42 months prior to a planned project release date.

* PlasmaPy should support at least the 3 latest minor versions of
  Python.

* PlasmaPy should support minor versions of NumPy initially released
  in the 24 months prior to a planned project release date or the
  oldest version that supports the minimum Python version (whichever is
  higher).

* PlasmaPy should support at least the 3 latest minor versions of
  NumPy and Astropy.

The required major and minor version numbers of upstream packages may
only change during major or minor releases of PlasmaPy, and never during
patch releases.

Exceptions to these guidelines should only be made when there are major
improvements or fixes to upstream functionality or when other required
packages have stricter requirements.

.. _ASCII: https://en.wikipedia.org/wiki/ASCII
.. _cognitive complexity: https://www.sonarsource.com/docs/CognitiveComplexity.pdf
.. _example notebook on particles: ../notebooks/getting_started/particles.ipynb
.. _example notebook on units: ../notebooks/getting_started/units.ipynb
.. _extract function refactoring pattern: https://refactoring.guru/extract-method
.. _NEP 29: https://numpy.org/neps/nep-0029-deprecation_policy.html
.. _not a number: https://en.wikipedia.org/wiki/NaN
.. _NumPy Enhancement Proposal 29: https://numpy.org/neps/nep-0029-deprecation_policy.html
.. _pyupgrade: https://github.com/asottile/pyupgrade
.. _rename refactoring in PyCharm: https://www.jetbrains.com/help/pycharm/rename-refactorings.html
