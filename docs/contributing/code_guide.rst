.. _coding guide:

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
important as consistency, readability, and maintainability.

This guide can (and should!) be regularly refined by the PlasmaPy
community as we collectively learn new practices and our shared coding
style changes. Please feel free to propose revisions to this guide by
:ref:`submitting a pull request <code-contribution>` or by bringing up
an idea at a community meeting.

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

Coding Style
============

PlasmaPy Code Style Guide, codified
-----------------------------------

* PlasmaPy generally follows the :pep:`8` style guide for Python code and the black_ code style.
  This helps ensure that the code will be consistent and
  readable.

  * Docstrings and comments should generally be limited to
    about 72 characters.

* During code development, use black_ to automatically format code and
  ensure a consistent code style throughout the package and isort_ to
  automatically sort imports.

* Follow the existing coding style within a subpackage.  This includes,
  for example, variable naming conventions.

* Use standard abbreviations for imported packages when possible, such
  as ``import numpy as np``, ``import matplotlib as mpl``, ``import
  matplotlib.pyplot as plt``, and ``import astropy.units as u``.

* ``__init__.py`` files for modules should not contain any significant
  implementation code, but it can contain a docstring describing the
  module and code related to importing the module.  Any substantial
  functionality should be put into a separate file.

* Use absolute imports, such as
  ``from plasmapy.particles import Particle``, rather than relative
  imports such as ``from ..particles import Particle``.

* Use ``Optional[type]`` for type hinted keyword arguments with a
  default value of `None`.

* There should be at least one pun per 1284 lines of code.

* Avoid using ``lambda`` to define functions, as this notation may be
  unfamiliar to newcomers to Python.

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

Commit Messages
---------------
Good commit messages communicate context and intention to other
developers and to our future selves.  They provide insight into why we
chose a particular implementation, and help us avoid past mistakes.

Suggestions on `how to write a git commit message
<https://cbea.ms/git-commit>`_:

* Separate subject from body with a blank line

* Limit the subject line to 50 characters

* Capitalize the subject line

* Do not end the subject line with a period

* Use the imperative mood in the subject line

* Wrap the body at 72 characters

* Use the body to explain what and why vs. how

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

>>> omega_ce = 1.76e7*(B/u.G)*u.rad/u.s   # doctest: +SKIP

  In contrast, the following line of code shows the exact formula
  which makes the code much more readable.

>>> omega_ce = (e * B) / (m_e * c)       # doctest: +SKIP

  The origins of numerical coefficients in formulae should be
  documented.

* Docstrings should describe the physics associated with these
  quantities in ways that are understandable to students who are
  taking their first course in plasma physics while still being useful
  to experienced plasma physicists.

* SI units that were named after a person should not be capitalized
  except at the beginning of a sentence.

* Some plasma parameters depend on more than one quantity with
  the same units.  In the following line, it is difficult to discern which
  is the electron temperature and which is the ion temperature.

  >>> ion_sound_speed(1e6*u.K, 2e6*u.K)  # doctest: +SKIP

  Remembering that "explicit is better than implicit", it is more
  readable and less prone to errors to write:

  >>> ion_sound_speed(T_i=1e6*u.K, T_e=2e6*u.K)    # doctest: +SKIP

* SI units that were named after a person should be lower case except at
  the beginning of a sentence, even if their symbol is capitalized. For
  example, kelvin is a unit while Kelvin was a scientist.

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
>>> f_ce = omega_ce.to(u.Hz, equivalencies=[(u.cy/u.s, u.Hz)])   # doctest: +SKIP

However, ``dimensionless_angles`` does work when dividing a velocity
by an angular frequency to get a length scale:

>>> d_i = (c/omega_pi).to(u.m, equivalencies=u.dimensionless_angles())    # doctest: +SKIP

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

.. _ASCII: https://en.wikipedia.org/wiki/ASCII
.. _rename refactoring in PyCharm: https://www.jetbrains.com/help/pycharm/rename-refactorings.html
