.. _testing-guidelines:

******************
Testing Guidelines
******************

.. _testing-guidelines-motivation:

Motivation
==========

Tests are vital for software reliability and maintainability.  Writing
tests requires additional effort now, but saves considerable time in the
long run.  Tests enable us to modify code and quickly discover when we
introduce errors [1]_.  Tests provide future contributors with
examples of how functions and classes were originally intended to be
used.

Tests should also be readable and maintainable.  Well-written tests are
easier to understand and modify when the behavior of a function or
method is intended to be changed.  Consequently, tests should be held to
the same coding standards as the rest of the package.

When bugs are discovered, they should be turned into test cases to
prevent the bug from emerging again in the future [2]_.

.. _testing-guidelines-overview:

Overview
========

Pull requests that create or change functionality must include tests
before being merged.

PlasmaPy uses `pytest <https://docs.pytest.org>`_ for software testing.
The test suite may be run locally or automatically via pull requests on
GitHub.  PlasmaPy undergoes continuous integration testing of the code
base by `Travis CI <https://travis-ci.org>`_ and `AppVeyor
<https://www.appveyor.com>`_, including code examples in docstrings.
`Codecov <https://codecov.io>`_ performs test coverage checks and shows
whether or not each line of code is run during the test suite.
`CircleCI <https://circleci.com/>`_ tests that the documentation can be
successfully built.  The results of the documentation test builds are
displayed using `Giles <https://github.com/apps/giles>`_.  PlasmaPy's
test suite is automatically run whenever a pull request to the main
repository is made or updated.

.. _testing-guidelines-running-tests:

Running Tests
=============

.. _testing-guidelines-running-tests-github:

Running tests on GitHub
-----------------------

The recommended way to run PlasmaPy's full test suite when contributing
code is to `create a pull request
<https://help.github.com/articles/creating-a-pull-request/>`_ from your
development branch to `PlasmaPy's GitHub repository
<https://github.com/PlasmaPy/PlasmaPy>`_.  The test suite will be run
when the pull request is created and every time your development branch
is subsequently updated.

`Travis CI <https://travis-ci.org>`_ and `AppVeyor
<https://www.appveyor.com>`_ run code tests and check that code examples
in docstrings produce the expected output.  Travis CI runs the tests in a
Linux/MacOS environment whereas AppVeyor runs the tests in a Windows
environment.

The results from Travis CI are used to generate test coverage reports
which are displayed by `Codecov <https://codecov.io>`_. These reports
show which lines of code are covered by tests and which are not, and
allow us to write targeted tests to fill in the gaps in test coverage.
The results displayed by Codecov will be marked as passing when the code
coverage is sufficiently high.

`Circle CI <https://circleci.com>`_ performs a test build of the
documentation in both HTML and LaTeX formats, and reports any errors
that arise.

If any inconsistencies with the `PEP 8 style guide
<https://www.python.org/dev/peps/pep-0008/?>`_ are found, then
`pep8speaks <https://pep8speaks.com/>`_ will comment on the pull request
and update that comment as the pull request is updated.

.. _testing-guidelines-running-tests-command-line:

Running tests from the command line
-----------------------------------

The recommended method for running the test suite locally on your
computer is running

.. code-block:: shell

  python setup.py test

in the repository's root directory.  This command will run all of the
tests and verify that examples in docstrings produce the expected
output.  This command (which was enabled by `integrating pytest with
setuptools
<https://docs.pytest.org/en/latest/goodpractices.html#integrating-with-setuptools-python-setup-py-test-pytest-runner>`_)
ensures that the package is set up and `Cython <http://cython.org>`_
code is compiled before the tests are run.  These tests should be run in
a Python environment in which PlasmaPy has not already been installed.

Command line options for pytest may be passed using the ``-a`` flag.
For example, if you want to stop pytest after two test failures, return
short traceback reports, and run tests only if the test path contains
``plasma`` and not ``blob``, then run

.. code-block:: shell

  python setup.py test -a "--maxfail=2 --tb=short -k 'plasma and not blob'"

One may also run ``pytest`` as a shortcut from the command line, though
this command may result in an error if pytest collects tests
of Cython code.

.. _testing-guidelines-running-tests-python:

Running tests within Python
---------------------------

After installing PlasmaPy by running ``pip install plasmapy`` or
``python setup.py install``, then PlasmaPy's test suite may be run
using

.. code-block:: python

  >>> import plasmapy
  >>> plasmapy.test()

.. _testing-guidelines-writing-tests:

Writing Tests
=============

Pull requests must include tests of new or changed functionality before
being merged.

.. _testing-guidelines-writing-tests-best-practices:

Best practices for writing tests
--------------------------------

The following guidelines are helpful ways for writing neat, readable,
and useful tests.

* Each function and method should have unit tests that check that it
  returns the expected results, issues the appropriate warnings, and
  raises the appropriate exceptions.

* Bugs should be turned into test cases.

* Tests are run frequently during code development, and slow tests may
  interrupt the flow of a contributor.  Tests should be minimal,
  sufficient enough to be complete, and as efficient as possible.

.. _testing-guidelines-writing-tests-organization:

Test organization and collection
--------------------------------

Pytest has certain `test discovery conventions
<https://docs.pytest.org/en/latest/goodpractices.html#conventions-for-python-test-discovery>`_
that are used to collect the tests to be run.

The tests for each subpackage are contained in a ``tests`` subfolder.
For example, the tests for `~plasmapy.atomic` are located in
``plasmapy/atomic/tests``.  Test files should begin with ``test_`` and
generally contain the name of the module or `object` that is being
tested.

The functions that are to be tested in each test file should likewise be
prepended with `test_` (e.g., ``test_atomic.py``).  Tests may also be
`grouped into classes
<https://docs.pytest.org/en/latest/getting-started.html#group-multiple-tests-in-a-class>`_.
In order for pytest to find tests in classes, the class name should
start with ``Test`` and the methods to be run as tests should start with
``test_``.  For example, ``test_particle_class.py`` could define the
``TestParticle`` class containing the method ``test_integer_charge``.

.. _testing-guidelines-writing-tests-asserts:

Assert statements
-----------------

* Pytest often runs tests by checking `assert` statements.

.. code-block:: python

  def test_addition():
      assert 2 + 2 == 4

When `assert` statements raise an `AssertionError`, pytest will display
the values of the expressions evaluated in the `assert` statement.  This
output is often sufficient for simple tests as above.  For more complex
tests

However, descriptive error messages that provide context help us
pinpoint the causes of test failures more quickly.

.. code-block:: python

  def test_addition():
      assert 2 + 2 == 4, "Addition is broken. Reinstall universe and reboot."

When `assert` statements raise an `AssertionError`, pytest will display
the value of the expressions in `assert` statements.  In cases where

Pytest display the value of each expression in `assert`
statements, so for simul.

To make the error statement easier to read, the values of
variables can be included in the error message by using `f-strings
<https://www.python.org/dev/peps/pep-0498/>`_.

.. code-block:: python

  def test_addition():
      result = 2 + 2
      expected = 4
      assert result == expected, f"2 + 2 returns {result} instead of {expected}."

.. _testing-guidelines-writing-tests-warnings:

Testing warnings and exceptions
-------------------------------

Robust testing frameworks should test that functions and methods return
the expected results, issue the expected warnings, and raise the
expected exceptions.  Pytest contains functionality to `test warnings
<https://docs.pytest.org/en/latest/warnings.html#warns>`_
and `test exceptions
<https://docs.pytest.org/en/latest/assert.html#assertions-about-expected-exceptions>`_.

To test that a function issues an appropriate warning, use
`pytest.warns`.

.. code-block:: python

  import pytest, warnings

  def issue_user_warning():
      warnings.warn("grumblemuffins", UserWarning)

  def test_issue_warning():
      with pytest.warns(UserWarning, message="UserWarning not issued."):
          issue_user_warning()

To test that a function raises an appropriate exception, use
`pytest.raises`.

.. code-block:: python

  def raise_value_error():
      raise ValueError

  def test_raise_value_error():
      with pytest.raises(ValueError, message="ValueError not raised."):
          raise_value_error()

.. _testing-guidelines-writing-tests-parametrize:

Test independence and parametrization
-------------------------------------

In this section, we'll discuss the issue of parametrization based on
an example of a `proof
<https://en.wikipedia.org/wiki/Riemann\_hypothesis#Excluded\_middle>`_
of Gauss's class number conjecture.

The proof goes along these lines:
* If the generalized Riemann hypothesis is true, the conjecture is true.
* If the generalized Riemann hypothesis is false, the conjecture is also true.
* Therefore, the conjecture is true.

One way to use pytest would be to write sequential test in a single
function.

.. code-block:: python

  def test_proof_by_riemann_hypothesis():
       assert proof_by_riemann(False)
       assert proof_by_riemann(True)  # only run if previous test passes

If the first test were to fail, then the second test will never be run.
We would therefore not know the potentially useful results of the second
test.  This drawback can be avoided by making independent tests that
will both be run.

.. code-block:: python

  def test_proof_if_riemann_false():
       assert proof_by_riemann(False)

  def test_proof_if_riemann_true():
       assert proof_by_riemann(True)

However, this approach can lead to cumbersome, repeated code if you are
calling the same function over and over.  If you wish to run multiple
tests for the same function, the preferred method is to use pytest's
`parametrization <https://docs.pytest.org/en/stable/parametrize.html>`_
capabilities.

.. code-block:: python

  @pytest.mark.parametrize("truth_value", [True, False])
  def test_proof_if_riemann(truth_value):
       assert proof_by_riemann(truth_value)

This code snippet will run ``proof_by_riemann(truth_value)`` for each
``truth_value`` in ``truth_values_to_test``.  Both of the above
tests will be run regardless of failures.  This approach is much cleaner
for long lists of arguments, and has the advantage that you would only
need to change the function call in one place if something changes.

With qualitatively different tests you would use either separate
functions or pass in tuples containing inputs and expected values.

.. code-block:: python

  @pytest.mark.parametrize("truth_value, expected", [(True, True), (False, True)])
  def test_proof_if_riemann(truth_value, expected):
       assert proof_by_riemann(truth_value) == expected

.. _testing-guidelines-writing-tests-helpers:

Pytest helpers
--------------

A robust testing framework should test not just that functions and
methods return the expected results, but also that they issue the
expected warnings and raise the expected exceptions. In PlasmaPy, tests
often need to compare a `float` against a `float`, an `~numpy.array`
against an `~numpy.array`, and `~astropy.units.Quantity` objects against
other `~astropy.units.Quantity` objects to within a certain tolerance.
Occasionally tests will be needed to make sure that a function will
return the same value for different arguments (e.g., due to symmetry
properties). PlasmaPy's `~plasmapy.utils` subpackage contains the
`~plasmapy.utils.run_test` and
`~plasmapy.utils.run_test_equivalent_calls` helper functions that can
generically perform many of these comparisons and checks.

The `~plasmapy.utils.run_test` function can be used to check that a
callable object returns the expected result, raises the expected
exception, or issues the expected warning for different positional and
keyword arguments.  This function is particularly useful when unit
testing straightforward functions when you have a bunch of inputs and
know the expected result.

Suppose that we want to test the trigonometric property that

.. math::

  \sin(\theta) = \cos(\theta + \frac{\pi}{2}).

We may use `~plasmapy.utils.run_test` as in the following example to
check the case of :math:`\theta \equiv 0`.

.. code-block:: python

  from numpy import sin, cos, pi
  from plasmapy.utils import run_test

  def test_trigonometric_properties():
      run_test(func=sin, args=0, expected_outcome=cos(pi/2), atol=1e-16)

We may use `pytest.mark.parametrize` with `~plasmapy.utils.run_test` to
check multiple cases.  If `~plasmapy.utils.run_test` only receives one
positional argument that is a `list` or `tuple`, then it will assume
that `list` or `tuple` contains the `callable`, the positional
arguments, the keyword arguments (which may be omitted), and the
expected outcome (which may be the returned `object`, a warning, or an
exception).

.. code-block:: python

  @pytest.mark.parametrize("input_tuple", [(sin, 0, cos(pi/2)), (sin, '.', TypeError)])
  def test_trigonometry(input_tuple):
      run_test(input_tuple, atol=1e-16)

This parametrized function will check that ``sin(0)`` is within
``1e-16`` of ``cos(pi/2)`` and that  ``sin('.')`` raises a `TypeError`.

We may use `~plasmapy.utils.run_test_equivalent_calls` to check symmetry
properties such as

.. math::

  \cos(\theta) = \cos(-\theta).

This property can be checked for :math:`\theta = 1` with the following
code.

.. code-block:: python

  def test_cosine_symmetry():
      """Test that cos(1) equals cos(-1)."""
      plasmapy.utils.run_test_equivalent_calls(cos, 1, -1)

We may also use `pytest.mark.parametrize` with
`~plasmapy.utils.run_test_equivalent_calls` to sequentially test
multiple symmetry properties.

.. code-block:: python

  @pytest.mark.parametrize('input_tuple', [(cos, 1, -1), ([cos, pi/2], [sin, 0])])
  def test_symmetry_properties(input_tuple):
      plasmapy.utils.run_test_equivalent_calls(input_tuple, atol=1e-16)

This parametrized function will check that ``cos(1)`` is withn ``1e-16``
of ``cos(-1)``, and that ``cos(pi/2)`` is within ``1e-16`` of ``sin(0)``.

Please refer to the documentation for `~plasmapy.utils.run_test` and
`~plasmapy.utils.run_test_equivalent_calls` to learn about the full
capabilities of these pytest helper functions.

.. _testing-guidelines-writing-tests-fixtures:

Fixtures
--------

`Fixtures <https://docs.pytest.org/en/stable/fixture.html>`_ provide a
way to set up well-defined states in order to have consistent tests.

.. TODO: finish writing this section!

We recommend using fixtures for complex tests that

.. code-block:: python

  import pytest

  @pytest.fixture
  def sample_fixture():
      return {'x': 1, 'y': 2}

  def test_fixture(sample_fixture):
      assert sample_fixture['x'] == 1
      assert sample_fixture['y'] == 2

Fixtures are recommended for complex tests...

.. _testing-guidelines-coverage:

Code Coverage
=============

PlasmaPy uses `Codecov <https://codecov.io>`_ to show what lines of code
are covered by the test suite and which lines are not.  At the end of
every Travis CI testing session, information on which lines were
executed is sent to Codecov.  Codecov comments on the pull request on
GitHub with a coverage report.

.. The following lines should be included if we end up using Numba JIT
   compiled functions:  "At the time of writing this, coverage.py has a
   known issue with being unable to check lines executed in Numba JIT
   compiled functions."

.. _testing-guidelines-coverage-testing:

Test coverage of contributed code
---------------------------------

Code contributions to PlasmaPy are required to be well-tested.  A good
practice is for new code to have a test coverage percentage of at least
about the current code coverage. Tests must be provided in the original
pull request, because often a delayed test ends up being a test not
written.  There is no strict cutoff percentage for how high the code
coverage must be in order to be acceptable, and it is not always
necessary to cover every line of code.  For example, it is often helpful
for methods that raise a `NotImplementedError` to be marked as untested
as a reminder of unfinished work.

Occasionally there will be some lines that do not require testing.
For example, testing exception handling for an `ImportError` when
importing an external package would usually be impractical.  In these
instances, we may end a line with ``# coverage: ignore`` to indicate
that these lines should be excluded from coverage reports (or add a
line to ``.coveragerc``).  This strategy should be used sparingly, since
it is often better to explicitly test exceptions and warnings and to
show the lines of code that are not tested.

.. _testing-guidelines-coverage-local:

Generating coverage reports locally
-----------------------------------

Coverage reports may be generated on your local computer by running

.. code-block:: shell

  python setup.py test --coverage
  coverage html

The coverage reports may be accessed by opening the newly generated
``htmlcov/index.html`` in your favorite web brower.  These commands
require the ``pytest`` and ``coverage`` packages to be installed.

.. _testing-guidelines-coverage-ignore:

Ignoring lines in coverage tests
--------------------------------

Occasionally there will be lines of code that do not require tests.  For
example, it would be impractical to test that an `ImportError` is raised
when running ``import plasmapy`` from Python 2.7.

To ignore a line of code in coverage tests, append it with
``# coverage: ignore``.  If this comment is used on a line with a
control flow structure (e.g., `if`, `for`, and `while`) that begins a
block of code, then all lines in that block of code will be ignored.  In
the following example, lines 3 and 4 will be ignored in coverage tests.

.. code-block:: python
  :linenos:
  :emphasize-lines: 3,4

  try:
      import numpy
  except ModuleNotFoundError as exc:  # coverage: ignore
      raise RuntimeError from exc

The ``.coveragerc`` file is used to specify lines of code and files that
should always be ignored in coverage tests.  For example, tests in
``astropy-helpers`` should not be run because those tests are performed
through the ``astropy-helpers`` repository.

.. note::

  In general, untested lines of code should remain marked as untested to
  give future developers a better idea of where tests should be added in
  the future and where potential bugs may exist.

Footnotes
=========

.. [1] In `Working Effectively With Legacy Code
   <https://www.oreilly.com/library/view/working-effectively-with/0131177052/>`__,
   Michael Feathers bluntly writes: "Code without tests is bad code.  It
   doesn't matter how well written it is; it doesn't matter how pretty
   or object-oriented or well-encapsulated it is.  With tests, we can
   change the behavior of our code quickly and verifiably.  Without
   them, we really don't know if our code is getting better or worse."

.. [2] In the chapter "Bugs Are Missing Tests" in `Beyond
   Legacy Code <https://pragprog.com/book/dblegacy/beyond-legacy-code>`__,
   David Bernstein writes: "Every bug exists because of a missing test
   in a system.  The way to fix bugs using TDD [test-driven development]
   is first write a failing test that represents the bug and then fix
   the bug and watch the failing test turn green.
