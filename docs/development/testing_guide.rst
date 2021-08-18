*************
Testing Guide
*************

Software testing is vital for software reliability and maintainability.
Software tests help us to:

* Find and fix bugs.
* Prevent old bugs from getting re-introduced.
* Provide confidence that our code is behaving correctly.
* Define what "behaving correctly" actually means.
* Speed up code development and refactoring.
* Show future contributors examples of how code was intended to be used.
* Confirm that our code works on different operating systems and
  with different versions of software dependencies.
* Enable us to change code with confidence that we are not unknowingly
  introducing bugs elsewhere in our program.

.. hint::

   Writing tests takes time, but debugging takes more time.

Every code contribution to PlasmaPy with new functionality must also
have corresponding tests. Creating or updating a pull request will
activate PlasmaPy's test suite to be run via `GitHub Actions`_, along
with some additional checks. The results of the test suite are shown at
the bottom of each pull request. Click on *Details* next to each test
run to find the reason for any test failures.

PlasmaPy's tests are set up using the pytest_ framework. The tests for
a subpackage are located in its :file:`tests` subdirectory in files with
names of the form :file:`test_*.py`. For example, tests for
`plasmapy.formulary.parameters` are located at
:file:`plasmapy/formulary/tests/test_parameters.py` relative to the top
of the package. Example code contained within docstrings is tested to
make sure that the actual printed output matches what is in the
docstring.

Running tests
=============

PlasmaPy's tests can be run in the following ways:

1. Creating and updating a pull request on GitHub_.
2. Running pytest_ from the command line.
3. Running tox_ from the command line.
4. Running tests from an `integrated development environment` (IDE).

We recommend that new contributors perform the tests via a pull request
on GitHub_. Creating a draft pull request and keeping it updated will
ensure that the necessary checks are run frequently. This approach is
also appropriate for pull requests with a limited scope. This advantage
of this approach is that the tests are run automatically and do not
require any extra work. The disadvantages are that running the tests on
GitHub_ is often slow and that navigating the test results is sometimes
difficult.

We recommend that experienced contributors run tests either by using
pytest_ from the command line or by using your preferred IDE.
Using tox_ is an alternative to pytest_, but running tests with tox_ is
typically slower than running tests with pytest_.

Using GitHub
------------

The recommended way for new contributors to run PlasmaPy's full test
suite is to `create a pull request`_ from your development branch to
`PlasmaPy's GitHub repository`_. The test suite will be run
automatically when the pull request is created and every time changes
are pushed to the development branch on GitHub_.

The following checks are performed with each pull request. The results
of the checks are found near the end of the *Conversation* tab in each
pull request. Most of these checks have been automated using `GitHub
Actions`_. Checks that pass are marked with ✔️, while tests that fail
are marked with ❌. Click on *Details* for information about why a
particular check failed.

* Checks with labels like **CI / Python 3.9 (pull request)** verify that
  PlasmaPy works with different versions of Python and other
  dependencies, and on different operating systems.

  * These tests are set up using tox_ and run with pytest_.

  * When multiple checks fail, investigate these tests first.

If multiple tests fail, investigate these tests first.

* Checks with labels like **CI / Python 3.9 with NumPy dev (pull
  request)** verify that PlasmaPy works the version of NumPy that is
  currently being developed on GitHub_. Occasionally these tests will
  fail for reasons not associated with a particular pull request.

* The **CI / Documentation (pull request)** check verifies that
  PlasmaPy's documentation is able to build correctly from the pull
  request. Warnings are treated as errors.

* The **docs/readthedocs.org:plasmapy** check allows us to preview
  how the documentation will appear if the pull request is merged.
  Click on *Details* link to access this preview.

* The check labeled **changelog: found** or **changelog: absent**
  indicates whether or not a changelog entry with the correct number
  is present, unless the pull request has been labeled with "No
  changelog entry needed".

  * The :file:`changelog/README.rst` file describes the process for
    adding a changelog entry to a pull request.

* The **codecov/patch** and **codecov/project** checks generate test
  coverage reports which are displayed by Codecov_. These reports show
  which lines of code are covered by tests and which are not. The
  Codecov_ checks will be marked as passing when the test coverage is
  satisfactorily high.

  * Test coverage reports help us write targeted tests to fill in gaps
    in test coverage and find unreachable blocks of code that can never
    be run.

* PlasmaPy uses black_ to format code and isort_ to sort `import`
  statements. The **CI / Linters (pull request)** and
  **pre-commit.ci - pr** checks verify that the pull request meets these
  style requirements. These checks will fail when inconsistencies with
  the output from black_ or isort_ are found or when there are syntax
  errors. These checks can usually be ignored until the pull request is
    nearing completion.

  .. tip::

     The required formatting fixes can be applied automatically by
     writing a comment with the message ``pre-commit.ci autofix`` to the
     *Conversation* tab on a pull request, as long as there are no
     syntax errors. This approach is much more efficient than making the
     style fixes manually. Remember to ``git pull`` afterwards!

* The **CI / Packaging (pull request)** check verifies that no errors
  arise that would prevent an official release of PlasmaPy from being
  made.

* The **Pull Request Labeler / triage (pull_request_target)** check
  applies appropriate GitHub_ labels to pull requests.

.. attention::

   The continuous integration checks performed for pull requests change
   frequently. If you notice that the above list has become out-of-date,
   please `submit an issue that this section needs updating
   <https://github.com/PlasmaPy/PlasmaPy/issues/new?title=Update%20information%20on%20GitHub%20checks%20in%20testing%20guide&labels=Documentation>`__.

Using pytest
------------

To install the packages necessary to run tests on your local computer,
run:

.. code-block:: shell

   pip install -r requirements.txt

To run PlasmaPy's tests from the command line, go to a directory within
PlasmaPy's repository and run:

.. code-block:: shell

   pytest

This command will run all of the tests found within your current
directory and all of its subdirectories. Because it takes time to run
PlasmaPy's tests, it is usually most convenient to specify that only a
subset of the tests be run. To run the tests contained within a
particular file or directory, include its name after ``pytest``. The
tests in :file:`test_atomic.py` can be run with:

.. code-block:: shell

   pytest test_atomic.py

The documentation for pytest_ describes `how to invoke pytest`_ and
specify which tests should or should not be run.

.. _`how to invoke pytest`: https://docs.pytest.org/en/latest/how-to/usage.html

Some tests in the test suite can take a long time to run, which can slow
down development of new features. These tests are decorated with
`pytest.mark.slow`. To skip the slow tests, run:

.. code-block:: shell

   pytest -m 'not slow'

To exclusively run the slow tests, run:

.. code-block:: shell

   pytest -m slow

Using tox
---------

PlasmaPy's continuous integration tests on GitHub_ are typically run
using tox_, a tool for automating Python testing. Using tox_ simplifies
testing PlasmaPy with different releases of Python, with different
versions of PlasmaPy's dependencies, and on different operating systems.
While testing with tox_ is more robust than testing with pytest_, using
tox_ to run tests is typically slower because tox_ creates its own
virtual environments.

The `tox environments`_ are found in :file:`tox.ini` in the
top-level directory of PlasmaPy's repository. To find a list of
the environments defined in :file:`tox.ini`, run:

.. code-block:: shell

   tox -a

The ``py39`` testing environment, for example, can be run with:

.. code-block:: shell

   tox -e py39

These commands can be run in any directory within PlasmaPy's repository
with the same effect.

Environments with names like ``py38``, ``py39``, and ``py310`` are
interpreted to mean that the tests should be performed with Python 3.8,
3.9, or 3.10, respectively. Running these tests requires that the
appropriate version of Python has been installed and can be found by
tox_.

Using an IDE
------------

Most `integrated development environments`_ (IDEs) have built-in tools
that simplify running tests. Setting up testing configurations within an
IDE generally saves considerable time

Testing techniques
==================

Assertions
----------

A software test runs a section of code and checks that a particular
condition is met.  If the condition is not met, then the test fails.
This check is most commonly made using an `assert` statement.

.. code-block:: python

  def test_addition():
      assert 2 + 2 == 4

If the condition is not met, then the `assert` statement will raise an
`AssertionError`.

When `assert` statements raise an `AssertionError`, pytest_ will display
the values of the expressions evaluated in the `assert` statement. The
automatic output from pytest is sufficient for simple tests as
above. For more complex tests, we can add a descriptive error message
to help us find the cause of a particular test failure.

.. code-block:: python

  def test_addition():
      result = 2 + 2
      expected = 4
      assert result == expected, f"2 + 2 returns {result} instead of {expected}."

.. TODO Python 3.8+: update this example to use the f"{result=}" syntax.

.. tip::

   Use `f-strings`_ to improve error message readability.

Floating point comparisons
--------------------------

Because of limited precision and rounding errors, comparisons between
floating point numbers with ``==`` is not recommended.  Additionally,
the values of fundamental constants in `astropy.constants` are
occasionally refined as improvements become available.

Using `numpy.isclose` when comparing floating point numbers and
`astropy.units.isclose` for |Quantity| instances lets us
avoid these difficulties.  The ``rtol`` keyword for each of these
functions allows us to set an acceptable relative tolerance.  Ideally,
``rtol`` should be set to be an order of magnitude or two greater than
the expected uncertainty.  For mathematical functions, a value of
``rtol=1e-14`` may be appropriate.  For quantities that depend on
physical constants, a value between ``rtol=1e-8`` and ``rtol=1e-5`` may
be required, depending on how much the accepted values for fundamental
constants are likely to change.  For comparing arrays, `numpy.allclose`
and `astropy.units.allclose` should be used instead.

Testing warnings and exceptions
-------------------------------

Robust testing frameworks should test that functions and methods return
the expected results, issue the expected warnings, and raise the
expected exceptions.  pytest_ contains functionality to `test warnings`_
and `test exceptions`_.

To test that a function issues an appropriate warning, use
`pytest.warns`.

.. code-block:: python

  import pytest
  import warnings

  def issue_warning():
      warnings.warn("Beware the ides of March", UserWarning)

  def test_issue_warning():
      with pytest.warns(UserWarning):
          issue_warning()

To test that a function raises an appropriate exception, use
`pytest.raises`.

.. code-block:: python

  def raise_exception():
      raise Exception

  def test_raise_exception():
      with pytest.raises(Exception):
          raise_exception()
          pytest.fail("Exception not raised.")

.. _testing-guidelines-writing-tests-parametrize:

Test independence and parametrization
-------------------------------------

In this section, we'll discuss the issue of parametrization based on
an example of a `proof
<https://en.wikipedia.org/wiki/Riemann\_hypothesis#Excluded\_middle>`_
of Gauss's class number conjecture.

The proof goes along these lines:

* If the generalized Riemann hypothesis is true, the conjecture is true.

* If the generalized Riemann hypothesis is false, the conjecture is also
  true.

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
tests for the same function, the preferred method is to use
`pytest.mark.parametrize`.

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



Writing Tests
=============

Pull requests must include tests of new or changed functionality before
being merged.

Best practices for writing tests
--------------------------------

The following guidelines are helpful suggestions for writing readable,
maintainable, and robust tests.

* Each function and method should have unit tests that check that it
  returns the expected results, issues the appropriate warnings, and
  raises the appropriate exceptions.

* Each unit test should test *one unit of behavior*, *quickly*, and
  *in isolation from other tests*.

.. add citation for above from the audiobook that I don't feel like
   looking up again

* Bugs should be turned into test cases.

* Tests are run frequently during code development, and slow tests may
  interrupt the flow of a contributor.  Tests should be minimal,
  sufficient enough to be complete, efficient.

* Decorate slow tests with `pytest.mark.slow`.

  .. code-block:: python

     import pytest

     @pytest.mark.slow
     def test_calculating_primes():
        calculate_all_primes()

* Write test code with the same quality as production code. Well-written
  tests are easier to modify when the tested behavior changes. Poorly
  written tests are difficult to change and slow down future development.


.. The following hint would be worth putting somewhere, at least after
   Python 3.10 is released, but maybe not here.

.. .. hint::
   Running tests in Python ≥3.10 will provide improved error messages
   compared to Python ≤3.9.

Test organization and collection
--------------------------------

Pytest has certain `test discovery conventions
<https://docs.pytest.org/en/latest/goodpractices.html#conventions-for-python-test-discovery>`_
that are used to collect the tests to be run.

The tests for each subpackage are contained in a :file:`tests/` subdirectory.
For example, the tests for `~plasmapy.particles` are located in
:file:`plasmapy/particles/tests`.  Test files should begin with :file:`test_`
and generally contain the name of the module or `object` that is being
tested.

The functions that are to be tested in each test file should likewise be
prepended with `test_` (e.g., :file:`test_atomic.py`).  Tests may also be
`grouped into classes
<https://docs.pytest.org/en/latest/getting-started.html#group-multiple-tests-in-a-class>`_.
In order for pytest to find tests in classes, the class name should
start with ``Test`` and the methods to be run as tests should start with
``test_``.  For example, :file:`test_particle_class.py` could define the
``TestParticle`` class containing the method ``test_charge_number``.

.. _testing-guidelines-writing-tests-fixtures:

Fixtures
--------

`Fixtures <https://docs.pytest.org/en/stable/fixture.html>`_ provide a
way to set up well-defined states in order to have consistent tests.
We recommend using fixtures for complex tests that would be unwieldy to
set up with parametrization as described above.

.. At some point in the future, we may wish to add more information
   and/or more references for pytest fixtures when we use them more
   frequently.

.. .. _testing-guidelines-coverage:

.. Code Coverage
.. =============

.. PlasmaPy uses `Codecov`_ to show what lines of code
are covered by the test suite and which lines are not.  At the end of
every testing session, information on which lines were
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
line to :file:`.coveragerc`).  This strategy should be used sparingly, since
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
:file:`htmlcov/index.html` in your favorite web brower.  These commands
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

The :file:`.coveragerc` file is used to specify lines of code and files that
should always be ignored in coverage tests.

.. note::

  In general, untested lines of code should remain marked as untested to
  give future developers a better idea of where tests should be added in
  the future and where potential bugs may exist.

.. Footnotes
   =========

.. .. [1] In `Working Effectively With Legacy Code
   <https://www.oreilly.com/library/view/working-effectively-with/0131177052/>`__,
   Michael Feathers bluntly writes: "Code without tests is bad code.  It
   doesn't matter how well written it is; it doesn't matter how pretty
   or object-oriented or well-encapsulated it is.  With tests, we can
   change the behavior of our code quickly and verifiably.  Without
   them, we really don't know if our code is getting better or worse."

.. .. [2] In the chapter "Bugs Are Missing Tests" in `Beyond
   Legacy Code <https://pragprog.com/book/dblegacy/beyond-legacy-code>`__,
   David Bernstein writes: "Every bug exists because of a missing test
   in a system.  The way to fix bugs using TDD [test-driven development]
   is first write a failing test that represents the bug and then fix
   the bug and watch the failing test turn green.


.. _Codecov: https://about.codecov.io/
.. _`create a pull request`: https://help.github.com/articles/creating-a-pull-request
.. _`f-strings`: https://docs.python.org/3/tutorial/inputoutput.html#tut-f-strings
.. _`integrated development environment`: https://en.wikipedia.org/wiki/Integrated_development_environment
.. _pytest: https://docs.pytest.org/
.. _`test warnings`: https://docs.pytest.org/en/latest/warnings.html#warns
.. _`test exceptions`: https://docs.pytest.org/en/latest/assert.html#assertions-about-expected-exceptions
.. _`tox environments`: https://tox.readthedocs.io/en/latest/config.html?highlight=py37#tox-environments
