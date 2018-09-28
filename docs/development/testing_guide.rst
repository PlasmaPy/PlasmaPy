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
method is intended to be changed. When bugs are discovered, they should
be turned into test cases to prevent the bug from emerging again in the
future [2]_.

.. _testing-guidelines-overview:

Overview
========

PlasmaPy uses the `pytest <https://docs.pytest.org>`_ framework for
software testing.  The test suite may be run locally or automatically
via pull requests on GitHub.  PlasmaPy undergoes continuous integration
testing of the code base by `Travis CI <https://travis-ci.org>`_ and
`AppVeyor <https://www.appveyor.com>`_, including code examples in
docstrings. `Codecov <https://codecov.io>`_ performs test coverage
checks and shows whether or not each line of code is run during the test
suite. `CircleCI <https://circleci.com/>`_ tests that the documentation
can be successfully built.  The results of the documentation test builds
are displayed using `Giles <https://github.com/apps/giles>`_.
PlasmaPy's test suite is automatically run whenever a pull request to
the main repository is made or updated.

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

The test suite may be performed locally on your computer by running

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



Tests are required when adding new functionality in code contribution.

Best practices
--------------

The following guidelines are helpful ways for writing neat, readable,
and useful tests.

* Unit tests should be provided for all methods when possible.

* Bugs should be turned into test cases.

..  The Travis CI, CircleCI, and AppVeyor integrations on GitHub run tests
  whenever pull requests to PlasmaPy are created or updated.  The pytest
  module may used on a local computer.

* Tests are run frequently during code development, and slow tests may
  interrupt the flow of a contributor.  Tests should be minimal,
  sufficient enough to be complete, and as efficient as possible.

Assert statements
-----------------

* Pytest runs tests by checking ``assert`` statements, so this is sufficient:

.. code-block:: python

  def test_universe_is_sane():
      assert 2 + 2 == 4

However, making assertions descriptive is better in most cases:

.. code-block:: python

  def test_universe_is_sane():
      assert 2 + 2 == 4, "Addition is broken. Reinstall the universe and reboot."

pytest should display the value of the ``2 + 2`` expression, but the value can be added to the thrown string:

.. code-block:: python

  def test_universe_is_sane():
      assert 2 + 2 == 4, f"Addition is broken, 2 + 2 giving {2 + 2}. Reinstall the universe and reboot."

A note on test independence and parametrization
-----------------------------------------------

In this section, we'll discuss the issue of parametrization based on
an example of a `proof
<https://en.wikipedia.org/wiki/Riemann\_hypothesis#Excluded\_middle>`_
of Gauss's class number conjecture.

The proof goes along these lines:
* If the generalized Riemann hypothesis is true, the conjecture is true.
* If the former is false, the latter is also true.
* Therefore, the latter is true.

One way to use pytest for testing is to write continuous assertions:

.. code-block:: python

  def test_proof_by_riemann_hypothesis():
       # if this step fails, the test stops
       assert proof_by_riemann(False)
       # and you have to run this again
       assert proof_by_riemann(True)

To do this the right way, what you technically should do to make the
tests independent:

.. code-block:: python

  def test_proof_if_riemann_false():
       assert proof_by_riemann(False)
  def test_proof_if_riemann_true():
       assert proof_by_riemann(True)

but that's a lot of typing so what you actually do is use pytest parametrization:

.. code-block:: python

  @pytest.mark.parametrize("truth", [True, False])
  def test_proof_if_riemann(truth):
       assert proof_by_riemann(truth)

And both of these are going to run regardless of failures, which is
awesome!

Of course, with qualitatively different tests you would use either
separate functions or you'd pass in pairs of inputs and expected
values:

.. code-block:: python

  @pytest.mark.parametrize("truth,expected", [(True, True), (False, True)])
  def test_proof_if_riemann(truth, expected):
       assert proof_by_riemann(truth) == expected

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

Generating coverage reports locally
-----------------------------------

Coverage reports may be generated on your local computer by running

.. code-block:: shell

  python setup.py test --coverage
  coverage html

The coverage reports may be accessed by opening the newly generated
``htmlcov/index.html`` in your favorite web brower.  These commands
require the ``pytest`` and ``coverage`` packages to be installed.

Ignoring lines in coverage tests
--------------------------------

Occasionally there will be lines of code that do not require tests

The ``.coveragerc`` file in the top level directory contains

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
