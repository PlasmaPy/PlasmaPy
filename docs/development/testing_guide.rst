******************
Testing Guidelines
******************

Testing
=======

PlasmaPy uses pytest for its testing needs.

The simplest way to run tests on your development branch, assuming you
submitted it as a pull request, is to just push your code to GitHub.
Pull requests are tested using Travis CI, CircleCI, and AppVeyor to make
sure that the code passes tests and the documentation is able to be built.

To run the test suite locally, use ``python setup.py test`` from the
repository's root directory.  This should detect our ``setup.cfg`` file and thus
also verify examples in docstrings, the narrative documentation. Doing it
this way also


A few guidelines to on writing neat, readable and useful tests:

* Unit tests should be provided for all methods when possible.

* Solved bugs should be turned into test cases.
  
* The Travis CI, CircleCI, and AppVeyor integrations on GitHub run
  tests whenever pull requests to PlasmaPy are updated.  The pytest
  module may used on a local computer.
  
* Tests are run frequently during code development, and slow tests may
  interrupt the flow of a contributor.  Tests should be minimal,
  sufficient enough to be complete, and as efficient as possible.


Assert statements
=================
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
===============================================

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

Code coverage
=============

PlasmaPy uses the coverage.py addon via Coveralls.io. At the end of
every Travis CI testing session, information on which lines were
executed in the test is sent to Coveralls.io. At the very least, try
to avoid test coverage decreasing if possible.

To run coverage.py locally, run ``python setup.py test --coverage``, then
generate a HTML description with ``coverage html``.

At the time of writing this, coverage.py has a known issue with being
unable to check lines executed in Numba JIT compiled functions.
  
Occasionally there will be some lines that do not require testing.
For example, testing exception handling for an `ImportError` when
importing an external package would usually be impractical.  In these
instances, we may end a line with `# coveralls: ignore` to indicate
that these lines should be excluded from coverage reports (or add a
line to `.coveragerc`).  This strategy should be used sparingly, since
it is often better to explicitly test exceptions and warnings and to
show the lines of code that are not tested.
