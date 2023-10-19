.. _pre-commit:

****************
Using pre-commit
****************

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Introduction
============

PlasmaPy uses |pre-commit| to automate code quality checks and perform
automated fixes. The configuration for pre-commit is in
|.pre-commit-config.yaml|_.

Troubleshooting pre-commit failures
===================================

Many common |pre-commit| test failures related to formatting can be
automatically fixed by adding a comment on a pull request that says
``pre-commit.ci autofix`` (like in
`this comment
<https://github.com/PlasmaPy/PlasmaPy/pull/1500#issuecomment-1216865989>`__).
This comment will lead to a new commit to the pull request branch that
applies the automatic fixes made by the different pre-commit hooks.

After doing this, please do a :bash:`git pull` in your clone of
PlasmaPy's repository to pull back the auto fixes to your computer.

The following sections contain suggestions for how to fix pre-commit
failures that were not corrected by commenting ``pre-commit.ci autofix``
on the issue.

.. tip::

   Make sure all other tests are passing before manually fixing
   |pre-commit| test failures. Failures from pre-commit should generally
   be fixed after making sure all other tests are passing.

ruff
----

PlasmaPy uses |ruff| as its primary linter and code quality tool. |ruff|
can quickly find code quality issues and is able to do many code quality
fixes.

Every issue detected by ruff corresponds to a specific lint rule. For
example, lint rule F401_ removes unused :py:`import` statements. If you
encounter a confusing ruff rule, search `ruff's documentation page on
rules`_ for the rule code and click on its name for more information.

Problems flagged by C901_ occur when a function is too complex (i.e.,
when it contains heavily nested control flow), which makes code much
more difficult to maintain.

.. tip::

   Reduce complexity by breaking up complicated functions into short
   functions that do exactly one thing with no side effects.

Disabling a ruff rule
~~~~~~~~~~~~~~~~~~~~~

While |ruff| usually suggests improvements, there will occasionally be
times where a departure from a |ruff| rule is (at least temporarily)
justified. In these cases, we can append a :samp:`# noqa {<rule-codes>}`
comment to the end of a line (where :samp:`{<rule-codes>}` is replaced
with the corresponding |ruff| rule codes, and ``noqa`` stands for "no
quality assurance") to tell |ruff| to ignore that error on that line.

For example, we can tell |ruff| to ignore a function with excessive
code complexity (C901_), too many branches (PLR0912_), and too many
statements (PLR0915_) by adding the following ``noqa`` comment:

.. code-block:: python

   def overly_complicated_function():  # noqa: C901, PLR0912, PLR0915
       """A function with 100+ lines of code and lots of if/else branches."""

.. important::

   When writing new code, it is almost always better to refactor the
   code to remove the error rather than add a ``noqa`` comment. In the
   above example, it would be better to refactor an overly complicated
   function into multiple short functions that do exactly one thing with
   no side effects so that the code is easier to understand, modify, and
   maintain. We should only add ``noqa`` statements when we have a good
   reason to.

codespell
---------

PlasmaPy uses codespell_ to find typos in source code. Rather than
checking if each word matches a dictionary entry, codespell tries to
match words to a set of common misspellings. This approach greatly
reduces the number of false positives, but will occasionally miss some
uncommon misspellings.

If you encounter a false positive with codespell, add it to
``ignore-words-list`` under ``[codespell]`` in :file:`pyproject.toml`.

Using pre-commit locally
========================

|pre-commit| checks are performed on GitHub for every pull request, but
it is also possible to set up pre-commit locally.

.. tip::

   We recommend enabling pre-commit for the clone of
   |PlasmaPy's GitHub repository| only *after* you have become
   comfortable with the |code contribution workflow|.

Enabling pre-commit
-------------------

To enable pre-commit on your computer:

#. |Open a terminal|.

#. If you use a Conda or virtual environment for developing PlasmaPy,
   activate it (i.e., with ``conda activate plasmapy-dev``).

#. Make sure that pre-commit is installed to your Python environment by
   running:

   .. tabs::

      .. group-tab:: Windows

         .. code-block:: bash

            py -m pip install pre-commit

      .. group-tab:: macOS

         .. code-block:: bash

            python -m pip install pre-commit

      .. group-tab:: Linux/WSL

         .. code-block:: bash

            python -m pip install pre-commit

#. Navigate to the :file:`PlasmaPy/` directory that contains your clone
   of PlasmaPy's repository. For example, if you cloned PlasmaPy into
   the :file:`~/repos/` directory, then run:

   .. code-block:: bash

      cd ~/repos/PlasmaPy

#. Enable pre-commit with:

   .. code-block:: bash

      pre-commit install

Changes to the workflow
-----------------------

Once |pre-commit| has been installed for a repository, pre-commit will
run every time you try to commit a change.

If any pre-commit checks fail, or if pre-commit changes any files, it
will be necessary to redo :bash:`git add` on the changed files and
:bash:`git commit` once again.

.. tip::

   To commit a change without running pre-commit, use the :bash:`-n`
   flag (short for :bash:`--no-verify`) with |git|.

.. tip::

   To run pre-commit on all files, use

   .. code-block:: bash

      pre-commit run --all-files

.. _C901: https://docs.astral.sh/ruff/rules/complex-structure
.. _codespell: https://github.com/codespell-project/codespell
.. _F401: https://docs.astral.sh/ruff/rules/unused-import
.. _PLR0912: https://docs.astral.sh/ruff/rules/too-many-branches
.. _PLR0915: https://docs.astral.sh/ruff/rules/too-many-statements
.. _ruff's documentation page on rules: https://docs.astral.sh/ruff/rules

.. _`.pre-commit-config.yaml`: https://github.com/PlasmaPy/PlasmaPy/blob/main/.pre-commit-config.yaml
.. |.pre-commit-config.yaml| replace:: :file:`.pre-commit-config.yaml`
