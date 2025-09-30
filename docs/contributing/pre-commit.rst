.. _using-pre-commit:

****************
Using pre-commit
****************

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Introduction
============

PlasmaPy uses |pre-commit| to run code quality checks (via _hooks_
defined in |.pre-commit-config.yaml|_) and apply automatic fixes, when
available. All |pre-commit| hooks are run for every pull request (PR).

.. important::

   Automatically fix most |pre-commit| failures by commenting

      pre-commit.ci autofix

   in the :guilabel:`Conversation` tab of a PR. This comment triggers a
   new commit that applies the available automatic fixes. Use
   :bash:`git pull` to bring changes on GitHub back to your computer.

   For confusing

Running pre-commit locally
==========================

After `installing uv`_, run all |pre-commit| hooks in your clone of the
repository with:

.. code-block:: bash

   uvx pre-commit run -a

The ``-a`` is short for ``--all-files``. This

.. tip::

   If you get tired of the verbose output of |pre-commit|, consider
   using |prek|: an emerging drop-in replacement to |pre-commit| with
   better performance and improved user experience. To only show failure
   and error messages, run this command instead:

   .. code-block:: bash

      uvx prek run -a -q

   The ``-q`` is short for ``--quiet``.

.. _pre-commit-troubleshooting:

Troubleshooting pre-commit failures
===================================

The following sections contain suggestions for how to fix pre-commit
failures that were not corrected by commenting ``pre-commit.ci autofix``
on the issue.

ruff
----

PlasmaPy uses |ruff| as its primary linter and code formatter. |ruff|
quickly finds and often fixes many code quality issues.

Every issue detected by |ruff| corresponds to a specific linter rule. For
example, lint rule F401_ removes unused :py:`import` statements.

.. tip::

   Find more information about issues flagged by |ruff| by searching
   `ruff's documentation page on rules`_ for the rule code and clicking
   on its name.

Disabling a ruff rule
~~~~~~~~~~~~~~~~~~~~~

While |ruff| usually suggests improvements, there are occasionally
times when a departure from a |ruff| rule is preferable.

To ignore a |ruff| rule on a specific line, append a comment of the form
:samp:`# noqa {<rule-codes>}`, where :samp:`<rule-codes>` is replaced by
the |ruff| rule code(s) to ignore.

For example, we can tell |ruff| to ignore a function with excessive code
complexity (C901_), too many branches (PLR0912_), and too many
statements (PLR0915_) by adding a ``noqa`` comment like in the following
example:

.. code-block:: python

   def overly_complicated_function():  # noqa: C901, PLR0912, PLR0915
       """A function with 97 lines and multiple nested if/else blocks."""

When writing new code, it is almost always preferable to refactor the
code to remove the error than to add a ``# noqa`` comment to ignore the
rule. Complex functions flagged by C901_ could be simplified by
extracting sections of code into separate functions that do exactly one
thing with no side effects.

.. important::

   Use ``# noqa`` comments sparingly, and only when you have a strong
   reason to do so.

Spellchecks
-----------

PlasmaPy uses codespell_ and typos_ to spellcheck source code. While
these tools generally work well, occasionally there will be false
positives.

.. tip::

   Add false positives found by codespell_ to ``ignore-words-list`` in
   the ``[tool.codespell]`` section of :file:`pyproject.toml`.

   Add false positives found by typos_ to the ``[default.extend-words]``
   section of :file:`_typos.toml`.

.. _C901: https://docs.astral.sh/ruff/rules/complex-structure
.. _codespell: https://github.com/codespell-project/codespell
.. _F401: https://docs.astral.sh/ruff/rules/unused-import
.. _installing pre-commit: https://pre-commit.com/#installation
.. _PLR0912: https://docs.astral.sh/ruff/rules/too-many-branches
.. _PLR0915: https://docs.astral.sh/ruff/rules/too-many-statements
.. _ruff's documentation page on rules: https://docs.astral.sh/ruff/rules
.. _typos: https://github.com/crate-ci/typos

.. _`.pre-commit-config.yaml`: https://github.com/PlasmaPy/PlasmaPy/blob/main/.pre-commit-config.yaml
.. |.pre-commit-config.yaml| replace:: :file:`.pre-commit-config.yaml`
