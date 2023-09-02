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

PlasmaPy uses pre-commit_ to automate code quality checks and perform
automated fixes.

The configuration for pre-commit is in |.pre-commit-config.yaml|_.

Troubleshooting pre-commit failures
===================================

Most common pre-commit_ failures can be automatically fixed by adding
a comment on a pull request that says ``pre-commit.ci autofix`` (like in
`this comment
<https://github.com/PlasmaPy/PlasmaPy/pull/1500#issuecomment-1216865989>`__).

After doing this, please do a :bash:`git pull` in your clone of
PlasmaPy's repository to pull back the auto fixes to your computer.

If these steps do not fix the pre-commit failure, then please check out
the sections below for more suggestions on how to make the fixes.

codespell
---------

PlasmaPy uses codespell_ to find typos in source code. Rather than
checking if each word matches a dictionary entry, codespell tries to
match words to a set of common misspellings. This approach greatly
reduces the number of false positives, but will occasionally miss some
less common misspellings.

If you encounter a false positive with codespell, add it to
``ignore-words-list`` under ``[codespell]`` in :file:`pyproject.toml`.

ruff
----

PlasmaPy uses |ruff| as its primary linter. Every issue detected by ruff
corresponds to a particular lint rule. For example, lint rule F401_
checks for and removes unused imports. If you encounter a confusing ruff
rule, try searching `ruff's documentation page on rules`_ for the rule
code and clicking on its name for more information.

Problems flagged by C901_ occur when a function is too complex (i.e.,
when it contains heavily nested control flow), which makes code much
more difficult to maintain. These problems can often be fixed by
breaking up complicated functions into short functions that do exactly
one thing.

.. _C901: https://beta.ruff.rs/docs/rules/complex-structure/
.. _F401: https://beta.ruff.rs/docs/rules/unused-import
.. _ruff's documentation page on rules: https://beta.ruff.rs/docs/rules/

Using pre-commit locally
========================

Installing pre-commit
---------------------

PlasmaPy uses pre-commit_ to automate code quality checks and perform
automated fixes. Because pre-commit checks are performed on GitHub, it
is optional to set up pre-commit locally.

.. tip::

   We recommend enabling pre-commit locally on your computer after you
   become comfortable with the :ref:`code contribution workflow
   <workflow>`.

To enable pre-commit on your computer:

#. `Open a terminal <opening-a-terminal>`_.

#. Navigate to the :file:`PlasmaPy/` directory that contains your clone
   of PlasmaPy's repository. For example, if you cloned PlasmaPy into
   the :file:`~/repos/` directory, then run:

   .. code-block:: bash

      cd ~/repos/PlasmaPy

#. If you created a Conda_ environment for PlasmaPy, activate it with:

   .. code-block:: bash

      conda activate plasmapy-dev

#. Make sure that pre-commit is installed by running:

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

#. Install pre-commit with:

   .. code-block:: bash

      pre-commit install

Using pre-commit locally
------------------------

Once pre-commit_ has been installed for a repository, pre-commit will
run every time you try to commit a change.

If any pre-commit checks fail, or if pre-commit changes any files, it
will be necessary to redo :bash:`git add` on the changed files and
:bash:`git commit` once again.

.. tip::

   To commit a change without running pre-commit, use the :bash:`-n` or
   :bash:`--no-verify` flag with |git|_.

.. tip::

   To run pre-commit on all files, use

   .. code-block:: bash

      pre-commit run --all-files

.. _codespell: https://github.com/codespell-project/codespell



.. _`.pre-commit-config.yaml`: https://github.com/PlasmaPy/PlasmaPy/blob/main/.pre-commit-config.yaml
.. |.pre-commit-config.yaml| replace:: :file:`.pre-commit-config.yaml`
