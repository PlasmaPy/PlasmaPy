.. _how-to-contribute:

==========================
Code Contribution Workflow
==========================

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Introduction
============

This page describes the workflow for contributing code, documentation,
and tests to PlasmaPy. This page assumes that you have already gone
through the steps described in the page on :ref:`install-plasmapy-dev`.

:ref:`

.. getting help: Element chat, OH, community meeting

.. tip::

   Use the :guilabel:`good first contribution` label to find good
   issues to contribute to.

.. hint::

   Making multiple focused pull requests usually works better than
   making a monolithic pull request.
  
.. danger::

  Avoid making pull requests from your ``main`` branch. (describe why)
  
Making a code contribution
==========================

Creating a new branch
---------------------

#. Open a terminal.

#. Navigate to the :file:`PlasmaPy` directory that contains the clone
   of your repository.

#. Download the current status of `PlasmaPy's GitHub repository`_ and
   your fork by running:
   
   .. code-block::

      git fetch --all

#. Create and switch to a new branch_ by running:

   .. code-block::

      git checkout -b new-branch-name upstream main

   where ``new-branch-name`` is changed to the name of the new branch.
      
   .. tip::

      Use descriptive branch names like ``refactor-particle-tests`` to
      make it easier to remember the purpose of each branch.

#. Connect your local branch to your fork_ of PlasmaPy on GitHub_ by
   running:

   .. code-block::

      git push --set-upstream origin new-branch-name

Add and commit changes
----------------------

After creating a branch, we may now go through the cycle of making
changes and saving them.

#. Edit a file and save the changes.

#. In a terminal, run:

   .. code-block:: bash

      git add filename

   where :samp:`{filename}` is replaced with the name of the file (or
   ``*`` if you want to add all files in the current directory and its
   subdirectories, excluding those specified in :file:`.gitignore`).

#. 

.. hint::

   The process of adding and commiting a change and pushing it to
   GitHub can be thought of as packing and mailing a box. Using ``git
   add`` is equivalent to putting an item into a box without closing
   it. The ``git commit`` command is akin to closing the box and
   putting it in the outgoing mail. The ``git push`` command is like
   mailing the packa



.. Branches, commits, and pull requests
   ====================================

.. Before making any changes, it is prudent to update your local
   repository with the most recent changes from the development
   repository:

.. ucode-block bash

..  git fetch upstream

.. Changes to PlasmaPy should be made using branches.  It is usually best
.. to avoid making changes on your main branch so that it can be kept
.. consistent with the upstream repository. Instead we can create a new
.. branch for the specific feature that you would like to work on:

.. .. code-block:: bash

..  git branch *your-new-feature*

.. Descriptive branch names such as ``grad-shafranov`` or
.. .. ``adding-eigenfunction-poetry`` are helpful, while vague names like
.. .. ``edits`` are considered harmful.  After creating your branch locally,
.. let your fork of PlasmaPy know about it by running:

.. .. code-block:: bash

..  git push --set-upstream origin *your-new-feature*

.. It is also useful to configure git so that only the branch you are
.. working on gets pushed to GitHub:

.. .. code-block:: bash

..  git config --global push.default simple

.. Once you have set up your fork and created a branch, you are ready to
   make edits to PlasmaPy.  Switch to your new branch by running:

.. .. code-block:: bash

..   git checkout *your-new-feature*

.. Go ahead and modify files with your favorite text editor.  Be sure to
   include tests and documentation with any new functionality.  We
   recommend reading about `best practices for scientific computing
   <https://doi.org/10.1371/journal.pbio.1001745>`_.  PlasmaPy uses the
   `PEP 8 style guide for Python code
   <https://www.python.org/dev/peps/pep-0008/>`_ and the `numpydoc format
   for docstrings
   <https://github.com/numpy/numpy/blob/main/doc/HOWTO_DOCUMENT.rst.txt>`_
   to maintain consistency and readability.  New contributors should not
   worry too much about precisely matching these styles when first
.. submitting a pull request, GitHub Actions will check pull requests
   for :pep:`8` compatibility, and further changes to the style can be
   suggested during code review.

.. You may periodically commit changes to your branch by running

.. .. code-block:: bash

..  git add filename.py
..  git commit -m "*brief description of changes*"

.. Committed changes may be pushed to the corresponding branch on your
.. GitHub fork of PlasmaPy using

.. .. code-block:: bash

..  git push origin *your-new-feature*

.. or, more simply,

.. .. code-block:: bash

..   git push

.. Once you have completed your changes and pushed them to the branch on
   GitHub, you are ready to make a pull request.  Go to your fork of
   PlasmaPy in GitHub.  Select "Compare and pull request".  Add a
   descriptive title and some details about your changes.  Then select
   "Create pull request".  Other contributors will then have a chance to
   review the code and offer constructive suggestions.  You can continue
   to edit the pull request by changing the corresponding branch on your
   PlasmaPy fork on GitHub.  After a pull request is merged into the 
   code, you may delete the branch you created for that pull request.


.. Beforehand
   ==========

.. 1. `Sign up for a free GitHub account <https://github.com/signup>`_
   2.


.. Create a GitHub account
   -----------------------

.. Install git
   -----------

.. Learning Python
 ---------------

.. Getting started
.. .. ===============

.. Fork the repository
   -------------------

.. Clone the repository
   --------------------

.. Set up remotes
   --------------

.. Workflow
   ========

.. Fetch recent changes
   --------------------

.. Create a new branch
   -------------------

.. Connect the branch to GitHub
   ----------------------------

.. Make changes
   ------------

.. Commit the changes
   ------------------

.. Push the changes to GitHub
   --------------------------

.. Create a pull request
   ---------------------

.. Add a changelog entry
   ---------------------

.. Code review
   -----------

.. Getting help
   ============

.. Many ways to contribute
   =======================

.. There are many ways to contribute to an open source project such as
   PlasmaPy beyond contributing code. You can create educational notebooks
   that introduce plasma concepts using PlasmaPy. You can

.. * `Request new features`_.
   * `Report bugs`_.
   * Write tutorials on how to use different PlasmaPy features.
   * Create educational notebooks that introduce plasma concepts using PlasmaPy.
   * Improve the project's documentation. 
   * Translate PlasmaPy's documentation into another language.
   * Organize events such as `Plasma Hack Week`_.


.. Resources
   ========

.. * `GitHub Documentation`_
     - `Collaborating with pull requests`_
   * `How to Contribute to Open Source`_

.. _`Collaborating with pull requests`: https://docs.github.com/en/github/collaborating-with-pull-requests
.. _`GitHub Documentation`: https://docs.github.com/
.. _`How to Contribute to Open Source`: https://opensource.guide/how-to-contribute/
.. _`Plasma Hack Week`: https://hack.plasmapy.org
.. _`Request new features`: https://github.com/PlasmaPy/PlasmaPy/issues/new?assignees=&labels=&template=Feature_request.md
.. _`Report bugs`: https://github.com/PlasmaPy/PlasmaPy/issues/new?assignees=&labels=&template=Bug_report.md
.. _real python: https://realpython.com/python-coding-setup-windows/
.. _Add a new SSH key to your GitHub account: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
.. _install git: https://github.com/git-guides/install-git
.. _sign up on GitHub: https://github.com/join
.. _opening a terminal on macOS: https://support.apple.com/guide/terminal/open-or-quit-terminal-apd5265185d-f365-44cb-8b09-71a064a42125/mac
.. _Powershell: https://learn.microsoft.com/en-us/powershell/
.. _Unix commands: https://www.unixtutorial.org/basic-unix-commands
.. _Unix shell: https://en.wikipedia.org/wiki/Unix_shell
.. _Windows Subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/install
