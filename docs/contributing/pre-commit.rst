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

Installing pre-commit
=====================

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

.. note::

   Once pre-commit has been installed for a repository, pre-commit will
   run every time you try to commit a change. If any pre-commit checks
   fail, or if pre-commit changes any files, it will be necessary to
   redo :bash:`git add` on the changed files and :bash:`git commit` once
   again.

.. tip::

   To commit a change without running pre-commit, use the :bash:`-n` or
   :bash:`--no-verify` flag with |git|_.

.. tip::

   To run :bash:`pre-commit` on all files, use

   .. code-block:: bash

      pre-commit run --all-files
