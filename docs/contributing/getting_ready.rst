.. _getting ready to contribute:

***************************
Getting Ready to Contribute
***************************

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Introduction
============

This page describes the tasks needed prior to contributing to PlasmaPy.
After finishing these tasks, the |code contribution workflow| describes
the steps necessary

Pre-requisites
==============

Using a terminal
----------------

The commands described below are intended for a terminal running the
`Unix shell`_. Here are some essential `Unix commands`_.

* For Windows users, we recommend `Windows Subsystem for Linux`_ (WSL).
  Another alternative is Powershell_.
* These instructions describe `opening a terminal on macOS`_.
* On Linux, a terminal can be opened using :kbd:`Ctrl + Alt + t`.

Installing Python
-----------------

PlasmaPy requires Python_ |minpython| or newer. These instructions
describe how to `download Python`_ and install it.

`Real Python`_ has a helpful guide on `setting up a Python coding
environment on Windows`_ which uses Powershell_.

Using git and GitHub
--------------------

Plasma code development is done using |git|_ and GitHub_. Before
contributing code to PlasmaPy, it is necessary to:

#. `Sign up on GitHub`_ for a free account.

#. `Install git`_, if necessary. The |git|_ installation can be verified
   using:

   .. code-block:: bash

      git --version

   .. note::

      |git|_ often comes pre-installed with WSL, in macOS, and on many
      Linux distributions.

#. Optionally configure |git|_ with your name and email using the
   following commands, where the name and email are replaced with your
   own.

   .. code-block:: bash

      git config --global user.name "Spacecat Q. Spacecat"
      git config --global user.email "spacecat@spacecats.com"

#. `Add a new SSH key to your GitHub account`_.

Initial setup
=============

#. Log in to GitHub_.

#. Go to `PlasmaPy's GitHub repository`_.

#. Create a fork_ of PlasmaPy by clicking on :guilabel:`Fork`, followed
   by :guilabel:`Create fork`.

#. Open a terminal, and create and/or navigate to the folder (e.g.,
   :file:`~/repos/`) in which you want to download PlasmaPy.

#. Clone_ PlasmaPy with the following command, replacing ``username``
   with your GitHub username. This will create a subdirectory called
   :file:`PlasmaPy/` containing the cloned repository.

   .. code-block:: bash

      git clone git@github.com:username/PlasmaPy.git

#. Enter the newly created directory with ``cd PlasmaPy``.

#. Add a remote_ called ``upstream`` for `PlasmaPy's GitHub repository`_
   by using the following command.

   .. code-block:: bash

      git remote add upstream git@github.com:PlasmaPy/PlasmaPy.git

   .. tip::

      The remote named ``origin`` refers to the

      The ``upstream

      .. code-block:: bash

         git remote rename origin username
         git remote rename upstream plasmapy

.. _clone: https://github.com/git-guides/git-clone
.. _fork: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks
.. _remote: https://github.com/git-guides/git-remote

#. Create a virtual environment and activate it.

#. Install PlasmaPy's requirements with:

   .. code-block:: bash

      pip install -r requirements.txt

#. Install your clone of `plasmapy` with:

   .. code-block:: bash

      pip install -e .

   The ``-e`` makes it an editable installation.

#. In the :file:`PlasmaPy/` directory, run:

   .. code-block:: bash

      pytest -m 'not slow'

pre-commit
----------

  Install pre-commit_ with:

   .. code-block:: bash

      pre-commit install

pandoc
------
