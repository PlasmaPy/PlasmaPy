.. _getting ready to contribute:

***************************
Getting Ready to Contribute
***************************

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

.. role:: bash(code)
   :language: bash

Introduction
============

Thank you for considering contributing to PlasmaPy — we really
appreciate it! |:seedling:| This page describes how to get set up to
contribute to PlasmaPy. After taking these steps, you'll be ready to go
through the :ref:`code contribution workflow <workflow>`. |:recycle:|

If you run into any problems, please feel free to reach out to us during
our |community meetings|.

Pre-requisites
==============

.. _opening-a-terminal:

Opening a terminal
------------------

The commands in this page are intended to be run in a :wikipedia:`Unix
<Unix shell>` terminal. If you are new to Unix, check out this `Unix
tutorial`_ and these `frequently used Unix commands`_.

These tabs describe how to open and use a terminal on different
operating systems.

.. tabs::

   .. group-tab:: Windows

      There are several options for terminals on Windows.

      * Powershell_ comes pre-installed with Windows. These instructions
        cover `opening Powershell`_. We recommend Powershell for a quick
        start, if Windows is the only operating system you use, or if
        you have not used Unix before.

      * We recommend `Windows Subsystem for Linux`_ (WSL) if you are
        familiar with Unix, you use macOS or Linux too, or you expect to
        contribute to PlasmaPy extensively. These instructions cover
        `installing WSL`_. If you choose WSL, follow the tabs for
        :guilabel:`Linux/WSL` below.

   .. group-tab:: macOS

      In the :guilabel:`Finder`, go to :guilabel:`Applications`. Enter
      the :guilabel:`Utilities` folder and double click on
      :guilabel:`Terminal`.

   .. group-tab:: Linux/WSL

      Open a terminal by using :kbd:`Ctrl + Alt + T`.

Using git and GitHub
--------------------

Code contributions to PlasmaPy are made using |git| and |GitHub|. Before
contributing code to PlasmaPy, please take the following steps:

#. `Sign up on GitHub`_ for a free account.

#. Verify that |git| is installed by
   :ref:`opening a terminal <opening-a-terminal>` and running:

   .. code-block:: bash

      git --version

   If there is an error, follow these instructions to `install git`_.

#. Optionally, configure |git| with your name with a command like:

   .. code-block:: bash

      git config --global user.name "Your Name"

   You can also configure |git| with your email with a command like:

   .. code-block:: bash

      git config --global user.email "your.email@example.com"

   You may also set your default editor with a command like the
   following, where ``notepad`` can be replaced with the name or path of
   your preferred editor:

   .. code-block:: bash

      git config --global core.editor notepad

   For different editor and configuration options, check out `git
   commands for setup and config`_.

#. `Add a new SSH key to your GitHub account`_. This step is needed for
   security and authentication.

.. _initial-setup:

Getting set up to contribute
============================

Clone the repository
--------------------

#. Log in to |GitHub|.

#. Go to |PlasmaPy's GitHub repository|.

#. Create a fork_ of PlasmaPy by clicking on :guilabel:`Fork`, followed
   by :guilabel:`Create fork`.

#. |Open a terminal|. Then create and/or navigate to the folder in which
   you want to download PlasmaPy. For example, to put PlasmaPy into a
   new directory called :file:`repos/` in your home directory (denoted
   by :bash:`~`), run:

   .. code-block:: bash

      mkdir ~/repos
      cd ~/repos

#. Clone_ the PlasmaPy repository with the following command, replacing
   ``YOUR-USERNAME`` with your GitHub username. This will create a
   subdirectory called :file:`PlasmaPy/` containing your local clone of
   the repository.

   .. code-block:: bash

      git clone git@github.com:YOUR-USERNAME/PlasmaPy.git

   .. important::

      If you have trouble connecting to GitHub, you may need to `add a
      new SSH key to your GitHub account`_.

#. Enter the newly created directory with:

   .. code-block:: bash

      cd PlasmaPy

#. Add a remote_ called ``upstream`` for |PlasmaPy's GitHub repository|
   by using the following command.

   .. code-block:: bash

      git remote add upstream git@github.com:PlasmaPy/PlasmaPy.git

   If you run :bash:`git remote -v`, you should see that :bash:`origin`
   corresponds to your fork_ and :bash:`upstream` corresponds to
   |PlasmaPy's GitHub repository|.

.. _install-uv:

Install uv
----------

|uv| is an extremely fast Python package and project manager used
ubiquitiously during PlasmaPy development.

|Open a terminal| and follow these instructions to |install uv|.

Install Nox
-----------

|Nox| is an automation tool used by PlasmaPy to run tests, build
documentation, and perform code quality checks. Install Nox with:

.. code-block:: bash

   uv tool install nox

Install pre-commit
------------------

|pre-commit| is a framework for running code quality checks and
performing automated fixes. Install pre-commit with:

.. code-block:: bash

   uv tool install pre-commit

Create a virtual environment
----------------------------

To create a virtual environment, run:

.. code-block:: bash

   uv venv

The output will provide a command to `activate the virtual environment`_.
This command is most often :bash:`source .venv/bin/activate` on
Linux, macOS, or WSL; and ``.venv\Scripts\activate`` when using
Windows PowerShell (which may need to be run as an administrator).

To sync the virtual environment with the development environment, run:

.. code-block:: bash

   uv sync

And with that, you're all set up to contribute to PlasmaPy!

.. tip::

   To install PlasmaPy and developer dependency groups into an activated
   virtual environment not managed by |uv|, enter the top-level directory
   of your local clone of PlasmaPy and run:

   .. code-block:: bash

      python -m pip install -e . --group dev

   The ``-e`` indicates that the installation is editable.

.. _Add a new SSH key to your GitHub account: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
.. _activate the virtual environment: https://docs.astral.sh/uv/pip/environments/#activating-a-virtual-environment
.. _clone: https://github.com/git-guides/git-clone
.. _creating an environment: https://www.anaconda.com/docs/tools/anaconda-navigator/tutorials/manage-environments#creating-a-new-environment
.. _download Python: https://www.python.org/downloads
.. _fork: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks
.. _frequently used Unix commands: https://www.geeksforgeeks.org/linux-unix/essential-linuxunix-commands/
.. _git commands for setup and config: https://git-scm.com/book/en/v2/Appendix-C%3A-Git-Commands-Setup-and-Config
.. _install git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
.. _install Graphviz: https://graphviz.org/download
.. _install pandoc: https://pandoc.org/installing.html
.. _installing Python: https://realpython.com/installing-python
.. _installing WSL: https://learn.microsoft.com/en-us/windows/wsl/install
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _opening Powershell: https://learn.microsoft.com/en-us/powershell/scripting/windows-powershell/starting-windows-powershell?view=powershell-7.4
.. _powershell: https://learn.microsoft.com/en-us/powershell
.. _Real Python: https://realpython.com
.. _remote: https://github.com/git-guides/git-remote
.. _sign up on GitHub: https://github.com/signup
.. _terminal user guide: https://support.apple.com/guide/terminal/welcome/mac
.. _this xkcd comic: https://xkcd.com/1987
.. _unix tutorial: https://www.hpc.iastate.edu/guides/unix-introduction/unix-tutorial-1
.. _using an environment: https://www.anaconda.com/docs/tools/anaconda-navigator/tutorials/manage-environments#using-an-environment
.. _venv: https://docs.python.org/3/library/venv.html
.. _virtual environment: https://realpython.com/python-virtual-environments-a-primer
.. _Windows Subsystem for Linux: https://learn.microsoft.com/en-us/windows/wsl
.. _WSL: https://learn.microsoft.com/en-us/windows/wsl
