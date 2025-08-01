.. _getting ready to contribute:

******************************
Getting Ready to Contribute 🎉
******************************

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

.. role:: bash(code)
   :language: bash

Introduction
============

Thank you for considering contributing to PlasmaPy — we really
appreciate it!  This page goes through the steps that first-time
contributors can take to get set up to contribute to PlasmaPy. After
taking these steps, you'll be ready to go through the :ref:`code
contribution workflow <workflow>`.

If you run into any problems, please feel free to reach out to us during
one of our |community meetings|.

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

.. _installing-python:

Installing Python
-----------------

.. note::

   PlasmaPy requires a version of Python between |minpython| and
   |maxpython|. We recommend using Python |maxpython|.

We suggest using Anaconda_ to install |Python|. Anaconda_ is a versatile
package and environment management system which is widely used in the
data science and scientific Python communities. Anaconda includes
`Anaconda Navigator`_ as its graphical user interface (GUI) and Conda_
as its command line interface (CLI).

* If you prefer a GUI, follow these instructions on `installing Anaconda
  Navigator`_.

* If you prefer a CLI, follow these instructions on `installing Conda`_.

.. note::

   There are many other equally good ways to install Python. Python's
   website describes how to `download Python`_. `Real Python`_ has
   instructions on `installing Python`_ for several different operating
   systems (if working within WSL_, follow the Linux instructions).

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
   authentication purposes.

.. _initial-setup:

Initial setup
=============

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

   .. tip::

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

Setting up a Python environment
===============================

If you plan to make multiple contributions, we recommend setting up a
Python environment specifically for PlasmaPy. This section describes how
to set up a Conda_ environment from the command line, which can be done
after installing Conda or `Anaconda Navigator`_ as described in the
section on :ref:`getting Python <installing-python>`. If you did not use
Conda or Anaconda to install Python, we suggest using a `virtual
environment`_ instead.

.. tip::

   Using Conda/virtual environments helps avoid situations as in `this
   xkcd comic`_.

#. |Open a terminal|.

#. Create a Conda environment named ``plasmapy-dev`` by running:

   .. code-block:: bash

      conda create -n plasmapy-dev python=3.12

   The :bash:`-n` flag is used to specify the name of the environment.
   The ``3.12`` can be replaced with any version of Python from
   |minpython| to |maxpython|.

#. Activate the environment with:

   .. code-block:: bash

      conda activate plasmapy-dev

   The :bash:`conda activate` command will need to be run every time you
   open a terminal, or can be added to the appropriate configuration
   file (i.e., :file:`.bashrc` for bash or :file:`.zshrc` for zsh).


Installing your clone of PlasmaPy 🦹
====================================

🏁 This section covers how to make an |editable installation| of your
clone of PlasmaPy. Making the PlasmaPy installation *editable* means
that if you modify the source code, then those changes will be included
when you :py:`import plasmapy`.

1. |Open a terminal|.

2. Navigate to the directory for your clone of PlasmaPy, which should be
   named :file:`PlasmaPy`. For example, if you ran the :bash:`git clone`
   command in the :file:`~/repos/` directory, then run:

   .. code-block:: bash

      cd ~/repos/PlasmaPy

   .. note::

      In Windows, the directory path will be :file:`C:\\Users\\<username>\\repos\\PlasmaPy`.

3. If you created a Conda_ environment for contributing to PlasmaPy,
   activate it with:

   .. code-block:: bash

      conda activate plasmapy-dev

4. Run 🏃 the command to install PlasmaPy for your operating system:

   .. tabs::

      .. group-tab:: Windows

         .. code-block:: PowerShell

            py -m pip install -e .[docs,tests]

      .. group-tab:: macOS

         .. code-block:: bash

            python -m pip install -e '.[docs,tests]'

      .. group-tab:: Linux/WSL

         .. code-block:: bash

            python -m pip install -e .[docs,tests]

   .. note::

      Replace ``py`` with ``python`` if you are not using conda.

   The :bash:`-e` flag specifies that this will be an
   |editable installation|.

   .. tip::

      If the above command does not work, try running

      .. code-block:: bash

         pip install -r requirements.txt

      This command will install that packages that PlasmaPy depends on,
      but not PlasmaPy itself.

.. hint::

   If you import a package after doing an editable installation, then
   changes made after the :py:`import` step will not be immediately
   available during a Python session. To re-import the package, use
   `importlib.reload`:

   .. code-block:: pycon

      >>> from importlib import reload
      >>> import plasmapy
      >>> # now change the source code
      >>> reload(plasmapy)

.. _Add a new SSH key to your GitHub account: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
.. _Anaconda Navigator: https://www.anaconda.com/docs/tools/anaconda-navigator/main
.. _Anaconda: https://www.anaconda.com/docs/main
.. _clone: https://github.com/git-guides/git-clone
.. _Conda: https://docs.conda.io
.. _creating an environment: https://www.anaconda.com/docs/tools/anaconda-navigator/tutorials/manage-environments#creating-a-new-environment
.. _download Python: https://www.python.org/downloads
.. _fork: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks
.. _frequently used Unix commands: https://www.geeksforgeeks.org/linux-unix/essential-linuxunix-commands/
.. _git commands for setup and config: https://git-scm.com/book/en/v2/Appendix-C%3A-Git-Commands-Setup-and-Config
.. _install git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
.. _install Graphviz: https://graphviz.org/download
.. _install pandoc: https://pandoc.org/installing.html
.. _installing Anaconda Navigator: https://www.anaconda.com/docs/tools/anaconda-navigator/install
.. _installing Conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
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
