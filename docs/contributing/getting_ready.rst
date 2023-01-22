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

Thank you for your interest in contributing to PlasmaPy! This page goes
through the steps that first-time contributors can take before
proceeding to the |code contribution workflow|.

Pre-requisites
==============

Opening a terminal
------------------

The commands described in the following sections are intended to be run
in a :wikipedia:`Unix <Unix shell>` terminal. If you are new to Unix,
check out this `Unix tutorial`_ and these `frequently used Unix
commands`_.

The following tabs describe how to use a terminal on different operating
systems.

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
        `installing WSL`_. If you choose WSL, please follow the tabs for
        :guilabel:`Linux/WSL` in the following sections.

   .. group-tab:: macOS

      In the :guilabel:`Finder`, go to :guilabel:`Applications`. Enter
      the :guilabel:`Utilities` folder and double click on
      :guilablel:`Terminal`.

   .. group-tab:: Linux/WSL

      Open a terminal by using :kbd:`Ctrl + Alt + T`.

Installing Python
-----------------

There multiple ways to install Python. PlasmaPy requires a version of
Python between |minpython| and |maxpython|. We recommend using Python
|maxpython|.

We recommend using Anaconda_, which is widely used as a package
management system and environment management system in the data science
and scientific Python communities. Anaconda includes `Anaconda
Navigator`_ as its graphical user interface (GUI) and Conda_ as its
command line interface (CLI). Anaconda allows us to install different
versions of Python into different Python, environments and can also be
used to distribute packages in languages besides Python.

* If you prefer a graphical user interface, please following these
  instructions on `installing Anaconda Navigator`_.
* If you prefer using the command line, please follow these instructions
  on `installing Conda`_.

Anaconda is not the only good way to install Python. Python's website
describes how to `download Python`_, and `Real Python`_ has instructions
on `installing Python`_ for several different operating systems (if
using WSL_, follow the Linux instructions).

Using git and GitHub
--------------------

Plasma code development is done using |git|_ and GitHub_. Before
contributing code to PlasmaPy:

#. `Sign up on GitHub`_ for a free account.

#. Verify that |git|_ is installed. Open a terminal and run:

   .. code-block:: bash

      git --version

   If there is a ``command not found`` error, then please follow these
   instructions to `install git`_.

#. Optionally configure |git|_ with your name and email using the
   following commands in the terminal, where the name and email are
   replaced with your own.

   .. code-block:: bash

      git config --global user.name "Your Name"
      git config --global user.email "your.email@example.com"

   You may also set your default editor with a command like one of the
   following:

   .. code-block:: bash

      git config --global core.editor emacs
      git config --global core.editor notepad

   For more editor options, see this page on `git commands for setup and
   config`_.

#. `Add a new SSH key to your GitHub account`_.

Initial setup
=============

#. Log in to GitHub_.

#. Go to `PlasmaPy's GitHub repository`_.

#. Create a fork_ of PlasmaPy by clicking on :guilabel:`Fork`, followed
   by :guilabel:`Create fork`.

#. Open a terminal, and navigate to or create the folder (e.g.,
   :file:`~/repos/`) in which you want to download PlasmaPy.

#. Clone_ PlasmaPy with the following command, replacing ``username``
   with your GitHub username. This will create a subdirectory called
   :file:`PlasmaPy/` containing your local clone of the repository.

   .. code-block:: bash

      git clone git@github.com:username/PlasmaPy.git

#. Enter the newly created directory with ``cd PlasmaPy``.

#. Add a remote_ called ``upstream`` for `PlasmaPy's GitHub repository`_
   by using the following command.

   .. code-block:: bash

      git remote add upstream git@github.com:PlasmaPy/PlasmaPy.git

   If you run ``git remote -v``, you should see that ``origin``
   corresponds to your fork_ and ``upstream`` corresponds to `PlasmaPy's
   GitHub repository`_.

Setting up a virtual environment
================================

If you plan to make multiple contributions to PlasmaPy, we recommend
setting up a Conda environment or a `virtual environment`_.

#. [Optional, but recommended] Create a `virtual environment`_ and
   activate it.

   .. tabs::

      .. tab:: Conda

         Create a conda_ environment named ``plasmapy`` by opening a
         terminal and running:

         .. code-block:: bash

            conda create -n plasmapy python=3.10

         The ``-n`` flag specifies the name of the environment. Activate
         this conda_ environment for your current terminal session with:

         .. code-block:: bash

            conda activate plasmapy

         This command will need to be run every time you open a
         terminal.

      .. tab:: Anaconda Navigator

         Please follow the instructions in Anaconda Navigator's
         documentation on `creating an environment`_ and then `using an
         environment`_.

      .. tab:: venv

         Python's documentation describes how to `create virtual
         environments`_.

Installing your clone of PlasmaPy
=================================

#. Use one of the following commands in the :file:`PlasmaPy/` directory
   to perform an editable (``-e``) installation of PlasmaPy, along with
   the Python packages needed to build documentation and run tests.

   .. tabs::

      .. group-tab:: Windows

         .. code-block:: bash

            py -m pip install -e .[docs,tests]

      .. group-tab:: macOS

         .. code-block:: bash

            python -m pip install -e .[docs,tests]

      .. group-tab:: Linux/WSL

         .. code-block:: bash

            python -m pip install -e .[docs,tests]

Optional tasks
==============

Requirements for building documentation
---------------------------------------

If you plan to build the documentation locally on your computer, you
might need to:

* `Install pandoc`_
* `Install Graphviz`_

These packages are not installed using the ``pip`` command above.

Installing pre-commit
---------------------

PlasmaPy uses pre-commit_ to automate code quality checks and
. Because the pre-commit checks are also performed on
GitHub, it is optional to set up pre-commit locally.

.. tip::

   We recommend installing pre-commit locally on your computer after you
   become comfortable with the |code contribution workflow|.

To enable pre-commit_ on your computer, enter the :file:`PlasmaPy/`
directory and run:

.. code-block:: bash

   pre-commit install

PlasmaPy uses |pre-commit|_ to perform code quality checks and apply
reformatting tools. The |pre-commit|_ checks are performed on every code
contribution made to GitHub,

manage and perform automated checks and
changes for code quality. The pre-commit_ checks are run on every code
contribution.

Install pre-commit_ with:

   .. code-block:: bash

      pre-commit install


Choosing a

.. _Add a new SSH key to your GitHub Account: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
.. _Anaconda Navigator: https://docs.anaconda.com/navigator/
.. _Anaconda: https://docs.anaconda.com/
.. _clone: https://github.com/git-guides/git-clone
.. _creating an environment: https://docs.anaconda.com/navigator/tutorials/manage-environments/#creating-a-new-environment
.. _download Python: https://www.python.org/downloads/
.. _fork: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks
.. _frequently used Unix commands: https://faculty.tru.ca/nmora/Frequently%20used%20UNIX%20commands.pdf
.. _git commands for setup and config: https://git-scm.com/book/en/v2/Appendix-C%3A-Git-Commands-Setup-and-Config
.. _install git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
.. _install Graphviz: https://graphviz.org/download/
.. _install pandoc: https://pandoc.org/installing.html
.. _installing Anaconda Navigator: https://docs.anaconda.com/navigator/install
.. _installing Conda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _installing Python: https://realpython.com/installing-python/
.. _installing WSL: https://learn.microsoft.com/en-us/windows/wsl/install
.. _powershell: https://learn.microsoft.com/en-us/powershell/
.. _Real Python: https://realpython.com/
.. _remote: https://github.com/git-guides/git-remote
.. _sign up on GitHub: https://github.com/join
.. _terminal user guide: https://support.apple.com/guide/terminal/welcome/mac
.. _using an environment: https://docs.anaconda.com/navigator/tutorials/manage-environments/#using-an-environment
.. _virtual environment: https://docs.python.org/3/library/venv.html
.. _Windows Subsystem for Linux: https://learn.microsoft.com/en-us/windows/wsl
.. _WSL: https://learn.microsoft.com/en-us/windows/wsl
