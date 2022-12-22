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

      Before opening a Unix terminal on Windows, it is necessary to
      install a version of it first. We recommend `Windows Subsystem for
      Linux`_ (WSL).  Please follow these instructions on `installing
      WSL`_.

   .. group-tab:: macOS

      In the :guilabel:`Finder`, go to :guilabel:`Applications`. Enter
      the :guilabel:`Utilities` folder and double click on
      :guilablel:`Terminal`.

   .. group-tab:: Linux

      Open a terminal by using :kbd:`Ctrl + Alt + T`.

Installing Python
-----------------

PlasmaPy requires Python_ |minpython| or newer. Python's website
contains instructions on how to `download Python`_.

`Real Python`_ has instructions on `installing Python`_ for several
operating systems. If using WSL, follow the Linux instructions.

Anaconda_ is widely used as a package management system and environment
management system in the data science and scientific Python communities.
Anaconda includes `Anaconda Navigator`_ as its graphical user interface
(GUI) and Conda_ as its command line interface (CLI). Anaconda can be
used to distribute packages in Python and many other languages. Anaconda
also makes it easy to install different versions of Python into
different environments. If you would like to use Anaconda, please follow
these instructions on `installing Anaconda Navigator`_ or on `installing
Conda`_.

.. tip::

   New versions of Python_ are released annually in October, and usually
   include improved error messages and performance improvements. Because
   it can take a few months for the scientific Python ecosystem to catch
   up, we recommend installing the most recent version of Python that is
   at least four months old.

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

   You may also set your default editor:

   .. code-block:: bash

      git config --global core.editor emacs

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

#. Create a virtual environment and activate it. If you installed Python
   by downloading the link from the website,

   .. tabs::

      .. tab:: venv

         Add instructions here...

      .. tab:: Anaconda Navigator

         Add instructions or links here...

      .. tab::

         Add instructions here

#. Install your clone of `plasmapy` with:

   .. tabs::

      .. group-tab:: Windows

         .. code-block:: bash

            py -m pip install -e .[docs,tests]

      .. group-tab:: macOS

         .. code-block:: bash

            python -m pip install -e .[docs,tests]

      .. group-tab:: Linux

         .. code-block:: bash

            python -m pip install -e .[docs,tests]

   The ``-e`` makes it an editable installation, the ``.`` refers to the
   current directory, and ``[docs,tests]`` indicates that `pip`

#. In the :file:`PlasmaPy/` directory, run:

   .. code-block:: bash

      pytest -m 'not slow'


Optional tasks
==============

Installing packages needed to build documentation

If you plan to build the documentation, it may be necessary to `install
pandoc`_ and `install Graphviz`_. This step may be skipped if you do not
plan to build the documentation locally.

Installing pre-commit
---------------------

  Install pre-commit_ with:

   .. code-block:: bash

      pre-commit install


Choosing a

.. _Add a new SSH key to your GitHub Account: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
.. _Anaconda Navigator: https://docs.anaconda.com/navigator/
.. _clone: https://github.com/git-guides/git-clone
.. _fork: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks
.. _frequently used Unix commands: https://faculty.tru.ca/nmora/Frequently%20used%20UNIX%20commands.pdf
.. _download Python: https://www.python.org/downloads/
.. _install git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
.. _install Graphviz: https://graphviz.org/download/
.. _install pandoc: https://pandoc.org/installing.html
.. _installing WSL: https://learn.microsoft.com/en-us/windows/wsl/install
.. _installing Anaconda Navigator: https://docs.anaconda.com/navigator/install
.. _installing Conda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _remote: https://github.com/git-guides/git-remote
.. _sign up on GitHub: https://github.com/join
.. _terminal user guide: https://support.apple.com/guide/terminal/welcome/mac
.. _Windows Subsystem for Linux: https://learn.microsoft.com/en-us/windows/wsl



.. _Anaconda: https://docs.anaconda.com/

.. _installing Python: https://realpython.com/installing-python/
.. _Real Python: https://realpython.com/
