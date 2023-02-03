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

.. _opening-a-terminal:

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


.. _installing-python:

Installing Python
-----------------

.. note::

   PlasmaPy requires a version of Python between |minpython| and
   |maxpython|. We recommend using Python |maxpython|.

We recommend using Anaconda_ to install Python_. Anaconda_ is a
versatile package and environment management system which is widely used
in the data science and scientific Python communities. Anaconda includes
`Anaconda Navigator`_ as its graphical user interface (GUI) and Conda_
as its command line interface (CLI).

* If you prefer a GUI, please following these instructions on
  `installing Anaconda Navigator`_.
* If you prefer using the command line, please follow these instructions
  on `installing Conda`_.

There are many other good ways to install Python. Python's website
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

   If there is an error, then please follow these instructions to
   `install git`_.

#. Optionally, configure |git|_ with your name with a command like:

   .. code-block:: bash

      git config --global user.name "Your Name"

   You can also configure |git|_ with your email with a command like:

   .. code-block:: bash

      git config --global user.email "your.email@example.com"

   You may also set your default editor with a command like the
   following, where ``notepad`` can be replaced with your preferred
   editor:

   .. code-block:: bash

      git config --global core.editor notepad

   For different editor options, check out `git commands for setup and
   config`_.

#. `Add a new SSH key to your GitHub account`_.

.. _initial-setup:

Initial setup
=============

#. Log in to GitHub_.

#. Go to `PlasmaPy's GitHub repository`_.

#. Create a fork_ of PlasmaPy by clicking on :guilabel:`Fork`, followed
   by :guilabel:`Create fork`.

#. :ref:`Open a terminal <opening-a-terminal>`. Then create and/or
   navigate to the folder in which you want to download PlasmaPy. For
   example, to put PlasmaPy into a new directory called :file:`repos/`
   in your home directory (denoted by ``~``), run:

   .. code-block::

      mkdir ~/repos
      cd ~/repos

#. Clone_ PlasmaPy with the following command, replacing
   ``YOUR-USERNAME`` with your GitHub username. This will create a
   subdirectory called :file:`PlasmaPy/` containing your local clone of
   the repository.

   .. code-block:: bash

      git clone git@github.com:YOUR-USERNAME/PlasmaPy.git

#. Enter the newly created directory with:

   .. code-block:: bash

      cd PlasmaPy

#. Add a remote_ called ``upstream`` for `PlasmaPy's GitHub repository`_
   by using the following command.

   .. code-block:: bash

      git remote add upstream git@github.com:PlasmaPy/PlasmaPy.git

   If you run ``git remote -v``, you should see that ``origin``
   corresponds to your fork_ and ``upstream`` corresponds to `PlasmaPy's
   GitHub repository`_.

Setting up a Python environment (optional)
==========================================

If you plan to make multiple contributions, we recommend setting up a
Python environment specifically for PlasmaPy. This section describes how
to set up a Conda_ environment from the command line, which can be done
after installing Conda_ or `Anaconda Navigator`_ as described in the
section on `getting Python <installing-python>`_.

1. `Open a terminal <opening-a-terminal>`_.

2. Create a Conda_ environment named ``plasmapy-dev`` by running

   .. code-block:: bash

      conda create -n plasmapy-dev python=3.10

   The ``-n`` flag is used to specify the name of the environment.

3. Activate the environment with

   .. code-block:: bash

      conda activate plasmapy-dev

   This command will need to be run every time you open a terminal, or
   can be added to the appropriate configuration file (i.e.,
   :file:`.bashrc` for bash or :file:`.zshrc` for zsh).

Installing your clone of PlasmaPy
=================================

.. fd#. Use one of the following commands in the :file:`PlasmaPy/` directory
   to perform an editable (``-e``) installation of PlasmaPy, along with
   the Python packages needed to build documentation and run tests.

1. `Open a terminal <opening-a-terminal>`_.

2. Navigate to the directory for your clone of PlasmaPy, which should be
   named :file:`PlasmaPy`. For example, if you ran the ``git clone``
   command in the :file:`~/repos/` directory, then run:

   .. code-block:: bash

      cd ~/repos/PlasmaPy

3. If you created a Conda_ environment for PlasmaPy, activate it with:

   .. code-block:: bash

      conda activate plasmapy-dev

4. Run the command to install PlasmaPy for your operating system:

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

   The ``-e`` specifies that this will be an `editable installation`_.

Optional tasks
==============

Requirements for building documentation
---------------------------------------

If you plan to build the documentation locally on your computer, you
might need to:

* `Install pandoc`_
* `Install Graphviz`_

These packages are not installed using the ``pip`` commands above.

Installing pre-commit
---------------------

PlasmaPy uses pre-commit_ to automate code quality checks and make
automatic fixes. Because pre-commit checks are performed on GitHub, it
is optional to set up pre-commit locally.

.. tip::

   We recommend installing pre-commit locally on your computer after you
   become comfortable with the |code contribution workflow|.

To enable pre-commit on your computer:

1. `Open a terminal <opening-a-terminal>`_.

2. Navigate to the :file:`PlasmaPy/` directory that contains your clone
   of PlasmaPy's repository. For example, if you cloned PlasmaPy into
   the :file:`~/repos/` directory, then run:

   .. code-block:: bash

      cd ~/repos/PlasmaPy

4. If you created a Conda_ environment for PlasmaPy, activate it with:

   .. code-block:: bash

      conda activate plasmapy-dev

5. Make sure that pre-commit is installed by running:

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

6. Install pre-commit with:

   .. code-block:: bash

      pre-commit install

Using pre-commit
~~~~~~~~~~~~~~~~

.. Probably need to simplify this!

Now suppose we added some trailing whitespace to :file:`some_file.py`
and attempted to commit it. If |pre-commit|_ has been installed, then
the ``trailing-whitespace`` hook will cause |pre-commit|_ to fail while
modifying :file:`some_file.py` to remove the trailing whitespace.

.. code-block:: console

   $ git add some_file.py
   $ git commit -m "Add trailing whitespace"
   Trim Trailing Whitespace.................................................Failed
   - hook id: trailing-whitespace
   - exit code: 1
   - files were modified by this hook

At this point it will be necessary to run these two commands again to
commit the changes. The changes made by |pre-commit|_ will be unstaged and
thus could be seen by running ``git diff``. Sometimes |pre-commit|_ will
not be able to automatically fix the files, such as when there are
syntax errors in Python code. In these cases, the files will need to be
changed manually before running the ``git add`` and ``git commit``
commands again. Alternatively, the |pre-commit|_ hooks can be skipped
using ``git commit -n`` instead.


After adding or updating |pre-commit|_ hooks, run the following command to
apply the changes to all files.

.. code-block:: bash

   pre-commit run --all-files



.. _add-ssh: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
.. _Anaconda Navigator: https://docs.anaconda.com/navigator/
.. _Anaconda: https://docs.anaconda.com/
.. _clone: https://github.com/git-guides/git-clone
.. _creating an environment: https://docs.anaconda.com/navigator/tutorials/manage-environments/#creating-a-new-environment
.. _download Python: https://www.python.org/downloads/
.. _editable installation: https://pip.pypa.io/en/latest/topics/local-project-installs/#editable-installs
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
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _powershell: https://learn.microsoft.com/en-us/powershell/
.. _Real Python: https://realpython.com/
.. _remote: https://github.com/git-guides/git-remote
.. _sign up on GitHub: https://github.com/join
.. _terminal user guide: https://support.apple.com/guide/terminal/welcome/mac
.. _using an environment: https://docs.anaconda.com/navigator/tutorials/manage-environments/#using-an-environment
.. _venv: https://docs.python.org/3/library/venv.html
.. _Windows Subsystem for Linux: https://learn.microsoft.com/en-us/windows/wsl
.. _WSL: https://learn.microsoft.com/en-us/windows/wsl
