.. _how-to-contribute:

=============================
How to contribute to PlasmaPy
=============================

This page describes the workflow for contributing code, documentation,
and tests to PlasmaPy.

.. getting help: Element chat, OH, community meeting

Preliminaries
=============

The steps described on this page are performed using a terminal running
the `Unix shell`_.  Here are guides on using `terminals on Linux`_ and
`terminals on macOS`_.



For
Windows users, we recommend installing Windows Subsystem for Linux.

The command on this are intended for the Unix shell.

.. How to open a terminal on macOS and Linux
.. How to install and use WSL



* `Create a GitHub account`_
* `Install git`_


Additionally, it will be necessary to `install Python`_.


.. _terminals on macOS: https://support.apple.com/guide/terminal/welcome/mac
.. _Unix shell: https://en.wikipedia.org/wiki/Unix_shell
.. _Windows Subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/install



Branches, commits, and pull requests
====================================

Before making any changes, it is prudent to update your local
repository with the most recent changes from the development
repository:

.. code-block:: bash

  git fetch upstream

Changes to PlasmaPy should be made using branches.  It is usually best
to avoid making changes on your main branch so that it can be kept
consistent with the upstream repository.  Instead we can create a new
branch for the specific feature that you would like to work on:

.. code-block:: bash

  git branch *your-new-feature*

Descriptive branch names such as ``grad-shafranov`` or
``adding-eigenfunction-poetry`` are helpful, while vague names like
``edits`` are considered harmful.  After creating your branch locally,
let your fork of PlasmaPy know about it by running:

.. code-block:: bash

  git push --set-upstream origin *your-new-feature*

It is also useful to configure git so that only the branch you are
working on gets pushed to GitHub:

.. code-block:: bash

  git config --global push.default simple

Once you have set up your fork and created a branch, you are ready to
make edits to PlasmaPy.  Switch to your new branch by running:

.. code-block:: bash

  git checkout *your-new-feature*

Go ahead and modify files with your favorite text editor.  Be sure to
include tests and documentation with any new functionality.  We
recommend reading about `best practices for scientific computing
<https://doi.org/10.1371/journal.pbio.1001745>`_.  PlasmaPy uses the
`PEP 8 style guide for Python code
<https://www.python.org/dev/peps/pep-0008/>`_ and the `numpydoc format
for docstrings
<https://github.com/numpy/numpy/blob/main/doc/HOWTO_DOCUMENT.rst.txt>`_
to maintain consistency and readability.  New contributors should not
worry too much about precisely matching these styles when first
submitting a pull request, GitHub Actions will check pull requests
for :pep:`8` compatibility, and further changes to the style can be
suggested during code review.

You may periodically commit changes to your branch by running

.. code-block:: bash

  git add filename.py
  git commit -m "*brief description of changes*"

Committed changes may be pushed to the corresponding branch on your
GitHub fork of PlasmaPy using

.. code-block:: bash

  git push origin *your-new-feature*

or, more simply,

.. code-block:: bash

  git push

Once you have completed your changes and pushed them to the branch on
GitHub, you are ready to make a pull request.  Go to your fork of
PlasmaPy in GitHub.  Select "Compare and pull request".  Add a
descriptive title and some details about your changes.  Then select
"Create pull request".  Other contributors will then have a chance to
review the code and offer constructive suggestions.  You can continue
to edit the pull request by changing the corresponding branch on your
PlasmaPy fork on GitHub.  After a pull request is merged into the
code, you may delete the branch you created for that pull request.


Beforehand
==========

1. `Sign up for a free GitHub account <https://github.com/signup>`_
2.


Create a GitHub account
-----------------------

Install git
-----------

Learning Python
---------------

Getting started
===============

Fork the repository
-------------------

Clone the repository
--------------------

Set up remotes
--------------

Workflow
========

Fetch recent changes
--------------------

Create a new branch
-------------------

Connect the branch to GitHub
----------------------------

Make changes
------------

Commit the changes
------------------

Push the changes to GitHub
--------------------------

Create a pull request
---------------------

Add a changelog entry
---------------------

Code review
-----------

Getting help
============



Many ways to contribute
=======================

There are many ways to contribute to an open source project such as
PlasmaPy beyond contributing code. You can create educational notebooks
that introduce plasma concepts using PlasmaPy. You can

* `Request new features`_.
* `Report bugs`_.
* Write tutorials on how to use different PlasmaPy features.
* Create educational notebooks that introduce plasma concepts using PlasmaPy.
* Improve the project's documentation.
* Translate PlasmaPy's documentation into another language.
* Organize events such as `Plasma Hack Week`_.


Resources
========

* `GitHub Documentation`_
  - `Collaborating with pull requests`_
* `How to Contribute to Open Source`_

.. _`Collaborating with pull requests`: https://docs.github.com/en/github/collaborating-with-pull-requests
.. _`GitHub Documentation`: https://docs.github.com/
.. _`How to Contribute to Open Source`: https://opensource.guide/how-to-contribute/
.. _`Plasma Hack Week`: https://hack.plasmapy.org
.. _`Request new features`: https://github.com/PlasmaPy/PlasmaPy/issues/new?assignees=&labels=&template=Feature_request.md
.. _`Report bugs`: https://github.com/PlasmaPy/PlasmaPy/issues/new?assignees=&labels=&template=Bug_report.md

.. _code-contribution:
