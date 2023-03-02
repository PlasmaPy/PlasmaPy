.. _workflow:

==========================
Code Contribution Workflow
==========================

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Introduction
============

This page describes the workflow for making a contribution to PlasmaPy.
This page assumes that you have finished the steps for
:ref:`getting ready to contribute`.

If you run into any problems, please feel free to reach out to us in
our `Matrix chat room`_ or during our weekly `office hours`_.

Making a code contribution
==========================

Create a new branch
-------------------

#. :ref:`Open a terminal <opening-a-terminal>`.

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
   Here ``upstream`` is the name of the remote_ and ``main`` is the name
   of the original branch.

   .. tip::

      Use descriptive branch names like ``update-contribution-workflow``
      to make it easier to remember the purpose of each branch.

#. Connect your local branch to your fork_ of PlasmaPy on GitHub_ by
   running:

   .. code-block::

      git push --set-upstream origin new-branch-name

Add and commit changes
----------------------

Next we can go through the cycle of making changes, which can be
repeated multiple times.

#. Edit a file and save the changes.

#. In a terminal, run:

   .. code-block:: bash

      git add filename

   where :samp:`{filename}` is replaced with the name of the edited
   file(s). Use ``git add *`` to add all files in the directory (except
   for files specified in :file:`.gitignore`. This step lets us line up
   the changes that we want to record as a snapshot in history.

#. To commit the changes, run:

   .. code-block:: bash

      git commit -m "<commit message>"

   where :samp:`{<commit message>}` is replaced with a descriptive
   commit message such as ``"Add gyroradius function"``.
   Committing a change is like preserving a snapshot of what each file
   looks like at this point in history.

   If it has been installed, pre-commit will perform automated checks
   and possibly make some automated changes. If pre-commit fails, then
   it'll be necessary to do the ``git add`` and ``git commit`` steps
   once more.

#. To push the changes to GitHub, run:

   .. code-block:: bash

      git push

.. tip::

   Try using the ``git status`` command after every step to get a better
   idea of what is happening.

.. note::

   The ``git`` workflow can be thought of as the process of mailing a
   package.

   * ``git add`` is like packing the contents of a package into a box.
     This step allows you to choose which changes to include in the next
     commit.

   * ``git commit`` is like sealing and labeling the package, and
     putting it in the outgoing mail.

   * ``git push`` is like sending the package off to its destination
     (i.e., GitHub).

Creating a pull request
-----------------------

#. Run ``git push`` to make sure that branch on GitHub is up-to-date.

#. Go to `PlasmaPy's GitHub repository`_.

#. If you recently pushed new changes, a pale yellow box will appear
   near the top of the screen. In that box, click
   :guilabel:`Compare & pull request`.

   .. note::

      If you did not recently push any new changes, click on
      :guilabel:`New pull request` and then the link saying "compare
      across forks." Select ``PlasmaPy/PlasmaPy`` for "base repository"
      and ``main`` for "base". Choose your fork of PlasmaPy for "head
      repository" and the name of the branch for "compare".  Then click
      on :guilabel:`Create pull request`.

#. Add a descriptive title, such as
   ``Add a function to calculate particle gyroradii``.

#. Write a description for the pull request. Describe the changes, and
   why they are being made. Include information that you think would be
   helpful for reviewers and for future reviewers and contributors.

   .. tip::

      If your pull request will close an issue, include
      :samp:`Fixes #{ISSUE-NUMBER}` in the pull request description,
      where :samp:`{ISSUE-NUMBER}` is replaced with the number of the
      issue.

#. Click on :guilabel:`Create pull request` or
   :guilabel:`Create draft pull request`.  ← update this!

#. Add a changelog entry...

#. Be sure to add/update tests and docs...

After the pull request had been created, running ``git push`` will
update the pull request.

At this stage, reviewers will perform a code review...

Note that there is a code review bottleneck...and ask people to consider
becoming code reviewers?

Code contribution tips
======================

.. danger::

   These tips are still being written!

* Choose a minor change for your very first pull request.

* Issues labeled as a `good first contribution`_ are a great place to
  get started contributing.

* In each pull request, include a set of closely related changes.

* Aim for pull requests that are ≲ 400 lines long. While longer pull
  requests are sometimes necessary, try to break up large pull requests
  into multiple pull requests, each with a set of closely related changes.

* If there is a part of code that you think could be improved, ask
  questions of reviewers.

* Avoid making pull requests from your ``main`` branch. (describe why)

.. Add something about the code review bottleneck?


.. Once you have completed your changes and pushed them to the branch on
   GitHub, you are ready to make a pull request.  Go to your fork of
   PlasmaPy in GitHub.  Select "Compare and pull request".  Add a
   descriptive title and some details about your changes.  Then select
   "Create pull request".  Other contributors will then have a chance to
   review the code and offer constructive suggestions.  You can continue
   to edit the pull request by changing the corresponding branch on your
   PlasmaPy fork on GitHub.  After a pull request is merged into the
   code, you may delete the branch you created for that pull request.


.. ... * `GitHub Documentation`_
   ...  - `Collaborating with pull requests`_
   ... * `How to Contribute to Open Source`_

.. _`Collaborating with pull requests`: https://docs.github.com/en/github/collaborating-with-pull-requests
.. _`GitHub Documentation`: https://docs.github.com/
.. _good first contribution: https://github.com/PlasmaPy/PlasmaPy/issues?q=is%3Aissue+is%3Aopen+label%3A%22Good+first+contribution%22
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
.. _remote: https://github.com/git-guides/git-remote
.. _branch: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches
.. _fork: https://docs.github.com/en/get-started/quickstart/fork-a-repo
