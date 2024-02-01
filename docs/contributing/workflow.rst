.. _workflow:

==========================
Code Contribution Workflow
==========================

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

.. role:: bash(code)
   :language: bash

Introduction
============

This page describes the workflow for making a contribution to PlasmaPy
via a `pull request`_ after having finished the steps for
|getting ready to contribute|.

If you run into any problems, please feel free to reach out to us in our
|Matrix chat room| or during our weekly |office hours|. Thank you for
contributing!

.. tip::

   Issues labeled as a `good first issue`_ are a great place to get
   started with contributing.

Making a code contribution
==========================

.. _create-branch:

Create a new branch
-------------------

#. |Open a terminal|.

#. Navigate to the :file:`PlasmaPy/` directory that contains the clone
   of your repository.

#. In the terminal, run:

   .. code-block:: bash

      git status

   If the output ends with ``nothing to commit, working tree clean``,
   then proceed to the next step.

   .. collapse:: Click here if there are uncommitted changes

      .. tip::

         If :bash:`git status` shows that any files are listed under
         ``Changes not staged for commit`` or ``Changes to be
         committed``, then do one of the following before proceeding
         to the next step:

         #. :ref:`Add and commit changes <commit-changes>`,

         #. Use `git stash`_ to temporarily file away the changes, or

         #. Use :bash:`git reset --hard` to **permanently** remove all
            changes to tracked files and return to the previous
            commit.

         If there are untracked files present, then you may delete the
         untracked files, :ref:`add and commit changes <commit-changes>`,
         or proceed to the next step.

#. Download the current status of |PlasmaPy's GitHub repository| and
   your fork by running:

   .. code-block:: bash

      git fetch --all

#. Create and switch to a new branch_ by running:

   .. code-block:: bash

      git checkout -b new-branch-name upstream/main

   where ``new-branch-name`` is changed to the name of the new branch.
   Here ``upstream`` is the name of the remote_ and ``main`` is the name
   of the original branch which the new branch will be based off of.

   .. tip::

      Use descriptive branch names like ``update-contribution-workflow``.

#. Connect your local branch to your fork_ of PlasmaPy on |GitHub| by
   running:

   .. code-block:: bash

      git push --set-upstream origin new-branch-name

.. _commit-changes:

Add and commit changes
----------------------

Next we can go through the cycle of making changes, which is usually
repeated multiple times. To get a better idea of what is being done in
each step, try running ``git status``.

#. Edit a file and save the changes.

#. In a terminal, navigate to the directory with the changed file and
   run:

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


   .. hint::

      If it has been installed, |pre-commit| will perform automated
      checks and possibly auto-fixes. If pre-commit fails, then
      it'll be necessary to fix any remaining problems and do the
      ``git add`` and ``git commit`` steps once more. Try using
      ``git diff`` and ``git diff --cached`` to view the changes, and
      :guilabel:`↑` and :guilabel:`↓` to scroll through previous
      commands in a terminal.

#. To push the changes to GitHub, run:

   .. code-block:: bash

      git push

.. tip::

   Try using the :bash:`git status` command after each step to get a
   better idea of what is happening.

.. note::

   The ``git`` workflow can be thought of as the process of mailing a
   package.

   * :bash:`git add` is like packing the contents of a package into a box.
     This step allows you to choose which changes to include in the next
     commit.

   * :bash:`git commit` is like sealing and labeling the package, and
     putting it in the outgoing mail.

   * :bash:`git push` is like sending the package off to its destination
     (i.e., GitHub).

.. _create-pr:

Creating a pull request
-----------------------

#. Run :bash:`git push` to make sure that branch on GitHub is up-to-date.

#. Go to |PlasmaPy's GitHub repository|.

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

#. Write a description for the pull request (PR). Describe the
   changes, and why they are being made. Include information that you
   think would be helpful for reviewers, future users, and future
   contributors..

   .. tip::

      If your pull request will resolve an issue, include
      :samp:`Closes #{ISSUE-NUMBER}` in the pull request description,
      where :samp:`{ISSUE-NUMBER}` is replaced with the number of the
      issue.

#. Select :guilabel:`Create pull request`.

   .. tip::

      If the pull request isn't ready for review, select the
      :guilabel:`▼` next to :guilabel:`Create pull request` to enable
      you to create a draft pull request instead.

#. :ref:`Add a changelog entry <add-changelog>`, except for minor
   changes like typo fixes.

At this stage, a reviewer will perform a code review, unless it has been
marked as a draft pull request. Thank you for contributing!

Pulling changes from GitHub
---------------------------

If your branch changes on GitHub, run

.. code-block:: bash

   git pull

to pull the changes from GitHub to your computer. If you'd like to pull
the changes from the ``main`` branch, instead run

.. code-block:: bash

   git pull upstream main

If any of the changes conflict with each other, it will be necessary to
`resolve the merge conflict`_.

.. note::

      After the pull request has been created, it can be updated by
      using :bash:`git push` to update the corresponding branch on
      GitHub.

.. important::

   If this is your first contribution, please add yourself to the author
   list in |CITATION.cff|_ (which uses |Citation File Format|) to make
   sure that you get credit for your contribution. The entry should be
   of the form:

   .. code-block:: yaml

      - given-names: <given names>
        family-names: <family names>
        affiliation: <affiliation>
        orcid: https://orcid.org/<ORCiD-iD>
        alias: <GitHub username>

   All fields are optional except ``alias``, which is your GitHub
   username. We encourage contributors to `sign up for an ORCID iD`_: a
   unique, persistent identifier used by researchers, authors, and open
   source contributors.

.. _Add a new SSH key to your GitHub account: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
.. _branch: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches
.. _fork: https://docs.github.com/en/get-started/quickstart/fork-a-repo
.. _GitHub Documentation: https://docs.github.com/
.. _git stash: https://git-scm.com/docs/git-stash
.. _good first issue: https://github.com/PlasmaPy/PlasmaPy/issues?q=is%3Aissue+is%3Aopen+label%3A%22Good+first+issue%22
.. _pull request: https://docs.github.com/en/github/collaborating-with-pull-requests
.. _remote: https://github.com/git-guides/git-remote
.. _resolve the merge conflict: https://www.atlassian.com/git/tutorials/using-branches/merge-conflicts
.. _sign up for an ORCID iD: https://orcid.org/register

.. _`CITATION.cff`: https://github.com/PlasmaPy/PlasmaPy/blob/main/CITATION.cff
.. |CITATION.cff| replace:: :file:`CITATION.cff`
