.. _release guide:

*************
Release Guide
*************

.. contents:: Table of Contents
   :depth: 2
   :local:
   :backlinks: none

Introduction
============

This document describes the procedure for making a release of PlasmaPy.
Developers should revise and expand these instructions while performing
each release, and may refer to `Astropy's release procedures`_ for
guidance.

Throughout this guide, ``0.9.0`` denotes the version you're releasing.

.. tip::

   Split up pre-release tasks into multiple focused pull requests to
   facilitate quicker code reviews.

.. When updating this guide, make sure that each bullet point is for
   doing exactly one task.

Announce the release timeline
=============================

* Create an issue on GitHub for the release with a checklist of tasks
  to be performed roughly one month before.

* About three weeks before a minor or major release, announce that a
  feature freeze will occur one week before the anticipated release
  date. Only pull requests with a limited scope that do not
  significantly change functionality should be merged during the feature
  freeze.

* Begin a code freeze about three weekdays before a release. Only
  bugfixes and pull requests that are directly related to the release
  should be merged during the code freeze.

Update metadata
===============

* Begin an upload to Zenodo_ for the new release using the
  ``team@plasmapy.org`` login, and reserve a digital object identifier
  (DOI).

* Open a pull request to update :file:`docs/about/citation.rst` to
  reflect the new version, and include the reserved DOI.

* Open a pull request to update and alphabetize the author list in
  :file:`docs/about/credits.rst`. Missing ORCID_ identifiers may be
  added.

* Open a pull request to update :file:`codemeta.json` and
  :file:`CITATION.cff`. Update the author list, version, and other
  metadata, as needed. Update the ``"identifier"`` tag with the DOI for
  the new release.

* Open a pull request to update :file:`.mailmap`.

.. Add a Python script here to update :file:`.mailmap`.

.. Use ``git shortlog -nse | cut -f 2 | vim -c "sort" -c "vsplit .mailmap" -c
   "windo diffthis"`` to compare the old and new :file:`.mailmap` version. Make sure
   the old addresses are preserved in the new version, then overwrite the
   existing :file:`.mailmap` file.
   This part may not be all that relevant anymore, except if we're using ``git
   shortlog``. ← put this in pre-release?

Perform code quality checks
===========================

* Create a pull request to revise changelog entries to make sure that
  they are understandable, necessary, and in the correct category.

  .. tip::

     Apply the :guilabel:`No changelog entry needed` label to pull
     requests that change multiple changelog entries in order to skip
     the changelog entry check.

* Open a pull request to ``main`` to update black_ and other pre-commit_
  hooks to their most recent versions in :file:`.pre-commit-config.yaml`
  for major and minor releases. Apply all changes by running:

  .. code-block:: bash

     pre-commit run --all-files

* Open a pull request to re-execute pre-executed notebooks, such as
  those for charged particle radiography.

* Run ``make linkcheck`` in :file:`docs/`, and if necessary, open a pull
  request to update redirects and fix broken links. The reserved DOI
  link in :file:`docs/about/citation.rst` will not work yet since the
  Zenodo_ record will not be published until after the official release.

  .. tip::

     Use ``linkcheck_allowed_redirects`` in :file:`docs/conf.py` to
     specify allowed redirects. For example, DOI links are always
     redirects, but are significantly more persistent than hyperlinks.

* Make sure that all tests are passing.

  - Go to the Actions_ page.
  - Click on the :guilabel:`CI` tab → :guilabel:`Run workflow`.
  - Click on the :guilabel:`fortnightly tests` tab →
    :guilabel:`Run workflow`.
  - Enjoy life for 15 minutes.
  - Fix any failures, and then repeat these steps.

Create the release branch
=========================

* Enter the :file:`PlasmaPy` directory and create a new branch for the
  release that is based off of the ``main`` branch. For a bugfix
  release, this branch should already exist.

  .. code-block:: bash

     git checkout -b v0.9.x upstream main

  The ``upstream`` remote corresponds to `PlasmaPy's GitHub repository`_.

* Push the branch to `PlasmaPy's GitHub repository`_.

  .. code-block:: bash

     git push -u upstream

* Open a pull request to transform the news fragments in
  :file:`changelog/` to a changelog page.

  - In the top-level directory, run:

    .. code-block:: bash

       towncrier build --version 0.9.0

    When asked about removing changelog entries, do so.

  - Copy the relevant parts of the generated :file:`CHANGELOG.rst` file
    into :file:`docs/whatsnew/0.9.0.rst`.

  - Add the entry for :file:`docs/whatsnew/0.9.0.rst` in the table of
    contents in :file:`docs/whatsnew/index.rst`.

    .. todo::

        Immediately following the ``v0.8.1`` release, we made (or
        planned to make) a few changes to the towncrier_ setup
        (:pr:`1623`, :pr:`1626`, :issue:`1627`). This guide may require
        some updates for the subsequent release.

    .. todo::

       We might be able to consolidate these steps into a single one.

* For major and minor releases, activate the new branch's version on
  `on Read the Docs <https://readthedocs.org/projects/plasmapy/versions>`_.

.. Use one of the following two methods to add the note on new
  contributors to :file:`docs/whatsnew/0.9.0.rst`.

..  If not done previously, add a `GitHub personal access token`_ and
    install Xonsh_. Download the `SunPy Xonsh script`_, and run:
    .. code-block::
       generate_releaserst.xsh \
           0.8.0 \
           --auth \
           --project-name=plasmapy \
           --pretty-project-name=PlasmaPy \
           --author-sort=alphabet
    Note that the argument is for the previous release. Double check
    that the above command works!!!!!!

Publish the release
===================

.. There used to be a step here to use the hub tool with `hub ci-status
   main -v [COMMIT]``. I kept getting a "Not Found" error when using the
   hub tool, and I'm not sure why.

.. Install `hub <https://hub.github.com/>`__ (if needed), and use it to
   check that the continuous integration is passing.
   ... code-block:: Shell
      hub ci-status main -v [COMMIT]
   Here, ``[COMMIT]`` is replaced by the hash from the latest commit on
   the `main <https://github.com/PlasmaPy/PlasmaPy/commits/main>`__
   branch of `PlasmaPy's GitHub repository`_.

* Go to the GitHub page to `draft a new release`_. We will perform a
  pre-release first.

  - Set the :guilabel:`Target` to ``v0.9.x``.
  - For :guilabel:`Choose a tag`, put ``0.9.0rc1``.
  - Under title, put ``v0.9.0rc1``.
  - Mark this as a pre-release.
  - Click on :guilabel:`Publish release`.

  The release is handled via |.github/workflows/python-publish.yml|_.

  In a few minutes, check `PlasmaPy releases on PyPI`_ to make sure that
  version ``0.9.0rc1`` has been released and is marked as pre-release.

  .. tip::

     If the release did not work, it may be necessary to create a new
     `API token for PyPI`_ and `update the secret on GitHub`_.

* Test that the new release is working. In a new virtual or conda
  environment, run

  .. code-block:: bash

     pip install plasmapy==0.9.0rc1

  to make sure that the new version installs correctly.

  - Open Python and run ``import plasmapy`` and ``dir(plasmapy)``.
  - Run ``plasma-calculator`` from the terminal to make sure that the
    plasma calculator is behaving correctly.

  Fix any errors that arise, and re-run the :guilabel:`CI` and
  :guilabel:`fortnightly tests` checks.

* Go to the GitHub page to `draft a new release`_. We will now perform
  the ``0.9.0`` release.

  - Set the :guilabel:`Target` to ``v0.9.x``.
  - For :guilabel:`Choose a tag`, put ``0.9.0``.
  - Under title, put ``v0.9.0``.
  - Copy the release notes from the changelog, using the beginning of
    :file:`docs/whatsnew/0.9.0.rst`
  - Click on :guilabel:`Publish release`.

  In a few minutes, check `PlasmaPy releases on PyPI`_ to make sure that
  the ``0.9.0`` release is present. If it is, congratulations!

  .. tip::

     To get the number of PRs merged and issues closed since the last
     release for the release notes, perform GitHub searches like
     ``is:pr merged:>=2021-11-19`` and ``is:issue closed:>=2021-11-19``,
     using the date one day after the last release.

* Merge the pull request from the ``v0.9.x`` branch to ``main``.

* In the ``v0.9.x`` branch, change the line in
  :file:`binder/requirements.txt` that has ``.`` to ``plasmapy == 0.9``.

  * Open one of the binder examples in the docs for ``v0.9.x``, and run
    the following commands to verify that the released version of
    PlasmaPy begins with ``0.9``.

    .. code-block:: python

       import plasmapy
       print(plasmapy.__version__)

* Merge the ``v0.9.x`` branch into the ``stable`` branch on GitHub:

  .. code-block:: bash

     git checkout v0.9.x
     git pull
     git checkout stable
     git merge v0.9.x
     git push

Post-release
============

* Make the release on conda-forge. The helpful conda-forge bots should
  automatically open up a PR on `conda-forge/plasmapy-feedstock
  <https://github.com/conda-forge/plasmapy-feedstock/pulls>`_. If
  nothing breaks, it'll even get auto-merged.

    * If tests fail, look at the :file:`recipe.yaml` file — usually it's
      either changed dependencies or the simple import tests there.

* Upload the release to the Zenodo_ record corresponding to the reserved
  DOI, making the metadata consistent with :file:`codemeta.json`.

Advertise the release
=====================

* Write a post on the release on `PlasmaPy's website`_.

* Notify plasma physics communities about the release on:

  * `PlasmaPy's Matrix chat room`_
  * PlasmaPy newsletter
  * Facebook_, LinkedIn_, and Twitter_
  * APS DPP Engage forum (requires login)

* Discuss how the release procedure went during the next community
  meeting.

* Update the release guide to reflect any changes.

.. |exclude bugfix| replace:: *Skip this step for bugfix releases.*

.. _Actions: https://github.com/PlasmaPy/PlasmaPy/actions
.. _API token for PyPI: https://pypi.org/help/#apitoken
.. _Astropy's release procedures: https://docs.astropy.org/en/stable/development/releasing.html
.. _Draft a new release: https://github.com/PlasmaPy/PlasmaPy/releases/new
.. _Facebook: https://www.facebook.com/people/PlasmaPy/100064083033291/
.. _LinkedIn: https://www.linkedin.com/company/plasmapy/
.. _ORCID: https://orcid.org
.. _PlasmaPy releases on PyPI: https://pypi.org/project/plasmapy/#history
.. _PlasmaPy's website: https://www.plasmapy.org
.. _SunPy Xonsh script: https://github.com/sunpy/sunpy/blob/v2.1dev/tools/generate_releaserst.xsh
.. _Twitter: https://twitter.com/plasmapy
.. _update the secret on GitHub: https://github.com/PlasmaPy/PlasmaPy/settings/secrets/actions

.. _`.github/workflows/python-publish.yml`: https://github.com/PlasmaPy/PlasmaPy/blob/main/.github/workflows/python-publish.yml
.. |.github/workflows/python-publish.yml| replace:: :file:`.github/workflows/python-publish.yml`
