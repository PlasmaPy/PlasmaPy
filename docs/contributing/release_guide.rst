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

This document describes the procedure for making a release of
PlasmaPy.  Developers should revise and expand these instructions
while performing each release, and may refer to `Astropy's release
procedures`_ for guidance.

Throughout this guide, ``2023.9.0`` denotes the version that is being
released.

.. tip::

   Split up pre-release tasks into multiple focused pull requests to
   facilitate quicker code reviews.

.. When updating this guide, make sure that each bullet point is for
   doing exactly one task.

Announce the release timeline
=============================

* Create an issue on GitHub for the release about one month ahead of
  time.

.. Automate the above step?

* Begin a code freeze ∼3 weekdays prior to a release. Only bugfixes
  and pull requests with a limited scope that do not significantly
  impact functionality should be merged during the code freeze.

Perform code quality checks
===========================

* Create a pull request to revise changelog entries to make sure that
  they are all understandable to users, necessary, and correctly
  categorized.

  .. tip::

     Apply the :guilabel:`No changelog entry needed` label to pull
     requests that change multiple changelog entries in order to skip
     the changelog entry check.

* Run ``make linkcheck`` in :file:`docs/`, and if necessary, open a
  pull request to update redirects and fix broken links. The reserved
  DOI link in :file:`docs/about/citation.rst` will not work yet since
  the Zenodo_ record will not be published until after the official
  release.

  .. tip::

     Use ``linkcheck_allowed_redirects`` in :file:`docs/conf.py` to
     specify allowed redirects. For example, DOI_ links are always
     redirects, but are significantly more persistent than hyperlinks.

* Make sure that all tests are passing just prior to the release.

  - Go to the Actions_ page.
  - Click on the :guilabel:`CI` tab → :guilabel:`Run workflow`.
  - Click on the :guilabel:`weekly tests` tab →
    :guilabel:`Run workflow`.
  - Enjoy life for 15 minutes.
  - Fix any failures, and then repeat these steps.

Update metadata
===============

* Look through |CITATION.cff|_ and make any necessary updates to the
  author list and other metadata. It's not necessary to update the
  DOI and version number yet.

* Reserve a DOI_ for Zenodo_
  - Log on to Zenodo under the ``team@plasmapy.org`` login (available
    among core contributors).
  - Navigate to the Zenodo record for the previous release.
  - Click on :guilabel:`New version`.
  - Under :guilabel:`Basic information`, click on :guilabel:`Reserve DOI`.
  - Copy the DOI (which will be needed below).

Publish the release
===================

* Go to the `Mint release 🍬
  <https://github.com/PlasmaPy/PlasmaPy/actions/workflows/release.yml>`_
  GitHub Action. Hit the :guilabel:`Run workflow` button, fill in the
  version number and DOI from Zenodo_, and hit :guilabel:`Run Workflow`.
  Refresh the page and make sure the new job goes through. This step
  will update the DOI and version number, build the changelog, tag the
  release, and create the ``v2023.9.x`` branch.

* Go to the GitHub page to `draft a new release`_.

  - Set the :guilabel:`Target` to ``v2023.9.x``.
  - For :guilabel:`Choose a tag`, put ``2023.9.0``.
  - Under title, put ``v2023.9.0``.
  - Click on :guilabel:`Publish release`.

  The release is handled via |.github/workflows/python-publish.yml|_.

Check the release
=================

* In a few minutes, check `PlasmaPy releases on PyPI`_ to make sure
  that version ``2023.9.0`` has been released and is marked as
  pre-release.

  .. tip::

     If the release did not work, it may be necessary to create a new
     `API token for PyPI`_ and `update the secret on GitHub`_.

* Test that the new release is working. In a new virtual or conda
  environment, run

  .. code-block:: bash

     pip install plasmapy==2023.9.0

  to make sure that the new version installs correctly.

  - Open Python and run ``import plasmapy`` and ``dir(plasmapy)``.
  - Run ``plasma-calculator`` from the terminal to make sure that the
    plasma calculator is behaving correctly.

  Fix any errors that arise, and re-run the :guilabel:`CI` and
  :guilabel:`weekly tests` checks.

  .. tip::

     To get the number of PRs merged and issues closed since the last
     release for the release notes, perform GitHub searches like
     ``is:pr merged:>=2021-11-19`` and ``is:issue
     closed:>=2021-11-19``, using the date one day after the last
     release.

* Merge the pull request from the ``v2023.9.x`` branch to ``main``.

* In the ``v2023.9.x`` branch, change the line in
  :file:`binder/requirements.txt` that has ``.`` to ``plasmapy ==
  2023.9``.

  * Open one of the binder examples in the docs for ``v2023.9.x``, and
    run the following commands to verify that the released version of
    PlasmaPy begins with ``2023.9``.

    .. code-block:: python

       import plasmapy

       print(plasmapy.__version__)

Post-release
============

* For major and minor releases, activate the new branch's version on
  `on Read the Docs
  <https://readthedocs.org/projects/plasmapy/versions>`_.


* Make the release on conda-forge. The helpful conda-forge bots should
  automatically open up a PR on `conda-forge/plasmapy-feedstock
  <https://github.com/conda-forge/plasmapy-feedstock/pulls>`_. If
  nothing breaks, it'll even get auto-merged.

    * If tests fail, look at the :file:`recipe.yaml` file — usually
      it's either changed dependencies or the simple import tests
      there.

* Upload the release to the Zenodo_ record corresponding to the
  reserved DOI, making the metadata consistent with |CITATION.cff|_.

Advertise the release
=====================

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
