.. _release guide:

*************
Release Guide
*************

This document describes the procedure for making a release of PlasmaPy.

The following is a partial list of tasks to be performed for each
release.  This list is currently under development.  Developers should
revise and expand the instructions while performing each release,
and may refer to `Astropy's release procedures
<https://docs.astropy.org/en/stable/development/releasing.html>`_ for
guidance.

Throughout this guide, ``0.9.0`` denotes the version you're releasing,
and ``0.8.0`` denotes the last released version.

Pre-pre-release
---------------

* Create an issue on GitHub for the release.

* Approximately ∼2–3 weeks before the release, announce that a feature
  freeze will occur one week before the anticipated release date. Only
  pull requests with a limited scope that do not significantly change
  functionality should be merged during the feature freeze.

* Announce a code freeze beginning ∼2–3 weekdays before the release.
  Only pull requests directly related to the release should be merged
  during the code freeze.

Pre-release
-----------

* Review and revise changelog entries to make sure that they are
  understandable and correctly categorized. For a pull request to revise
  multiple changelog entries, apply the :guilabel:`No changelog entry
  needed` label.

* Update hooks in :file:`.pre-commit-config.yaml` to the most recent
  versions, and run ``pre-commit run --all-files`` to apply all changes.
  Skip this step for bugfix releases.

* Re-run the pre-executed notebooks, including those for charged
  particle radiography.

* Reserve a digital object identifier (DOI) on Zenodo_ for the new
  release using the ``team@plasmapy.org`` login.

* Update :file:`docs/about/citation.rst` for the new version, and
  include the reserved DOI.

* Update and alphabetize the author list in
  :file:`docs/about/credits.rst`, with ORCID_ numbers when possible.

* Update the author list, version, and other metadata in
  :file:`codemeta.json`.  Update the ``"identifier"`` tag with the DOI
  for the new release.

* Update :file:`.mailmap`.  (Add bash command for this?)

* Build the documentation using ``make linkcheck`` and fix broken links,
  except for the reserved DOI link in :file:`docs/about/citation.rst`
  which will not work until after the Zenodo_ record has been published.

Release
-------

.. I kept getting a "Not Found" error when using the hub tool, and I'm
   not sure why.

.. Install `hub <https://hub.github.com/>`__ (if needed), and use it to
   check that the continuous integration is passing.
   ... code-block:: Shell
      hub ci-status main -v [COMMIT]
   Here, ``[COMMIT]`` is replaced by the hash from the latest commit on
   the `main <https://github.com/PlasmaPy/PlasmaPy/commits/main>`__
   branch of `PlasmaPy's GitHub repository`_.

* Enter the :file:`PlasmaPy` directory and create a new branch for the
  release that is based off of the ``main`` branch. For a bugfix
  release, this branch should already exist.

  .. code-block:: Shell

     git checkout -b v0.9.x upstream main

  The ``upstream`` remote corresponds to `PlasmaPy's GitHub repository`_.

* Push the branch to `PlasmaPy's GitHub repository`_.

  .. code-block:: Shell

     git push -u upstream

* Turn changelog entries into a :file:`CHANGELOG.rst` file.

  .. code-block::

     towncrier build --version 0.9.0

  When asked about removing changelog entries, do so.

.. Turn changelog entries into a :file:`CHANGELOG.rst` file via ``towncrier --version
   v0.9.0``. When asked about removing changelog entries, do so.

* Copy the relevant part of the generated :file:`CHANGELOG.rst` file
  into :file:`docs/whatsnew/0.9.0.rst`.

* Add the entry for :file:`docs/whatsnew/0.9.0.rst` in the table of
  contents in :file:`docs/whatsnew/index.rst`.

* Use one of the following two methods to add the note on new
  contributors to :file:`docs/whatsnew/0.9.0.rst`.

  * If not done previously, add a `GitHub personal access token`_ and
    install Xonsh_. Download the `SunPy Xonsh script`_, and run:

    .. code-block::

       generate_releaserst.xsh \
           0.8.0 \
           --auth \
           --project-name=plasmapy \
           --pretty-project-name=PlasmaPy \
           --author-sort=alphabet

    Note that the argument is for the previous release.  Double check
    that the above command works!!!!!!

.. double check this ↑

.. Use ``git shortlog -nse | cut -f 2 | vim -c "sort" -c "vsplit .mailmap" -c
   "windo diffthis"`` to compare the old and new :file:`.mailmap` version. Make sure
   the old addresses are preserved in the new version, then overwrite the
   existing :file:`.mailmap` file.
   This part may not be all that relevant anymore, except if we're using ``git
   shortlog``.  ← put this in pre-release?

* Commit and push your changes up until now.

* Open a pull request from the ``0.9.x`` branch to the ``main`` branch.

* Go to `Actions <https://github.com/PlasmaPy/PlasmaPy/actions>`__, and
  click on :guilabel:`Run workflow` under both the :guilabel:`CI` and
  :guilabel:`fortnightly tests`. Verify that all continuous integration
  checks are passing.

.. Make sure that tests pass and that documentation builds without issue.

.. No, really, check twice. Let the tests do their thing. You want things tip
    top, and by now, you want that cuppa tea anyway. Treat yourself! Celebrate
    the new release and let the darn tests pass.

..    If you want to do any rebase to clean up the commit history on your ``0.6.x``
    branch, now is the time to do that. Ensure that no tests broke.

* Create a GPG key, if not done previously.

* After verifying that all continuous integration checks are passing for
  a second time, tag the new version with

  .. code-block:: Shell

     git tag -s v0.9.0 -m "Version v0.9.0"

  The ``-s`` signs the commit with your GPG key.

* After verifying that all continuous integration checks are passing for
  a third time, push the tagged commit to the ``0.9.x`` branch on GitHub.

  .. code-block:: Shell

     git push --force --follow-tags upstream v0.9.x

  The ``--force`` is necessary to trigger a rebuild with the tagged
  version. Be careful during this step, as tags cannot be deleted once
  they have been pushed to GitHub.

At this point, the GitHub Actions packaging workflow should do most of
the work for you! `Ensure that the pipeline goes through.
<https://dev.azure.com/plasmapy/PlasmaPy/_build>`_. When ``sdist`` and
``wheels_universal`` finish, check PyPI_ for the new version!

Post-release
------------

* Merge the pull request from the ``0.9.x`` branch to ``main``.

* For major and minor releases, activate the new branch's version on
  `on Read the Docs <https://readthedocs.org/projects/plasmapy/versions>`_.

* In the ``0.9.x`` branch, change the line in
  :file:`binder/requirements.txt` that has ``.`` to ``plasmapy == 0.9``.

  * Open one of the binder examples in the docs for ``0.9.x``, and run
    the following commands to verify that the released version of
    PlasmaPy begins with ``0.9``.

    .. code-block:: python

       import plasmapy
       print(plasmapy.__version__)

* Merge the ``v0.9.x`` branch into the ``stable`` branch on GitHub:

  .. code-block::

     git checkout v0.9.x
     git pull
     git checkout stable
     git merge v0.9.x
     git push

* Make the release on conda-forge. The helpful conda-forge bots should
  automatically open up a PR on `conda-forge/plasmapy-feedstock
  <https://github.com/conda-forge/plasmapy-feedstock/pulls>`_. If nothing
  breaks, it'll even get automerged.

    * If tests fail, look at the :file:`recipe.yaml` file — usually it's
      either changed dependencies or the simple import tests there.

* Upload the release to the Zenodo_ record corresponding to the reserved
  DOI, making the metadata consistent with :file:`codemeta.json`.

.. As of July 2022, Zenodo doesn't have CodeMeta support but does have
   Citation File Format (CFF) support. Should we switch to CFF?

* Write a short post on the PlasmaPy release on PlasmaPy's website.

* Notify plasma physics communities about the release.

  * Post the release announcement in PlasmaPy's chat room.

  * Post the release announcement on social media sites (Twitter,
    Facebook).

  * Send the release announcement to the mailing list.

  * Post on the APS DPP Engage forum.

* Discuss how the release procedure went during the next community
  meeting.

* Update this very release guide to reflect any changes.

* Drop support for the versions of Python_ that will have been released
  more than 42 months prior to the next expected PlasmaPy release, as
  per the drop schedule in `NEP 29`_. Consider bumping the minimum
  supported versions of NumPy_ and Astropy_ too.

Compatibility with Prior Versions of Python, NumPy, and Astropy
===============================================================

PlasmaPy releases will generally abide by the following standards,
which are adapted from `NEP 29`_ for the support of old versions of
Python_, NumPy_, and Astropy_.

* PlasmaPy should support at least the minor versions of Python
  initially released 42 months prior to a planned project release date.

* PlasmaPy should support at least the 3 latest minor versions of
  Python.

* PlasmaPy should support minor versions of NumPy initially released
  in the 24 months prior to a planned project release date or the
  oldest version that supports the minimum Python version (whichever is
  higher).

* PlasmaPy should support at least the 3 latest minor versions of
  NumPy and Astropy.

The required major and minor version numbers of upstream packages may
only change during major or minor releases of PlasmaPy, and never during
patch releases.

Exceptions to these guidelines should only be made when there are major
improvements or fixes to upstream functionality or when other required
packages have stricter requirements.

.. _GitHub personal access token:
.. _`NEP 29`: https://numpy.org/neps/nep-0029-deprecation_policy.html
.. _ORCID: https://orcid.org
.. _SunPy Xonsh script: https://github.com/sunpy/sunpy/blob/v2.1dev/tools/generate_releaserst.xsh
