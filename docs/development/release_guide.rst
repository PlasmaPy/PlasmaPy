*************
Release Guide
*************

This document describes the procedure for making a release of
PlasmaPy.

The following is a partial list of tasks to be performed for each
release.  This list is currently under development.  Developers should
revise and expand the instructions while performing each release,
and may refer to `Astropy's release procedures
<http://docs.astropy.org/en/stable/development/releasing.html>`_ for
guidance.

Throughout this guide, ``0.6.0`` denotes the version you're releasing,
and ``0.5.0`` denotes the last released version.

Release
-------

* Reserve a digital object identifier (DOI) on `Zenodo <https://zenodo.org>`_
  for version ``0.6.0``.

* Update ``docs/about/citation.rst`` with the DOI for version ``0.6.0``.

* Update version metadata in ``codemeta.json``.  In particular, update the
  ``"identifier"`` tag with the DOI for version ``0.6.0``.

* Update the author list (with affiliations and ORCIDs, when possible) to be
  consistent with the Zenodo record.  Update any other tags if necessary. Check
  ``.mailmap``, ``codemeta.json``, and ``docs/about/credits.rst``

* ``hub ci-status master -v`` — Check that the Continuous Integration is passing for the correct
  version `(see the latest commit on master)
  <https://github.com/PlasmaPy/PlasmaPy/commits/master>`_. You can use the handy `hub <https://github.com/github/hub>`_ command line interface (CLI) tool.

* ``git checkout -b v0.6.x`` — create a new branch for the release that is
  separate from the master branch, with the bugfix version replaced by ``x``, for
  example, ``v0.6.x``. This is the branch for the entire series of releases — if
  you're releasing, say, ``0.6.1``, the main repository should already have a
  branch for that.

* ``git push -u upstream`` to create the branch on the main repository.

* Turn changelog entries into a ``CHANGELOG.rst`` file via ``towncrier --version
  v0.6.0``. When asked about removing changelog entries, do so. Ensure
  the entries are in proper categories.

* Copy the relevant part of the generated ``CHANGELOG.rst`` file into
  ``docs/whatsnew/0.6.0.rst``. Add the corresponding entry in the
  table of contents in ``docs/whatsnew/index.rst``.

* Add the note on new contributors to ``docs/whatsnew/{version_number}.rst``. To
  do this efficiently, borrow the `SunPy Xonsh script
  <https://github.com/sunpy/sunpy/blob/v2.1dev/tools/generate_releaserst.xsh>`_
  ``generate_releaserst.xsh 0.5.0 --auth --project-name=plasmapy
  --pretty-project-name=PlasmaPy``.

    * Note that you'll need `a GitHub personal access token
      <https://github.com/settings/tokens>`_ for that.

* Use ``git shortlog -nse | cut -f 2 | vim -c "sort" -c "vsplit .mailmap" -c
  "windo diffthis"`` to compare the old and new ``.mailmap`` version. Make sure
  the old addresses are preserved in the new version, then overwrite the
  existing ``.mailmap`` file

  .. note::

     This part may not be all that relevant anymore, except if we're using ``git
     shortlog``.

* Commit and push your changes up until now.

* Open them up as a Pull Request from the ``0.6.x`` branch to the master branch.

* Make sure that tests pass and that documentation builds without issue.

  * No, really, check twice. Let the tests do their thing. You want things tip
    top, and by now, you want that cuppa tea anyway. Treat yourself! Celebrate
    the new release and let the darn tests pass.

  * If you want to do any rebase to clean up the commit history on your ``0.6.x``
    branch, now is the time to do that. Ensure that no tests broke.

* Tag the new version with ``git tag -s v<version> -m "Version v<version>"``

  * Note that ``-s`` signs the commit with your GPG key.

* Push the tagged commit to the version's branch on GitHub: ``git push --force
  --follow-tags upstream v0.6.x``. Note that ``--force`` is necessary to trigger
  a rebuild with the tagged version. This kicked us in the posterior for ``0.4.0``.

At this point, the GitHub Actions packaging workflow should do most of the work
for you! `Ensure that the pipeline goes through.
<https://dev.azure.com/plasmapy/PlasmaPy/_build>`_. When ``sdist`` and
``wheels_universal`` finish, check `PyPI <https://pypi.org/project/plasmapy/>`_
for the new version!

Post-release
------------

* Merge the pull request from the version branch to master.

* If necessary (for MINOR+ and not for BUGFIX versions) activate the new
  branch's version `on Read the Docs
  <https://readthedocs.org/projects/plasmapy/versions/>`_.

* Update the ``stable`` branch on GitHub: ``git checkout v0.6.x; git pull; git
  checkout stable; git merge v0.6.x; git push``.

* Make the release on conda-forge. The helpful conda-forge bots should
  automatically open up a PR on `conda-forge/plasmapy-feedstock
  <https://github.com/conda-forge/plasmapy-feedstock/pulls>`_. If nothing
  breaks, it'll even get automerged.

    * If tests fail, look at the ``recipe.yaml`` file - usually it's either
      changed dependencies or the simple import tests they've got there.

* Upload the release to the Zenodo record corresponding to the reserved
  DOI

* Notify plasma physics communities about the release

  * Post release announcement on social media sites (Twitter, Facebook)

  * Send release announcement to mailing list

* Update this very release guide to reflect any changes

Compatibility with Prior Versions of Python, NumPy, and Astropy
===============================================================

PlasmaPy releases will generally abide by the following standards,
which are adapted from `NumPy Enhancement Proposal 29
<https://numpy.org/neps/nep-0029-deprecation_policy.html>`_ for the
support of old versions of Python, NumPy, and Astropy.

* PlasmaPy should support at least the minor versions of Python
  initially released 42 months prior to a planned project release date.
* PlasmaPy should support at least the 2 latest minor versions of
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
