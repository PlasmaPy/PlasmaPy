*************
Release Guide
*************

This document describes the procedure for making a release of
PlasmaPy.  This document is under development and should be updated
during all releases.

The following is a partial list of tasks to be performed for each
release.  This list is currently under development.  Developers should
expand the instructions while performing each release, and may refer to
`Astropy's release procedures
<http://docs.astropy.org/en/stable/development/releasing.html>` for
guidance.

Release
-------

* Create a new branch for the release that is separate from the master
  branch, e.g. `v0.3.x`
  
* Check that the Continuous Integration is passing for the correct
  version `(see the latest commit on master)
  <https://github.com/PlasmaPy/PlasmaPy/commits/master>`_. You can use `hub
  ci-status master` with the `hub` CLI tool.

* Turn changelog entries into a `CHANGELOG.rst` file via `towncrier --version
  v0.3.0` or equivalent. When asked about removing changelog entries, do so. Ensure
  the entries are in proper categories.

* Move the generated `CHANGELOG.rst` file into
  `docs/whatsnew/{version_number}.rst`. Add the corresponding entry in the
  table of contents in `docs/whatsnew/index.rst`. 

* Add the note on include new contributors. To do this efficiently, borrow the
  SunPy Xonsh script `generate_releaserst.xsh 0.2.0 --auth
  --project-name=plasmapy --pretty-project-name=PlasmaPy`.

* Use ``git shortlog -nse | cut -f 2 | vim`` for ``.mailmap``

* Use ``astropy-tools/author_lists.py`` for ``credits.rst``

.. note:

I would think about limiting this to the credits in new release entries in
`docs/whatsnew` due to maintenance burden. ~Dominik

* Update ``setup.cfg``

  * Remove ``.dev`` from ``version = x.y.z.dev``
  * Update minimum versions of required packages, including
    ``python_requires``, ``install_requires``, and other variables

* Commit your changes up until now

* Make sure that tests pass and that documentation builds without issue (``tox``)

* Tag the new version with ``git tag -s v<version> -m "Version v<version>"``

  * Note that ``-s`` signs the commit with a GPG key

* Push the tagged commit to the version's branch on GitHub: `git push --force --follow-tags upstream v0.3.x`

At this point, `the OpenAstronomy Azure Pipelines
<https://openastronomy-azure-pipelines.readthedocs.io/en/latest/publish.html>`
infrastructure should do most of the work for you! `Ensure that the pipeline
goes through. <https://dev.azure.com/plasmapy/PlasmaPy/_build>`

* Update ``setup.cfg`` (increment version and add ``.dev`` suffix)

Post-release
------------

* If necessary (for MINOR+ and not for BUGFIX versions) activate the new
  branch's version `on RTD
  <https://readthedocs.org/projects/plasmapy/versions/>`.

* Update the `stable` branch on GitHub.

* Make the release on conda-forge

* Reserve a digital object identifier on Zenodo

* Update code metadata in ``codemeta.json``

  * The `Codemeta standard <https://codemeta.github.io/>`_ is
    relatively new, so check the standard for terms that have changed
    and new terms that may apply

* Upload the release to the Zenodo record corresponding to the reserved
  DOI

* Notify plasma physics communities about the release

* Post release announcement on social media sites

* Send release announcement to mailing list

* Update the release guide to reflect any changes

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
