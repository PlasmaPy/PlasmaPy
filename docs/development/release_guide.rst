*************
Release Guide
*************

This document describes the procedure for making a release of
PlasmaPy.  This document is under development and should be updated
during all releases.

Release Tasks
=============

The following is a partial list of tasks to be performed for each
release.  This list is currently under development.  Developers should
expand the instructions while performing each release, and may use
`Astropy's release procedures`<https://github.com/astropy/astropy/blob/master/docs/development/releasing.rst>
for guidance.

* Update ``CHANGE_LOG.rst``

* Update ``RELEASE_NOTES.rst``

* Update ``docs/credits.rst`` to include new contributors

* Edit ``plasmapy/_metadata.py`` to remove ``.dev`` from ``version``

* Update ``astropy_helpers`` to the most recent release

* Create a new branch for the release that is separate from the master
  branch

* Make sure all tests pass

* Make the release on PyPI

* Make the release on conda-forge

* Mint a release on Zenodo and get a digital object identifier (DOI)

* Alert plasma physics communities about the release
