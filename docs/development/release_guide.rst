*************
Release Guide
*************

This document describes the procedure for making a release of
PlasmaPy.  This document is under development and should be updated
during all releases.

The following is a partial list of tasks to be performed for each
release.  This list is currently under development.  Developers should
expand the instructions while performing each release, and may use
`Astropy's release procedures <http://docs.astropy.org/en/stable/development/releasing.html>`
for guidance.

Pre-release
-----------

* Check that the Continuous Integration is passing for the correct version `(see the latest commit on master) <https://github.com/PlasmaPy/PlasmaPy/commits/master>`_

* Update ``docs/about/change_log.rst``

* Update ``docs/about/release_notes.rst``

* Update ``.mailmap`` and ``docs/about/credits.rst`` to include new contributors

  * use ``git shortlog -n -s -e`` for ``.mailmap``
  * use ``astropy-tools/author_lists.py`` for ``credits.rst``

* Edit ``setup.cfg`` to remove ``.dev`` from ``version = 0.2.0.dev``

* Update ``astropy_helpers`` to the most recent release

* Make sure tests pass (``python setup.py test``) and documentation builds without issue (``python setup.py build_docs -W``)

* commit your changes up until now

  * double-check CI on pushing this to a branch

* tag the new version with ``git tag -s v<version> -m "Tagging v<version>"``

  * note that ``-s`` signs the commit with a GPG key

* test that RTD is building the documentation correctly on release branch (and the version is correct)

* perform a source distribution release

  * get a clean copy of the repository (``git clean -fxd`` or ``git clone``)
  * use ``python setup.py build sdist``
  * test that the sdist installs in a clean environment::

       $ conda create -n plasmapy_release_test_v<version> numpy
       $ conda activate plasmapy_release_test_v<version>
       $ pip install dist/plasmapy-<version>.tar.gz
       $ python -c 'import plasmapy; plasmapy.test()'
       $ conda deactivate

  * Check that the `plasmapy.__version__` number is correct (``python -c 'import plasmapy; print(plasmapy.__version__)'``)

Release
-------

* Create a new branch for the release that is separate from the master
  branch
  
* Merge (via fast-forward merge with `git merge --ff-only`) changes from `master` into `stable`

* Make sure all tests pass

* Make the release on PyPI::
    
    twine upload dist/plasmapy-X.Y.Z.tar.gz dist/plasmapy-X.Y.Z.tar.gz.asc

* Make the release on conda-forge


Post-release
------------

* Mint a release on Zenodo and get a digital object identifier (DOI)

* Notify plasma physics communities about the release

* Post release announcement on social media sites

* Send release announcement to mailing list

* Update ``CHANGELOG.rst`` (Add a new heading for the next release)
