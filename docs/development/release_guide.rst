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

Pre-release
-----------

* Check that the Continuous Integration is passing for the correct
  version `(see the latest commit on master)
  <https://github.com/PlasmaPy/PlasmaPy/commits/master>`_

* Update ``docs/about/change_log.rst``

* Update ``docs/about/release_notes.rst``

* Update ``.mailmap`` and ``docs/about/credits.rst`` to include new
  contributors

  * Use ``git shortlog -n -s -e`` for ``.mailmap``
  * Use ``astropy-tools/author_lists.py`` for ``credits.rst``

* Update ``setup.cfg``

  * Remove ``.dev`` from ``version = x.y.z.dev``
  * Update minimum versions of required packages, including
    ``python_requires``, ``install_requires``, and other variables

* Reserve a digital object identifier on Zenodo, and update citation
  information (e.g., in ``plasmapy.__citation__`` and ``README.md``)

* Update code metadata in ``codemeta.json``

  * The `Codemeta standard <https://codemeta.github.io/>`_ is
    relatively new, so check the standard for terms that have changed
    and new terms that may apply

* Make sure that tests pass  and that
  documentation builds without issue (``tox``)

* Commit your changes up until now

  * Double-check CI on pushing this to a branch

* Tag the new version with ``git tag -s v<version> -m "Tagging v<version>"``

  * Note that ``-s`` signs the commit with a GPG key

* Test that RTD is building the documentation correctly on release
  branch (and the version is correct)

* Perform a source distribution release

  * Get a clean copy of the repository (``git clean -fxd`` or ``git clone``)
  * Use ``python setup.py build sdist``
  * Test that the sdist installs in a clean environment::

       $ conda create -n plasmapy_release_test_v<version> numpy
       $ conda activate plasmapy_release_test_v<version>
       $ pip install dist/plasmapy-<version>.tar.gz
       $ python -c 'import plasmapy; plasmapy.test()'
       $ conda deactivate

  * Check that the `plasmapy.__version__` number is correct
    (``python -c 'import plasmapy; print(plasmapy.__version__)'``)

Release
-------

* Create a new branch for the release that is separate from the master
  branch
  
* Merge (via fast-forward merge with `git merge --ff-only`) changes
  from ``master`` into ``stable``

* Make sure all tests pass

* Make the release on PyPI::
    
    twine upload dist/plasmapy-X.Y.Z.tar.gz dist/plasmapy-X.Y.Z.tar.gz.asc

* Make the release on conda-forge

Post-release
------------

* Update ``docs/about/change_log.rst`` (Add a new heading for the next
  release)

* Update ``setup.cfg`` (increment version and add ``dev`` suffix)

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
