PlasmaPy v2026.2.0 (2026-02-20)
===============================

Documentation Improvements
--------------------------

- The ``mamba-org/setup-micromamba`` action is now used to install graphviz
  and pandoc during documentation build workflows on GitHub. Because this
  action caches installed packages between runs, this change speeds up
  documentation builds on GitHub by ∼1–2 minutes. (:pr:`3140`)
- Updated the acknowledgements page to include a new NASA award supporting
  work on
  `PlasmaPy infrastructure and maintenance
  <https://doi.org/10.5281/zenodo.17148678>`__. (:pr:`3148`)
- The |getting ready to contribute| page of the |contributor guide| has been
  updated to show how to set up a virtual environment for development purposes
  using |uv| and install the dependency groups used for development.
  (:pr:`3157`)
- Bumped documentation builds from using Python 3.13 to using Python 3.14.
  (:pr:`3159`)
- Made minor documentation updates and fixed a configuration issue with the
  Read the Docs build. (:pr:`3188`)
- Updated subpackage headers and docstrings in `plasmapy.formulary`.
  (:pr:`3190`)
- Updated submodule headings in `plasmapy.formulary.collisions`. (:pr:`3191`)
- Switched the documentation builds from using the ``sphinx-collapse``
  extension to using ``sphinx_toolbox.collapse``. (:pr:`3203`)


Backwards Incompatible Changes
------------------------------

- Development dependencies are now defined using dependency groups (see
  :pep:`735`) rather than via optional
  dependencies. Commands like ``pip install plasmapy[docs,tests]`` will no
  longer install the dependencies used for documentation and test builds.
  (:pr:`3157`)


Bug Fixes
---------

- Dropped usage of the ``newshape`` parameter from `numpy.reshape` that was
  deprecated in the ``v2.1.0`` and removed in the ``v2.4.0`` releases of NumPy.
  Usage of this parameter led to an error
  in `~plasmapy.analysis.swept_langmuir.helpers.merge_voltage_clusters`.
  (:pr:`3172`)
- Fixed bug in tests for `~plasmapy.diagnostics.thomson` with NumPy 2.4 where
  `numpy.random.uniform` returned
  an array that needed to be cast as a float. (:pr:`3186`)


Internal Changes and Refactorings
---------------------------------

- Updated the release checklist following the ``v2025.10.0`` release.
  (:pr:`3131`)
- Updated the workflow to verify that a PlasmaPy release can be installed via
  ``conda-forge``. (:pr:`3134`)
- Renamed :file:`.github/workflows/conda.yml` to
  :file:`.github/workflows/installability.yml`,
  and added workflows to test that a new release of PlasmaPy can be installed
  from the Python Package Index
  with pip and uv, and from conda-forge with miniconda. (:pr:`3136`)
- Updated the GitHub workflow for labeling pull requests that do not need
  changelog entries. (:pr:`3146`)
- Updated the Nox session that upgrades :file:`uv.lock` to only delete
  :file:`uv.lock` beforehand if the file is invalid. (:pr:`3150`)
- Simplified the Nox session to regenerate the :file:`uv.lock` lockfile.
  (:pr:`3174`)
- Renamed and revised |Nox| sessions for updating :file:`uv.lock` and
  validating :file:`uv.lock` against the requirements defined
  in :file:`pyproject.toml`. (:pr:`3179`)
- Adopted ``uv build`` as the build frontend. (:pr:`3202`)


Updates to Software Testing
---------------------------

- Switched the pytest configuration in :file:`pyproject.toml` to use the
  modern TOML configuration (available for ``pytest>=9``) rather than the
  legacy
  INI compatibility mode. (:pr:`3147`)
- Began using ``nox-uv`` in :file:`noxfile.py` to simplify installation of
  developer dependency groups. (:pr:`3157`)
- Removed several tests that were checking that passing a wrong type into
  certain functions raised a `TypeError`. These tests had been issuing a
  deprecation warning from Astropy related to multiplication of a unit and
  `str`
  instance. (:pr:`3187`)
- Updated the pytest configuration in :file:`pyproject.toml` and
  :file:`noxfile.py`. (:pr:`3194`)
- Separated the continuous integration checks for running tests and building
  documentation against unreleased versions of upstream dependencies into
  separate GitHub workflows. (:pr:`3195`)


Additional Changes
------------------

- Added CI workflows to check PHEP 3 compliance and PyHC Environment
  compatibility. (:pr:`3211`)
