PlasmaPy v2025.10.0 (2025-10-25)
================================

Documentation Improvements
--------------------------

- Set documentation builds in continuous integration tests to use Python 3.13.
  (:pr:`3091`)


Deprecations and Planned Removals
---------------------------------

- Removed the ``binding_energy`` attribute of |Particle|, which had
  previously been deprecated in favor
  of `~plasmapy.particles.particle_class.Particle.nuclear_binding_energy`.
  (:pr:`3074`)
- Dropped support for Python 3.11. (:pr:`3110`)


Bug Fixes
---------

- Fixed an `AttributeError` that was raised when attempting to set
  the ``__doc__`` attribute of the type aliases |ParticleLike| and
  |ParticleListLike| in
  Python 3.14. The ``__doc__`` attribute is no longer set directly because
  attributes of `typing.Union` objects may no longer be set in Python 3.14.
  This error prevents previous versions of PlasmaPy from being imported in
  Python 3.14 (see :issue:`3123`). Calling `help` on |ParticleLike| or
  |ParticleListLike|
  now pulls up the docstring for a generic `typing.Union` object
  rather than the docstrings defined in the source code. The docstrings
  for these widely used type aliases are still available in the online
  documentation for
  PlasmaPy, complemented by the glossary definition for |particle-like|.
  (:pr:`3110`)


Internal Changes and Refactorings
---------------------------------

- Reduced the number of warnings reported during tests to zero when tests
  are run against the most recent versions of dependencies. (:pr:`3074`)
- Warnings issued during tests are now converted to errors, except when
  testing against the lowest allowed versions of direct dependencies.
  (:pr:`3074`)
- Updated the checklist for performing a release. (:pr:`3089`)
- Removed obsolete configuration settings from :file:`mypy.ini`. (:pr:`3094`)
- Added Python 3.14 to the continuous integration suite. (:pr:`3110`)
- Split up the weekly continuous integration tests previously
  in :file:`.github/workflows/weekly.yml` into :file:`ci-comprehensive.yml`
  for comprehensive tests, :file:`ci-upstream.yml` for tests
  and documentation builds performed with the latest unreleased versions
  of upstream dependencies, and :file:`conda.yml` for attempting an
  installation of PlasmaPy using miniconda. (:pr:`3110`)
- Removed ``pre-commit-search-and-replace`` as a |pre-commit| hook.
  (:pr:`3113`)
- Replaced the |Nox| session and GitHub check to validate :file:`CITATION.cff`
  with a |pre-commit| hook. (:pr:`3122`)
- Added a |pre-commit| hook to update :file:`uv.lock` if any changes to
  requirements were made in :file:`pyproject.toml`. This hook is not run
  on ``pre-commit.ci`` because it requires network access. (:pr:`3122`)
- Removed some |pre-commit| hooks that are covered by |ruff| rules, such
  as ``python-check-blanket-noqa``. (:pr:`3122`)
