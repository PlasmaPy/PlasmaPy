PlasmaPy v2024.10.0 (2024-10-30)
================================

New Features
------------

- Added the option to pass |inf| to the
  `~plasmapy.particles.particle_class.Particle.ionize` method of |Particle| to
  return the nucleus of the particle. (:pr:`2800`)
- Implemented
  `~plasmapy.simulation.resolution_constraints.CFL_limit_electromagnetic_yee`
  to calculate the CFL condition for electromagnetic simulations. (:pr:`2832`)
-
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.synthetic_radiograph`
  now accepts a file path to an
  HDF5 file saved by
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.Tracker`
  as input to create
  a synthetic radiograph. (:pr:`2868`)
- Save routines in `~plasmapy.simulation.particle_tracker.save_routines` now
  take an optional keyword argument ``output_basename``
  that sets the basename of the saved file(s).
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.Tracker`
  now
  also accepts an ``output_basename`` keyword, which is passed to the save
  routine. (:pr:`2868`)


Documentation Improvements
--------------------------

- Updated the docstring of `~plasmapy.formulary.dimensionless.beta` to
  explicitly
  use definitions of :math:`p_{th}` and :math:`p_{mag}`, and added links
  to `~plasmapy.formulary.misc.thermal_pressure`
  and `~plasmapy.formulary.misc.magnetic_pressure`. (:pr:`2822`)
- Restructured sections in the `~plasmapy.formulary.mathematics.rot_a_to_b`
  docstring to be consistent with other docstrings, and added 'Raises' and
  'Examples' sections. (:pr:`2824`)
- Added the :term:`force-free` definition to the |glossary|. (:pr:`2830`)
- Updated docstrings in `plasmapy.formulary` to follow the numpydoc standard.
  (:pr:`2831`)
- Added examples to the |ParticleTracker| docstring. (:pr:`2833`)
- Added a section in the installation instructions for installing PlasmaPy with
  |uv|. (:pr:`2861`)


Backwards Incompatible Changes
------------------------------

- Removed the ``optical_density`` keyword argument in
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography`.
  (:pr:`2843`)
- The property `~plasmapy.plasma.grids.AbstractGrid.recognized_quantities` of
  `~plasmapy.plasma.grids.AbstractGrid` is now a class method
  instead of a class property, as using `classmethod` and `property` decorators
  together is no longer
  allowed in Python 3.13. The syntax for accessing this dictionary has
  therefore changed
  from :py:`AbstractGrid.recognized_quantities` to
  :py:`AbstractGrid.recognized_quantities()`. (:pr:`2871`)


Bug Fixes
---------

- Patched a bug in
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography` in
  which particles stopped
  before the detector were still included in synthetic radiographs.
  (:pr:`2843`)


Internal Changes and Refactorings
---------------------------------

- Updated the release checklist. (:pr:`2784`)
- Added a |Nox| session to verify that the pinned requirements files used in
  continuous integration tests are consistent with the requirements in
  :file:`pyproject.toml`. (:pr:`2794`)
- Fixed a bug in the |Nox| session for running tests that prevented
  doctests from being run, and fixed doctest errors that were introduced
  while doctests were not enabled. (:pr:`2834`)
- Removed Numba as a project dependency. Consequently,
  `~plasmapy.formulary.frequencies.plasma_frequency_lite` and
  `~plasmapy.formulary.speeds.thermal_speed_lite` are no longer just-in-time
  compiled by Numba. (:pr:`2841`)
- Adjusted the Sphinx configuration to account for recent deprecations in Read
  the Docs. (:pr:`2857`)
- Added testing support for Python 3.13. (:pr:`2869`)
- Updated the versions of Python used in continuous integration workflows.
  (:pr:`2879`)


Additional Changes
------------------

- Tentatively reverted :pr:`2715` because it introduced doctest errors during a
  time when doctests were not enabled. (:pr:`2834`)


PlasmaPy v2024.7.0 (2024-07-21)
===============================

New Features
------------

- Implemented `~plasmapy.particles.atomic.stopping_power` to calculate stopping
  powers using the NIST's ASTAR and PSTAR data. (:pr:`2555`)
- Added ionization energy data from NIST to the |Particle| class.
  This can now be accessed using the
  `~plasmapy.particles.particle_class.Particle.ionization_energy` attribute
  from the |Particle| class. (:pr:`2657`)
- Renamed the `~plasmapy.particles.particle_class.Particle.binding_energy`
  attribute of |Particle| to
  `~plasmapy.particles.particle_class.Particle.nuclear_binding_energy` to avoid
  confusion with
  `~plasmapy.particles.particle_class.Particle.electron_binding_energy`.
  (:pr:`2693`)
- Added electron binding energy data, relying on ionization energy data from
  NIST, to the |Particle| class.
  This can now be accessed using the
  `~plasmapy.particles.particle_class.Particle.electron_binding_energy`
  attribute
  from the |Particle| class. (:pr:`2693`)
- Added a ``return_interpolator`` keyword to
  `~plasmapy.particles.atomic.stopping_power` to allow the user to specify the
  return of an interpolator function (`~scipy.interpolate.CubicSpline` under
  the hood). (:pr:`2712`)
- Added the `~plasmapy.formulary.collisions.misc.Bethe_stopping` function to
  the `~plasmapy.formulary.collisions` subpackage. (:pr:`2712`)
- Added the ability to enable particle stopping in the |ParticleTracker|.
  (:pr:`2712`)


Documentation Improvements
--------------------------

- Updated the |coding guide| with information on |static type checking|
  with |mypy|. (:pr:`2454`)
- Updated the section in the |coding guide| about requirements and
  dependencies. (:pr:`2720`)
- Updated docstrings in `plasmapy.dispersion`. (:pr:`2735`)
- Updated docstrings in `plasmapy.formulary.collisions`. (:pr:`2736`)
- Updated docstrings in `plasmapy.formulary`. (:pr:`2737`)
- Updated docstrings for `plasmapy.diagnostics` and `plasmapy.plasma.grids`.
  (:pr:`2738`)
- Added :file:`README.md` files in the :file:`src/plasmapy`, :file:`tests`,
  :file:`docs`, :file:`type_stubs`, :file:`.github/content`,
  :file:`.github/scripts`, and :file:`.github/workflows` directories. The
  contents of these files now appear as local documentation for each of these
  directories in |PlasmaPy's GitHub repository|. (:pr:`2742`)
- Automated creation of the index file for the release changelogs. The page for
  unreleased changes is included in the table of contents only if there are
  unreleased changes. (:pr:`2754`)
- Re-wrote the "Test independence and parametrization" section of the |testing
  guide| to use extremely simple math. (:pr:`2763`)
- Added functionality to generate a table of global substitutions in the
  |documentation guide|. (:pr:`2766`)
- Renamed :file:`docs/_cff_to_rst.py` to :file:`docs/_author_list_from_cff.py`.
  (:pr:`2766`)
- Based the version of PlasmaPy that gets included in development documentation
  builds on the current date and most recent git hash. (:pr:`2775`)
- Merged the release guide into the |coding guide|. (:pr:`2777`)
- Added a new page to the |contributor guide| on |many ways to contribute| to
  an open source project. (:pr:`2777`)
- Updated the |coding guide|, |testing guide|, and |documentation guide|
  within the |contributor guide|. (:pr:`2777`)
- Moved the |contributor guide| section on example Jupyter notebooks from the
  |coding guide| to the |documentation guide|. (:pr:`2777`)
- Added ``sphinxemoji`` as a |Sphinx| extension. (:pr:`2781`)

Backwards Incompatible Changes
------------------------------

- Added a ``__str__`` method to the |CustomParticle|
  class that returns the symbol of the particle if provided, and
  otherwise falls back to using ``__repr__``. (:pr:`2702`)
- Changed default keyword argument for the ``fraction_exited_threshold`` in
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.Tracker`
  and
  `~plasmapy.simulation.particle_tracker.termination_conditions.AllParticlesOffGridTerminationCondition`
  to correspond with the fraction of particles that have entered and
  subsequently exited the grids. Previously this keyword was a misnomer,
  causing the simulation to instead terminate when the specified fraction of
  particles remaining on the grids was less than or equal to the provided
  ``fraction_exited_threshold``. (:pr:`2712`)
- Convert ``particle`` to a required argument of the
  `~plasmapy.simulation.particle_tracker.particle_tracker.ParticleTracker.load_particles`
  method of |ParticleTracker|. (:pr:`2746`)


Bug Fixes
---------

- - Enabled |validate_quantities| to be compatible with postponed evaluation of
    annotations (see :pep:`563`). (:pr:`2479`) (:pr:`2506`)
- Changed the |charge number| (:math:`Z`) dependence of the ion contribution to
  the optical Thomson scattering
  spectral density function in
  `~plasmapy.diagnostics.thomson.spectral_density_lite` from :math:`Z`
  to :math:`z^2 / \bar{z}` to match Eq. 5.1.2 and following equations in
  :cite:t:`sheffield:2011`.
  The result is a small change in the ion acoustic wave spectrum for plasmas
  with multiple ion species. (:pr:`2699`)
- Add axes removed by `numpy.squeeze` to arrays in
  `~plasmapy.dispersion.analytical.mhd_waves_` (:pr:`2715`)


Internal Changes and Refactorings
---------------------------------

- Converted the tox environment for regenerating the requirements files
  used in continuous integration checks to |Nox|. (:pr:`2664`)
- Created a parametrized |Nox| session to run tests. (:pr:`2681`)
- Added |Nox| sessions to test importing PlasmaPy, validating
  :file:`CITATION.cff`,
  and building a source distribution and wheel. (:pr:`2682`)
- Switched the GitHub workflows for running tests from using tox environments
  to using |Nox| sessions. (:pr:`2685`)
- Added ``pytest-filter-subpackage`` to the ``tests`` dependency set. This
  dependency enables
  us to run, for example, ``pytest -P particles`` to invoke tests for
  `plasmapy.particles`. (:pr:`2688`)
- Added |Nox| sessions to run tests and build documentation against unreleased
  versions
  of major dependencies. (:pr:`2694`)
- Deleted :file:`tox.ini`, since all tox environments defined therein
  have been converted to |Nox| sessions. (:pr:`2694`)
- Removed :file:`requirements.txt`, along with the requirements files
  in :file:`ci_requirements/` that were used in tox environments
  that have since been replaced with |Nox| sessions. (:pr:`2694`)
- Switched over weekly tests to use |Nox| sessions rather than tox
  environments. (:pr:`2694`)
- Added the ``lint`` and ``manifest`` sessions for |Nox| to run |pre-commit| on
  all files
  and verify :file:`MANIFEST.in` with ``check-manifest``, respectively.
  (:pr:`2695`)
- Added a |Nox| session that invokes ``autotyping`` to automatically
  add |type hint annotations|, using either the ``--safe`` or
  ``--aggressive`` options. (:pr:`2696`)
- Added ``typos`` as a |pre-commit| hook to perform spellchecking. (:pr:`2700`)
- Added a condition to check if the GitHub API can be reached to be used by the
  `~plasmapy.utils.data.downloader.Downloader` object. (:pr:`2710`)
- Applied |type hint annotations| using ``autotyping``, and made other updates
  to type
  hint annotations and docstrings. (:pr:`2728`)
- Added |type hint annotations| to `plasmapy.utils.roman`. (:pr:`2733`)
- Added |type hint annotations| to ``plasmapy.utils._units_helpers``.
  (:pr:`2734`)
- Added a |Nox| session for building the changelog. (:pr:`2744`)
- Added an experimental |Nox| session for adding |type hint annotations| using
  `MonkeyType <https://github.com/Instagram/MonkeyType>`__.
  This session creates a database of variable types from running pytest, and
  then applies the observed types to a particular module. (:pr:`2747`)
- Updated |Nox| sessions, including docstrings and troubleshooting messages.
  (:pr:`2750`)
- Enabled tests to pass with ``numpy == 2.0.0``. (:pr:`2772`)


Additional Changes
------------------

- Refactored
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.Tracker`
  to use |ParticleTracker|. (:pr:`2704`)
- Included :file:`src/plasmapy/_version.py` in :file:`MANIFEST.in`. This file
  is automatically generated using ``setuptools_scm``, but is necessary for the
  version to be correct in the titles of pages in development documentation
  builds. (:pr:`2756`)
- Updated the comment that gets posted to new pull requests via a GitHub
  workflow. (:pr:`2765`)
