PlasmaPy v2024.2.0 (2024-02-06)
===============================

Backwards Incompatible Changes
------------------------------

- Created new `~plasmapy.simulation.particle_tracker.ParticleTracker`
  object for general particle pushing simulations involving
  electromagnetic fields. The tracker replaces the old
  ``~plasmapy.simulation.particletracker.ParticleTracker``. (:pr:`2245`)
- Imports of classes and functions from `plasmapy.utils` must now be made
  directly from the subpackages and modules, rather than from `plasmapy.utils`
  itself. (:pr:`2403`)
- Moved ``mass_density`` and its alias ``rho_`` from `plasmapy.formulary.misc`
  to `plasmapy.formulary.densities`. (:pr:`2410`)


Features
--------

- Added a ``notch`` argument to
  `~plasmapy.diagnostics.thomson.spectral_density`,
  which allows users to output spectrum over one or multiple wavelength
  ranges to correspond to a notch filter commonly applied to experimentally
  measured Thomson scattering spectra. Changed the
  `~plasmapy.diagnostics.thomson.spectral_density_model` function to allow
  ``notch`` to be applied during fitting. (:pr:`2058`)
- Changed unit annotations to |Quantity| type hint annotations. (:pr:`2421`)
- Enabled |particle_input| to decorate functions which have variadic
  positional arguments followed by keyword arguments. See :issue:`2150`.
  (:pr:`2428`)
- Added a ``random_seed`` keyword in
  the
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.Tracker.create_particles`
  method of the
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.Tracker`
  class
  to make the function deterministic to resolve intermittent test failures.
  (:pr:`2487`)


Bug Fixes
---------

- Enabled |particle_input| to be compatible with postponed evaluation of
  annotations (see :pep:`563`). (:pr:`2479`)


Improved Documentation
----------------------

- Switched from ``sphinx_toolbox.collapse`` to ``sphinx_collapse`` for
  including collapsible content in our documentation. (:pr:`2387`)
- Added a discussion of |type hint annotations| in the |testing guide|
  within the |contributor guide|. (:pr:`2440`)


Trivial/Internal Changes
------------------------

- Added Python 3.12 to the suite of continuous integration tests. (:pr:`2368`)
- Added a GitHub Action that will create an issue containing a release
  checklist, which will largely supersede the release process outlined
  in the |contributor guide|. (:pr:`2376`)
- Separated the GitHub Action for checking hyperlinks in the
  documentation into its own GitHub Action. (:pr:`2392`)
- Replaced |black| with |ruff| in the |pre-commit| configuration. (:pr:`2394`)
- Modified the ``__exit__`` method of ``HDF5Reader`` for context management.
  (:pr:`2402`)
- Added an initial configuration for |mypy| that temporarily ignores existing
  errors. (:pr:`2424`)
- Added a |tox| environment for running |mypy|. (:pr:`2431`)
- Added |mypy| to the suite of continuous integration checks. (:pr:`2432`)
- Used ``autotyping`` to implement |type hint annotations| for special
  methods like ``__init__`` and ``__str__``, and changed ``-> typing.NoReturn``
  annotations to ``-> None``. (:pr:`2437`)
- Used ``autotyping`` to add :py:`-> None` return annotations to functions
  and methods with no :py:`return` statement. (:pr:`2439`)
- Added a stub file containing |type hint annotations| for
  ``@wrapt.decorator``. (:pr:`2442`)
- Improved |type hint annotations| for `plasmapy.particles.decorators`,
  which includes |particle_input|, and the corresponding tests. (:pr:`2443`)
- Dropped the |pre-commit| hook for ``isort`` and enabled all ``isort``
  rules in |ruff|. (:pr:`2453`)
- Added a :file:`py.typed` marker to indicate that PlasmaPy contains
  type hint annotations as per :pep:`561`. (:pr:`2473`)
- Changed ``_nearest_neighbor_interpolator`` method in `~plasmapy.plasma.grids`
  to interpolate quantities array instead of producing an intermediate index.
  (:pr:`2475`)
- Enabled the ``sphinx`` linkchecker in quiet mode to make it easier to find
  problem links from the console output. (:pr:`2476`)
- Bumped the minimum versions of dependencies to drop support for
  minor releases older than two years old. In particular, the minimum
  version of NumPy was bumped to ``1.23.0``. (:pr:`2488`)


PlasmaPy v2023.10.0 (2023-10-20)
================================

Backwards Incompatible Changes
------------------------------

- Renamed the ``plasmapy.dispersion.dispersionfunction`` module
  to `plasmapy.dispersion.dispersion_functions`. Both
  of ``plasma_dispersion_func`` and ``plasma_dispersion_func_deriv``
  are temporarily still available, but will issue
  a `~plasmapy.utils.exceptions.PlasmaPyFutureWarning` and will be
  removed in a subsequent release. (:pr:`2271`)
- Removed the |lite-functions|
  for `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func`
  and `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func_deriv`.
  Instead, the performance of the original functions has been improved
  by using a :py:`try` and :py:`except` block instead of having multiple
  :py:`if` statements to check for preconditions. (:pr:`2361`)
- Providing a real number to the ``charge`` parameter in |CustomParticle|
  will now result in a |InvalidParticleError| instead of a deprecation warning.
  Now,
  ``charge`` must be a |Quantity| with units of electrical charge. To
  express the charge as a multiple of the elementary charge, provide a
  real number to the ``Z`` parameter instead. (:pr:`2369`)


Features
--------

- Added the `~plasmapy.plasma.equilibria1d.HarrisSheet` class to calculate
  magnetic field, current density, and plasma pressure for 1D Harris sheets.
  (:pr:`2068`)
- Added module `plasmapy.dispersion.analytical.mhd_waves_`
  with classes for storing and calculating parameters of
  magnetohydrodynamic waves. (:pr:`2206`)
- Added the `plasmapy.analysis.time_series.conditional_averaging` module
  including the
  `~plasmapy.analysis.time_series.conditional_averaging.ConditionalEvents`
  class for calculating the conditional average and variance of time series.
  (:pr:`2275`)
- Added the `~plasmapy.plasma.cylindrical_equilibria.ForceFreeFluxRope`
  class to calculate magnetic field for the Lundquist solution for
  force-free cylindrical equilibria. (:pr:`2289`)


Bug Fixes
---------

- Fixed a bug that had been causing incorrect results
  in `~plasmapy.formulary.collisions.helio.collisional_analysis.temp_ratio`.
  (:pr:`2248`)
- Enabled the ``time_step`` parameter in
  `~plasmapy.analysis.time_series.excess_statistics.ExcessStatistics` class to
  be a |Quantity| with a unit. (:pr:`2300`)


Improved Documentation
----------------------

- Updated the |code contribution workflow| in the |contributor guide| to
  describe how to use ``git pull``. (:pr:`2193`)
- Expanded the troubleshooting section of the |documentation guide| to
  describe how to resolve warnings related to documents not being
  included in any toctrees. (:pr:`2257`)
- Added a step to the |code contribution workflow| about using ``git status``
  to verify that there have been no changes to tracked files before
  creating and switching to a new branch. (:pr:`2263`)
- Added a page to the |contributor guide| about |pre-commit|, including
  how to troubleshoot test failures. (:pr:`2265`)
- Added :file:`CONTRIBUTING.md` to PlasmaPy's GitHub repository.  This
  page refers contributors to PlasmaPy's |contributor guide|. (:pr:`2266`)
- Enabled the `sphinx.ext.duration` extension to show the times required
  to process different pages during documentation builds. (:pr:`2268`)
- Enabled the `sphinx.ext.viewcode` extension for adding links in the
  documentation to pages containing the source code. (:pr:`2269`)
- Moved definitions of certain |reStructuredText| substitutions
  from :file:`docs/common_links.rst` to the
  file :file:`docs/contributing/doc_guide.rst` in order to speed up the
  documentation build (see :issue:`2277`\ ). (:pr:`2272`)
- Implemented ``sphinxcontrib-globalsubs`` to enable global
  |reStructuredText| substitutions to be used throughout the documentation,
  and moved the definition of substitutions from :file:`docs/common_links.rst`
  to the ``global_substitutions`` `dict` in
  :file:`docs/_global_substitutions.py`. (:pr:`2281`)
- Changed :py:`from astropy import units as u` to :py:`import astropy.units as
  u`
  and :py:`from astropy import constants as const`
  to :py:`import astropy.constants as const` throughout the code in order to
  increase consistency of import statements. (:pr:`2282`)
- Added and applied ``nbqa-ruff`` to our suite of |pre-commit| hooks so
  that |ruff| can perform code quality checks on our example notebooks.
  (:pr:`2302`)
- Renamed :file:`docs/cff_to_rst.py` to :file:`docs/_cff_to_rst.py`, and
  updated the functionality contained within that file for converting
  author information in :file:`CITATION.cff` into a |reStructuredText|
  author list to be included in the documentation. (:pr:`2307`)
- Fixed broken hyperlinks and |reStructuredText| references. (:pr:`2308`)
- Replaced :py:`from plasmapy.particles import *` in
  :file:`docs/notebooks/getting_started/particles.ipynb` with imports of
  the actual functions and classes that were used. (:pr:`2311`)
- Applied minor refactorings and formatting improvements to
  :file:`docs/notebooks/dispersion/stix_dispersion.ipynb`. (:pr:`2312`)
- Updated the |coding guide| by discussing when to use aliases and
  applied the ``:py:`` role so that in-line code gets formatted the same
  as Python code blocks. (:pr:`2324`)
- Updated the docstrings and type hint annotations in
  `plasmapy.formulary.lengths`. (:pr:`2356`)
- Refactored :file:`docs/conf.py` to improve organization. (:pr:`2363`)
- Updated the narrative documentation on particle objects to
  include |CustomParticle|, |DimensionlessParticle|, and |ParticleList|
  objects. (:pr:`2377`)


Trivial/Internal Changes
------------------------

- Modernized :file:`MANIFEST.in`. (:pr:`2189`)
- Applied automated refactorings from Sourcery. (:pr:`2219`)
- Distributions defined in the `~plasmapy.formulary.distribution` module
  will now raise a `ValueError` for an improper ``units`` parameter.
  (:pr:`2229`)
- Added "decorators" section to the |coding guide|. (:pr:`2231`)
- Improved the error message issued by
  `~plasmapy.formulary.speeds.Alfven_speed`
  when the argument provided to ``density`` has a physical type of
  number density and ``ion`` is not provided. (:pr:`2262`)
- Exposed `plasmapy.dispersion.analytical` and `plasmapy.dispersion.numerical`
  to the `plasmapy.dispersion` namespace. (:pr:`2271`)
- Expanded the |ruff| settings to include more linter rules. (:pr:`2295`)
- Add |ruff| linter rules that check for `print` and :py:`pprint`, as
  the `logging` library is generally preferred for production code.
  (:pr:`2296`)
- Updated and corrected author information in :file:`CITATION.cff`.
  (:pr:`2307`)
- Reduced the number of warnings emitted by `plasmapy.particles` during
  tests by decorating test functions with `pytest.mark.filterwarnings`.
  (:pr:`2314`)
- Fixed a `pytest` deprecation warning that had been issued
  by ``plasmapy.utils._pytest_helpers/pytest_helpers.run_test``
  so that `None` is no longer passed to the `pytest.warns` context
  manager. (:pr:`2314`)
- Changed the default configuration for `pytest` so that if a test is marked
  as expected to fail actually passes, then that test will issue an error
  to indicate that the `pytest.mark.xfail` mark can be removed. (:pr:`2315`)
- Added a weekly linkcheck test that verifies that hyperlinks in the
  documentation are up-to-date. (:pr:`2328`)
- Enabled |validate_quantities| to accept annotations of the form
  :py:`u.Quantity[u.m]`, where we have previously run :py:`import astropy.units
  as u`. (:pr:`2346`)
- Both `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func`
  and `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func_deriv`
  now allow |inf| and |nan| arguments without raising a `ValueError`.
  (:pr:`2361`)
- Modified the
  `~plasmapy.particles.particle_class.CustomParticle.charge_number`
  attribute of |CustomParticle| to return a real number rather than a
  dimensionless |Quantity|. (:pr:`2377`)
- Made minor updates to :file:`plasmapy/__init__.py`, including to the
  top-level package docstring. (:pr:`2378`)
- Improved the consistency and specificity of the names of various
  GitHub Actions. (:pr:`2379`)
- Added a |pre-commit| hook to validate |GitHub Actions|. (:pr:`2380`)


PlasmaPy v2023.5.1 (2023-06-07)
===============================

Trivial/Internal Changes
------------------------

- Loosened the requirement on ``ipykernel`` for compatibility with
  Google Colab. (:pr:`2202`)


PlasmaPy v2023.5.0 (2023-05-31)
===============================

Backwards Incompatible Changes
------------------------------

- The signature of `~plasmapy.formulary.relativity.relativistic_energy`
  has changed. The parameter ``m`` has been replaced with ``particle``,
  which now accepts a broader variety of |particle-like| arguments,
  including but not limited to a |Quantity| representing mass. The
  parameter ``v`` has been replaced with ``V`` for consistency with other
  functionality. (:pr:`1871`)
- Changed the minimum required version of Python from 3.8 to 3.9.
  Accordingly, increased the minimum versions of ``numpy`` to ``1.21.0``,
  ``pandas`` to
  ``1.2.0``, ``h5py`` to ``3.1.0``, ``scipy`` to ``1.6.0``, ``voila`` to
  ``0.3.0``, and ``xarray`` to ``0.17.0``. (:pr:`1885`)
- Made |ParticleList| raise a `TypeError` when provided with a string.
  This change was made to avoid potentially ambiguous situations like
  :py:`ParticleList("He")` which was previously equivalent to
  :py:`ParticleList(["H", "e"])` instead of the possibly expected value of
  :py:`ParticleList(["He"])`. (:pr:`1892`)
- In `~plasmapy.dispersion.analytical.two_fluid_.two_fluid`,
  `~plasmapy.dispersion.numerical.hollweg_.hollweg`, and
  `~plasmapy.dispersion.numerical.kinetic_alfven_.kinetic_alfven`
  in `plasmapy.dispersion`, providing the |charge number| as a keyword
  argument (now ``Z``, formerly ``z_mean``) will no longer override the
  charge number provided in ``ion``. (:pr:`2022`, :pr:`2181`, :pr:`2182`)
- |particle_input| no longer enforces that |parameters| named
  ``ionic_level`` are ions or neutral atoms. For equivalent behavior,
  name the parameter ``ion`` instead. (:pr:`2034`)
- Removed ``plasmapy.utils.pytest_helpers`` from PlasmaPy's public API. It is
  still available as ``plasmapy.utils._pytest_helpers``, but might be removed
  in the future. (:pr:`2114`)
- Removed ``plasmapy.tests.helpers`` from PlasmaPy's public API. It is
  still available as ``plasmapy.tests._helpers``, but might be removed in
  the future. (:pr:`2114`)
- The ``ion_species`` |parameter| to
  `~plasmapy.formulary.radiation.thermal_bremsstrahlung` has been renamed to
  ``ion`` in order to provide a more consistent API to functions that accept
  ions as arguments. (:pr:`2135`)


Deprecations and Removals
-------------------------

- In `plasmapy.dispersion`, the ``z_mean`` parameter
  to `~plasmapy.dispersion.analytical.two_fluid_.two_fluid`,
  `~plasmapy.dispersion.numerical.hollweg_.hollweg`, and
  `~plasmapy.dispersion.numerical.kinetic_alfven_.kinetic_alfven` has
  been deprecated. Provide the |charge number| to ``Z`` instead.
  (:pr:`2022`, :pr:`2181`, :pr:`2182`)
- When a function decorated with |particle_input| is provided with
  ``z_mean`` as a keyword |argument|, it will change ``z_mean`` to ``Z``
  and issue a `~plasmapy.utils.exceptions.PlasmaPyDeprecationWarning` if
  the decorated function accepts ``Z`` as a parameter. This capability
  is intended to temporarily preserve the current behavior of several
  functions in `plasmapy.dispersion` and `plasmapy.formulary` as they get
  decorated with |particle_input| over the next few releases. (:pr:`2027`)
- The ``z_mean`` parameter to `~plasmapy.formulary.speeds.ion_sound_speed`
  and `~plasmapy.formulary.speeds.Alfven_speed` has been deprecated and
  may be removed in a future release. Use ``Z`` instead. (:pr:`2134`, :pr:`2179`)


Features
--------

- Added `~plasmapy.dispersion.numerical.kinetic_alfven_.kinetic_alfven`,
  which numerically solves dispersion relations for kinetic Alfvén waves.
  (:pr:`1665`)
- Added the :file:`stix_dispersion.ipynb` notebook
  which contains Stix cold-plasma dispersion examples. (:pr:`1693`)
- Added the `~plasmapy.formulary.frequencies.Buchsbaum_frequency` function.
  (:pr:`1828`)
- Decorated `~plasmapy.formulary.frequencies.gyrofrequency` with
  |particle_input| so that it can accept a broader variety of
  |particle-like| arguments. (:pr:`1869`)
- After having been decorated with |particle_input|, the
  `~plasmapy.formulary.relativity.relativistic_energy` function now
  accepts a broader variety of |particle-like| objects rather than only
  |Quantity| objects representing mass. (:pr:`1871`)
- After having been decorated with |particle_input|, |RelativisticBody|
  now accepts a broader variety of |particle-like| objects. (:pr:`1871`)
- Enabled |particle_input| to accept values of the |charge number| that
  are real numbers but not integers. This capability can now be used by
  many of the functions in `plasmapy.formulary` and elsewhere that are
  decorated with |particle_input|. (:pr:`1884`)
- Decorated `~plasmapy.particles.atomic.reduced_mass` with
  |particle_input| so that it can now accept a broader variety of
  |particle-like| arguments. (:pr:`1921`)
- Added the `plasmapy.analysis.time_series.excess_statistics` module
  including the
  `~plasmapy.analysis.time_series.excess_statistics.ExcessStatistics`
  class for calculating excess statistics of time series. (:pr:`1984`)
- Added `plasmapy.formulary.collisions.helio.collisional_analysis`.
  (:pr:`1986`)
- Enabled |ParticleList| to accept |Quantity| objects of physical type
  mass or electrical charge. (:pr:`1987`)
- The following functions have been decorated with |particle_input| and
  now accept a broader variety of |particle-like| arguments (see also
  :issue:`341`):

  - `~plasmapy.dispersion.analytical.two_fluid_.two_fluid` (:pr:`2022`)
  - `~plasmapy.formulary.frequencies.plasma_frequency` (:pr:`2026`)
  - `~plasmapy.formulary.speeds.ion_sound_speed` (:pr:`2134`)
  - `~plasmapy.formulary.speeds.kappa_thermal_speed` (:pr:`2136`)
  - `~plasmapy.formulary.speeds.Alfven_speed` (:pr:`2179`)
  - `~plasmapy.dispersion.numerical.hollweg_.hollweg` (:pr:`2181`)
  - `~plasmapy.dispersion.numerical.kinetic_alfven_.kinetic_alfven` (:pr:`2182`)

- Refactored `~plasmapy.formulary.lengths.gyroradius` to reduce cognitive
  complexity and increase readability. (:pr:`2031`)
- Added ``mass_numb`` and ``Z`` as parameters to functions decorated
  with |particle_input| in `plasmapy.formulary.lengths`
  and `plasmapy.formulary.distribution`. (:pr:`2140`)


Bug Fixes
---------

- When attempting to create a |Particle| object representing a proton,
  calls like :py:`Particle("H", Z=1, mass_numb=1)` no longer incorrectly
  issue a |ParticleWarning| for redundant particle information. (:pr:`1992`)
- Updated the docstring of
  `~plasmapy.dispersion.numerical.kinetic_alfven_.kinetic_alfven`. (:pr:`2016`)
- Fixed a slight error in `~plasmapy.formulary.frequencies.plasma_frequency`
  and `~plasmapy.formulary.speeds.Alfven_speed` when the charge number was
  provided via ``z_mean`` (or now ``Z``) and inconsistent with the
  charge number provided to ``particle`` (or zero, if ``particle``
  represented an element or isotope with no charge
  information. Previously, if we represented a proton with
  :py:`particle="H-1"` and :py:`z_mean=1`, then the mass used to
  calculate the plasma frequency would have been the mass of a neutral
  hydrogen atom rather than the mass of a proton. However, using
  :py:`particle="p+"` would have produced the correct mass. This
  behavior has been corrected by decorating this function with
  |particle_input|. See also :issue:`2178` and :pr:`2179`. (:pr:`2026`)
- The ``plasmapy.analysis.nullpoint._vector_space`` function now returns a
  list for its delta values instead of an array. (:pr:`2133`)


Improved Documentation
----------------------

- Enabled `sphinx-codeautolink
  <https://sphinx-codeautolink.readthedocs.io/en/latest/>`_ to make code
  examples clickable and give quick access to API documentation. (:pr:`1410`)
- Added an example notebook on ionization states in the solar wind.
  (:pr:`1513`)
- Moved the location of the changelog pages for past releases from
  :file:`docs/whatsnew/` to :file:`docs/changelog/`, and set up
  appropriate redirects. (:pr:`1639`)
- Removed outdated instructions on installing the development version
  of PlasmaPy contained in :file:`docs/contributing/install_dev.rst`.
  (:pr:`1656`)
- Converted :file:`docs/CONTRIBUTING.rst` to :file:`.github/contributing.md`.
  (:pr:`1656`)
- Added a new page to the |contributor guide| on the
  |code contribution workflow|, replacing content previously contained in
  the |coding guide|. (:pr:`1656`)
- Added a page to the |contributor guide| on |getting ready to contribute|.
  (:pr:`1656`)
- Updated docstrings in `plasmapy.formulary.collisions.frequencies`.
  (:pr:`1793`)
- Updated the docstring for |particle_input|. (:pr:`1883`)
- Updated the introductory paragraphs to the |contributor guide|. (:pr:`2014`)
- Moved PlasmaPy's `vision statement
  <https://doi.org/10.5281/zenodo.7734998>`__
  from the online documentation to a Zenodo record. (:pr:`2017`)
- Restructured the |documentation guide| by putting information on writing
  documentation prior to instructions for building documentation. (:pr:`2038`)
- Restructured the |testing guide| by putting information on writing
  tests prior to instructions for running tests. (:pr:`2041`)
- Updated the introduction on the documentation landing page and the
  citation instructions. (:pr:`2055`)
- Updated the |changelog guide|. (:pr:`2059`)
- Added admonitions for functionality that is under development and for
  which backwards incompatible changes might occur in the future. (:pr:`2112`)
- Updated the code contribution workflow instructions in the |contributor
  guide|
  to reflect that first-time contributors should add themselves to the author
  list in :file:`CITATION.cff` instead of in |docs/about/credits.rst|_. (:pr:`2155`)
- Added functionality to automatically generate the author list included
  in |docs/about/credits.rst|_ directly from :file:`CITATION.cff`. The script
  is located at :file:`docs/cff_to_rst.py`. (:pr:`2156`)


Trivial/Internal Changes
------------------------

- Included Python 3.11 in continuous integration tests. (:pr:`1775`)
- Turned the root-level :file:`requirements.txt` into a lockfile for continuous
  integration purposes. (:pr:`1864`)
- Enabled the particle creation factory in
  ``plasmapy.particles._factory`` used by |particle_input| to create
  |CustomParticle| instances of an element or isotope with a
  |charge number| that is a real number but not an integer. (:pr:`1884`)
- Implemented the new private |CustomParticle| constructor from
  :pr:`1881` into the private particle creation factory used by
  |particle_input|. (:pr:`1884`)
- Dropped ``dlint`` from the tests requirements, as it is no longer
  being maintained. (:pr:`1906`)
- Modified |particle_input| to allow |CustomParticle|\ -like objects with
  a defined charge to be passed through to decorated functions when a
  |parameter| to that function annotated with |ParticleLike| is named
  ``ion``. Previously, only |Particle| objects representing ions or
  neutral atoms were allowed to pass through when the parameter was named
  ``ion``. (:pr:`2034`)
- Updated package metadata in :file:`pyproject.toml`. (:pr:`2075`)
- Set minimum versions for all explicitly listed dependencies. (:pr:`2075`)
- Enabled and applied changes for additional rule sets for ``ruff``, and
  removed corresponding ``flake8`` extensions. (:pr:`2080`)
- Changed from ``indexserver`` to ``PIP_INDEX_URL`` to index nightly `numpy`
  builds (:pr:`2138`)
- Updated the function and docstring of
  `~plasmapy.formulary.collisions.helio.collisional_analysis`. (:pr:`2151`)
- Dropped ``flake8`` and its extensions as linters. Instead, |ruff| is now used as
  the primary linter. (:pr:`2170`)
- Expanded the variety of arguments that could be provided to a function
  decorated by `~plasmapy.utils.decorators.converter.angular_freq_to_hz`,
  and refactored this decorator to use ``wrapt``. (:pr:`2175`)


PlasmaPy v2023.1.0 (2023-01-13)
===============================

Backwards Incompatible Changes
------------------------------

- Moved the charged particle radiography analysis codes into a new module
  `~plasmapy.diagnostics.charged_particle_radiography` containing synthetic
  radiography tools in
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography`
  and detector stack calculation tools in
  `~plasmapy.diagnostics.charged_particle_radiography.detector_stacks`.
  (:pr:`1274`)
- Changed the `~plasmapy.formulary.lengths.gyroradius` function so that it
  takes relativistic effects into account by default. (:pr:`1813`)


Deprecations and Removals
-------------------------

- Changed the `~plasmapy.formulary.lengths.gyroradius` function so it no
  longer accepts deprecated ``T_i``. (:pr:`1824`)
- Removed ``plasmapy.formulary.parameters``, which was deprecated in
  the ``0.7.0`` release.  The functionality in that module had previously
  been migrated to modules that are broken down by physical type, such as:
  `plasmapy.formulary.densities`, `plasmapy.formulary.dimensionless`,
  `plasmapy.formulary.frequencies`, `plasmapy.formulary.lengths`,
  `plasmapy.formulary.misc`, and `plasmapy.formulary.speeds`. (:pr:`1833`)
- Deprecated providing a real number to the ``charge`` parameter of
  |CustomParticle| to represent the |charge number|. Use ``Z`` instead.
  (:pr:`1866`)


Features
--------

- Added the
  `~plasmapy.diagnostics.charged_particle_radiography.detector_stacks.Stack`
  and |Layer| objects to the
  `~plasmapy.diagnostics.charged_particle_radiography` module, which represent
  a stack of detector media layers. The
  `~plasmapy.diagnostics.charged_particle_radiography.detector_stacks.Stack.deposition_curves`
  and
  `~plasmapy.diagnostics.charged_particle_radiography.detector_stacks.Stack.energy_bands`
  methods of
  `~plasmapy.diagnostics.charged_particle_radiography.detector_stacks.Stack`
  calculate the particle energies deposited in each detector layer.
  (:pr:`1274`)
-
  `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.Tracker`
  now supports multiple field grids, provided as an iterable. (:pr:`1799`)
- Added the `plasmapy.analysis.time_series.running_moments` module
  including two functions for calculating running moments of time
  series. (:pr:`1803`)
- Added ``lorentzfactor`` as an optional keyword-only argument
  to `~plasmapy.formulary.lengths.gyroradius`. Also added ``relativistic``
  as an optional keyword-only argument which can be set to `False` for the
  non-relativistic approximation. (:pr:`1813`)
- Modified |Particle| attributes to return |nan| in the appropriate units
  when undefined rather than raising exceptions. (:pr:`1825`)
- Added the `~plasmapy.particles.particle_class.CustomParticle.charge_number`
  attribute to |CustomParticle|. (:pr:`1866`)
- Added ``Z`` as a |keyword-only| |parameter| representing the
  |charge number| to |CustomParticle|. (:pr:`1866`)


Improved Documentation
----------------------

- Updated docstrings and annotations in `plasmapy.diagnostics.thomson`.
  (:pr:`1756`)
- Updated the discussion on type descriptions and parameter descriptions
  for docstrings in the |documentation guide|. (:pr:`1757`)
- Updated troubleshooting sections of the |documentation guide|. (:pr:`1817`)
- Added a summary section to the |testing guide|. (:pr:`1823`)
- Updated the |changelog guide|. (:pr:`1826`)
- Reorganized the |coding guide|. (:pr:`1856`)
- Added a documentation page on performance tips. (:pr:`1887`)


Trivial/Internal Changes
------------------------

- Updated warning messages in
  `~plasmapy.formulary.collisions.coulomb.Coulomb_logarithm`. (:pr:`1586`)
- Transferred most of the contents of :file:`setup.py` and :file:`setup.cfg` to
  :file:`pyproject.toml` (see :pep:`518` and :pep:`621`). Simplified ``extras``
  requirements (``pip install plasmapy[all]`` and ``[extras]`` are gone).
  (:pr:`1758`)
- Added `blacken-docs <https://github.com/adamchainz/blacken-docs>`__ to
  the |pre-commit| configuration. (:pr:`1807`)
- Removed ``pytest-xdist`` from the testing requirements (see also
  :issue:`750`). (:pr:`1822`)
- Refactored tests of `~plasmapy.formulary.relativity.Lorentz_factor`
  and `~plasmapy.formulary.relativity.relativistic_energy`. (:pr:`1844`)
- Applied refactorings from ``ruff`` and ``refurb`` to `plasmapy.utils`.
  (:pr:`1845`)
- Applied changes from ``refurb`` to `plasmapy.particles`. (:pr:`1846`)
- Applied changes from ``refurb`` to `plasmapy.formulary`. (:pr:`1847`)
- Apply changes from ``ruff`` and ``refurb`` to `plasmapy.analysis`,
  `plasmapy.diagnostics`, `plasmapy.dispersion`, and `plasmapy.plasma`.
  (:pr:`1853`)
- Added `ruff <https://github.com/charliermarsh/ruff>`__ to the ``pre-commit``
  configuration. (:pr:`1854`)
- Added the ``strict`` and ``allowed_physical_types`` parameters to
  ``plasmapy.utils._units_helpers._get_physical_type_dict``. (:pr:`1880`)
- Added a private constructor method to |CustomParticle| with an API
  that is better suited for use in |ParticleList| and the particle
  factory function used by |particle_input|. (:pr:`1881`)
- Dropped the dependency on ``cached-property`` in favor of
  `functools.cached_property`. (:pr:`1886`)


PlasmaPy v0.9.0 (2022-11-11)
============================

Backwards Incompatible Changes
------------------------------

- Removed the ``none_shall_pass`` parameter from
  |particle_input|. Instead, `typing.Optional` should be used to create
  the annotation (e.g., ``Optional[ParticleLike]``). (:pr:`1057`)
- Renamed the first |parameter| of |particle_input| from
  ``wrapped_function`` to ``wrapped``. (:pr:`1057`)
- Refactored the return pattern
  of
  `~plasmapy.analysis.swept_langmuir.floating_potential.find_floating_potential`
  to follow the ``(vf, extras)`` pattern, where ``vf`` is the computed
  floating potential and ``extras`` is an instance of
  `~plasmapy.analysis.swept_langmuir.floating_potential.VFExtras`
  containing extra parameters from the computation. (:pr:`1565`)
- Moved ``plasmapy.particles.particle_collections.ionic_levels`` to
  `plasmapy.particles.atomic.ionic_levels`. (:pr:`1697`)
- Deprecated
  ``plasmapy.formulary.collisions.fundamental_electron_collision_freq``.
  The
  `~plasmapy.formulary.collisions.frequencies.MaxwellianCollisionFrequencies.Maxwellian_avg_ei_collision_freq`
  attribute of
  `~plasmapy.formulary.collisions.frequencies.MaxwellianCollisionFrequencies`
  should be used instead. (:pr:`1703`)
- Deprecated ``plasmapy.formulary.collisions.fundamental_ion_collision_freq``.
  The
  `~plasmapy.formulary.collisions.frequencies.MaxwellianCollisionFrequencies.Maxwellian_avg_ii_collision_freq`
  attribute of
  `~plasmapy.formulary.collisions.frequencies.MaxwellianCollisionFrequencies`
  should be used instead. (:pr:`1703`)
- The parameters ``Z`` and ``mass_numb`` are now |keyword-only| for
  `~plasmapy.particles.symbols.ionic_symbol` and
  `~plasmapy.particles.symbols.particle_symbol`. (:pr:`1718`)
- Moved the ``valid_categories`` attribute of
  `~plasmapy.particles.particle_class.AbstractPhysicalParticle.is_category`
  to `plasmapy.particles.particle_class.valid_categories`. (:pr:`1720`)
- Changed the behavior of |IonicLevel|, |IonizationState|, and
  |IonizationStateCollection| so that an equality comparison with an
  `object` of a different type returns `False` instead of raising a
  `TypeError`. (:pr:`1721`)
- When the argument provided to
  `~plasmapy.formulary.magnetostatics.GeneralWire` is not callable, a
  `TypeError` will be raised instead of a `ValueError`. (:pr:`1782`)
- In `~plasmapy.diagnostics.thomson.spectral_density` and
  `~plasmapy.diagnostics.thomson.spectral_density_model`, a `TypeError` is
  now raised when ``ions`` is an unsupported type. (:pr:`1782`)
- In `~plasmapy.plasma.grids.AbstractGrid`, a `TypeError` is now raised
  instead of a `ValueError` when ``start`` or ``stop`` are not of the
  appropriate type. (:pr:`1783`)


Deprecations and Removals
-------------------------

- The capability of |particle_input| to process arguments annotated with
  |Particle| or ``(Particle, Particle)`` is now deprecated and will be
  removed in a future release. Use |ParticleLike| as an annotation
  instead. (:pr:`1057`)
- The ``integer_charges`` attribute of |IonizationState| has been
  removed after being deprecated in ``v0.7.0``. Use
  `~plasmapy.particles.ionization_state.IonizationState.charge_numbers`
  instead. (:pr:`1633`)
- The ``integer_charge`` attributes of |Particle| and |IonicLevel| have
  been removed after being deprecated in ``v0.7.0``. Use the
  ``charge_number`` attribute instead. (:pr:`1633`)
- The ``plasmapy.particles.atomic.integer_charge`` function has been
  removed after being deprecated in ``v0.7.0``. Use
  `~plasmapy.particles.atomic.charge_number` instead. (:pr:`1633`)
- Deprecated `plasmapy.formulary.collisions.frequencies.collision_frequency`
  in favor of collision frequency classes in `plasmapy.formulary.collisions`.
  See also :pr:`1703`. (:pr:`1676`)


Features
--------

- Expanded the functionality of the |particle_input| decorator to convert
  |particle-like| and |particle-list-like| arguments into |Particle|,
  |CustomParticle|, and |ParticleList| objects. This change is part of an
  ongoing effort to improve compatibility of functions in subpackages like
  `plasmapy.particles` and `plasmapy.formulary` with |CustomParticle| and
  |ParticleList| objects. (:pr:`1057`)
- Added the
  `~plasmapy.analysis.swept_langmuir.ion_saturation_current.find_ion_saturation_current`
  function to the `~plasmapy.analysis.swept_langmuir` module.  The function
  fits the tail of a swept Langmuir probe trace and returns the linear
  fit corresponding to the ion-saturation current. (:pr:`1469`)
- Created `plasmapy.utils.data` to contain functionality for downloading data
  from
  |PlasmaPy's data repository|. This module contains a new prototype function
  `plasmapy.utils.data.downloader.get_file` which downloads a file from the
  repository. (:pr:`1486`)
- Added the |RelativisticBody| class to facilitate calculation of the
  relativistic properties of a body in motion. (:pr:`1540`)
- Added ``inplace`` as an optional argument to
  `~plasmapy.simulation.particle_integrators.boris_push`. (:pr:`1556`)
- Added a function to calculate the dimensionless Lundquist number.
  (:pr:`1642`)
- Created the `plasmapy.formulary.densities` module. (:pr:`1664`)
- Added `~plasmapy.formulary.densities.critical_density` to calculate the
  critical density of a plasma for a given frequency of radiation. (:pr:`1664`)
- Added the ``plasmapy.formulary.collisions.CollisionFrequencies`` class.
  This class can be used to calculate collision frequencies for two
  interacting species in a plasma. Superseded by :pr:`1703`. (:pr:`1676`)
- Reimplemented `~plasmapy.formulary.quantum.chemical_potential`. (:pr:`1678`)
- Allowed `~plasmapy.formulary.relativity.Lorentz_factor` to accept and
  return |nan| values. (:pr:`1681`)
- Added a test for `~plasmapy.formulary.dimensionless.Hall_parameter` in
  :file:`plasmapy/formulary/test/test_dimensionless.py`\ . (:pr:`1689`)
- Replaced usage of `os.path` with the more modern `pathlib`. (:pr:`1690`)
- Replaced ``pkg_resources`` with the more modern `importlib.metadata`.
  (:pr:`1692`)
- Added the `~plasmapy.particles.particle_class.CustomParticle.categories`
  attribute to |CustomParticle|, and added the ``"custom"`` particle
  category. (:pr:`1700`)
- Moved the ``is_category`` method of |Particle| to
  `~plasmapy.particles.particle_class.AbstractPhysicalParticle`. This
  method is now inherited by both |Particle| and |CustomParticle|. (:pr:`1700`)
- Added
  `~plasmapy.formulary.collisions.frequencies.MaxwellianCollisionFrequencies`
  for
  calculating relevant collision frequencies for Maxwellian populations.
  (:pr:`1703`)
- Refactored `~plasmapy.formulary.collisions`.  The file
  :file:`collisions.py` was converted into a subpackage (directory) and
  it's contents was split into appropriately categorized and named
  sub-modules (files). (:pr:`1769`)


Bug Fixes
---------

- Modified tests in the class ``TestSyntheticRadiograph`` to try to fix an
  intermittent failure of ``test_optical_density_histogram``. (:pr:`1685`)


Improved Documentation
----------------------

- Added the Hollweg dispersion notebook. (:pr:`1392`)
- Creates an example notebook for fitting Thomson scattering spectra using the
  `~plasmapy.diagnostics.thomson.spectral_density_model` function. (:pr:`1520`)
- Updated the |release guide| following the ``0.8.1`` release. (:pr:`1615`)
- Added :file:`docs/whatsnew/dev.rst` as a stub file for the changelogs
  between releases. (:pr:`1623`)
- Added customizations for |towncrier| in :file:`pyproject.toml`. (:pr:`1626`)
- Updated the introductory paragraphs of the |coding guide|. (:pr:`1649`)
- Added a section to the |coding guide| on best practices for naming
  variables. (:pr:`1650`)
- Updated the section of the contributor guide on |pre-commit|, and
  moved it to :file:`docs/contributing/install_dev.rst`. (:pr:`1651`)
- Added sections to the |coding guide| on units and particles. (:pr:`1655`)
- Updated the section of the |coding guide| on code style. (:pr:`1657`)
- Added sections to the |coding guide| on :term:`lite-functions` and
  :term:`aliases`. (:pr:`1658`)
- Added sections to the |coding guide| on imports and requirements.
  (:pr:`1659`)
- Added sections on best practices for comments and error messages to the
  |coding guide|. (:pr:`1660`)
- Updated the section of the |documentation guide| with more detail on the
  "Parameters", "Raises", and "Warns" sections of docstrings. (:pr:`1667`)
- Added a guideline to the |coding guide| specifying how |nan| values
  should be treated in functions that accept |array_like| or |Quantity|
  inputs. (:pr:`1673`)
- Added an admonition to the |changelog guide| that describes how to
  change |reStructuredText| links for removed code objects into inline literals in old
  changelog entries. (:pr:`1674`)
- Split the patent clause from the license file (:file:`LICENSE.md`) into
  its own file (:file:`PATENT.md`). (:pr:`1686`)
- Added explanatory text to the "Notes" sections in the docstrings for
  functions within `~plasmapy.formulary.magnetostatics`. (:pr:`1695`)
- Enabled ``:py:`` as a |reStructuredText| role for inline code formatting in the
  documentation. (:pr:`1698`)
- Increased the minimum version of |Sphinx| to 5.0.0. (:pr:`1698`)
- Updated docstrings and annotations for |ParticleList| and its methods.
  (:pr:`1713`)
- Updated docstrings and annotations in `plasmapy.particles`, including by
  marking parameters as |particle-like| or |atom-like|. (:pr:`1718`)
- Added a section to the |documentation guide| on troubleshooting. (:pr:`1752`)


Trivial/Internal Changes
------------------------

- Moved the functionality responsible for converting |particle-like|
  arguments to particle objects from the |particle_input| decorator into a
  separate class that is now used inside |particle_input|. (:pr:`1057`)
- Added `wrapt <https://wrapt.readthedocs.io>`__ as a dependency. (:pr:`1057`)
- The |particle_input| decorator now processes arguments annotated with
  |ParticleLike|. (:pr:`1057`)
- Added ``tomli`` to the ``tests`` category of requirements. (:pr:`1500`)
- Added tests to verify that the requirements given in the :file:`.txt`
  files in the :file:`requirements` directory are consistent with the
  requirements given in :file:`setup.cfg` and :file:`pyproject.toml`.
  (:pr:`1500`)
- Restricted the required version of
  `sphinx-gallery <https://sphinx-gallery.github.io/stable/index.html>`__
  to ``< 0.11.0``, since
  ``sphinx-gallery`` changed their thumbnail containers to flex containers.
  See pull request
  `sphinx-gallery/#906
  <https://github.com/sphinx-gallery/sphinx-gallery/pull/906>`__
  and issue
  `sphinx-gallery/#905
  <https://github.com/sphinx-gallery/sphinx-gallery/issues/905>`__ for more
  detail. (:pr:`1654`)
- Moved the ``plasmapy.formulary.dimensionless.quantum_theta`` function to
  `plasmapy.formulary.quantum.quantum_theta`. This function can still be
  called from the `plasmapy.formulary.dimensionless` module without issue.
  (:pr:`1671`)
- Reimplemented ``plasmapy.formulary.quantum._chemical_potential_interp``.
  (:pr:`1678`)
- Re-enabled value testing for the ``quantum`` keyword argument in
  `~plasmapy.formulary.collisions.dimensionless.coupling_parameter`.
  (:pr:`1678`)
- Increased the minimum version of |NumPy| to 1.20.0. (:pr:`1694`)
- Added the `~plasmapy.utils.decorators.validators.validate_class_attributes`
  decorator to the `~plasmapy.utils.decorators` module.
  This decorator is useful for class methods that require optional parameters
  to be specified during class instantiation. (:pr:`1703`)
- Made minor improvements to
  ``plasmapy.formulary.collisions.CollisionFrequencies``. (:pr:`1705`)
- Changed the |towncrier| requirement to ``>= 19.2.0, < 22.8.0``.
  Superseded by :pr:`1717`\ . (:pr:`1710`)
- Applied automated refactorings from `Sourcery
  <https://sourcery.ai/>`__. (:pr:`1714`)
- Changed the minimum version of |towncrier| to 22.8.0 and the minimum
  version of |sphinx_changelog|_ to 1.2.0. (:pr:`1717`)
- Changed `~plasmapy.formulary.quantum.chemical_potential` to use the
  :wikipedia:`Broyden-Fletcher-Goldfarb-Shanno algorithm` to implicitly
  solve for the ideal chemical potential. (:pr:`1726`)
- Increased the minimum version of |Astropy| to 5.0.1. (:pr:`1727`)
- Simplified the pull request template. (:pr:`1729`)
- Added a GitHub Action to automatically comment on pull requests with a
  code review checklist. (:pr:`1729`)
- The following functions are now decorated by |particle_input|\ :
  `~plasmapy.formulary.dimensionless.Hall_parameter`,
  `~plasmapy.formulary.distribution.kappa_velocity_1D`,
  `~plasmapy.formulary.distribution.kappa_velocity_3D`,
  `~plasmapy.formulary.distribution.Maxwellian_1D`,
  `~plasmapy.formulary.distribution.Maxwellian_velocity_2D`,
  `~plasmapy.formulary.distribution.Maxwellian_velocity_3D`,
  `~plasmapy.formulary.distribution.Maxwellian_speed_1D`,
  `~plasmapy.formulary.distribution.Maxwellian_speed_2D`,
  `~plasmapy.formulary.distribution.Maxwellian_speed_3D`,
  `~plasmapy.formulary.lengths.gyroradius`, and
  `~plasmapy.formulary.quantum.deBroglie_wavelength`. (:pr:`1732`)
- Changed |particle_input| to raise a `~astropy.units.UnitConversionError`
  when the annotated argument has a physical type other than mass or
  electrical charge. (:pr:`1732`)
- Set up issue forms on |PlasmaPy's GitHub repository| to replace
  issue templates. (:pr:`1733`)
- Made ``pytest`` an ``install`` requirement instead of a ``testing``
  requirement. (:pr:`1749`)
- Added a step to validate :file:`CITATION.cff` as part of the ``linters``
  |tox| testing environment. (:pr:`1771`)
- Added ``cffconvert`` to the ``testing`` requirements. (:pr:`1771`)
- Deleted :file:`codemeta.json`, which recorded project metadata using
  the `CodeMeta <https://codemeta.github.io>`__ metadata
  schema. Instead, project metadata is now stored in :file:`CITATION.cff` which
  uses the `Citation File Format`_ and was created in :pr:`1640`. See also
  :pr:`676` and :issue:`794`. (:pr:`1772`)
- Added the ``flake8`` extensions ``flake8-use-pathlib``,
  ``flake8-builtins``, and ``flake8-comments`` to the testing
  requirements. (:pr:`1777`)
- Added ``tryceratops`` as a ``flake8`` extension. (:pr:`1782`)


PlasmaPy 0.8.1 (2022-07-05)
===========================

Backwards Incompatible Changes
------------------------------

- In `~plasmapy.diagnostics.thomson.spectral_density`, the arguments ``Te`` and ``Ti`` have been renamed ``T_e`` and ``T_i``
  and are now required :term:`keyword-only` arguments. (`#974 <https://github.com/plasmapy/plasmapy/pull/974>`__)
- Moved the ``grid_resolution`` attribute from `~plasmapy.plasma.grids.AbstractGrid`
  to `~plasmapy.plasma.grids.CartesianGrid` and `~plasmapy.plasma.grids.NonUniformCartesianGrid`
  separately. This fixes a potential future bug, because this attribute is only valid as written
  when all axes share the same units. (`#1295 <https://github.com/plasmapy/plasmapy/pull/1295>`__)
- Changed the behavior of the ``__repr__`` method of |CustomParticle| to
  display the symbol as well if it was provided. (`#1397 <https://github.com/plasmapy/plasmapy/pull/1397>`__)
- Removed a block of code that printed out special particle properties
  when ``plasmapy.particles.special_particles`` (renamed to
  ``plasmapy.particles._special_particles``) was executed. (`#1440 <https://github.com/plasmapy/plasmapy/pull/1440>`__)
- Renamed ``plasmapy.particles.elements`` to ``plasmapy.particles._elements``,
  ``plasmapy.particles.isotopes`` to ``plasmapy.particles._isotopes``,
  ``plasmapy.particles.parsing`` to ``plasmapy.particles._parsing``, and
  ``plasmapy.particles.special_particles`` to
  ``plasmapy.particles._special_particles``. Consequently, these modules
  are no longer part of PlasmaPy's public API. Most of these modules did
  not contain any public objects, except for
  ``plasmapy.particles.special_particles.ParticleZoo`` which was renamed
  to ``plasmapy.particles._special_particles.particle_zoo`` and removed
  from the public API. (`#1440 <https://github.com/plasmapy/plasmapy/pull/1440>`__)
- The parameters ``Z`` and ``mass_numb`` to |Particle| are now
  :term:`keyword-only`. (`#1456 <https://github.com/plasmapy/plasmapy/pull/1456>`__)


Deprecations and Removals
-------------------------

- Officially deprecated `plasmapy.formulary.parameters` and scheduled its
  permanent removal for the ``v0.9.0`` release. (`#1453 <https://github.com/plasmapy/plasmapy/pull/1453>`__)
- Dropped support for Python 3.7 in accordance with the deprecation policy
  laid out in `NumPy Enhancement Proposal 29
  <https://numpy.org/neps/nep-0029-deprecation_policy.html>`__. (`#1465 <https://github.com/plasmapy/plasmapy/pull/1465>`__)
- The ``[all]`` option when using |pip| to install `plasmapy` is now
  deprecated and may be removed in a future release. Packages that were
  previously optional (|h5py|, |lmfit|, |mpmath|, and |Numba|) are now
  installed by default when running ``pip install plasmapy``. To install
  all packages required for code development of PlasmaPy, instead run
  ``pip install plasmapy[developer]``. (`#1482 <https://github.com/plasmapy/plasmapy/pull/1482>`__)
- Removed ``plasmapy.optional_deps``. (`#1482 <https://github.com/plasmapy/plasmapy/pull/1482>`__)


Features
--------

- `~plasmapy.diagnostics.thomson.spectral_density` and `~plasmapy.diagnostics.thomson.spectral_density_model`
  now support `~plasmapy.particles.particle_collections` objects as input to the ``ions`` keywords. (`#974 <https://github.com/plasmapy/plasmapy/pull/974>`__)
- Created a :term:`lite-function` for `~plasmapy.diagnostics.thomson.spectral_density`, `~plasmapy.diagnostics.thomson.spectral_density_lite`. (`#974 <https://github.com/plasmapy/plasmapy/pull/974>`__)
- Added a fitting function for 1D spectra, `~plasmapy.diagnostics.thomson.spectral_density_model`, to the Thomson scattering diagnostic module. (`#974 <https://github.com/plasmapy/plasmapy/pull/974>`__)
- Created function ``plasmapy.formulary.parameters.thermal_speed_coefficients``
  to support ``plasmapy.formulary.parameters.thermal_speed_lite`` usage by
  calculating the various thermal speed coefficients.
  ``plasmapy.formulary.parameters.thermal_speed_coefficients`` is also bound
  to ``plasmapy.formulary.parameters.thermal_speed`` as the ``coefficients``
  attribute. (`#1145 <https://github.com/plasmapy/plasmapy/pull/1145>`__)
- Created decorator `~plasmapy.utils.decorators.lite_func.bind_lite_func`
  for handling the binding of :term:`lite-functions` and any supporting
  functions to a parent function. (`#1145 <https://github.com/plasmapy/plasmapy/pull/1145>`__)
- Introduced the concept of :term:`lite-functions`, by creating the lite-function
  ``plasmapy.formulary.parameters.thermal_speed_lite``, which is a simplified
  and |Numba| jitted version of ``plasmapy.formulary.parameters.thermal_speed``.
  These functions are intended for computational use and as such have no
  validation of input or output values.
  ``plasmapy.formulary.parameters.thermal_speed_lite`` is also bound to
  ``plasmapy.formulary.parameters.thermal_speed`` as the ``lite`` attribute. (`#1145 <https://github.com/plasmapy/plasmapy/pull/1145>`__)
- Added the :file:`hollweg_.py` module to the `~plasmapy.dispersion.numerical`
  subpackage to numerically solve the dispersion relation using Hollweg's method
  :cite:p:`hollweg:1999,bellan:2012`. (`#1189 <https://github.com/plasmapy/plasmapy/pull/1189>`__)
- Implemented non-breaking speed improvements on the methods
  `~plasmapy.plasma.grids.CartesianGrid.nearest_neighbor_interpolator`
  and `~plasmapy.plasma.grids.CartesianGrid.volume_averaged_interpolator`
  for `~plasmapy.plasma.grids.CartesianGrid`. The new interpolators now
  require that the grid axes be sorted, which is always true for uniform
  grids. Added a new test to ensure this stays true. (`#1295 <https://github.com/plasmapy/plasmapy/pull/1295>`__)
- Refactored the interpolator methods on objects defined in `~plasmapy.plasma.grids`.
  All interpolators are now defined in the subclasses of `~plasmapy.plasma.grids.AbstractGrid`.
  Calling the interpolator methods on `~plasmapy.plasma.grids.AbstractGrid`
  raises a `NotImplementedError` exception. (`#1295 <https://github.com/plasmapy/plasmapy/pull/1295>`__)
- Created :term:`lite-function` ``plasmapy.formulary.parameters.plasma_frequency_lite``. (`#1308 <https://github.com/plasmapy/plasmapy/pull/1308>`__)
- Added the `~plasmapy.particles.particle_class.molecule` function to build
  |CustomParticle| objects from a `str` representing a molecule symbol. (`#1309 <https://github.com/plasmapy/plasmapy/pull/1309>`__)
- Added the `~plasmapy.particles.particle_collections.ParticleList.is_category`
  method for |ParticleList| objects.  This method is analogous to the
  `~plasmapy.particles.particle_class.Particle.is_category` method for
  |Particle| objects. (`#1378 <https://github.com/plasmapy/plasmapy/pull/1378>`__)
- Created the prototype analysis tool `plasmapy.analysis.nullpoint` for finding the
  null points in a vector space using the trilinear interpolation method of
  :cite:t:`haynes:2007`. (`#1383 <https://github.com/plasmapy/plasmapy/pull/1383>`__)
- Created `plasmapy.formulary.lengths` to contain length related plasma
  parameters, and migrated `~plasmapy.formulary.lengths.Debye_length`,
  `~plasmapy.formulary.lengths.gyroradius`, and
  `~plasmapy.formulary.lengths.inertial_length` from
  `plasmapy.formulary.parameters` to the new module.  Related aliases were
  also migrated. (`#1434 <https://github.com/plasmapy/plasmapy/pull/1434>`__)
- Created `plasmapy.formulary.frequencies` to contain frequency related
  plasma parameters, and migrated
  `~plasmapy.formulary.frequencies.gyrofrequency`,
  `~plasmapy.formulary.frequencies.plasma_frequency`,
  `~plasmapy.formulary.frequencies.plasma_frequency_lite`,
  `~plasmapy.formulary.frequencies.lower_hybrid_frequency`, and
  `~plasmapy.formulary.frequencies.upper_hybrid_frequency` from
  `plasmapy.formulary.parameters` to the new module.  Related aliases were
  also migrated. (`#1439 <https://github.com/plasmapy/plasmapy/pull/1439>`__)
- Migrated
  `~plasmapy.formulary.dimensionless.Debye_number`, and
  `~plasmapy.formulary.dimensionless.Hall_parameter` from
  `plasmapy.formulary.parameters` to `plasmapy.formulary.dimensionless`.
  Related aliases were also migrated. (`#1444 <https://github.com/plasmapy/plasmapy/pull/1444>`__)
- Created `plasmapy.formulary.speeds` to contain frequency related
  plasma parameters, and migrated
  `~plasmapy.formulary.speeds.Alfven_speed`,
  `~plasmapy.formulary.speeds.ion_sound_speed`,
  `~plasmapy.formulary.speeds.kappa_thermal_speed`,
  `~plasmapy.formulary.speeds.thermal_speed`,
  `~plasmapy.formulary.speeds.thermal_speed_coefficients`, and
  `~plasmapy.formulary.speeds.thermal_speed_lite` from
  `plasmapy.formulary.parameters` to the new module.  Related aliases were
  also migrated. (`#1448 <https://github.com/plasmapy/plasmapy/pull/1448>`__)
- Created `plasmapy.formulary.misc` to contain functionality for
  miscellaneous plasma parameters, and migrated
  ``~plasmapy.formulary.misc._grab_charge``,
  `~plasmapy.formulary.misc.Bohm_diffusion`,
  `~plasmapy.formulary.misc.magnetic_energy_density`,
  `~plasmapy.formulary.misc.magnetic_pressure`,
  `~plasmapy.formulary.misc.mass_density`, and
  `~plasmapy.formulary.misc.thermal_pressure` from
  `plasmapy.formulary.parameters` to the new module.  Related aliases were
  also migrated. (`#1453 <https://github.com/plasmapy/plasmapy/pull/1453>`__)
- Created :term:`lite-functions`
  `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func_lite` and
  `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func_deriv_lite`
  for `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func`
  and `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func_deriv`
  respectively. (`#1473 <https://github.com/plasmapy/plasmapy/pull/1473>`__)
- Created :term:`lite-function`
  `plasmapy.formulary.dielectric.permittivity_1D_Maxwellian_lite` for
  `plasmapy.formulary.dielectric.permittivity_1D_Maxwellian`. (`#1476 <https://github.com/plasmapy/plasmapy/pull/1476>`__)
- Added the :file:`stix_.py` module to the `~plasmapy.dispersion.analytical`
  subpackage which contains the Stix cold-plasma dispersion solution
  :func:`~plasmapy.dispersion.analytical.stix_.stix`,
  :cite:p:`stix:1992,bellan:2012`. (`#1511 <https://github.com/plasmapy/plasmapy/pull/1511>`__)
- ``Particle("Li").ionize()`` no longer results in a `~plasmapy.particles.exceptions.ChargeError`. Instead, ionization of a neutral atom is assumed. (`#1514 <https://github.com/plasmapy/plasmapy/pull/1514>`__)
- Created the |ParticleListLike| typing construct and added
  :term:`particle-list-like` to the |glossary|. (`#1528 <https://github.com/plasmapy/plasmapy/pull/1528>`__)
- Added a null point classifier function which determines the
  type of a given 3D magnetic null point. (`#1554 <https://github.com/plasmapy/plasmapy/pull/1554>`__)
- Added support for arbitrarily shaped input arrays to the function `plasmapy.formulary.collisions.lengths.impact_parameter`. (`#1604 <https://github.com/plasmapy/plasmapy/pull/1604>`__)


Bug Fixes
---------

- Fixed a bug in the ``_make_grid`` method of `~plasmapy.plasma.grids.AbstractGrid`
  that would fail to smoothly handle invalid user input if the ``start``,
  ``stop``, or ``num`` keywords were not the correct type. (`#1295 <https://github.com/plasmapy/plasmapy/pull/1295>`__)
- Fixed a bug with |Particle| where ``Particle("p+") == Particle("H", Z=1,
  mass_numb=1)`` led to a |ParticleError|. (`#1366 <https://github.com/plasmapy/plasmapy/pull/1366>`__)
- For ``plasmapy.formulary.parameters.gyroradius``, updated the default
  keyword arguments and conditional for issuing the
  `~plasmapy.utils.exceptions.PlasmaPyFutureWarning`.  This addresses the
  incorrect behavior where a `ValueError` is raised if an array is passed
  to the deprecated keyword ``T_i``. (`#1430 <https://github.com/plasmapy/plasmapy/pull/1430>`__)
- Exposed `plasmapy.formulary.misc` to the `plasmapy.formulary` namespace. (`#1471 <https://github.com/plasmapy/plasmapy/pull/1471>`__)
- Replaced misuse of ``max_exp_bias - max_exp_bias`` with ``max_exp_bias - min_exp_bias``
  when creating seed parameters for the bimaxwellian fit function inside
  :func:`~plasmapy.diagnostics.langmuir.get_electron_temperature`. (`#1487 <https://github.com/plasmapy/plasmapy/pull/1487>`__)
- Corrected the improper inversion of the electron temperature for the
  non-bimaxwellian case for
  :func:`~plasmapy.diagnostics.langmuir.get_electron_temperature`.
  The electron temperature, and not the slope, is a fit parameter of the
  curve used by
  :func:`~plasmapy.diagnostics.langmuir.get_electron_temperature`,
  so there is no need for the inversion.  The returned value is now the
  electron temperature and not its reciprocal. (`#1487 <https://github.com/plasmapy/plasmapy/pull/1487>`__)
- Exposed the `~plasmapy.analysis` and `~plasmapy.dispersion` subpackages
  to the `plasmapy` namespace. (`#1512 <https://github.com/plasmapy/plasmapy/pull/1512>`__)
- Changed the :meth:`~plasmapy.analysis.fit_functions.Linear.curve_fit`
  method on `plasmapy.analysis.fit_functions.Linear` so that the
  arbitrary keyword arguments get passed to `scipy.stats.linregress`.
  Previously, :meth:`~plasmapy.analysis.fit_functions.Linear.curve_fit`
  had accepted arbitrary keyword arguments but did not pass them along to
  `~scipy.stats.linregress`. (`#1518 <https://github.com/plasmapy/plasmapy/pull/1518>`__)
- Fixed a bug in :func:`~plasmapy.dispersion.numerical.hollweg_.hollweg`
  that did not allow for arguments ``theta`` and ``k`` to simultaneously
  be arrays. (`#1529 <https://github.com/plasmapy/plasmapy/pull/1529>`__)
- Fixed the ``Z`` dependence in
  `~plasmapy.formulary.collisions.frequencies.fundamental_electron_collision_freq`,
  by replacing ``n_e`` with ``n_i`` while calling
  `~plasmapy.formulary.collisions.frequencies.collision_frequency`. (`#1546 <https://github.com/plasmapy/plasmapy/pull/1546>`__)
- Updated the regular expression matching used by
  `~plasmapy.particles.particle_class.Particle` to parse and identify a
  |particle-like| string.  This fixes the bug where a string with
  a trailing space (e.g. ``"Ar "``) was converted into a negatively charged
  ion (e.g. ``"Ar -1"``). (`#1555 <https://github.com/plasmapy/plasmapy/pull/1555>`__)
- Exposed `plasmapy.formulary.radiation` and functions therein to the
  `plasmapy.formulary` namespace. (`#1572 <https://github.com/plasmapy/plasmapy/pull/1572>`__)


Improved Documentation
----------------------

- Added a :term:`lite-function` group to the configuration value
  :confval:`automodapi_custom_groups` that introduces the
  ``__lite_funcs__`` dunder for listing the lite-functions in a module
  (akin to the ``__all__`` dunder). (`#1145 <https://github.com/plasmapy/plasmapy/pull/1145>`__)
- Added a page in the |contributor guide| that describes how to add
  changelog entries. (`#1198 <https://github.com/plasmapy/plasmapy/pull/1198>`__)
- Created an example notebook that lets users input plasma properties and get plasma parameters. (`#1229 <https://github.com/plasmapy/plasmapy/pull/1229>`__)
- The file
  :file:`docs/_static/css/admonition_color_contrast.css` was added to
  include color customizations for |Sphinx| admonitions that originally
  came from
  `sphinx_rtd_theme_ext_color_contrast
  <https://github.com/AaltoSciComp/sphinx_rtd_theme_ext_color_contrast>`_. (`#1287 <https://github.com/plasmapy/plasmapy/pull/1287>`__)
- Changed the color contrast of links and admonitions to be consistent
  with the `Web Content Accessibility Guidelines 2 Level AA Conformance
  <https://www.w3.org/TR/2021/WD-WCAG22-20210521/#contrast-minimum>`__
  for contrast. (`#1287 <https://github.com/plasmapy/plasmapy/pull/1287>`__)
- Re-organized CSS_ files for the online documentation. The file
  :file:`docs/_static/rtd_theme_overrides.css` was re-organized,
  renamed to :file:`docs/_static/css/plasmapy.css`, and updated with
  comments to help someone unfamiliar with CSS_ to understand the file and
  syntax. (`#1287 <https://github.com/plasmapy/plasmapy/pull/1287>`__)
- Put references from `plasmapy.formulary` into :file:`docs/bibliography.bib`
  in BibTeX format. (`#1299 <https://github.com/plasmapy/plasmapy/pull/1299>`__)
- Added a discussion of test parametrization with argument unpacking to
  the |testing guide| in the |contributor guide|. (`#1316 <https://github.com/plasmapy/plasmapy/pull/1316>`__)
- Adopted the `Contributor Covenant Code of Conduct version 2.1
  <https://www.contributor-covenant.org/version/2/1/code_of_conduct/>`__
  and updated the
  :ref:`Contributor Covenant Code of Conduct <plasmapy-code-of-conduct>`
  page accordingly. (`#1324 <https://github.com/plasmapy/plasmapy/pull/1324>`__)
- Updated deprecated meeting and calendar links in :file:`README.md`. (`#1327 <https://github.com/plasmapy/plasmapy/pull/1327>`__)
- Enabled the `sphinx-hoverxref <https://sphinx-hoverxref.readthedocs.io>`_
  extension to |Sphinx|. (`#1353 <https://github.com/plasmapy/plasmapy/pull/1353>`__)
- Added bullet points on module level docstrings and ``__all__`` to the
  documentation guide. (`#1359 <https://github.com/plasmapy/plasmapy/pull/1359>`__)
- Reverted the code syntax highlighting style back to the ``pygments``
  default. The minimum version of ``pygments`` was set to ``2.11.0`` because
  the default style was changed to meet accessibility guidelines for
  contrast in this release. (`#1361 <https://github.com/plasmapy/plasmapy/pull/1361>`__)
- Described additional environments for building the documentation with make_
  in the |documentation guide|. (`#1373 <https://github.com/plasmapy/plasmapy/pull/1373>`__)
- Moved references from individual docstrings to the |bibliography|. (`#1374 <https://github.com/plasmapy/plasmapy/pull/1374>`__)
- Fixed the docstring of `~plasmapy.formulary.collisions.dimensionless.coupling_parameter`. (`#1379 <https://github.com/plasmapy/plasmapy/pull/1379>`__)
- Added an example notebook that introduces how to use `astropy.units`. (`#1380 <https://github.com/plasmapy/plasmapy/pull/1380>`__)
- Added a "Getting Started" page to the documentation sidebar and a "Getting
  Started" section to the examples gallery. (`#1380 <https://github.com/plasmapy/plasmapy/pull/1380>`__)
- Added an example notebook that introduces how to use `plasmapy.particles`. (`#1382 <https://github.com/plasmapy/plasmapy/pull/1382>`__)
- Described the |plasma-calculator| in the narrative documentation. (`#1390 <https://github.com/plasmapy/plasmapy/pull/1390>`__)
- Updated the cold magnetized plasma dielectric permittivity tensor
  notebook. (`#1396 <https://github.com/plasmapy/plasmapy/pull/1396>`__)
- Configured the |Sphinx| extension `sphinx-hoverxref`. (`#1437 <https://github.com/plasmapy/plasmapy/pull/1437>`__)
- Removed the following files from :file:`docs/api_static`\ :
  ``plasmapy.particles.elements.rst``,
  ``plasmapy.particles.isotopes.rst``,
  ``plasmapy.particles.parsing.rst``, and
  ``plasmapy.particles.special_particles.rst``. These files corresponded
  to modules that were renamed with a leading underscore to indicate that
  they are no longer part of the public API. (`#1440 <https://github.com/plasmapy/plasmapy/pull/1440>`__)
- Updated the docstring for `plasmapy.particles.particle_class.molecule`. (`#1455 <https://github.com/plasmapy/plasmapy/pull/1455>`__)
- Hid the documentation page that contained the subpackage stability
  matrix. (`#1466 <https://github.com/plasmapy/plasmapy/pull/1466>`__)
- Added a discussion of doctests to the |documentation guide|. (`#1478 <https://github.com/plasmapy/plasmapy/pull/1478>`__)
- Removed the section on package requirements from the instructions on how
  to install `plasmapy`. (`#1482 <https://github.com/plasmapy/plasmapy/pull/1482>`__)
- Updated the instructions on how to install `plasmapy`. (`#1482 <https://github.com/plasmapy/plasmapy/pull/1482>`__)
- Defined ``autodoc_typehints_format="short"`` so signature type hints
  are displayed in short form, i.e. without the leading module names. (`#1488 <https://github.com/plasmapy/plasmapy/pull/1488>`__)
- Set minimum version of `sphinx` to ``v4.4``. (`#1488 <https://github.com/plasmapy/plasmapy/pull/1488>`__)
- Defined the :confval:`nitpick_ignore_regex` configuration variable in
  :file:`docs/conf.py` to specify regular expressions for objects to
  ignore in nitpicky documentation builds. (`#1509 <https://github.com/plasmapy/plasmapy/pull/1509>`__)
- Made numerous minor updates and fixes to |reStructuredText| links in docstrings and
  the narrative documentation. (`#1509 <https://github.com/plasmapy/plasmapy/pull/1509>`__)
- Described the GitHub Action for `codespell <https://github.com/codespell-project/codespell>`__
  in the |testing guide|. (`#1530 <https://github.com/plasmapy/plasmapy/pull/1530>`__)
- Added the |sphinx-issues|_ extension to |Sphinx| to simplify linking to
  GitHub issues, pull requests, users, and commits. (`#1532 <https://github.com/plasmapy/plasmapy/pull/1532>`__)
- Added the `sphinx.ext.extlinks` extension to |Sphinx| to simplify adding
  links to external domains which have a common base URL. (`#1532 <https://github.com/plasmapy/plasmapy/pull/1532>`__)
- Added the |sphinx-notfound-page|_ extension to |Sphinx| so that the
  documentation now has a :wikipedia:`404 <HTTP_404>` page in the same
  style as the rest of the documentation. (`#1532 <https://github.com/plasmapy/plasmapy/pull/1532>`__)
- Added a notebook on using `~plasmapy.formulary.dimensionless.beta`
  from the `plasmapy.formulary` module to calculate plasma β in
  different parts of the solar atmosphere. (`#1552 <https://github.com/plasmapy/plasmapy/pull/1552>`__)
- Added an example notebook for the null point finder module. (`#1554 <https://github.com/plasmapy/plasmapy/pull/1554>`__)
- Added an example notebook that calculates plasma parameters associated
  with the Magnetospheric Multiscale Mission (MMS). (`#1568 <https://github.com/plasmapy/plasmapy/pull/1568>`__)
- Added an example notebook that discusses Coulomb collisions. (`#1569 <https://github.com/plasmapy/plasmapy/pull/1569>`__)
- Increased the strictness of the ``build_docs`` |tox| environment so that
  broken |reStructuredText| links now emit warnings which are then treated as errors,
  fixed the new errors, removed the ``build_docs_nitpicky`` |tox|
  environment, and updated the |documentation guide| accordingly. (`#1587 <https://github.com/plasmapy/plasmapy/pull/1587>`__)
- Renamed the :file:`magnetic_statics.ipynb` notebook to
  :file:`magnetostatics.ipynb`, and made some minor edits to its text
  and plotting code. (`#1588 <https://github.com/plasmapy/plasmapy/pull/1588>`__)
- Added examples sections to the documentation pages for several modules
  within `plasmapy.formulary`. (`#1590 <https://github.com/plasmapy/plasmapy/pull/1590>`__)
- Re-organized the directory structure for example notebooks. (`#1590 <https://github.com/plasmapy/plasmapy/pull/1590>`__)
- Alphabetized the author list in :file:`docs/about/credits.rst`, and
  added missing authors from using ``git log`` and the pull request
  history. (`#1599 <https://github.com/plasmapy/plasmapy/pull/1599>`__)
- Renamed :file:`docs/development` → :file:`docs/contributing`, and set up
  redirects from the original hyperlinks to the new ones for the
  contributor guide. (`#1605 <https://github.com/plasmapy/plasmapy/pull/1605>`__)
- Added |sphinx-reredirects|_ as a |Sphinx| extension to allow website redirects. (`#1605 <https://github.com/plasmapy/plasmapy/pull/1605>`__)
- Added a :file:`robots.txt` file to the online documentation to tell web
  crawlers to ignore all but ``stable`` and ``latest`` documentation
  builds when indexing for search engines. (`#1607 <https://github.com/plasmapy/plasmapy/pull/1607>`__)


Trivial/Internal Changes
------------------------

- Streamlined `~plasmapy.utils.decorators.helpers.preserve_signature` such that it only
  binds ``__signature__`` to the wrapped function, i.e. it no longer touches
  any other attribute of the wrapped function. (`#1145 <https://github.com/plasmapy/plasmapy/pull/1145>`__)
- Moved all tests associated with calculating the thermal speed from test
  file :file:`plasmapy/formulary/tests/test_parameters.py` to
  :file:`plasmapy/formulary/tests/test_thermal_speed.py`. (`#1145 <https://github.com/plasmapy/plasmapy/pull/1145>`__)
- Applied |reStructuredText| substitutions for `plasmapy.particles` and
  |ParticleTracker| in the narrative documentation. (`#1158 <https://github.com/plasmapy/plasmapy/pull/1158>`__)
- Added `csslint <https://github.com/CSSLint/csslint>`_ to the
  |pre-commit| configuration to check the formatting and style of CSS_
  files. (`#1287 <https://github.com/plasmapy/plasmapy/pull/1287>`__)
- Added Python 3.10 to the `GitHub Actions`_ test suite. (`#1292 <https://github.com/plasmapy/plasmapy/pull/1292>`__)
- Parametrized tests for ``plasmapy.formulary.parameters.ion_sound_speed``. (`#1313 <https://github.com/plasmapy/plasmapy/pull/1313>`__)
- Added cron tests of the development versions of |matplotlib| and |SciPy|,
  while changing the cadence of cron tests to be run approximately
  fortnightly. (`#1333 <https://github.com/plasmapy/plasmapy/pull/1333>`__)
- Applied `pytest.warns` in several tests to catch warnings that are being
  issued during execution of the test suite. (`#1345 <https://github.com/plasmapy/plasmapy/pull/1345>`__)
- Split the tests running on pull requests into multiple stages. The
  various `pytest` test environments, including code coverage, now run
  conditionally given successful execution of a basic test environment and
  the linter checks. This change also prevents code coverage prompts from
  appearing twice, with incomplete information on the first time. (`#1350 <https://github.com/plasmapy/plasmapy/pull/1350>`__)
- Added a helper function that takes an iterable and creates a `dict` with
  physical types as keys and the corresponding objects from that iterable
  as values. This change updates the minimum required version of |Astropy|
  to 4.3.1. (`#1360 <https://github.com/plasmapy/plasmapy/pull/1360>`__)
- Added the module ``plasmapy.particles._factory`` which contains a
  private function that accepts arguments that can be provided to
  |Particle|, |CustomParticle|, or |ParticleList| and returns the
  appropriate instance of one of those three classes. (`#1365 <https://github.com/plasmapy/plasmapy/pull/1365>`__)
- Used the extract method refactoring pattern on the initialization of
  |Particle| objects. (`#1366 <https://github.com/plasmapy/plasmapy/pull/1366>`__, `#1368 <https://github.com/plasmapy/plasmapy/pull/1368>`__)
- Refactored tests in `plasmapy.particles`. (`#1369 <https://github.com/plasmapy/plasmapy/pull/1369>`__)
- |CustomParticle| and |DimensionlessParticle| no longer emit a warning
  when the charge and/or mass is not provided and got assigned a value of
  |nan| in the appropriate units. (`#1399 <https://github.com/plasmapy/plasmapy/pull/1399>`__)
- Added unit test cases for manual entry of vector values in order to improve code coverage
  in the null point finder. (`#1427 <https://github.com/plasmapy/plasmapy/pull/1427>`__)
- Consolidated and parametrized tests associated with
  ``plasmapy.formulary.parameters.gyroradius``. (`#1430 <https://github.com/plasmapy/plasmapy/pull/1430>`__)
- Within `plasmapy.particles` modules, the ``_elements``, ``_isotopes``,
  ``_parsing``, and ``_special_particles`` modules are now imported
  directly. Before this, objects within these modules were typically
  imported. (`#1440 <https://github.com/plasmapy/plasmapy/pull/1440>`__)
- Renamed objects within the source code for `plasmapy.particles` to
  conform with :pep:`8` naming conventions (e.g., ``ParticleZooClass``
  → ``ParticleZoo``, ``ParticleZoo`` → ``particle_zoo``, and ``Particles``
  → ``particles``). (`#1440 <https://github.com/plasmapy/plasmapy/pull/1440>`__)
- Applied automated refactorings from `Sourcery <https://sourcery.ai/>`__
  to `plasmapy.utils`. (`#1463 <https://github.com/plasmapy/plasmapy/pull/1463>`__)
- Applied automated refactorings from
  `Sourcery <https://sourcery.ai/>`__ to `plasmapy.plasma`. (`#1464 <https://github.com/plasmapy/plasmapy/pull/1464>`__)
- Bumped the minimum version of `h5py` to ``3.0.0``. (`#1465 <https://github.com/plasmapy/plasmapy/pull/1465>`__)
- Changed the raised exception to `ImportError` (from a general `Exception`)
  when attempting to import `plasmapy` from a Python version below the
  minimum supported version. (`#1465 <https://github.com/plasmapy/plasmapy/pull/1465>`__)
- Added a workflow to label pull requests based on size. (`#1467 <https://github.com/plasmapy/plasmapy/pull/1467>`__, `#1492 <https://github.com/plasmapy/plasmapy/pull/1492>`__)
- Separated ``plasmapy.analysis.nullpoint.null_point_find`` into
  two functions named
  `~plasmapy.analysis.nullpoint.null_point_find` and
  `plasmapy.analysis.nullpoint.uniform_null_point_find`.
  `~plasmapy.analysis.nullpoint.null_point_find` finds
  the null points of a vector space whose values are manually
  entered.  `plasmapy.analysis.nullpoint.uniform_null_point_find`
  finds the null points of a uniform vector space whose values
  are generated by a function provided by the user. (`#1477 <https://github.com/plasmapy/plasmapy/pull/1477>`__)
- Applied automated refactorings from
  `Sourcery <https://sourcery.ai/>`__ to `plasmapy.particles`. (`#1479 <https://github.com/plasmapy/plasmapy/pull/1479>`__)
- Applied automated refactorings from
  `Sourcery <https://sourcery.ai/>`__ to `plasmapy.formulary`. (`#1480 <https://github.com/plasmapy/plasmapy/pull/1480>`__)
- Bumped the minimum versions of |mpmath| to ``1.2.1``, `numpy` to
  ``1.19.0``, `pandas` to ``1.0.0``, `pytest` to ``5.4.0``, `scipy` to
  ``1.5.0``, and |xarray| to ``0.15.0``. (`#1482 <https://github.com/plasmapy/plasmapy/pull/1482>`__)
- Moved |h5py|, |lmfit|, |mpmath|, and |Numba| out of the ``extras``
  requirements category and into the ``install`` requirements category.
  These packages are now installed when running ``pip install plasmapy``. (`#1482 <https://github.com/plasmapy/plasmapy/pull/1482>`__)
- Added ``dlint``, ``flake8``, ``flake8-absolute-import``,
  ``flake8-rst-docstrings``, ``flake8-use-fstring``,
  |pydocstyle|, and ``pygments`` into the ``tests`` requirements category and
  |pre-commit| into the ``extras`` requirements category. These
  dependencies are not required for basic installation with |pip|. (`#1482 <https://github.com/plasmapy/plasmapy/pull/1482>`__)
- Updated :file:`docs/environment.yml` to use |pip| to install all
  requirements specified by :file:`requirements.txt` when creating a
  Conda environment. (`#1482 <https://github.com/plasmapy/plasmapy/pull/1482>`__)
- Used `codespell <https://github.com/codespell-project/codespell>`__
  to fix typos. (`#1493 <https://github.com/plasmapy/plasmapy/pull/1493>`__)
- Used `contextlib.suppress` to suppress exceptions, instead of ``try`` &
  ``except`` blocks. (`#1494 <https://github.com/plasmapy/plasmapy/pull/1494>`__)
- Added a |pre-commit| hook that transforms relative imports to absolute
  imports, except in :file:`docs/plasmapy_sphinx`. (`#1499 <https://github.com/plasmapy/plasmapy/pull/1499>`__)
- Added a test that ``import plasmapy`` does not raise an exception. (`#1501 <https://github.com/plasmapy/plasmapy/pull/1501>`__)
- Added a GitHub Action for `codespell
  <https://github.com/codespell-project/codespell>`__, and updated the
  corresponding |tox| environment to print out contextual information. (`#1530 <https://github.com/plasmapy/plasmapy/pull/1530>`__)
- Added :file:`plasmapy/utils/units_definitions.py` to precompute units
  which were applied to optimize functionality in
  :file:`plasmapy/formulary/distribution.py`. (`#1531 <https://github.com/plasmapy/plasmapy/pull/1531>`__)
- Replaced ``except Exception`` clauses in ``formulary``, ``particles``, and ``utils`` with specific exception statements. (`#1541 <https://github.com/plasmapy/plasmapy/pull/1541>`__)
- Added tests for passing array valued ``k`` and ``theta`` arguments
  to :func:`~plasmapy.dispersion.numerical.hollweg_.hollweg`, which was
  an added feature in :pr:`1529`. (`#1549 <https://github.com/plasmapy/plasmapy/pull/1549>`__)
- Added `flake8-implicit-str-concat
  <https://github.com/flake8-implicit-str-concat/flake8-implicit-str-concat>`__
  and `flake8-mutable <https://github.com/ebeweber/flake8-mutable>`__
  as extensions for ``flake8``. (`#1557 <https://github.com/plasmapy/plasmapy/pull/1557>`__)
- Added `flake8-simplify <https://github.com/MartinThoma/flake8-simplify>`__
  as an extension for ``flake8``. (`#1558 <https://github.com/plasmapy/plasmapy/pull/1558>`__)
- Applied automated refactorings from
  `Sourcery <https://sourcery.ai/>`__ to `plasmapy.dispersion`. (`#1562 <https://github.com/plasmapy/plasmapy/pull/1562>`__)
- Applied automated refactorings from
  `Sourcery <https://sourcery.ai/>`__ to `plasmapy.diagnostics`. (`#1563 <https://github.com/plasmapy/plasmapy/pull/1563>`__)
- Applied automated refactorings from
  `Sourcery <https://sourcery.ai/>`__ to `plasmapy.analysis`. (`#1564 <https://github.com/plasmapy/plasmapy/pull/1564>`__)
- Removed an extraneous `print` statement from
  `~plasmapy.formulary.collisions.frequencies.frequencies.collision_frequency` that
  activated when the colliding particles were both electrons. (`#1570 <https://github.com/plasmapy/plasmapy/pull/1570>`__)
- Changed the type hints for ``z_mean`` in `plasmapy.formulary.collisions`
  functions from ``astropy.units.dimensionless_unscaled`` to
  `~numbers.Real`. Consequently, ``z_mean`` will no longer be processed by
  `~plasmapy.utils.decorators.validators.validate_quantities`. Previously,
  ``z_mean`` issued a warning when a real number was provided instead of a
  dimensionless |Quantity|. (`#1570 <https://github.com/plasmapy/plasmapy/pull/1570>`__)
- Updated the version of black to 22.3.0 in PlasmaPy's |pre-commit|
  configuration. This update included a formatting change where spaces
  around power operators were removed for sufficiently simple operands
  (e.g., ``a ** b`` → ``a**b``). (`#1582 <https://github.com/plasmapy/plasmapy/pull/1582>`__)
- Renamed ``units_definitions`` to ``_units_definitions`` and
  ``units_helpers`` to ``_units_helpers`` in `plasmapy.utils` to mark
  these modules as private. (`#1587 <https://github.com/plasmapy/plasmapy/pull/1587>`__)
- Updated the :file:`codemeta.json` file with metadata for the
  version ``0.8.1`` release. (`#1606 <https://github.com/plasmapy/plasmapy/pull/1606>`__)


PlasmaPy v0.7.0 (2021-11-18)
============================

Backwards Incompatible Changes
------------------------------

- Removed alias ``tfds_`` to
  ``plasmapy.dispersion.two_fluid_dispersion.two_fluid_dispersion_solution``,
  with the reasoning behind the removal outlined in the pull request. (`#1101 <https://github.com/plasmapy/plasmapy/pull/1101>`__)
- Removed the ``Tracker.synthetic_radiograph()`` method and created the
  standalone function
  :func:`~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiograph`
  in its place.  This new function takes either a
  `~plasmapy.diagnostics.charged_particle_radiography.Tracker` object or
  a dictionary equivalent to
  `~plasmapy.diagnostics.charged_particle_radiography.Tracker.results_dict`. (`#1134 <https://github.com/plasmapy/plasmapy/pull/1134>`__)
- Renamed subpackage ``plasmapy.diagnostics.proton_radiography`` to
  `plasmapy.diagnostics.charged_particle_radiography`, and renamed the
  ``SyntheticProtonRadiograph`` class within that module to
  `~plasmapy.diagnostics.charged_particle_radiography.Tracker`. (`#1134 <https://github.com/plasmapy/plasmapy/pull/1134>`__)
- `~plasmapy.diagnostics.charged_particle_radiography.Tracker` no longer
  supports making changes to an instantiated object and
  re-running the simulation.  Subsequent simulations should be performed
  by instantiating a new
  `~plasmapy.diagnostics.charged_particle_radiography.Tracker` object and
  running its simulation. (`#1134 <https://github.com/plasmapy/plasmapy/pull/1134>`__)
- For `~plasmapy.plasma.grids.CartesianGrid` the
  `~plasmapy.plasma.grids.CartesianGrid.volume_averaged_interpolator`
  now returns `numpy.nan` values for any interpolation not bounded by
  the grid points. (`#1173 <https://github.com/plasmapy/plasmapy/pull/1173>`__)
- Renamed file :file:`two_fluid_dispersion.py` to :file:`two_fluid_.py`
  and moved it into the `plasmapy.dispersion.analytical` subpackage.  The
  function ``two_fluid_dispersion_solution()`` contained within that file
  was renamed to `~plasmapy.dispersion.analytical.two_fluid_.two_fluid`. (`#1208 <https://github.com/plasmapy/plasmapy/pull/1208>`__)
- Changed |ParticleList| so that if it is provided with no arguments, then it creates
  an empty |ParticleList|.  This behavior is analogous to how `list` and `tuple` work. (`#1223 <https://github.com/plasmapy/plasmapy/pull/1223>`__)
- Changed the behavior of |Particle| in equality comparisons. Comparing a
  |Particle| with an object that is not |particle-like| will now
  return `False` instead of raising a `TypeError`. (`#1225 <https://github.com/plasmapy/plasmapy/pull/1225>`__)
- Changed the behavior of `~plasmapy.particles.particle_class.CustomParticle`
  so that it returns `False` when compared for equality with another type.
  Previously, a `TypeError` was raised. (`#1315 <https://github.com/plasmapy/plasmapy/pull/1315>`__)


Deprecations and Removals
-------------------------

- In `plasmapy.particles`, use of the term "integer charge" has
  been deprecated in favor of the term "charge number". The
  `~plasmapy.particles.particle_class.Particle.integer_charge` attribute
  of |Particle| has been deprecated in favor of
  `~plasmapy.particles.particle_class.Particle.charge_number`. The
  `~plasmapy.particles.ionization_state.IonicLevel.integer_charge`
  attribute of |IonicLevel| (formerly ``IonicFraction``) has been
  deprecated in favor of
  `~plasmapy.particles.ionization_state.IonicLevel.charge_number`. The
  `~plasmapy.particles.ionization_state.IonizationState.integer_charges`
  attribute of |IonizationState| has been deprecated in favor of
  `~plasmapy.particles.ionization_state.IonizationState.charge_numbers`. (`#1136 <https://github.com/plasmapy/plasmapy/pull/1136>`__)
- The ``particle`` attribute of |Particle|
  has been removed after having been deprecated in 0.6.0. (`#1146 <https://github.com/plasmapy/plasmapy/pull/1146>`__)
- Use more generalized keyword argument ``T`` instead of ``T_i`` in ``plasmapy.formulary.parameters.gyroradius``.
  The ``T_i`` argument has been deprecated and will be removed in a subsequent release. (`#1210 <https://github.com/plasmapy/plasmapy/pull/1210>`__)


Features
--------

- Add the `~plasmapy.particles.ionization_state.IonizationState.average_ion`
  method to |IonizationState|. (`#1028 <https://github.com/plasmapy/plasmapy/pull/1028>`__)
- Added the
  `~plasmapy.particles.ionization_state_collection.IonizationStateCollection.average_ion`
  method to |IonizationStateCollection|. (`#1028 <https://github.com/plasmapy/plasmapy/pull/1028>`__)
- Added the ``plasmapy.formulary.mathematics.Chandrasekhar_G`` function, which is
  helpful in neoclassical transport theory. This change was
  reverted in `#1233 <https://github.com/plasmapy/plasmapy/pull/1233>`__. (`#1084 <https://github.com/plasmapy/plasmapy/pull/1084>`__)
- Enabled slicing of |IonizationState| instances to return a list of
  |IonicLevel| instances. (`#1130 <https://github.com/plasmapy/plasmapy/pull/1130>`__)
- |IonizationState| instances can now be compared to an |IonizationState|
  of a different element without raising an exception. (`#1130 <https://github.com/plasmapy/plasmapy/pull/1130>`__)
- Allowed `len` to be used on |IonizationState| instances. (`#1130 <https://github.com/plasmapy/plasmapy/pull/1130>`__)
- |IonicLevel| and |IonizationState| now accept an additional, optional ion
  temperature argument for each of the ionic levels. (`#1130 <https://github.com/plasmapy/plasmapy/pull/1130>`__)
- Added the
  :meth:`~plasmapy.diagnostics.charged_particle_radiography.Tracker.save_results`
  method to `~plasmapy.diagnostics.charged_particle_radiography.Tracker`
  for saving results to the :file:`.npz` file format (see `numpy.lib.format` for
  details on the file format). (`#1134 <https://github.com/plasmapy/plasmapy/pull/1134>`__)
- Added the `plasmapy.utils.decorators.deprecation` module. The module includes
  `~plasmapy.utils.decorators.deprecation.deprecated`, which is a decorator that
  is based on `astropy.utils.decorators.deprecated`. (`#1136 <https://github.com/plasmapy/plasmapy/pull/1136>`__)
- Created the `~plasmapy.particles.ionization_state.IonizationState.to_list`
  method of |IonizationState| to provide a |ParticleList| instance that
  contains the different ionic levels. (`#1154 <https://github.com/plasmapy/plasmapy/pull/1154>`__)
- The behavior of the function ``plasmapy.formulary.parameters.gyroradius`` has
  been changed. If `numpy.nan` values are provided for ``T_i`` or ``Vperp``,
  then instead of raising a slightly misleading error, `numpy.nan` in the
  appropriate units is returned. (`#1187 <https://github.com/plasmapy/plasmapy/pull/1187>`__)
- Added the `~plasmapy.particles.particle_collections.ParticleList.average_particle`
  method to |ParticleList|. This method returns a particle with the mean mass and
  charge of the |ParticleList|. The ``use_rms_charge`` and ``use_rms_mass`` keyword
  arguments make this method calculate the root mean square charge and mass, respectively.
  The ``abundances`` keyword argument allows the calculation of the mean or root
  mean square to be weighted. (`#1204 <https://github.com/plasmapy/plasmapy/pull/1204>`__)
- Restructured the `plasmapy.dispersion` subpackage by creating the
  `~plasmapy.dispersion.analytical` subpackage to contain functionality
  related to analytical dispersion solutions. (`#1208 <https://github.com/plasmapy/plasmapy/pull/1208>`__)
- Implemented ``__eq__``, ``__ne__`` and ``__hash__`` to allow
  |CustomParticle| instances to be used as `dict` keys. (`#1216 <https://github.com/plasmapy/plasmapy/pull/1216>`__)
- Added the `~plasmapy.particles.particle_collections.ionic_levels` function to create a
  |ParticleList| initialized with different ionic levels of an element or isotope. (`#1223 <https://github.com/plasmapy/plasmapy/pull/1223>`__)


Bug Fixes
---------

- Made |Particle| instances pickleable. (`#1122 <https://github.com/plasmapy/plasmapy/pull/1122>`__)
- Fixed the behavior of ``plasmapy.formulary.mathematics.Chandrasekhar_G``
  at very small and very large argument values. This change was reverted
  in `#1233 <https://github.com/plasmapy/plasmapy/pull/1233>`__. (`#1125 <https://github.com/plasmapy/plasmapy/pull/1125>`__)
- Running `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiograph`
  with the keyword ``optical_density=True`` will now return `numpy.inf`
  where the source profile intensity is zero. Previously, an incorrect value
  was returned since zero entries were replaced with values of ``1`` before
  taking the logarithm. (`#1134 <https://github.com/plasmapy/plasmapy/pull/1134>`__)
- Fixed a bug in the volume-averaged interpolator for
  `~plasmapy.plasma.grids.CartesianGrid`
  (`~plasmapy.plasma.grids.CartesianGrid.volume_averaged_interpolator`).
  The old method miss interpreted where the interpolation point was
  inside the nearest neighbor cell volume. So, if an interpolation point
  was at the lower bounds of the nearest neighbor cell volume, then the
  position was flipped and interpreted as being at the upper bounds of the
  cell volume, and visa-versa. (`#1173 <https://github.com/plasmapy/plasmapy/pull/1173>`__)
- Fixed the normalization of the wavevector in the Thomson spectral
  density function,
  :func:`~plasmapy.diagnostics.thomson.spectral_density`. The previous
  version was not properly normalizing the wavevector to unity. (`#1190 <https://github.com/plasmapy/plasmapy/pull/1190>`__)
- Reverted most of
  `#1084 <https://github.com/plasmapy/plasmapy/pull/1084>`__ and
  `#1125 <https://github.com/plasmapy/plasmapy/pull/1125>`__,
  removing our implementation of the
  Chandrasekhar G function (for now!). This function may get brought
  back at a later date, once we have an implementation we numerically
  trust. (`#1233 <https://github.com/plasmapy/plasmapy/pull/1233>`__)


Improved Documentation
----------------------

- Improved consistency of documentation style and made
  |reStructuredText| fixes in several subpackages. (`#1073 <https://github.com/plasmapy/plasmapy/pull/1073>`__)
- Added a pre-release section to the release guide.
  This section now includes steps for having a feature freeze about a week before the release,
  followed by a code freeze about two days before the release. (`#1081 <https://github.com/plasmapy/plasmapy/pull/1081>`__)
- Created the |Sphinx| extension package `plasmapy_sphinx` and used it to replace
  ``sphinx_automodapi``.  `plasmapy_sphinx` creates directives :rst:dir:`automodapi`
  and :rst:dir:`automodsumm` to replace the same directives defined by
  ``sphinx_automodapi``.  The documentation was updated so the slight syntax differences
  in the newly defined directives will still render the same as before. (`#1105 <https://github.com/plasmapy/plasmapy/pull/1105>`__)
- The term "integer charge" has been replaced in the documentation with
  the term "charge number". (`#1136 <https://github.com/plasmapy/plasmapy/pull/1136>`__)
- Implemented a framework to define and use common `Sphinx substitutions
  <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
  #substitutions>`__ across the narrative documentation and docstrings.
  These substitutions are defined in :file:`docs/common_links.rst`. (`#1147 <https://github.com/plasmapy/plasmapy/pull/1147>`__)
- Began a project glossary at :file:`docs/glossary.rst`. (`#1149 <https://github.com/plasmapy/plasmapy/pull/1149>`__)
- Changed the default branch name to ``main``.  Locations in the code
  and documentation that referred to the default branch of PlasmaPy (and
  certain other packages) were changed to reflect the new name (including,
  for example, in the development guide in the documentation). (`#1150 <https://github.com/plasmapy/plasmapy/pull/1150>`__)
- Updated information on how to write and build documentation in the
  development guide. (`#1156 <https://github.com/plasmapy/plasmapy/pull/1156>`__)
- Updated information on how to write and run tests in the contributor
  guide. (`#1163 <https://github.com/plasmapy/plasmapy/pull/1163>`__)
- Created an outline of a page in the development guide to describe the workflow
  required to contribute to PlasmaPy. (`#1178 <https://github.com/plasmapy/plasmapy/pull/1178>`__)
- Added brief description about the physics of the upper-hybrid resonance
  to the docstring of the function ``plasmapy.formulary.parameters.upper_hybrid_frequency``. (`#1180 <https://github.com/plasmapy/plasmapy/pull/1180>`__)
- Added a brief description about the physics of the lower-hybrid resonance
  to the docstring of the function ``plasmapy.formulary.parameters.lower_hybrid_frequency``. (`#1181 <https://github.com/plasmapy/plasmapy/pull/1181>`__)
- Made the function ``plasmapy.formulary.parameters.gyrofrequency`` more general
  by removing the indications that it might only work for ions. (`#1183 <https://github.com/plasmapy/plasmapy/pull/1183>`__)
- Make `plasmapy.analysis.fit_functions.AbstractFitFunction.FitParamTuple` a
  property to fix the documentation build warning caused by the release
  of |Sphinx| ``v4.1.0``. (`#1199 <https://github.com/plasmapy/plasmapy/pull/1199>`__)
- Included a step in the release guide to update Binder requirements
  so that the release of PlasmaPy on |PyPI| gets installed when opening
  example notebooks from the stable and release branches of the online
  documentation. (`#1205 <https://github.com/plasmapy/plasmapy/pull/1205>`__)
- Updated the documentation guide to include updates to |tox| environments
  for building the documentation. (`#1206 <https://github.com/plasmapy/plasmapy/pull/1206>`__)
- Fixed numerous broken |reStructuredText| links in prior changelogs. (`#1207 <https://github.com/plasmapy/plasmapy/pull/1207>`__)
- Improve the docstring for `plasmapy.online_help`. (`#1213 <https://github.com/plasmapy/plasmapy/pull/1213>`__)
- Renamed "Development Guide" to "Contributor Guide", and temporarily removed
  the incomplete :file:`docs/development/workflow.rst` from the ``toctree`` of the
  Contributor Guide. (`#1217 <https://github.com/plasmapy/plasmapy/pull/1217>`__)
- Fixed a typo in the docstring of ``plasmapy.formulary.parameters.Alfven_speed``. (`#1218 <https://github.com/plasmapy/plasmapy/pull/1218>`__)
- Fixed broken |reStructuredText| links in docstrings for aliases in `plasmapy.formulary`. (`#1238 <https://github.com/plasmapy/plasmapy/pull/1238>`__)
- Fixed multiple broken and redirected links. (`#1257 <https://github.com/plasmapy/plasmapy/pull/1257>`__)
- Updated the documentation guide to include a description on how to
  add and cite references to PlasmaPy's global bibliography BibTeX file,
  :file:`docs/bibliography.bib`. (`#1263 <https://github.com/plasmapy/plasmapy/pull/1263>`__)
- Added sphinxcontrib-bibtex_ as a |Sphinx| extension to enable references
  to be stored in a BibTeX file. (`#1263 <https://github.com/plasmapy/plasmapy/pull/1263>`__)
- Began a documentation-wide bibliography page. (`#1263 <https://github.com/plasmapy/plasmapy/pull/1263>`__)
- Updated documentation guide to describe where formulae should go in
  docstrings and how to use glossary entries. (`#1264 <https://github.com/plasmapy/plasmapy/pull/1264>`__)
- Updated and fixed hyperlinks in the documentation. (`#1267 <https://github.com/plasmapy/plasmapy/pull/1267>`__)
- Adopted the ``"xcode"`` code highlighting style for
  ``pygments`` to increase color contrast and improve web accessibility. (`#1268 <https://github.com/plasmapy/plasmapy/pull/1268>`__)
- Updated the feedback and communication page. (`#1272 <https://github.com/plasmapy/plasmapy/pull/1272>`__)
- Updated the requirements for the documentation build to include no
  restrictions on ``docutils`` and ``sphinx_rtd_theme >= 1.0.0``.
  ``docutils == 0.17`` is not compatible with ``sphinx_rtd_theme < 1.0``
  (see `#1107 <https://github.com/PlasmaPy/PlasmaPy/pull/1107>`__ and
  `#1230 <https://github.com/PlasmaPy/PlasmaPy/issues/1230>`__). (`#1275 <https://github.com/plasmapy/plasmapy/pull/1275>`__)
- Added a screenshot of the link for the |Read the Docs| preview of the
  documentation for a pull request. (`#1298 <https://github.com/plasmapy/plasmapy/pull/1298>`__)
- Incorporated citations in the
  `~plasmapy.dispersion.analytical.two_fluid_.two_fluid` docstring into
  the PlasmaPy bibliography framework. (`#1301 <https://github.com/plasmapy/plasmapy/pull/1301>`__)


Trivial/Internal Changes
------------------------

- Simplified handling of package dependencies.  Removed duplicated
  requirements files and centralized them instead. Developer dependencies
  can now be installed with either ``pip install plasmapy[developer]`` or
  ``pip install -r requirements.txt``. (`#789 <https://github.com/plasmapy/plasmapy/pull/789>`__)
- Reconfigured ``flake8`` settings in CI. (`#1062 <https://github.com/plasmapy/plasmapy/pull/1062>`__)
- Added |pydocstyle| to continuous integration (CI), to hopefully make
  writing prettier docstrings easier. (`#1062 <https://github.com/plasmapy/plasmapy/pull/1062>`__)
- Added ``flake8-rst-docstrings`` to catch |reStructuredText| formatting
  errors in documentation in the linter stage of
  CI. (`#1062 <https://github.com/plasmapy/plasmapy/pull/1062>`__)
- Added `pytest-regressions
  <https://pytest-regressions.readthedocs.io/en/latest/>`__ to testing
  dependencies, to make regression tests a little easier to write. (`#1084 <https://github.com/plasmapy/plasmapy/pull/1084>`__)
- Fixed a minor error in the :math:`\mathbf{E} × \mathbf{B}` drift
  notebook. (`#1088 <https://github.com/plasmapy/plasmapy/pull/1088>`__)
- Upgrade ``nbqa`` to latest available version (0.6.0). (`#1104 <https://github.com/plasmapy/plasmapy/pull/1104>`__)
- Moved our custom |pre-commit| style testing suite to ``pre-commit.ci``,
  taking advantage of the new ``pre-commit.ci autofix`` command that
  allows manually calling for pre-commit to be run by typing
  that command as a comment to a pull request. (`#1106 <https://github.com/plasmapy/plasmapy/pull/1106>`__)
- Added tests using hypothesis_. (`#1125 <https://github.com/plasmapy/plasmapy/pull/1125>`__)
- Added to :file:`setup.cfg` the configuration
  ``flake8.per-file-ignores=plasmapy/formulary/__init__.py:F403`` to
  ignore warnings resulting from imports like ``from xx import *``. (`#1127 <https://github.com/plasmapy/plasmapy/pull/1127>`__)
- Re-enabled several ``flake8`` checks by removing the following codes from
  the ``flake8.extend-ignore`` configuration in :file:`setup.cfg`: ``D100``, ``D102``,
  ``D103``, ``D104``, ``D200``, ``D210``, ``D301``, ``D401``, ``D407``,
  ``D409``, ``D412``, ``E712``, ``E713``, ``F403``, ``F541``, ``RST213``,
  ``RST306``, and ``RST902``. Addressed any failed linter checks from this
  modification. (`#1127 <https://github.com/plasmapy/plasmapy/pull/1127>`__)
- `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiograph`
  now determines the default detector size to be the smallest detector
  plane centered on the origin that includes all particles. (`#1134 <https://github.com/plasmapy/plasmapy/pull/1134>`__)
- Added ion velocity input to the :file:`thomson.ipynb` diagnostics notebook. (`#1171 <https://github.com/plasmapy/plasmapy/pull/1171>`__)
- Added |tox| and removed `pytest` as extra requirements. (`#1195 <https://github.com/plasmapy/plasmapy/pull/1195>`__)
- Updated |tox| test environments for building the documentation. Added the
  ``build_docs_nitpicky`` environment to check for broken |reStructuredText| links. (`#1206 <https://github.com/plasmapy/plasmapy/pull/1206>`__)
- Added the ``--keep-going`` flag to the ``build_docs*`` |tox| environments with
  the ``-W`` option so that test failures will not stop after the first warning
  (that is treated as an error). (`#1206 <https://github.com/plasmapy/plasmapy/pull/1206>`__)
- Make queries to `plasmapy.online_help` for ``"quantity"`` or ``"quantities"`` redirect to the
  help page for `astropy.units` (which was already the case for ``"unit"`` and ``"units"``). (`#1213 <https://github.com/plasmapy/plasmapy/pull/1213>`__)
- Bumped the |Python| version for |Read the Docs| builds from ``3.7`` to ``3.8``. (`#1248 <https://github.com/plasmapy/plasmapy/pull/1248>`__)
- Refactored :file:`plasmapy/dispersion/tests/test_dispersion.py` to use
  hypothesis_ for property based testing. (`#1249 <https://github.com/plasmapy/plasmapy/pull/1249>`__)
- Defined redirects to allow and anchors to avoid checking when using |Sphinx|
  to verify that hyperlinks are correct via ``make linkcheck``. (`#1267 <https://github.com/plasmapy/plasmapy/pull/1267>`__)
- Replaced usage of `eval` inside |IonizationStateCollection| with `getattr`. (`#1280 <https://github.com/plasmapy/plasmapy/pull/1280>`__)
- Added using `dlint <https://github.com/dlint-py/dlint>`__
  to the ``linters`` testing environment in :file:`tox.ini`
  as a static analysis tool to search for security issues. (`#1280 <https://github.com/plasmapy/plasmapy/pull/1280>`__)
- Enabled using
  `flake8-use-fstring <https://github.com/MichaelKim0407/flake8-use-fstring>`__
  in the ``linters`` testing environment in :file:`tox.ini` to enforce
  usage of formatted string literals (f-strings). (`#1281 <https://github.com/plasmapy/plasmapy/pull/1281>`__)
- Switched usage of `str.format` to formatted string literals (f-strings)
  in several files. (`#1281 <https://github.com/plasmapy/plasmapy/pull/1281>`__)
- Added `flake8-absolute-import <https://github.com/bskinn/flake8-absolute-import>`_
  to the ``linters`` |tox| environment. (`#1283 <https://github.com/plasmapy/plasmapy/pull/1283>`__)
- Removed unused imports, and changed several imports from relative to absolute. (`#1283 <https://github.com/plasmapy/plasmapy/pull/1283>`__)
- Added |pre-commit| hooks to auto-format :file:`.ini`,
  :file:`.toml`, and :file:`.yaml` files, and applied changes from
  those hooks to existing files. (`#1284 <https://github.com/plasmapy/plasmapy/pull/1284>`__)
- Changed the validated units for the ``theta`` input argument of
  `~plasmapy.dispersion.analytical.two_fluid_.two_fluid` from degrees to
  radians. (`#1301 <https://github.com/plasmapy/plasmapy/pull/1301>`__)
- Replaced usage of ``distutils.version.StrictVersion`` with
  ``packaging.version.Version`` because ``distutils`` has been deprecated.
  As part of this change, `packaging <https://packaging.pypa.io>`__ has been
  added as a dependency. (`#1306 <https://github.com/plasmapy/plasmapy/pull/1306>`__)
- Increased the minimum version of matplotlib to 3.3.0 and updated
  `plasmapy.diagnostics.langmuir.swept_probe_analysis` to be compatible
  with matplotlib 3.5.0. (`#1334 <https://github.com/plasmapy/plasmapy/pull/1334>`__)


PlasmaPy v0.6.0 (2021-03-14)
============================

Backwards Incompatible Changes
------------------------------

- The ``State`` namedtuple was changed to the `~plasmapy.particles.IonicFraction`
  class. (Note: #1046 subsequently changed that to
  `~plasmapy.particles.IonicLevel`). (`#796 <https://github.com/plasmapy/plasmapy/pull/796>`__)
- Now, when the `~plasmapy.particles.IonizationState` class is provided with an ion,
  the ionic fraction for that ion is set to 100% for the corresponding element or
  isotope. (`#796 <https://github.com/plasmapy/plasmapy/pull/796>`__)
- ``AtomicError`` was renamed to `~plasmapy.particles.exceptions.ParticleError`
  and ``MissingAtomicDataError`` was renamed to
  `~plasmapy.particles.exceptions.MissingParticleDataError`. (`#796 <https://github.com/plasmapy/plasmapy/pull/796>`__)
- In `plasmapy.particles`, the ``IonizationStates`` class was renamed to
  `~plasmapy.particles.IonizationStateCollection`.  Argument ``n`` of
  ``IonizationStates`` was changed to ``n0`` in
  `~plasmapy.particles.IonizationStateCollection`. (`#796 <https://github.com/plasmapy/plasmapy/pull/796>`__)
- Moved and refactored error message formatting functionality from
  ``plasmapy.utils.error_messages`` to `plasmapy.utils.code_repr`. (`#920 <https://github.com/plasmapy/plasmapy/pull/920>`__)
- Renamed the available "methods" for computing the Coulomb logarithm in an attempt
  to make the names more explicit.  This is implemented using the ``method`` keyword
  for functions `~plasmapy.formulary.collisions.coulomb.Coulomb_logarithm` and
  `~plasmapy.formulary.collisions.lengths.impact_parameter`, and then propagated throughout
  the functionality in `plasmapy.formulary.collisions`. (`#962 <https://github.com/plasmapy/plasmapy/pull/962>`__)
- Add dependency ``pandas >= 1.0.0``.  Modify |xarray| dependency to be
  ``xarray >= 0.14.0``. (`#963 <https://github.com/plasmapy/plasmapy/pull/963>`__)
- The `~plasmapy.plasma.grids.AbstractGrid` property
  `~plasmapy.plasma.grids.AbstractGrid.grid` is now dimensioned (has units) and
  cannot be accessed if all dimensions do not share the same units. (`#981 <https://github.com/plasmapy/plasmapy/pull/981>`__)
- Renamed attribute ``is_uniform_grid`` on `~plasmapy.plasma.grids.AbstractGrid`
  to ``is_uniform``. (`#981 <https://github.com/plasmapy/plasmapy/pull/981>`__)
- Drop Python 3.6 support. (`#987 <https://github.com/plasmapy/plasmapy/pull/987>`__)
- The ``__getitem__`` method of `~plasmapy.plasma.grids.AbstractGrid` now returns
  a `~astropy.units.Quantity` array instead of a reference to a `xarray.DataArray`. (`#1027 <https://github.com/plasmapy/plasmapy/pull/1027>`__)
- Renamed `IonicFraction` to `~plasmapy.particles.ionization_state.IonicLevel`.
  This lays groundwork for future changes, where that class is going to become
  more than a fraction. (`#1046 <https://github.com/plasmapy/plasmapy/pull/1046>`__)


Deprecations and Removals
-------------------------

- The ``particle`` attribute of `~plasmapy.particles.particle_class.Particle`
  has been deprecated in favor of the new ``symbol`` attribute.  The ``particle``
  attribute now issues a `FutureWarning` to indicate that it will be removed in
  a future release. (`#984 <https://github.com/plasmapy/plasmapy/pull/984>`__)


Features
--------

- Created the `~plasmapy.simulation.abstractions.AbstractNormalizations` class
  to serve as an abstract interface for future classes that represent normalizations. (`#859 <https://github.com/plasmapy/plasmapy/pull/859>`__)
- Create the analysis sub-package `plasmapy.analysis.swept_langmuir` for analysis
  code related to analyzing swept Langmuir traces.  Sub-package is initiated with
  functionality for calculating the floating potential,
  `~plasmapy.analysis.swept_langmuir.floating_potential.find_floating_potential`. (`#889 <https://github.com/plasmapy/plasmapy/pull/889>`__)
- Added a proton radiography diagnostic module containing a tool for generating synthetic proton radiographs from simulated or calculated fields using a particle tracking algorithm. (`#895 <https://github.com/plasmapy/plasmapy/pull/895>`__)
- Created new grid objects for representing plasma quantities as functions of space. (`#909 <https://github.com/plasmapy/plasmapy/pull/909>`__)
- Added functions in `plasmapy.utils.code_repr` to reproduce strings
  that represent a call to a method or attribute of an object. These
  functions are used, for example, in error messages. (`#920 <https://github.com/plasmapy/plasmapy/pull/920>`__)
- Add the function
  :func:`~plasmapy.dispersion.two_fluid_dispersion.two_fluid_dispersion_solution` to
  `plasmapy.dispersion`, which gives an analytical solution to the dispersion relation as
  derived by P. M. Bellan 2012 (DOI: `10.1029/2012JA017856
  <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2012JA017856>`_). (`#932 <https://github.com/plasmapy/plasmapy/pull/932>`__)
- Refactor out the `~plasmapy.simulation.particle_integrators.boris_push` tracking
  integrator algorithm from ``plasmapy.simulation.particletracker.ParticleTracker``. (`#953 <https://github.com/plasmapy/plasmapy/pull/953>`__)
- For `plasmapy.plasma.grids` functionality, add better support for recognizing and
  handling physical quantities (e.g. spatial position, magnetic field, etc.) added
  to a grid object. (`#963 <https://github.com/plasmapy/plasmapy/pull/963>`__)
- For `plasmapy.plasma.grids` functionality, improve interpolation performance on
  non-uniform grids. (`#963 <https://github.com/plasmapy/plasmapy/pull/963>`__)
- Added the `~plasmapy.formulary.drifts.diamagnetic_drift` function to
  `~plasmapy.formulary.drifts`. (`#966 <https://github.com/plasmapy/plasmapy/pull/966>`__)
- Add properties `~plasmapy.plasma.grids.AbstractGrid.grid_resolution` and
  `~plasmapy.plasma.grids.AbstractGrid.quantities` to
  `~plasmapy.plasma.grids.AbstractGrid`. (`#981 <https://github.com/plasmapy/plasmapy/pull/981>`__)
- Make several speed improvements to the functionality in `~plasmapy.plasma.grids`,
  including the addition of keyword ``persistent`` to
  `~plasmapy.plasma.grids.AbstractGrid` (and child class) methods
  `~plasmapy.plasma.grids.AbstractGrid.nearest_neighbor_interpolator` and
  `~plasmapy.plasma.grids.AbstractGrid.volume_averaged_interpolator`.  This keyword
  allows the interpolators to assume the last grid setup and contents if input
  arguments have not changed. (`#981 <https://github.com/plasmapy/plasmapy/pull/981>`__)
- Add methods `~plasmapy.plasma.grids.AbstractGrid.on_grid` and
  `~plasmapy.plasma.grids.AbstractGrid.vector_intersects` to
  `~plasmapy.plasma.grids.AbstractGrid`. (`#981 <https://github.com/plasmapy/plasmapy/pull/981>`__)
- The `~plasmapy.particles.particle_class.Particle` class now contains an
  attribute named ``symbol`` that is intended to replace ``particle``. The
  ``symbol`` attribute has been added as a property to
  `~plasmapy.particles.particle_class.AbstractParticle`,
  `~plasmapy.particles.particle_class.CustomParticle`, and
  `~plasmapy.particles.particle_class.DimensionlessParticle`. (`#984 <https://github.com/plasmapy/plasmapy/pull/984>`__)
- Added new ``can_be_zero`` check parameter to
  `~plasmapy.utils.decorators.checks.CheckValues` and its dependents (
  `~plasmapy.utils.decorators.checks.check_values`,
  `~plasmapy.utils.decorators.validators.ValidateQuantities`,
  `~plasmapy.utils.decorators.validators.validate_quantities`). (`#999 <https://github.com/plasmapy/plasmapy/pull/999>`__)
- Both `plasmapy.particles.CustomParticle` and `plasmapy.particles.DimensionlessParticle`
  now allow users to define a custom symbol via the ``symbol`` keyword argument, which
  can then be accessed by the ``symbol`` attribute in each of these classes. (`#1015 <https://github.com/plasmapy/plasmapy/pull/1015>`__)
- The greater than (``>``) operator can now be used between
  `~plasmapy.particles.Particle` and/or `~plasmapy.particles.ParticleList`
  instances to get the nuclear reaction energy. (`#1017 <https://github.com/plasmapy/plasmapy/pull/1017>`__)
- Create `plasmapy.particles.ParticleList` as a list-like collection for
  instances of `plasmapy.particles.Particle` and
  `plasmapy.particles.CustomParticle`.  Adding `~plasmapy.particles.Particle`
  and/or `~plasmapy.particles.CustomParticle` instances will now create a
  `~plasmapy.particles.ParticleList`. (`#1017 <https://github.com/plasmapy/plasmapy/pull/1017>`__)
- Added method `~plasmapy.plasma.grids.AbstractGrid.require_quantities` to
  `~plasmapy.plasma.grids.AbstractGrid` that verifies a list of quantities is present
  on the grid.  Method is also incorporated into
  `~plasmapy.diagnostics.proton_radiography.SyntheticProtonRadiograph`. (`#1027 <https://github.com/plasmapy/plasmapy/pull/1027>`__)
- Added the
  `~plasmapy.diagnostics.proton_radiography.SyntheticProtonRadiograph.add_wire_mesh()`
  method to `~plasmapy.diagnostics.proton_radiography.SyntheticProtonRadiograph`
  to allow the creation of synthetic proton radiographs that include a wire mesh
  reference grid. (`#1049 <https://github.com/plasmapy/plasmapy/pull/1049>`__)
- Created a function, `~plasmapy.formulary.mathematics.rot_a_to_b`, that calculates
  the rotation matrix that will rotate one 3D vector onto another. (`#1054 <https://github.com/plasmapy/plasmapy/pull/1054>`__)
- Made `~plasmapy.plasma.grids.AbstractGrid.is_uniform` a properly-documented
  public attribute of `~plasmapy.plasma.grids.AbstractGrid`. (`#1072 <https://github.com/plasmapy/plasmapy/pull/1072>`__)


Bug Fixes
---------

- Fixed a minus sign bug in the Particle Tracker simulation that caused the
  E×B drift to go in the incorrect direction. (`#953 <https://github.com/plasmapy/plasmapy/pull/953>`__)
- Bugfix :meth:`plasmapy.analysis.fit_functions.Linear.root_solve` to handle the case
  where the slope is zero and no finite roots exist. (`#959 <https://github.com/plasmapy/plasmapy/pull/959>`__)
- Fixed a bug that prevented nested iterations of a single
  `~plasmapy.particles.IonizationState` or
  `~plasmapy.particles.IonizationStateCollection` instance. (`#1025 <https://github.com/plasmapy/plasmapy/pull/1025>`__)
- Fixed a bug in :file:`grids.py` for non-uniform grids that arose when |xarray| upgraded
  to `v0.17.0` (`#1027 <https://github.com/plasmapy/plasmapy/pull/1027>`__)
- In `~plasmapy.diagnostics.proton_radiography.SyntheticProtonRadiograph`,
  adaptive ``dt`` now calculates the cyclotron period using the provided particle
  charge and mass (previously assumed protons). (`#1035 <https://github.com/plasmapy/plasmapy/pull/1035>`__)
- In `~plasmapy.diagnostics.proton_radiography.SyntheticProtonRadiograph`,
  the adaptive timestep algorithm now works when particles are provided using
  `~plasmapy.diagnostics.proton_radiography.SyntheticProtonRadiograph.load_particles`. (`#1035 <https://github.com/plasmapy/plasmapy/pull/1035>`__)
- In `~plasmapy.diagnostics.proton_radiography.SyntheticProtonRadiograph`, removed
  highly deflected particles so the call of
  `~plasmapy.diagnostics.proton_radiography.SyntheticProtonRadiograph.max_deflection`
  does not raise an exception. (`#1035 <https://github.com/plasmapy/plasmapy/pull/1035>`__)


Improved Documentation
----------------------

- Add narrative documentation on ionization state functionality. (`#796 <https://github.com/plasmapy/plasmapy/pull/796>`__)
- Added description to ``plasmapy.formulary.parameters.Hall_parameter``
  signature and equation in docstrings. (`#934 <https://github.com/plasmapy/plasmapy/pull/934>`__)
- Updated documentation for the `plasmapy.particles` and `plasmapy.utils` subpackages. (`#942 <https://github.com/plasmapy/plasmapy/pull/942>`__)
- Improves documentation of `plasmapy/formulary/quantum.py` by cleaning up docstrings of contained functionality. (`#951 <https://github.com/plasmapy/plasmapy/pull/951>`__)
- Update all docstrings associated with computing the Coulomb logarithm and the
  possible methods of calculation. (`#962 <https://github.com/plasmapy/plasmapy/pull/962>`__)
- Add two Jupyter notebooks for functionality contained in `plasmapy.plasma.grids`:
  `grids_cartesian.ipynb` and `grids_nonuniform.ipynb`. (`#963 <https://github.com/plasmapy/plasmapy/pull/963>`__)
- Added the ExB drift notebook, which demonstrates the analytical solution for the
  drift and the implementation of the corresponding formulary drift functions,
  `~plasmapy.formulary.drifts.force_drift` and `~plasmapy.formulary.drifts.ExB_drift`. (`#971 <https://github.com/plasmapy/plasmapy/pull/971>`__)
- Describe what constitutes a valid representation of a particle in the docstring
  for the `plasmapy.particles.particle_class.ParticleLike` typing construct. (`#985 <https://github.com/plasmapy/plasmapy/pull/985>`__)
- Put the docstring for `plasmapy.particles.Particle.is_category` into
  `numpydoc` format. (`#1039 <https://github.com/plasmapy/plasmapy/pull/1039>`__)
- Adds formulas (which were missing) to the docstrings of
  ``plasmapy.formulary.dimensionless.quantum_theta`` and
  `~plasmapy.formulary.dimensionless.beta`. (`#1041 <https://github.com/plasmapy/plasmapy/pull/1041>`__)
- Add live rendering of changelog entries on documentation builds, based on
  `sphinx-changelog <https://github.com/OpenAstronomy/sphinx-changelog>`_. (`#1052 <https://github.com/plasmapy/plasmapy/pull/1052>`__)
- Created an example notebook demonstrating how the
  `~plasmapy.diagnostics.proton_radiography.SyntheticProtonRadiograph` class can be
  used to generate synthetic proton radiographs with arbitrary source profiles.  Add
  code documentation links to all proton radiograph notebooks. (`#1054 <https://github.com/plasmapy/plasmapy/pull/1054>`__)
- Update formatting and broken `sphinx.ext.intersphinx` links in `plasmapy.formulary` docstrings. (`#1058 <https://github.com/plasmapy/plasmapy/pull/1058>`__)
- Make minor fixes in `plasmapy.particles` docstrings. (`#1064 <https://github.com/plasmapy/plasmapy/pull/1064>`__)
- Organize the layout of the example Jupyter notebooks on the Read the Docs
  example page. (`#1066 <https://github.com/plasmapy/plasmapy/pull/1066>`__)
- Fix formatting and broken `sphinx.ext.intersphinx` links in docstrings in
  various places in the code base. Improve installation instructions in the docs;
  the subpackage stability matrix, and funding acknowledgments. (`#1076 <https://github.com/plasmapy/plasmapy/pull/1076>`__)


Trivial/Internal Changes
------------------------

- Removed `colorama` as a dependency. (`#920 <https://github.com/plasmapy/plasmapy/pull/920>`__)
- Moved remaining CI from CircleCI to GitHub Actions. (`#996 <https://github.com/plasmapy/plasmapy/pull/996>`__)
- Add notebook CI through `nbqa`. (`#997 <https://github.com/plasmapy/plasmapy/pull/997>`__)
- Remove `lambda` expressions from `plasmapy.particles` and `plasmapy.utils`. (`#1013 <https://github.com/plasmapy/plasmapy/pull/1013>`__)
- Add unicode particle aliases for electrons (``"β-"``, ``"β⁻"``), muons
  (``"μ-"``, ``"μ⁻"``), anti-muons (``"μ+"``, ``"μ⁺"``), tau particles
  (``"τ"``, ``"τ-"``, ``"τ⁻"``), anti-tau particles (``"τ+"``, ``"τ⁺"``)
  electron neutrinos (``"ν_e"``), muon neutrinos (``"ν_μ"``), tau neutrinos
  (``"ν_τ"``), and alpha particles (``"α"``). (`#1036 <https://github.com/plasmapy/plasmapy/pull/1036>`__)
- A set containing all valid particle categories may now be accessed via
  `plasmapy.particles.Particle.is_category.valid_categories`. (`#1039 <https://github.com/plasmapy/plasmapy/pull/1039>`__)
- Properly handled warnings in `test_proton_radiography.py` (`#1050 <https://github.com/plasmapy/plasmapy/pull/1050>`__)


PlasmaPy v0.5.0 (2020-12-09)
============================

Backwards Incompatible Changes
------------------------------

- Created `plasmapy.dispersion` in accordance with PlasmaPy Enhancement Proposal 7
  (`PLEP 7 <https://github.com/PlasmaPy/PlasmaPy-PLEPs/blob/main/PLEP-0007.rst>`_)
  and migrated the dispersion functionality (:file:`dispersionfunction.py`) from
  `plasmapy.formulary` to `plasmapy.dispersion`. (`#910 <https://github.com/plasmapy/plasmapy/pull/910>`__)
- Removed default values for the `ion` and `particle` arguments of functions contained in ``plasmapy.formulary.parameters``, in accordance with issue [#453](https://github.com/PlasmaPy/PlasmaPy/issues/453), and updated all relevant calls to modified functionality. (`#911 <https://github.com/plasmapy/plasmapy/pull/911>`__)
- Moved test helper exceptions from ``plasmapy.utils.pytest_helpers`` to ``plasmapy.tests.helpers``. (`#919 <https://github.com/plasmapy/plasmapy/pull/919>`__)
- Update ``plasmapy.formulary.parameters.mass_density`` so it calculates the mass
  density for a specific particle from a given number density.  Original function
  calculated the total mass density (ion + electron). (`#957 <https://github.com/plasmapy/plasmapy/pull/957>`__)


Features
--------

- Added a function to calculate the power spectrum of thermal bremsstrahlung emitted by a Maxwellian plasma. (`#892 <https://github.com/plasmapy/plasmapy/pull/892>`__)
- Added support for multiple electron components to diagnostics.thomson.spectral_density. Also fixed a bug for multiple ion populations. (`#893 <https://github.com/plasmapy/plasmapy/pull/893>`__)
- Add dependency `pygments >= 2.4.1`. (`#898 <https://github.com/plasmapy/plasmapy/pull/898>`__)
- Create the `plasmapy.analysis` package as per
  `PLEP-7 <https://github.com/PlasmaPy/PlasmaPy-PLEPs/blob/main/PLEP-0007.rst>`_ and
  initialize the package with the `~plasmapy.analysis.fit_functions` module.  Fit
  functions are designed to wrap together an analytical function, a curve fitter,
  uncertainty propagation, and a root solver to make curve fitting a little less
  painful. (`#908 <https://github.com/plasmapy/plasmapy/pull/908>`__)
- Created a new subpackage, `plasmapy.tests.helpers`, to contain test helper functionality. (`#919 <https://github.com/plasmapy/plasmapy/pull/919>`__)
- Create decorator `~plasmapy.utils.decorators.helpers.modify_docstring`, which allows
  for programmatically prepending and/or appending a docstring. (`#943 <https://github.com/plasmapy/plasmapy/pull/943>`__)


Bug Fixes
---------

- Allowed implicit conversions of AstroPy units in inputs and outputs of validated functions to happen without warnings. Most notably, this removes warnings on eV inputs to temperature fields. (`#886 <https://github.com/plasmapy/plasmapy/pull/886>`__)
- Update ``plasmapy.formulary.parameters.Alfven_speed`` to properly use the updated
  ``plasmapy.formulary.parameters.mass_density`` and maintain the same behavior.
  Also add handling of the ``ion`` input keyword, so |Particle| and
  the |Particle| convertible representations can be used as inputs. (`#957 <https://github.com/plasmapy/plasmapy/pull/957>`__)


Improved Documentation
----------------------

- Improved the release guide after the release of 0.4.0. (`#872 <https://github.com/plasmapy/plasmapy/pull/872>`__)
- Add various improvements to the documentation.
      * Replace home link with the plasmapy logo.
      * Add module and index navigation links to sidebar header.
      * Replace raw html on the main page that simulates a `nbgallery` with a real
        `nbgallery` directive.
      * Move link to view page source code from the header to footer.
      * Add link to footer the jumps the user back to the top of the page.
      * Create and add custom CSS stylesheet.
      * Create `_templates` directory and templates to customize page elements. (`#875 <https://github.com/plasmapy/plasmapy/pull/875>`__)
- Add static stub files to `docs/api_static` so all modules of `plasmapy` are indexed.
  This is necessary to expose all of `plasmapy` since not all modules are indexed in
  the narrative documentation. (`#878 <https://github.com/plasmapy/plasmapy/pull/878>`__)
- Decompose sub-package `plasmapy/utils/roman/` into the `plasmapy/utils/roman.py`
  file.  Move definition of `roman` specific `Exceptions` into
  `plasmapy.utils.exceptions`. (`#883 <https://github.com/plasmapy/plasmapy/pull/883>`__)
- Replaced references to Riot.im with references to Element.io or Matrix, as appropriate, following their recent rebranding. (`#891 <https://github.com/plasmapy/plasmapy/pull/891>`__)
- Update the information on how to cite PlasmaPy, including in the release guide. (`#900 <https://github.com/plasmapy/plasmapy/pull/900>`__)


Trivial/Internal Changes
------------------------

- Apply isort to entire codebase, bringing it back to the pre-commit hook suite. (`#857 <https://github.com/plasmapy/plasmapy/pull/857>`__)
- Expand package metadata contained in ``codemeta.json``, following the CodeMeta standard. (`#902 <https://github.com/plasmapy/plasmapy/pull/902>`__)
- Changed remaining instances of @u.quantity_input to @validate_quantities in response to issue #880. (`#905 <https://github.com/plasmapy/plasmapy/pull/905>`__)
- Switched from Azure Pipelines to GitHub Actions for PR tests to make things
  easier for contributors. Moved away from Travis CI for test cron jobs. (`#952 <https://github.com/plasmapy/plasmapy/pull/952>`__)


PlasmaPy v0.4.0 (2020-07-20)
============================

Backwards Incompatible Changes
------------------------------

- Rename ``plasmapy.atomic`` to `~plasmapy.particles`.  In
  `~plasmapy.formulary.collisions` and `~plasmapy.formulary.braginskii`,
  change arguments named particles to ``species`` and arguments named
  ``ion_particle`` to ``ion`` for multiple functions. (`#742 <https://github.com/plasmapy/plasmapy/pull/742>`__)
- Officially delete :mod:`plasmapy.examples`. (`#822 <https://github.com/plasmapy/plasmapy/pull/822>`__)
- Move :mod:`plasmapy.data` to :mod:`plasmapy.particle.data`. (`#823 <https://github.com/plasmapy/plasmapy/pull/823>`__)
- Renamed the `plasmapy.classes` subpackage to `plasmapy.plasma`. (`#842 <https://github.com/plasmapy/plasmapy/pull/842>`__)


Features
--------

- Added units to reprs of .formulary.magnetostatics classes. (`#743 <https://github.com/plasmapy/plasmapy/pull/743>`__)
- Create prototype abstract interfaces for plasma simulations (`#753 <https://github.com/plasmapy/plasmapy/pull/753>`__)
- Created classes to represent custom and dimensionless particles in ``plasmapy.particles``. (`#755 <https://github.com/plasmapy/plasmapy/pull/755>`__)
- Create :func:`~plasmapy.formulary.relativity.relativistic_energy` function, which uses the established :func:`~plasmapy.formulary.relativity.Lorentz_factor` function to aid in the calculation of the relativistic energy of an object. (`#805 <https://github.com/plasmapy/plasmapy/pull/805>`__)
- Create :func:`~plasmapy.formulary.dimensionless.Reynolds_number` function. (`#815 <https://github.com/plasmapy/plasmapy/pull/815>`__)
- Create :func:`~plasmapy.formulary.dimensionless.Mag_Reynolds` function. (`#820 <https://github.com/plasmapy/plasmapy/pull/820>`__)
- Create ``plasmapy.formulary.parameters.Bohm_diffusion`` function. (`#830 <https://github.com/plasmapy/plasmapy/pull/830>`__)
- Added a new diagnostics module `thomson` containing a function
  `spectral_density` that calculates Thomson scattering spectra for
  Maxwellian plasmas in both the collective and non-collective regimes. As
  a followup to PR #835, set the minimal required Numpy version to 1.18.1 to
  finally fix unit dropping bugs. (`#831 <https://github.com/plasmapy/plasmapy/pull/831>`__)
- Revised parameters.thermal_speed to support 1D and 2D distributions as well as 3D, and added an example notebook for this function. (`#850 <https://github.com/plasmapy/plasmapy/pull/850>`__)
- Create `plasmapy/formulary/ionization.py`
  Create :func:`~plasmapy.formulary.ionization.Z_bal` function. (`#851 <https://github.com/plasmapy/plasmapy/pull/851>`__)
- Create :func:`~plasmapy.formulary.ionization.Saha` function. (`#860 <https://github.com/plasmapy/plasmapy/pull/860>`__)
- Added aliases (with trailing underscores) for parameters in the formulary:

      * `plasmapy.formulary.dimensionless.Reynolds_number` -> `~plasmapy.formulary.dimensionless.Re_`
      * `plasmapy.formulary.dimensionless.Mag_Reynolds` -> `~plasmapy.formulary.dimensionless.Rm_`
      * `plasmapy.formulary.drifts.ExB_drift` -> `~plasmapy.formulary.drifts.veb_`
      * `plasmapy.formulary.drifts.force_drift` -> `~plasmapy.formulary.drifts.vfd_`
      * ``plasmapy.formulary.parameters.mass_density`` -> ``plasmapy.formulary.parameters.rho_``
      * ``plasmapy.formulary.parameters.Afven_speed`` -> ``plasmapy.formulary.parameters.va_``
      * ``plasmapy.formulary.parameters.ion_sound_speed`` -> ``plasmapy.formulary.parameters.cs_``
      * ``plasmapy.formulary.parameters.thermal_speed`` -> ``plasmapy.formulary.parameters.vth_``
      * ``plasmapy.formulary.parameters.thermal_pressure`` -> ``plasmapy.formulary.parameters.pth_``
      * ``plasmapy.formulary.parameters.kappa_thermal_speed`` -> ``plasmapy.formulary.parameters.vth_kappa_``
      * ``plasmapy.formulary.parameters.inertial_length`` -> ``plasmapy.formulary.parameters.cwp_``
      * ``plasmapy.formulary.parameters.Hall_parameter`` -> ``plasmapy.formulary.parameters.betaH_``
      * ``plasmapy.formulary.parameters.gyrofrequency`` -> ``plasmapy.formulary.parameters.oc_``, ``plasmapy.formulary.parameters.wc_``
      * ``plasmapy.formulary.parameters.gyroradius`` -> ``plasmapy.formulary.parameters.rc_``, ``plasmapy.formulary.parameters.rhoc_``
      * ``plasmapy.formulary.parameters.plasma_frequency`` -> ``plasmapy.formulary.parameters.wp_``
      * ``plasmapy.formulary.parameters.Debye_length`` -> ``plasmapy.formulary.parameters.lambdaD_``
      * ``plasmapy.formulary.parameters.Debye_number`` -> ``plasmapy.formulary.parameters.nD_``
      * ``plasmapy.formulary.parameters.magnetic_pressure`` -> ``plasmapy.formulary.parameters.pmag_``
      * ``plasmapy.formulary.parameters.magnetic_energy_density`` -> ``plasmapy.formulary.parameters.ub_``
      * ``plasmapy.formulary.parameters.upper_hybrid_frequency`` -> ``plasmapy.formulary.parameters.wuh_``
      * ``plasmapy.formulary.parameters.lower_hybrid_frequency`` -> ``plasmapy.formulary.parameters.wlh_``
      * ``plasmapy.formulary.parameters.Bohm_diffusion`` -> ``plasmapy.formulary.parameters.DB_``
      * `plasmapy.formulary.quantum.deBroglie_wavelength` -> `~plasmapy.formulary.quantum.lambdaDB_`
      * `plasmapy.formulary.quantum.thermal_deBroglie_wavelength` -> `~plasmapy.formulary.quantum.lambdaDB_th_`
      * `plasmapy.formulary.quantum.Fermi_energy` -> `~plasmapy.formulary.quantum.Ef_` (`#865 <https://github.com/plasmapy/plasmapy/pull/865>`__)
- Add `json_dumps` method to `~plasmapy.particles.particle_class.AbstractParticle` to
  convert a particle object into a JSON string. Add `json_dump` method to
  `~plasmapy.particles.particle_class.AbstractParticle` to serialize a particle
  object and writes it to a file.  Add JSON decoder
  `~plasmapy.particles.serialization.ParticleJSONDecoder` to deserialize JSON objects
  into particle objects.  Add `plasmapy.particles.serialization.json_loads_particle`
  function to convert JSON strings to particle objects (using
  `~plasmapy.particles.serialization.ParticleJSONDecoder`). Add
  `plasmapy.particles.json_load_particle` function to deserialize a JSON file into a
  particle object (using `~plasmapy.particles.serialization.ParticleJSONDecoder`).
  (`#836 <https://github.com/plasmapy/plasmapy/pull/836>`__)


Bug Fixes
---------

- Fix incorrect use of `pkg.resources` when defining `plasmapy.__version__`.  Add
  `setuptools` to package dependencies.  Add a definition of `__version__` for
  developers using source files. (`#774 <https://github.com/plasmapy/plasmapy/pull/774>`__)
- Repair notebook links that are defined in the `nbsphinx_prolog` sphinx configuration
  variable. (`#828 <https://github.com/plasmapy/plasmapy/pull/828>`__)
- Increase the required Astropy version from 3.1 to 4.0, Numpy from 1.14 to 1.16.6, Scipy from 0.19 to 1.2 and lmfit from 0.9.7 to 1.0.1. This fixes long-standing issues with Numpy operations dropping units from AstroPy quantities. (`#835 <https://github.com/plasmapy/plasmapy/pull/835>`__)


Improved Documentation
----------------------

- - Added documentation to file test_converters (`#756 <https://github.com/plasmapy/plasmapy/pull/756>`__)
- - Updated installation instructions. (`#772 <https://github.com/plasmapy/plasmapy/pull/772>`__)
- Reorder documentation page (`#777 <https://github.com/plasmapy/plasmapy/pull/777>`__)
- Fix failing documentation build due to duplicate docstrings for
  `ParticleTracker.kinetic_energy_history` and incompatibility of `sphinx-automodapi`
  with `sphinx` `v3.0.0`. (`#780 <https://github.com/plasmapy/plasmapy/pull/780>`__)
- Automate definition of documentation `release` and `version` in `docs/conf.py` with
  `plasmapy.__version__`. (`#781 <https://github.com/plasmapy/plasmapy/pull/781>`__)
- Add a docstring to ``__init__.py`` in `plasmapy.formulary`. (`#788 <https://github.com/plasmapy/plasmapy/pull/788>`__)
- Replaced sphinx-gallery with nbsphinx, turning `.py` example files into `.ipynb` files and allowing for easier example submission. (`#792 <https://github.com/plasmapy/plasmapy/pull/792>`__)
- Linked various instances of classes and functions in the `.ipynb` examples in `docs/notebooks/` to the respective API docs. (`#825 <https://github.com/plasmapy/plasmapy/pull/825>`__)
- Fixed a few documentation formatting errors. (`#827 <https://github.com/plasmapy/plasmapy/pull/827>`__)
- Add notes on the PlasmaPy benchmarks repository to documentation. (`#841 <https://github.com/plasmapy/plasmapy/pull/841>`__)
- Improve readability of the `plasmapy.formulary` page by replacing the `toctree`
  list with a cleaner reStructuredText table. (`#867 <https://github.com/plasmapy/plasmapy/pull/867>`__)


Trivial/Internal Changes
------------------------

- Remove mutable arguments from `Particle.is_category` method. (`#751 <https://github.com/plasmapy/plasmapy/pull/751>`__)
- Remove all occurrences of default mutable arguments (`#754 <https://github.com/plasmapy/plasmapy/pull/754>`__)
- Handle `ModuleNotFoundError` when trying to import `__version__` but `setuptools_scm` has not
  generated the `version.py` file.  This commonly happens during development when `plasmapy` is
  not installed in the python environment. (`#763 <https://github.com/plasmapy/plasmapy/pull/763>`__)
- Updated pep8speaks/flake8 configuration and added `.pre-commit-config.yaml` to simplify automated style checks during development. (`#770 <https://github.com/plasmapy/plasmapy/pull/770>`__)
- Removes some lint from setup.py and setup.cfg. Use pkg_resources for version
  checking in code. Remove version.py file in favor of pkg_resources. (`#771 <https://github.com/plasmapy/plasmapy/pull/771>`__)
- Default settings for isort were set to be consistent with default settings for black. (`#773 <https://github.com/plasmapy/plasmapy/pull/773>`__)
- Update community meeting and funding information in docs. (`#784 <https://github.com/plasmapy/plasmapy/pull/784>`__)
- Improved pull request template to include more information about changelog entries. (`#843 <https://github.com/plasmapy/plasmapy/pull/843>`__)
- Added GitHub actions that apply pre-commit and flake8 (separately) to incoming pull requests. (`#845 <https://github.com/plasmapy/plasmapy/pull/845>`__)
- Apply pre-commit hooks to entire repository, so that GitHub actions do not shout at contributors needlessly. (`#846 <https://github.com/plasmapy/plasmapy/pull/846>`__)
- Update :class:`~plasmapy.particles.particle_class.CustomParticle` so input parameters
  `mass` and `charge` can accept string representations of astropy `Quantities`. (`#862 <https://github.com/plasmapy/plasmapy/pull/862>`__)


PlasmaPy v0.3.0 (2020-01-25)
============================

Backwards Incompatible Changes
------------------------------

- Create simulation subpackage; move Species particle tracker there; rename to particletracker (`#665 <https://github.com/plasmapy/plasmapy/pull/665>`__)
- Changed `plasmapy.classes.Species` to `plasmapy.simulation.ParticleTracker` (`#668 <https://github.com/plasmapy/plasmapy/pull/668>`__)
- Move pytest helper functionality from `plasmapy.utils` to
  ``plasmapy.utils.pytest_helpers`` (`#674 <https://github.com/plasmapy/plasmapy/pull/674>`__)
- Move `plasmapy.physics`, `plasmapy.mathematics` and `plasmapy.transport` into the common `plasmapy.formulary` subpackage (`#692 <https://github.com/plasmapy/plasmapy/pull/692>`__)
- Change `ClassicalTransport` methods into attributes (`#705 <https://github.com/plasmapy/plasmapy/pull/705>`__)

Deprecations and Removals
-------------------------

- Remove `parameters_cython.pyx`, switching to Numba for the future of computationally intensive code in PlasmaPy (`#650 <https://github.com/plasmapy/plasmapy/pull/650>`__)
- Remove plasmapy.constants, which was a thin wrapper around astropy.constants
  with no added value (`#651 <https://github.com/plasmapy/plasmapy/pull/651>`__)

Features
--------

- Generalize `ion_sound_speed` function to work for all values of :math:`k^2 \lambda_{D}^2` (i.e. not just in the non-dispersive limit). (`#700 <https://github.com/plasmapy/plasmapy/pull/700>`__)
- Optimize `add__magnetostatics` for a 16x speedup in tests! (`#703 <https://github.com/plasmapy/plasmapy/pull/703>`__)

Bug Fixes
---------

- Define `preserve_signature` decorator to help IDEs parse signatures of decorated functions. (`#640 <https://github.com/plasmapy/plasmapy/pull/640>`__)
- Fix Pytest deprecations of `message` argument to `raise` and `warn` functions. (`#666 <https://github.com/plasmapy/plasmapy/pull/666>`__)
- Fix `h5py` warning in OpenPMD module, opening files in read mode by default (`#717 <https://github.com/plasmapy/plasmapy/pull/717>`__)


Improved Documentation
----------------------

- Added real-world examples to examples/plot_physics.py and adjusted the plots to be more human-friendly. (`#448 <https://github.com/plasmapy/plasmapy/pull/448>`__)
- Add examples images to the top of the main doc page in `docs\index.rst` (`#655 <https://github.com/plasmapy/plasmapy/pull/655>`__)
- Added exampes to the documentation to mass_density
   and Hall_parameter functions (`#709 <https://github.com/plasmapy/plasmapy/pull/709>`__)
- Add docstrings to decorator :func:`plasmapy.utils.decorators.converter.angular_freq_to_hz`. (`#729 <https://github.com/plasmapy/plasmapy/pull/729>`__)


Trivial/Internal Changes
------------------------

- Replace decorator :func:`plasmapy.utils.decorators.checks.check_quantity` with decorator
  :func:`plasmapy.utils.decorators.validators.validate_quantities`.  Permanently delete decorator
  :func:`~plasmapy.utils.decorators.checks.check_quantity` and its supporting code.  For functions
  :func:`plasmapy.formulary.quantum.chemical_potential` and
  :func:`plasmapy.formulary.quantum._chemical_potential_interp`, add a `RaiseNotImplementedError` due
  to bug outlined in issue `<https://github.com/PlasmaPy/PlasmaPy/issues/726>`_.  Associated pytests
  are marked with `pytest.mark.xfails` and doctests are marked with `doctests: +SKIP`. (`#722 <https://github.com/plasmapy/plasmapy/pull/722>`__)
- Add `Towncrier <https://github.com/hawkowl/towncrier>`_ automated changelog creation support (`#643 <https://github.com/plasmapy/plasmapy/pull/643>`__)
- Move existing "check" decorators to new ``plasmapy.utils.decorators`` module (`#647 <https://github.com/plasmapy/plasmapy/pull/647>`__)
- Allow running our sphinx-gallery examples as Jupyter notebooks via Binder (`#656 <https://github.com/plasmapy/plasmapy/pull/656>`__)
- Overhaul CI setup, following the example of SunPy (`#657 <https://github.com/plasmapy/plasmapy/pull/657>`__)
- Patch `sphinx_gallery.binder` to output custom links to Binder instance (`#658 <https://github.com/plasmapy/plasmapy/pull/658>`__)
- Remove the now unnecessary `astropy_helpers` submodule (`#663 <https://github.com/plasmapy/plasmapy/pull/663>`__)
- Followup PR to CI overhaul (`#664 <https://github.com/plasmapy/plasmapy/pull/664>`__)
- Add a Codemeta file (``codemeta.json``) (`#676 <https://github.com/plasmapy/plasmapy/pull/676>`__)
- Overhaul and simplify CI, add Python 3.8 to tests, bump minimal required package versions, fix docs. (`#712 <https://github.com/plasmapy/plasmapy/pull/712>`__)
- Update communication channels in docs (`#715 <https://github.com/plasmapy/plasmapy/pull/715>`__)
- Code style fixes to the `atomic` subpackage (`#716 <https://github.com/plasmapy/plasmapy/pull/716>`__)
- Clean up main package namespace, removing `plasmapy.test` (`#718 <https://github.com/plasmapy/plasmapy/pull/718>`__)
- Reduce precision of tests and doctests to allow for refinements of
  fundamental constants. (`#731 <https://github.com/plasmapy/plasmapy/pull/731>`__)
- Create decorators for checking/validating values and units of function/method input
  and return arguments.  Defined decorators include
  :func:`~plasmapy.utils.decorators.checks.check_values`,
  :func:`~plasmapy.utils.decorators.checks.check_units`, and
  :func:`~plasmapy.utils.decorators.validators.validate_quantities`.  These decorators are
  fully defined by "decorator classes" :class:`~plasmapy.utils.decorators.checks.CheckBase`,
  :class:`~plasmapy.utils.decorators.checks.CheckValues`,
  :class:`~plasmapy.utils.decorators.checks.CheckUnits`, and
  :class:`~plasmapy.utils.decorators.validators.ValidateQuantities`. (`#648 <https://github.com/plasmapy/plasmapy/pull/648>`__)
- Create a decorator to change output of physics functions from "radians/s" to "hz" (`#667 <https://github.com/plasmapy/plasmapy/pull/667>`__)
- Added pytest.mark.slow to pytest markers.
  Updated documentation to notify developers of functionality. (`#677 <https://github.com/plasmapy/plasmapy/pull/677>`__)
