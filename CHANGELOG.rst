Plasmapy v0.6.0 (2021-03-14)
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
  for functions `~plasmapy.formulary.collisions.Coulomb_logarithm` and
  `~plasmapy.formulary.collisions.impact_parameter`, and then propogated throughout
  the functionality in `plasmapy.formulary.collisions`. (`#962 <https://github.com/plasmapy/plasmapy/pull/962>`__)
- Add dependency ``pandas >= 1.0.0``.  Modify `xarray` dependency to be
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
  integrator algorithm from `~plasmapy.simulation.particletracker.ParticleTracker`. (`#953 <https://github.com/plasmapy/plasmapy/pull/953>`__)
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
  on the grid.  Method is also incorported into
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
- Fixed a bug in `grids.py` for non-uniform grids that arose when `xarray` upgraded
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
- Added description to :func:`~plasmapy.formulary.parameters.Hall_parameter`
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
  `~plasmapy.formulary.dimensionless.quantum_theta` and
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


Plasmapy v0.5.0 (2020-12-09)
============================

Backwards Incompatible Changes
------------------------------

- Created `plasmapy.dispersion` in accordance with PlasmaPy Enhancement Proposal 7
  (`PLEP 7 <https://github.com/PlasmaPy/PlasmaPy-PLEPs/blob/master/PLEP-0007.rst>`_)
  and migrated the dispersion functionality (`dispersionfunction.py`) from
  `plasmapy.formulary` to `plasmapy.dispersion`. (`#910 <https://github.com/plasmapy/plasmapy/pull/910>`__)
- Removed default values for the `ion` and `particle` arguments of functions contained in `plasmapy.formulary.parameters`, in accordance with issue [#453](https://github.com/PlasmaPy/PlasmaPy/issues/453), and updated all relevant calls to modified functionality. (`#911 <https://github.com/plasmapy/plasmapy/pull/911>`__)
- Moved test helper exceptions from `plasmapy.utils.pytest_helpers` to `plasmapy.tests.helpers`. (`#919 <https://github.com/plasmapy/plasmapy/pull/919>`__)
- Update :func:`plasmapy.formulary.parameters.mass_density` so it calculates the mass
  density for a specific particle from a given number density.  Original function
  calculated the total mass density (ion + electron). (`#957 <https://github.com/plasmapy/plasmapy/pull/957>`__)


Features
--------

- Added a function to calculate the power spectrum of thermal bremsstrahlung emitted by a Maxwellian plasma. (`#892 <https://github.com/plasmapy/plasmapy/pull/892>`__)
- Added support for multiple electron components to diagnostics.thomson.spectral_density. Also fixed a bug for multiple ion populations. (`#893 <https://github.com/plasmapy/plasmapy/pull/893>`__)
- Add dependency `pygments >= 2.4.1`. (`#898 <https://github.com/plasmapy/plasmapy/pull/898>`__)
- Create the `plasmapy.analysis` package as per
  `PLEP-7 <https://github.com/PlasmaPy/PlasmaPy-PLEPs/blob/master/PLEP-0007.rst>`_ and
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
- Update :func:`plasmapy.formulary.parameters.Alfven_speed` to properly use the updated
  :func:`~plasmapy.formulary.parameters.mass_density` and maintain the same behavior.
  Also add handling of the ``ion`` input keyword, so `~plasmapy.particles.Particle` and
  the `~plasmapy.particles.Particle` convertible representations can be used as inputs. (`#957 <https://github.com/plasmapy/plasmapy/pull/957>`__)


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


Plasmapy v0.4.0 (2020-07-20)
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
- Create :func:`~plasmapy.formulary.relativity.relativistic_energy` function, which uses the established :func:`~plamsapy.formulary.relativity.Lorentz_factor` function to aid in the calculation of the relativistic energy of an object. (`#805 <https://github.com/plasmapy/plasmapy/pull/805>`__)
- Create :func:`~plasmapy.formulary.dimensionless.Reynolds_number` function. (`#815 <https://github.com/plasmapy/plasmapy/pull/815>`__)
- Create :func:`~plasmapy.formulary.dimensionless.Mag_Reynolds` function. (`#820 <https://github.com/plasmapy/plasmapy/pull/820>`__)
- Create :func:`~plasmapy.formulary.parameters.Bohm_diffusion` function. (`#830 <https://github.com/plasmapy/plasmapy/pull/830>`__)
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
      * `plasmapy.formulary.parameters.mass_density` -> `~plasmapy.formulary.parameters.rho_`
      * `plasmapy.formulary.parameters.Afven_speed` -> `~plasmapy.formulary.parameters.va_`
      * `plasmapy.formulary.parameters.ion_sound_speed` -> `~plasmapy.formulary.parameters.cs_`
      * `plasmapy.formulary.parameters.thermal_speed` -> `~plasmapy.formulary.parameters.vth_`
      * `plasmapy.formulary.parameters.thermal_pressure` -> `~plasmapy.formulary.parameters.pth_`
      * `plasmapy.formulary.parameters.kappa_thermal_speed` -> `~plasmapy.formulary.parameters.vth_kappa_`
      * `plasmapy.formulary.parameters.inertial_length` -> `~plasmapy.formulary.parameters.cwp_`
      * `plasmapy.formulary.parameters.Hall_parameter` -> `~plasmapy.formulary.parameters.betaH_`
      * `plasmapy.formulary.parameters.gyrofrequency` -> `~plasmapy.formulary.parameters.oc_`, `~plasmapy.formulary.parameters.wc_`
      * `plasmapy.formulary.parameters.gyroradius` -> `~plasmapy.formulary.parameters.rc_`, `~plasmapy.formulary.parameters.rhoc_`
      * `plasmapy.formulary.parameters.plasma_frequency` -> `~plasmapy.formulary.parameters.wp_`
      * `plasmapy.formulary.parameters.Debye_length` -> `~plasmapy.formulary.parameters.lambdaD_`
      * `plasmapy.formulary.parameters.Debye_number` -> `~plasmapy.formulary.parameters.nD_`
      * `plasmapy.formulary.parameters.magnetic_pressure` -> `~plasmapy.formulary.parameters.pmag_`
      * `plasmapy.formulary.parameters.magnetic_energy_density` -> `~plasmapy.formulary.parameters.ub_`
      * `plasmapy.formulary.parameters.upper_hybrid_frequency` -> `~plasmapy.formulary.parameters.wuh_`
      * `plasmapy.formulary.parameters.lower_hybrid_frequency` -> `~plasmapy.formulary.parameters.wlh_`
      * `plasmapy.formulary.parameters.Bohm_diffusion` -> `~plasmapy.formulary.parameters.DB_`
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
  list with a cleaner reST table. (`#867 <https://github.com/plasmapy/plasmapy/pull/867>`__)


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


Plasmapy v0.3.0 (2020-01-25)
============================

Backwards Incompatible Changes
------------------------------

- Create simulation subpackage; move Species particle tracker there; rename to particletracker (`#665 <https://github.com/plasmapy/plasmapy/pull/665>`__)
- Changed `plasmapy.classes.Species` to `plasmapy.simulation.ParticleTracker` (`#668 <https://github.com/plasmapy/plasmapy/pull/668>`__)
- Move pytest helper functionality from `plasmapy.utils` to
  `~plasmapy.utils.pytest_helpers` (`#674 <https://github.com/plasmapy/plasmapy/pull/674>`__)
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
