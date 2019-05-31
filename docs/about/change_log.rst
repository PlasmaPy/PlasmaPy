.. _change-log:

##########
Change Log
##########

This document lists the changes made during each release of PlasmaPy,
including bug fixes and changes to the application programming interface
(API).  The :ref:`release-notes` summarize the changes for each version.

.. _change-log-0.2.0:

Version 0.2.0
-------------

Version 0.2.0 is the second development release of PlasmaPy. Alongside a few
new features, it brings plentiful refactoring, documentation and back stage
improvements.

.. _change-log-0.2.0-new:

New Features
~~~~~~~~~~~~

- Implement machinery for a ``Plasma`` class factory based on
  `PLEP 6 <http://doi.org/10.5281/zenodo.1460977>`__
- Create an openPMD ``Plasma`` subclass
- Create classes to represent ionization state distributions for one
  or more elements or isotopes.
- Add basic particle drifts to `plasmapy.physics.drifts`
- Turn most dependencies into optional, subpackage-specific ones

.. _change-log-0.2.0-bugfix:

Bug Fixes
~~~~~~~~~

- Improve handling of NumPy arrays for plasma parameter and transport functions.
- Vendor the `roman` package so as to allow installation via Conda
- Decrease strictness of `check_quantity` to allow `nan` and `inf` by default

.. _change-log-0.2.0-api:

Changes to API
~~~~~~~~~~~~~~

- Move `~plasmapy.transport` from `~plasmapy.physics` to its own
  subpackage.

.. _change-log-0.1.1:

Version 0.1.1
-------------

Version 0.1.1 is a bugfix patch release correcting a number of issues
that arose during the release process and adding two minor convenience
features.

.. _change-log-0.1.1-new:

New Features
~~~~~~~~~~~~

- Add `plasmapy.online_help()`
- Add `plasmapy.__citation__` containing a BibTeX reference.

.. _change-log-0.1.1-bugfix:

Bug Fixes
~~~~~~~~~

- Bring back mistakenly removed Cython versions of plasma parameters.
- Optimize `check_relativistic`.
- Correct a failing import statement.
- Fix a number of issues with the Maxwellian distribution in `physics.distribution`.

.. _change-log-0.1.0:

Version 0.1.0
-------------

Version 0.1.0 is the initial development release of PlasmaPy.  This
version is a prototype and a preview, and is not feature complete.
Significant changes to the API are expected to occur between versions
0.1.0 and 0.2.0, including backward incompatible changes.

.. _change-log-0.1.0-new:

New Features
~~~~~~~~~~~~

* Composed :ref:`plasmapy-vision-statement`.

* Adopted the :ref:`plasmapy-code-of-conduct`.

* Created a guide on :ref:`contributing-to-plasmapy`.

* Adopted a permissive BSD 3-clause `license
  <https://github.com/PlasmaPy/PlasmaPy/blob/master/LICENSE.md>`_ with
  protections against software patents.

* Set up continuous integration testing with `Travis CI
  <https://travis-ci.org/>`_, `CircleCI <https://circleci.com/>`_, and
  `AppVeyor <https://www.appveyor.com/>`_, along with test coverage
  checks with `Coveralls <https://coveralls.io/>`_.

* Decided upon code and docstring style conventions and set up
  automated code style checks with `pep8speaks
  <https://pep8speaks.com/>`_.

* Developed `online documentation for PlasmaPy
  <http://docs.plasmapy.org>`_ that is hosted by `Read the Docs
  <https://readthedocs.org/>`_.

  - Automated documentation builds with `Sphinx
    <http://www.sphinx-doc.org/>`_.

  - Wrote narrative documentation for each subpackage.

* Adopted use of `~astropy.units` as a units package.

* Created the `~plasmapy.atomic` subpackage to provide easy access to
  commonly used atomic data.

  - Created a functional interface to access particle properties and
    find the energy released from nuclear reactions.

  - Created the `~plasmapy.atomic.Particle` class as an object-oriented
    interface to the `~plasmapy.atomic` subpackage.

  - Created the `~plasmapy.atomic.particle_input` decorator.

* Created the `~plasmapy.classes` subpackage that includes the prototype
  `~plasmapy.classes.Plasma3D`, `~plasmapy.classes.PlasmaBlob`, and
  `~plasmapy.classes.Species` classes.

* Created the `~plasmapy.constants` subpackage.

* Created the `~plasmapy.mathematics` subpackage that contains
  analytical functions commonly used in plasma physics.

* Created the `~plasmapy.physics` subpackage with its
  `~plasmapy.physics.transport` module to calculate plasma parameters,
  transport coefficients, dielectric tensor elements, collision rates,
  and relativity/quantum physics parameters used in plasma physics.

* Created the `~plasmapy.utils` subpackage.

  - Created `~plasmapy.utils.check_quantity` and
    `~plasmapy.utils.check_relativistic` decorators.

  - Created custom exceptions.

  - Added import helper and test helper functionality.

* Began development of the `~plasmapy.diagnostics` subpackage.

  - Created a module to interpret Langmuir probe data.

* Created a repository for `PlasmaPy Enhancement Proposals
  <https://github.com/PlasmaPy/PlasmaPy-PLEPs>`_.

* Began using `type hint annotations
  <https://docs.python.org/3/library/typing.html>`_.

* Set up architecture to incorporate `Cython <http://cython.org/>`_ into
  performance-critical sections of code.

* Incorporated import and setup tools from the `~astropy_helpers`
  package.

* Set up a page describing the :ref:`subpackage-stability`.

.. _change-log-0.1.0-api:

Changes to API
~~~~~~~~~~~~~~

- PlasmaPy now has an API.

.. _change-log-0.1.0-bugfix:

Bug Fixes
~~~~~~~~~

- Fixed bug in universe that cause solar neutrinos to oscillate
  between different flavors.

.. I went to a talk on neutrinos once, but it all just went in one ear
   and out the other.
