.. _change-log:

===================
PlasmaPy Change Log
===================

This document provides a detailed list of changes associated with each
release of PlasmaPy.  The release_notes_ contain a narrative description
for each version.

.. _change-log-0.1.0:

Version `0.1.0` (2018-04-27)
--------------------------

Version `0.1.0` is the initial development release of PlasmaPy.  This
version is a prototype and a preview, and is not feature complete.
Significant changes to the API are expected to occur between versions
`0.1.0` and `0.2.0`.

.. _change-log-0.1.0-new:

New Features
~~~~~~~~~~~~

* The vision statement describes the motivation of early PlasmaPy
  developers in creating a fully open source Python package for plasma
  physics.

* Added a code of conduct.

* Added a contribution guide.

* Adopted a BSD 3-clause license and added protections against
  software patents.

* Set up continuous integration testing with Travis CI, CircleCI, and
  AppVeyor, along with test coverage checks with Coveralls.

* Decided upon code and docstring style conventions and set up
  automated code style checks with pep8speaks.

* Set up automated documentation builds with Sphinx that are hosted on
  Read the Docs.

* Adopted use of `~astropy.units` as a units package.  Allowed import of
  units as ``import plasmapy.units as u``.

* Created the `~plasmapy.atomic` subpackage to provide easy access to commonly
  used atomic data.

  - Created functional interface to access particle properties and find
    the energy released from nuclear reactions.

  - Created `~plasmapy.atomic.Particle` class.

  - Created `~particle_input` decorator.

* Created `~plasmapy.classes` subpackage with the
  `~plasmapy.classes.Plasma` and `~plasmapy.classes.Species` classes.

* Created `~plasmapy.constants` subpackage.

* Created `~plasmapy.mathematics` subpackage with functionality commonly
  used in plasma physics.

* Created `~plasmapy.physics` subpackage with a
  `~plasmapy.physics.transport` module to calculate plasma parameters,
  transport coefficients, collision rates, and relativity/quantum
  physics parameters used in plasma physics.

* Created `~plasmapy.utils` subpackage.

  - Created `~plasmapy.utils.check_quantity` and
    `~plasmapy.utils.check_relativistic` decorators.

  - Added custom exceptions, including several for the
    `~plasmapy.atomic` subpackage.

  - Added import helper and test helper functionality.

* Created basic framework for `~plasmapy.diagnostics` subpackage.

  - Added initial support for the analysis of Langmuir probes.

* Created a repository for `PlasmaPy Enhancement Proposals
  <https://github.com/PlasmaPy/PlasmaPy-PLEPs>`_.

* Included astropy-helpers as a submodule.

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