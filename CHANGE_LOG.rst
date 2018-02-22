===================
PlasmaPy Change Log
===================

This document provides a detailed list of changes associated with each
release of PlasmaPy.  The [release
notes](https://github.com/PlasmaPy/PlasmaPy/blob/master/RELEASE_NOTES.md)
contain a narrative description for each version.

Version 0.1.0 (2018-03-xx)
--------------------------

Version 0.1.0 is the initial development release of PlasmaPy.  This
version is a prototype and a preview, and is not feature complete.
Significant changes to the API are expected to occur between versions
0.1.0 and 0.2.0.

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
  AppVeyor.  Set up test coverage checks with Coveralls.  Automated
  code style checks with pep8speaks.  Decided upon code and docstring
  style conventions.

* Set up automated documentation builds with sphinx that are hosted on
  Read the Docs.

* Adopted use of `astropy.units` as a units package.  Allow import of
  units as ``import plasmapy.units`` as u.

* Created the `atomic` subpackage to provide easy access to commonly
  used atomic data.

  - Created ``Particle`` class

  - Created ``@particle_input`` decorator

* Created `classes` subpackage.

  - ``Plasma`` class

* Created `constants` subpackage.

* Created `mathematics` subpackage with functionality commonly used in
  plasma physics.

* Created `physics` subpackage 

  - Plasma parameters

  - Transport parameters

  - Relativity and quantum physics parameters related to plasma
    physics


* Created `utils` subpackage.

  - Created ``@check_quantity`` and ``@check_relativistic``
    decorators.

  - Added custom exceptions, including several for the `atomic`
    subpackage.

* Created basic framework for `diagnostics` subpackage.

* Created a repository for 

* Included astropy-helpers as a submodule.

Changes to API
~~~~~~~~~~~~~~

- PlasmaPy now has an API.

Bug Fixes
~~~~~~~~~

- Fixed bug in universe that cause solar neutrinos to oscillate
  between different flavors.
