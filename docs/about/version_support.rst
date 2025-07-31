.. _policies:

********
Policies
********

.. _version-support:

Version support policy for dependencies
=======================================

.. image:: https://img.shields.io/badge/SPEC-0-green?labelColor=%23004811&color=%235CA038
   :target: https://scientific-python.org/specs/spec-0000/

PlasmaPy generally follows the time-based policy for dropping
dependencies as recommended in |SPEC 0|. SPEC 0 recommends a common
dependency support window across the scientific Python ecosystem, while
limiting the maintenance burden associated with supporting older
versions of dependencies. The specific recommendations are that:

1. Support for Python versions be dropped **3 years** after their
   initial release.
2. Support for core package dependencies be dropped **2 years** after
   their initial release.

Support for dependencies may be dropped sooner than what SPEC 0
recommends if a more recent release of a dependency contains critical
bug fixes are important new features.

.. note::

   Dependency groups used primarily for code development activities such
   as running tests and building documentation are excluded from this
   policy because they are not intended for use by end users.
