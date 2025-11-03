.. _policies:

**********************
Policies and Practices
**********************

.. contents:: Contents
   :local:

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
bug fixes or includes important new features.

.. note::

   Dependency groups used primarily for code development activities such
   as running tests and building documentation are excluded from this
   policy because they are not intended for use by end users.

.. _security-policy:

Security policy
===============

PlasmaPy's `security policy`_ is maintained within its GitHub repository.
Please use this link to `privately report a security vulnerability`_.

.. _privately report a security vulnerability: https://github.com/plasmapy/plasmapy/security/advisories/new
.. _security policy: https://github.com/PlasmaPy/PlasmaPy?tab=security-ov-file
