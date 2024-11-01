.. py:module:: plasmapy.utils

.. _plasmapy-utils:

************************************
Package utilities (`plasmapy.utils`)
************************************

.. currentmodule:: plasmapy.utils

Introduction
============

The `~plasmapy.utils` subpackage contains functionality that is needed
across multiple subpackages or does not fit nicely in any other
subpackage. Functionality contained in `~plasmapy.utils` includes:

 * Warnings and exceptions used in PlasmaPy, such as
   `~plasmapy.utils.exceptions.RelativityWarning` or
   `~plasmapy.utils.exceptions.PhysicsError`.
 * Decorators we use for reusable physical `~astropy.units.Quantity`
   computation and checking, such as
   `~plasmapy.utils.decorators.validators.validate_quantities`
   and `~plasmapy.utils.decorators.checks.check_relativistic`.
 * Helper utilities for importing and testing packages.
 * Functionality for downloading files from |PlasmaPy's data repository|.

API
===

.. automodapi:: plasmapy.utils.decorators
   :include-heading:

.. automodapi:: plasmapy.utils.exceptions
   :include-heading:

.. automodapi:: plasmapy.utils.code_repr
   :include-heading:

.. automodapi:: plasmapy.utils.calculator
   :include-heading:

.. automodapi:: plasmapy.utils.data
   :include-heading:
