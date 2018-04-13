.. py:module:: utils

.. _plasmapy-utils:

************************************
The Utils Package (`plasmapy.utils`)
************************************

.. currentmodule:: plasmapy.utils

Introduction
============

`PlasmaPy.utils` is where we store functionality that helps us write (what we
try to think of as) clean, readable and informative code. This means:

 * the many kinds of warnings and exceptions you may (hopefully not!) encounter
   while working with PlasmaPy, such as `plasmapy.utils.RelativityWarning` or
   `plasmapy.utils.PhysicsError`.
 * decorators we use for reusable physical quantity computation and checking,
   such as `plasmapy.utils.check_quantity` and `plasmapy.utils.check_relativistic`
 * Some helper utilities for importing and testing packages such as
   `plasmapy.utils.call_string`

Submodules
=======================

.. toctree::
   :maxdepth: 1

   checks
   exceptions
   pytest_helpers
   import_helpers

