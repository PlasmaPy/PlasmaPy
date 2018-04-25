.. _plasmapy-constants:

********************************
Constants (`plasmapy.constants`)
********************************

The `~plasmapy.constants` module contains many physical and mathematical
constants that are commonly used in the plasma sciences.

Mathematical constants (such as `pi`) are provided as `float` objects.

Physical constants are imported directly from `~astropy.constants` as
instances of the `~astropy.constants.Constant` class (which is a
subclass of `~astropy.units.Quantity`).  We recommend reviewing the
`documentation for Astropy's constants subpackage
<http://docs.astropy.org/en/stable/constants/index.html>`_ for more
details (including Astropy's strategy for versioning constants).
Physical constants by default are provided in SI units.

.. warning::
    The values of constants that are not defined exactly will change
    slightly as improved measurements become available.

Reference/API
=============

.. automodapi:: astropy.constants
