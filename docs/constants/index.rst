.. _plasmapy-constants:

********************************
Constants (`plasmapy.constants`)
********************************

The `~plasmapy.constants` module contains many physical and mathematical
constants that are commonly used in the plasma sciences.

Mathematical constants (such as `pi`) are provided as `float` objects.

Physical constants are imported directly from Astropy's
`~astropy.constants` subpackage as instances of the
`~astropy.constants.Constant` class (which is a subclass of
`~astropy.units.Quantity`).  We recommend reviewing the
`documentation for Astropy's constants subpackage
<http://docs.astropy.org/en/stable/constants/index.html>`_ for more
details (including Astropy's strategy for versioning constants).
Physical constants by default are provided in SI units.

Getting Started
===============

Constants in SI units may be imported directly from PlasmaPy's
`~plasmapy.constants` subpackage.  These objects may be treated very
much like `~astropy.units.Quantity` instances.

A sample use case would be to calculate the fine structure constant
.. math::
    α = \frac{1}{4 \pi \epsilon_0} \frac{e^2}{ħ c}

using `~plasmapy.constants`.

    >>> from plasmapy.constants import eps0, e, hbar, c, pi as π

    >>> α = e ** 2 / (4 * π * eps0 * hbar * c)
    >>> (1 / α).to('')
    <Quantity 137.03599911>

.. warning::
    The values of constants that are not defined exactly will change
    slightly as improved measurements become available.

Reference/API
=============

.. automodapi:: plasmapy.constants
