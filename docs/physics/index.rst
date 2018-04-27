.. py:module:: plasmapy.physics

.. _plasmapy-physics:

********************************************
Plasma physics formulas (`plasmapy.physics`)
********************************************

.. currentmodule:: plasmapy.physics

Introduction
============

`plasmapy.physics` provides theoretical formulas for calculation of physical quantities helpful for plasma physics.
The layout of the subpackage is still in flux.

The subpackage makes heavy use of `astropy.units.Quantity` for handling conversions between different unit systems.
This is especially important for electron volts, commonly used in plasma physics to denote temperature, although
it is technically a unit of energy.

Most functions expect `astropy.units.Quantity` as input, however some will use the `plasmapy.utils.check_quantity` decorator
to automatically cast arguments to Quantities. If that happens, you will be notified via an `astropy.units.UnitsWarning`.

For a general overview of how unit-based input works, take a look at the following example:

.. topic:: Examples:

   * :ref:`sphx_glr_auto_examples_plot_physics.py`
   * :ref:`sphx_glr_auto_examples_plot_cold_plasma_tensor_elements.py`
   * :ref:`sphx_glr_auto_examples_plot_distribution.py`

Reference/API
=============

.. automodapi:: plasmapy.physics.parameters
.. automodapi:: plasmapy.physics.dielectric
.. automodapi:: plasmapy.physics.dimensionless
.. automodapi:: plasmapy.physics.distribution
.. automodapi:: plasmapy.physics.quantum
.. automodapi:: plasmapy.physics.relativity

Notes for developers
====================
Values should be returned as an Astropy Quantity in SI units.

If a quantity has several names, then the function name should be the
one that provides the most physical insight into what the quantity
represents.  For example, 'gyrofrequency' indicates gyration, while
Larmor frequency indicates that this frequency is somehow related to a
human (or perhaps a cat?) named Larmor.  Similarly, using omega_ce as
a function name for this quantity will make the code less readable to
people who are unfamiliar with the notation or use a different symbol.

The docstrings for plasma parameter methods should describe the
physics associated with these quantities in ways that are
understandable to students who are taking their first course in plasma
physics while still being useful to experienced plasma physicists.

