.. py:module:: physics

.. _plasmapy-physics:

*****************************************
Plasma parameters (`plasmapy.parameters`)
*****************************************

.. currentmodule:: plasmapy.physics

Introduction
============

`plasmapy.physics` provides theoretical formulas for calculation of physical quantities helpful for plasma physics.
The layout of the subpackage is still in flux, but for now we have settled on providing a single `plasmapy.physics` namespace
for some of the most common functions. The actual functions are located in modules, subjectively grouped by topic,
with sub-namespaces for larger sub-packages.

We thus have:

* `plasmapy.physics.parameters` for plasma parameters such as the plasma frequency or Debye length
* `plasmapy.physics.dielectric` deals with tensor dielectric functions
* `plasmapy.physics.distribution` handles distribution functions commonly encountered in plasma physics, such as Maxwellian and Kappa
* `plasmapy.physics.quantum` contains functionality for degenerate, cold and dense plasmas for which quantum effects are important
* `plasmapy.physics.relativity` provides ample room for the Lorentz factor.

We also have a `plasmapy.physics.transport` subpackage, as transport and collision theory turn out to be large enough and sufficiently 
self-contained to warrant their own modules.

The subpackage makes heavy use of `astropy.units.Quantity` for handling conversions between different unit systems.
This is especially important for electron volts, commonly used in plasma physics to denote temperature, although
it is technically a unit of energy.

Most functions expect `astropy.units.Quantity` as input, however some will use the `plasmapy.utils.check_quantity` decorator
to automatically cast arguments to Quantities. If that happens, you will be notified via an `astropy.units.UnitsWarning`.

.. topic:: Examples:

   * :ref:`sphx_glr_auto_examples_plot_physics.py`
   * :ref:`sphx_glr_auto_examples_plot_distribution.py`
   * :ref:`sphx_glr_auto_examples_plot_cold_plasma_tensor_elements.py`

Reference/API
=============

.. automodapi:: plasmapy.physics
   :no-inheritance-diagram:
   :no-heading:
   :no-main-docstr:

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

