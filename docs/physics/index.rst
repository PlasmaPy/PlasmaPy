.. py:module:: physics

.. _plasmapy-physics:

****************************************
The Physics Package (`plasmapy.physics`)
****************************************

.. currentmodule:: plasmapy.physics

Introduction
============

`plasmapy.physics` provides theoretical formulas for calculation of physical quantities helpful for plasma physics.
The layout of the subpackage is still in flux, but for now we have settled on providinga single `plasmapy.physics` namespace
for some of the most common functions. The actual functions are located in modules, subjectively grouped by topic.

We thus have:

* `plasmapy.physics.parameters` for plasma parameteres such as the plasma frequency or Debye length
* `plasmapy.physics.dielectric` deals with tensor dielectric functions
* `plasmapy.physics.distribution` handles distribution functions commonly encountered in plasma physics, such as Maxwellian and Kappa
* `plasmapy.physics.quantum` contains functionality for degenerate, cold and dense plasmas for which quantum effects are important
* `plasmapy.physics.relativity` provides ample room for the Lorentz factor.

We also have a `plasmapy.physics.transport` subpackage, as transport and collision theory turns out to be 

The subpackage makes heavy use of `astropy.units.Quantity` for handling conversions between different unit systems.
This is especially important for electron volts, commonly used in plasma physics to denote temperature, although
it is technically a unit of energy.

Most functions expect `astropy.units.Quantity` as input, however some will use the `plasmapy.utils.check_quantity` decorator
to automatically cast arguments to Quantities. If that happens, you will be notified via an `astropy.units.UnitsWarning`.

Using `plasmapy.physics`
========================


.. toctree::
   :maxdepth: 2

   transport/index

Reference/API
=============

.. automodapi:: plasmapy.physics
   :no-inheritance-diagram:
