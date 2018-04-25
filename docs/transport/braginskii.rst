.. _braginskii
.. py:module:: physics.transport.braginskii

Braginskii
**********

The `~plasmapy.physics.transport.braginskii` subpackage includes a 
class interface 
`~plasmapy.physics.transport.braginskii.ClassicalTransport` with 
multiple methods which are involved in modeling classical transport.

.. topic:: Examples:

   * :ref:`sphx_glr_auto_examples_plot_braginskii.py`

.. _Classical transport models

Classical transport models
==========================

Broad overview of classical transport models implemented within 
`braginskii.py`. These include:
Braginskii
Spitzer-Harm
Epperlein-Haines (under construction)
Ji-Held

.. _Classical-transport-class

Classical transport class
=========================

Explain how to use the 
`~plasmapy.physics.transport.braginskii.ClassicalTransport` interface.

.. _Conductivities-resistivities

Conductivities and Resistivities
================================

Explain:
resistivity
thermoelectric_conductivity
ion_thermal_conductivity
electron_thermal_conductivity

.. _Viscosities

Viscosities
======================

These include:
ion_viscosity
electron_viscosity

Reference/API
=============

.. automodapi:: plasmapy.physics.transport.braginskii
    :skip: Coulomb_logarithm
    :skip: Hall_parameter
    :skip: collision_rate_electron_ion
    :skip: collision_rate_ion_ion
    :no-heading:
    :no-main-docstr:
