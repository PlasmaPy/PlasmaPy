******
Plasma
******

.. currentmodule:: plasmapy.classes.plasma

Introduction
============

The :class:`Plasma3D` class is a basic structure to contain spatial
information about a plasma.  To initialize a Plasma3D system, first
create an instance of the :class:`Plasma3D` class and then set the
:attr:`Plasma3D.density`, :attr:`Plasma3D.momentum`,
:attr:`Plasma3D.pressure` and the :attr:`Plasma3D.magnetic_field`.

This feature is currently under development.

The :class:`PlasmaBlob` class is a basic structure to contain spatial
information about a plasma.  To initialize a PlasmaBlob system, call
it with arguments: electron temperature :attr:`PlasmaBlob.T_e`, 
and electron density :attr:`PlasmaBlob.n_e`. You may also optionally
define the ionization, :attr:`PlasmaBlob.Z`, and relevant plasma 
particle, :attr:`PlasmaBlob.particle`.

This feature is currently under development.

.. automodapi:: plasmapy.classes.plasma
   :no-inheritance-diagram:
