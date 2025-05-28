.. _performance-tips:

****************
Performance Tips
****************

Most of the time, readability is more important than the performance of
scientific software. This page contains tips for improving the
performance of PlasmaPy for situations where performance becomes a
bottleneck 🐢.

Python versions
===============

Upgrade to the newest version of Python to take advantage of ongoing
performance improvements from the Faster CPython project. New versions
of Python also have improved error message, which can speed up the
debugging process too.

A new version of Python is released in October of each year, and can be
used with PlasmaPy a few months later.

Using Astropy units
===================

PlasmaPy makes heavy use of `astropy.units` and |Quantity| operations.
See Astropy's documentation for `performance tips for Quantity
operations`_.

Because PlasmaPy uses SI units internally, performance can be improved
slightly by providing |Quantity| objects in SI units to functions in
`plasmapy.formulary`. Unit conversions done by the |validate_quantities|
decorator then do not need to be performed.

Particles
=========

Many of the functions in `plasmapy.formulary` accept |particle-like|
|arguments|. Arguments that are not already a |Particle|,
|CustomParticle|, or |ParticleList| are converted into one (usually via
the |particle_input| decorator. When a formulary function is repeatedly
called, performance can be improved by creating the particle object
ahead of time.

For example, suppose we are calculating the
`~plasmapy.formulary.frequencies.gyrofrequency` of a proton. If we
represent the particle as a string, then the function will need to
create a |Particle| each time the function is called.

.. code-block:: python

   from plasmapy.formulary import gyrofrequency
   import astropy.units as u

   for i in range(1000):
       gyrofrequency(B=0.02 * u.T, particle="p+")  # create the Particle repeatedly

If we create the |Particle| and |Quantity| ahead of time, they will need
to be create only once instead of repeatedly.

.. code-block:: python

   from plasmapy.particles import Particle

   B = 0.02 * u.T
   proton = Particle("p+")  # create the Particle once

   for i in range(1000):
       gyrofrequency(B=B, particle=proton)

Lite-functions
==============

PlasmaPy includes :term:`lite-functions` for some `plasmapy.formulary`
functions for situations when performance matters. For example,
`~plasmapy.formulary.frequencies.plasma_frequency_lite` is the
lite-function for `~plasmapy.formulary.frequencies.plasma_frequency`.

Lite-functions accept and return NumPy arrays (assumed to be in SI
units) instead of |Quantity| objects. Lite-functions have improved
performance because they do not perform any checks on the inputs, and
thus should be used cautiously.

If you need a lite-function version of a `plasmapy.formulary` function
that has not already been implemented, please `raise an issue`_.

.. _performance tips for Quantity operations: https://docs.astropy.org/en/stable/units/index.html#astropy-units-performance
.. _raise an issue: https://github.com/PlasmaPy/PlasmaPy/issues/new
