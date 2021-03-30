.. _particles-functions:

Functions
*********

In addition to the `~plasmapy.particles.Particle` class, the
`~plasmapy.particles` subpackage has a functional interface.

.. _particles-func-symbols:

Symbols and names
=================

Several functions in `~plasmapy.particles` return string representations
of particles, including `~plasmapy.particles.atomic_symbol`,
`~plasmapy.particles.isotope_symbol`, `~plasmapy.particles.ionic_symbol`,
and `~plasmapy.particles.element_name`.

>>> from plasmapy.particles import *
>>> atomic_symbol('alpha')
'He'
>>> isotope_symbol('alpha')
'He-4'
>>> ionic_symbol('alpha')
'He-4 2+'
>>> particle_symbol('alpha')
'He-4 2+'
>>> element_name('alpha')
'helium'

The full symbol of the particle can be found using
`~plasmapy.particles.particle_symbol`.

>>> particle_symbol('electron')
'e-'

.. _particles-func-properties:

Particle properties
===================

The `~plasmapy.particles.atomic_number` and `~plasmapy.particles.mass_number`
functions are analogous to the corresponding attributes in the
`~plasmapy.particles.Particle` class.

>>> atomic_number('iron')
26
>>> mass_number('T+')
3

Charge information may be found using `~plasmapy.particles.integer_charge`
and `~plasmapy.particles.electric_charge`.

>>> integer_charge('H-')
-1
>>> electric_charge('muon antineutrino')
<Quantity 0. C>

These functions will raise a `~plasmapy.particles.exceptions.ChargeError` for
elements and isotopes that lack explicit charge information.

>>> electric_charge('H')
Traceback (most recent call last):
  ...
plasmapy.particles.exceptions.ChargeError: Charge information is required for electric_charge.

The standard atomic weight for the terrestrial environment may be
accessed using `~plasmapy.particles.standard_atomic_weight`.

>>> standard_atomic_weight('Pb').to('u')
<Quantity 207.2 u>

The mass of a particle may be accessed through the
`~plasmapy.particles.particle_mass` function.

>>> particle_mass('deuteron')
<Quantity 3.34358372e-27 kg>

.. _particles-isotopes:

Isotopes
========

The relative isotopic abundance of each isotope in the terrestrial
environment may be found using `~plasmapy.particles.isotopic_abundance`.

>>> isotopic_abundance('H-1')
0.999885
>>> isotopic_abundance('D')
0.000115

A list of all discovered isotopes in order of increasing mass number
can be found with `~plasmapy.particles.known_isotopes`.

>>> known_isotopes('H')
['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7']

The isotopes of an element with a non-zero isotopic abundance may be
found with `~plasmapy.particles.common_isotopes`.

>>> common_isotopes('Fe')
['Fe-56', 'Fe-54', 'Fe-57', 'Fe-58']

All stable isotopes of an element may be found with
`~plasmapy.particles.stable_isotopes`.

>>> stable_isotopes('Pb')
['Pb-204', 'Pb-206', 'Pb-207', 'Pb-208']

.. _particles-func-stability:

Stability
=========

The `~plasmapy.particles.is_stable` function returns `True` for stable
particles and `False` for unstable particles.

>>> is_stable('e-')
True
>>> is_stable('T')
False

The `~plasmapy.particles.half_life` function returns the particle's
half-life as a `~astropy.units.Quantity` in units of seconds, if known.

>>> half_life('n')
<Quantity 881.5 s>

For stable particles (or particles that have not been discovered to be
unstable), `~plasmapy.particles.half_life` returns `~numpy.inf` seconds.

>>> half_life('p+')
<Quantity inf s>

If the particle's half-life is not known to sufficient precision, then
`~plasmapy.particles.half_life` returns a `str` with the estimated value
while issuing a `~plasmapy.particles.exceptions.MissingParticleDataWarning`.

Reduced mass
============

The `~plasmapy.particles.reduced_mass` function is useful in cases of
two-body collisions.

>>> reduced_mass('e-', 'p+')
<Quantity 9.10442514e-31 kg>
>>> reduced_mass('D+', 'T+')
<Quantity 2.00486597e-27 kg>
