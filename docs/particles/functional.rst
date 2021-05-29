.. _particles-functions:

Functions
*********

In addition to the |Particle| class, the `plasmapy.particles`
subpackage has a functional interface.

.. _particles-func-symbols:

Symbols and names
=================

Several functions in `plasmapy.particles` return string representations
of particles, including |atomic_symbol|, |isotope_symbol|, |ionic_symbol|,
and |element_name|.

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

The full symbol of the particle can be found using |particle_symbol|.

>>> particle_symbol('electron')
'e-'

.. _particles-func-properties:

Particle properties
===================

The |atomic_number| and |mass_number| functions are analogous to the
corresponding attributes in the |Particle| class.

>>> atomic_number('iron')
26
>>> mass_number('T+')
3

Charge information may be found using |charge_number| and
|electric_charge|.

>>> charge_number('H-')
-1
>>> electric_charge('muon antineutrino')
<Quantity 0. C>

These functions will raise a |ChargeError| for elements and isotopes
that lack explicit charge information.

>>> electric_charge('H')
Traceback (most recent call last):
  ...
plasmapy.particles.exceptions.ChargeError: Charge information is required for electric_charge.

The standard atomic weight for the terrestrial environment may be
accessed using |standard_atomic_weight|.

>>> standard_atomic_weight('Pb').to('u')
<Quantity 207.2 u>

The mass of a particle may be accessed through the |particle_mass|
function.

>>> particle_mass('deuteron')
<Quantity 3.34358372e-27 kg>

.. _particles-isotopes:

Isotopes
========

The relative isotopic abundance of each isotope in the terrestrial
environment may be found using |isotopic_abundance|.

>>> isotopic_abundance('H-1')
0.999885
>>> isotopic_abundance('D')
0.000115

A list of all discovered isotopes in order of increasing mass number
can be found with |known_isotopes|.

>>> known_isotopes('H')
['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7']

The isotopes of an element with a non-zero isotopic abundance may be
found with |common_isotopes|.

>>> common_isotopes('Fe')
['Fe-56', 'Fe-54', 'Fe-57', 'Fe-58']

All stable isotopes of an element may be found with |stable_isotopes|.

>>> stable_isotopes('Pb')
['Pb-204', 'Pb-206', 'Pb-207', 'Pb-208']

.. _particles-func-stability:

Stability
=========

The |is_stable| function returns `True` for stable particles and
`False` for unstable particles.

>>> is_stable('e-')
True
>>> is_stable('T')
False

The |half_life| function returns the particle's half-life as a
|Quantity| in units of seconds, if known.

>>> half_life('n')
<Quantity 881.5 s>

For stable particles (or particles that have not been discovered to be
unstable), |half_life| returns |inf| seconds.

>>> half_life('p+')
<Quantity inf s>

If the particle's half-life is not known to sufficient precision, then
|half_life| returns a `str` with the estimated value while issuing a
|MissingParticleDataWarning|.

Reduced mass
============

The |reduced_mass| function is useful in cases of two-body collisions.

>>> reduced_mass('e-', 'p+')
<Quantity 9.10442514e-31 kg>
>>> reduced_mass('D+', 'T+')
<Quantity 2.00486597e-27 kg>

.. |atomic_number| replace:: :func:`~plasmapy.particles.atomic.atomic_number`
.. |atomic_symbol| replace:: :func:`~plasmapy.particles.atomic.atomic_symbol`
.. |element_name| replace:: :func:`~plasmapy.particles.atomic.element_name`
.. |half_life| replace:: :func:`~plasmapy.particles.atomic.half_life`
.. |ionic_symbol| replace:: :func:`~plasmapy.particles.atomic.is_stable`
.. |is_stable| replace:: :func:`~plasmapy.particles.atomic.is_stable`
.. |isotope_symbol| replace:: :func:`~plasmapy.particles.atomic.isotope_symbol`
.. |isotopic_abundance| replace:: :func:`~plasmapy.particles.atomic.isotopic_abundance`
.. |mass_number| replace:: :func:`~plasmapy.particles.atomic.mass_number`
.. |charge_number| replace:: :func:`~plasmapy.particles.atomic.charge_number`
.. |electric_charge| replace:: :func:`~plasmapy.particles.atomic.electric_charge`
.. |standard_atomic_weight| replace:: :func:`~plasmapy.particles.atomic.standard_atomic_weight`
.. |particle_mass| replace:: :func:`~plasmapy.particles.atomic.particle_mass`
.. |particle_symbol| replace:: :func:`~plasmapy.particles.atomic.particle_symbol`
.. |known_isotopes| replace:: :func:`~plasmapy.particles.atomic.known_isotopes`
.. |common_isotopes| replace:: :func:`~plasmapy.particles.atomic.common_isotopes`
.. |reduced_mass| replace:: :func:`~plasmapy.particles.atomic.reduced_mass`
.. |stable_isotopes| replace:: :func:`~plasmapy.particles.atomic.stable_isotopes`
