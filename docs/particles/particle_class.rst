.. _particle-class:

`~plasmapy.particles.Particle` Class
************************************

The `~plasmapy.particles.Particle` class provides an object-oriented
interface to access and represent particle information.

.. _particle-class-instantiation:

Creating a `~plasmapy.particles.Particle` object
================================================

The simplest way to create a `~plasmapy.particles.Particle` object
is to pass it a `str` representing a particle.

>>> from plasmapy.particles import Particle
>>> electron = Particle('e-')

The `~plasmapy.particles.Particle` class accepts a variety of different
`str` formats to represent particles. Atomic symbols are case-sensitive,
but element names and many aliases are not.

>>> alpha = Particle('alpha')
>>> deuteron = Particle('D+')
>>> triton = Particle('tritium 1+')
>>> iron56 = Particle('Fe-56')
>>> helium = Particle('helium')
>>> muon = Particle('mu-')
>>> antimuon = Particle('antimuon')
>>> hydride = Particle('H-')

An `int` may be used as the first positional argument to
`~plasmapy.particles.Particle` to represent an atomic number.  For isotopes
and ions, the mass number may be represented with the ``mass_numb``
keyword and the integer charge may be represented with the ``Z``
keyword.

>>> proton = Particle(1, mass_numb=1, Z=1)

The most frequently used `~plasmapy.particles.Particle` objects may be
imported directly from the atomic subpackage.

>>> from plasmapy.particles import proton, electron

The `~plasmapy.particles.Particle` objects that may be imported
directly are: `~plasmapy.particles.proton`,
`~plasmapy.particles.electron`, `~plasmapy.particles.neutron`,
`~plasmapy.particles.positron`, `~plasmapy.particles.deuteron`,
`~plasmapy.particles.triton`, and `~plasmapy.particles.alpha`.

.. _particle-class-properties:

Accessing particle properties
=============================

The properties of each particle may be accessed using the attributes of
the corresponding `~plasmapy.particles.Particle` object.

>>> proton.atomic_number
1
>>> electron.integer_charge
-1
>>> triton.mass_number
3

Some of these properties are returned as a `~astropy.units.Quantity` in
SI units.

>>> alpha.charge
<Quantity 3.20435324e-19 C>
>>> deuteron.mass
<Quantity 3.34358372e-27 kg>
>>> triton.half_life
<Quantity 3.888e+08 s>
>>> iron56.binding_energy.to('GeV')
<Quantity 0.49225958 GeV>

Strings representing particles may be accessed using the
`~plasmapy.particles.Particle.symbol`,
`~plasmapy.particles.Particle.element`,
`~plasmapy.particles.Particle.isotope`, and
`~plasmapy.particles.Particle.ionic_symbol` attributes.

>>> antimuon.symbol
'mu+'
>>> triton.element
'H'
>>> alpha.isotope
'He-4'
>>> deuteron.ionic_symbol
'D 1+'

.. _particle-class-categories:

Categories
==========

The `~plasmapy.particles.Particle.categories` attribute returns a `set`
with the classification categories corresponding to the particle.

>>> sorted(electron.categories)
['charged', 'electron', 'fermion', 'lepton', 'matter', 'stable']

Membership of a particle within a category may be checked using the
`~plasmapy.particles.Particle.is_category` method.

>>> alpha.is_category('lepton')
False
>>> electron.is_category('fermion', 'lepton', 'charged')
True
>>> iron56.is_category(['element', 'isotope'])
True

The particle must be in all of the categories in the ``require``
keyword, at least one of the categories in the ``any_of`` keyword, and
none of the categories in the ``exclude`` in order for it to return
`True`.

>>> deuteron.is_category(require={'element', 'isotope', 'ion'})
True
>>> iron56.is_category(any_of=['charged', 'uncharged'])
False
>>> alpha.is_category(exclude='lepton')
True

Calling the `~plasmapy.particles.Particle.is_category` method with no
arguments returns a set containing all of the valid categories for any
particle.  Valid categories include: ``'actinide'``, ``'alkali metal'``,
``'alkaline earth metal'``, ``'antibaryon'``, ``'antilepton'``,
``'antimatter'``, ``'antineutrino'``, ``'baryon'``, ``'boson'``,
``'charged'``, ``'electron'``, ``'element'``, ``'fermion'``,
``'halogen'``, ``'ion'``, ``'isotope'``, ``'lanthanide'``, ``'lepton'``,
``'matter'``, ``'metal'``, ``'metalloid'``, ``'neutrino'``,
``'neutron'``, ``'noble gas'``, ``'nonmetal'``, ``'positron'``,
``'post-transition metal'``, ``'proton'``, ``'stable'``,
``'transition metal'``, ``'uncharged'``, and ``'unstable'``.

.. _particle-class-conditionals-equality:

Conditionals and equality properties
====================================

Equality between particles may be tested either between two
`~plasmapy.particles.Particle` objects, or between a
`~plasmapy.particles.Particle` object and a `str`.

>>> Particle('H-1') == Particle('protium 1+')
False
>>> alpha == 'He-4 2+'
True

The `~plasmapy.particles.Particle.is_electron` and
`~plasmapy.particles.Particle.is_ion` attributes provide a quick way to
check whether or not a particle is an electron or ion, respectively.

>>> electron.is_electron
True
>>> hydride.is_electron
False
>>> deuteron.is_ion
True

The `~plasmapy.particles.Particle.element` and
`~plasmapy.particles.Particle.isotope` attributes return `None` when the
particle does not correspond to an element or isotope.  Because
non-empty strings evaluate to `True` and `None` evaluates to `False`
when converted to a `bool`, these attributes may be used in conditional
statements to test whether or not a particle is in one of these
categories.

.. code-block:: python

    particles = [Particle('e-'), Particle('Fe-56'), Particle('alpha')]

    for particle in particles:
        if particle.element:
            print(f"{particle} corresponds to element {particle.element}")
        if particle.isotope:
            print(f"{particle} corresponds to isotope {particle.isotope}")

.. _particle-class-antiparticles:

Returning antiparticles
=======================

The antiparticle of an elementary particle or antiparticle may be found
by either using Python's unary invert operator (``~``) or the
`~plasmapy.particles.Particle.antiparticle` attribute of a
`~plasmapy.particles.Particle` object.

>>> ~electron
Particle("e+")
>>> antimuon.antiparticle
Particle("mu-")
