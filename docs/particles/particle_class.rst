.. _particle-class:

Particle objects
****************

PlasmaPy contains several classes to represent particles, including
|Particle|, |CustomParticle|, |ParticleList|, and
|DimensionlessParticle|.

.. _particle-class-instantiation:

Particles
=========

To create a |Particle| object, pass it a |particle-like| string that
represents a particle.

>>> from plasmapy.particles import Particle
>>> electron = Particle('e-')

The |Particle| class accepts a variety of different `str` formats to
represent particles. Atomic symbols are case-sensitive, but element
names and many aliases are not.

>>> alpha = Particle('alpha')
>>> deuteron = Particle('D+')
>>> triton = Particle('tritium 1+')
>>> iron56 = Particle('Fe-56')
>>> helium = Particle('helium')
>>> muon = Particle('mu-')
>>> antimuon = Particle('antimuon')
>>> hydride = Particle('H-')

An `int` may be used as the first positional argument to |Particle| to
represent an atomic number. For isotopes and ions, the mass number may
be represented with the ``mass_numb`` keyword and the |charge number|
may be represented with the ``Z`` keyword.

>>> proton = Particle(1, mass_numb=1, Z=1)

The most frequently used |Particle| objects may be imported directly
from `plasmapy.particles`.

>>> from plasmapy.particles import proton, electron

The |Particle| objects that may be imported directly are:
`~plasmapy.particles.proton`, `~plasmapy.particles.electron`,
`~plasmapy.particles.neutron`, `~plasmapy.particles.positron`,
`~plasmapy.particles.deuteron`, `~plasmapy.particles.triton`, and
`~plasmapy.particles.alpha`.

.. _particle-class-properties:

Accessing particle properties
-----------------------------

The properties of each particle may be accessed using the attributes of
the corresponding |Particle| object.

>>> proton.atomic_number
1
>>> electron.charge_number
-1
>>> triton.mass_number
3

These properties are often returned as a |Quantity| in SI units.

>>> alpha.charge
<Quantity 3.20435324e-19 C>
>>> deuteron.mass
<Quantity 3.34358372e-27 kg>
>>> triton.half_life
<Quantity 3.888e+08 s>
>>> iron56.binding_energy.to('GeV')
<Quantity 0.49225958 GeV>
>>> hydrogen.ionization_energy
<Quantity 2.17870942e-18 J>

Strings representing particles may be accessed using the
`~plasmapy.particles.particle_class.Particle.symbol`,
`~plasmapy.particles.particle_class.Particle.element`,
`~plasmapy.particles.particle_class.Particle.isotope`, and
`~plasmapy.particles.particle_class.Particle.ionic_symbol` attributes.

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
----------

The `~plasmapy.particles.particle_class.Particle.categories` attribute
returns a `set` with the classification categories corresponding to the
particle.

>>> sorted(electron.categories)
['charged', 'electron', 'fermion', 'lepton', 'matter', 'stable']

Membership of a particle within a category may be checked using
|is_category|.

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

Valid particle categories are listed in the docstring for |is_category|.

.. _particle-class-conditionals-equality:

Conditionals and equality properties
------------------------------------

Equality between particles may be tested either between two |Particle|
objects, or between a |Particle| object and a `str`.

>>> Particle('H-1') == Particle('protium 1+')
False
>>> alpha == 'He-4 2+'
True

The `~plasmapy.particles.particle_class.Particle.is_electron` and
`~plasmapy.particles.particle_class.Particle.is_ion` attributes
provide a quick way to check whether or not a particle is an electron or
ion, respectively.

>>> electron.is_electron
True
>>> hydride.is_electron
False
>>> deuteron.is_ion
True
.. _particle-class-antiparticles:

Antiparticles
-------------

The antiparticle of an elementary particle or antiparticle may be found
by either using Python's unary invert operator (``~``) or the
`~plasmapy.particles.particle_class.Particle.antiparticle` attribute
of a |Particle| object.

>>> ~electron
Particle("e+")
>>> antimuon.antiparticle
Particle("mu-")

Custom particles
================

We can use |CustomParticle| to create particle objects with a mass,
charge, and/or symbol that we provide. The mass and charge must be
|Quantity| objects from `astropy.units`.

>>> import astropy.units as u
>>> from plasmapy.particles import CustomParticle
>>> cp = CustomParticle(mass = 9.3e-26 * u.kg, charge = 1.5e-18 * u.C, symbol = "Fe 9.5+")

|CustomParticle| has many of the same attributes and methods as
|Particle|, and can often be used interchangeably.

>>> cp.charge
<Quantity 1.52e-18 C>
>>> cp.mass
<Quantity 9.3e-26 kg>
>>> cp.symbol
'Fe 9.5+'

If the charge and/or mass is not provided, the attribute will return
|nan| in the appropriate units.

Molecules
---------

We can use `~plasmapy.particles.particle_class.molecule` to convert a
chemical symbol into a |CustomParticle| object with the appropriate
mass, charge, and symbol.

>>> from plasmapy.particles import molecule
>>> molecule("CO2 1+")  # carbon dioxide cation
CustomParticle(mass=7.30786637819994e-26 kg, charge=1.602176634e-19 C, symbol=CO2 1+)

Particle lists
==============

|ParticleList| lets us work with multiple particles at once. A
|ParticleList| can contain |Particle| and/or |CustomParticle| objects.

We can create a |ParticleList| by providing it with a
|particle-list-like| object (i.e., a `list` containing |particle-like|
objects). For example, we could provide |ParticleList| with a `list` of
strings that represent individual particles.

>>> from plasmapy.particles import ParticleList
>>> helium_ions = ParticleList(["He-4 0+", "He-4 1+"])

|ParticleList| objects behave similarly to `list` objects, but convert
its contents into the appropriate |Particle| or |CustomParticle|
objects.

>>> helium_ions.append("alpha")
>>> print(helium_ions)
ParticleList(['He-4 0+', 'He-4 1+', 'He-4 2+'])
>>> helium_ions[1]
Particle("He-4 1+")

|ParticleList| shares many of the same attributes as |Particle| and
|CustomParticle|. Attributes of |Particle| and |CustomParticle| that
provide a scalar |Quantity| will provide a |Quantity| array from
|ParticleList|.

>>> helium_ions.charge
<Quantity [0.00000000e+00, 1.60217663e-19, 3.20435327e-19] C>
>>> helium_ions.mass
<Quantity [6.64647907e-27, 6.64556813e-27, 6.64465719e-27] kg>

If we provide a |Quantity| with units of mass or charge, it will get
converted into a |CustomParticle|.

>>> cp_list = ParticleList([1 * u.kg, 1 * u.C])
>>> cp_list[0]
CustomParticle(mass=1.0 kg, charge=nan C)
>>> cp_list.charge
<Quantity [nan,  1.] C>
>>> cp_list.mass
<Quantity [ 1., nan] kg>

We can create a |CustomParticle| with the mean mass and charge of the
particles in a |ParticleList| with its
`~plasmapy.particles.particle_collections.ParticleList.average_particle`
method.

>>> helium_ions.average_particle()
CustomParticle(mass=6.645568133213004e-27 kg, charge=1.602176634e-19 C)

We can create a |ParticleList| by adding |Particle|, |CustomParticle|,
and/or |ParticleList| objects together.

>>> helium_ions + cp + proton
ParticleList(['He-4 0+', 'He-4 1+', 'He-4 2+', 'Fe 9.5+', 'p+'])

As with an individual |Particle| and |CustomParticle|, we can check whether
all the particles in a list fall within a category using |is_category|.

>>> helium_ions.is_category("ion")
False

We may also check each particle in the list individually by setting
the keyword ``particlewise`` to `True`.

>>> helium_ions.is_category("ion", particlewise=True)
[False, True, True]

The machinery contained with |ParticleList| lets us calculate plasma
parameters from `plasmapy.formulary` for multiple particles at once.

>>> from plasmapy.formulary import gyroradius
>>> gyroradius(B = 5 * u.nT, particle=["e-", "p+"], Vperp = 100 * u.km/u.s)
<Quantity [1.13712608e+02, 2.08793710e+05] m>

Dimensionless particles
=======================

We can use |DimensionlessParticle| to represent particles that have been
normalized (i.e., both the mass and charge are dimensionless).

>>> dp = DimensionlessParticle(mass=1, charge=-1)
>>> dp.charge
-1.0
>>> dp.mass
1.0

Because |DimensionlessParticle| objects do not directly represent
physical particles without normalization information, they cannot be
contained within a |ParticleList| or used in `plasmapy.formulary`.

.. |is_category| replace:: `~plasmapy.particles.particle_class.AbstractPhysicalParticle.is_category`
