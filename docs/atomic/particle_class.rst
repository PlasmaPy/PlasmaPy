=================
Atomic Subpackage
=================

The physics of plasmas depends heavily on the properties of the
particles that make up the plasma.  The `~plasmapy.atomic` subpackage
allows us to access necessary particle information.

--------------
Particle Class
--------------

The `~plasmapy.atomic.Particle` class provides an object-oriented
interface to particle information. The simplest way to create an
instance of the `~plasmapy.atomic.Particle` class is to pass it a `str`
representing a particle.

>>> from plasmapy.atomic import Particle
>>> electron = Particle('e-')

The `~plasmapy.atomic.Particle` class accepts a variety of different
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
`plasmapy.atomic.Particle` to represent an atomic number.  For isotopes
and ions, the mass number may be represented with the `mass_numb`
keyword and the integer charge may be represented with the `Z` keyword.

>>> proton = Particle(1, mass_numb=1, Z=1)

The properties of each particle may be accessed using attributes of the
`~plasmapy.atomic.Particle` instance.

>>> proton.atomic_number
1
>>> electron.integer_charge
1
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

Strings representing particles may be accessed using the `particle`,
`element`, `isotope`, and `ion` attributes.

>>> antimuon.particle
'mu-'
>>> triton.element
'H'
>>> alpha.isotope
'He-4'
>>> deuteron.ion
'D 1+'

Equality between particles may be tested either between two
`~plasmapy.atomic.Particle` instances, or between a
`~plasmapy.atomic.Particle` instance and a `str`.

>>> Particle('H-1') == Particle('protium 1+')
True
>>> alpha == 'He-4 1+'
False

The `is_electron` attribute provides a quick way to check whether or not
a particle is an electron.

>>> electron.is_electron
True
>>> hydride.is_electron
False

The `element`, `isotope`, and `ion` return `None` when the particle is
not the respective category.  Because non-empty strings evaluate to
`True` and `None` evaluates to `False` when converted to a `bool`, these
attributes may be used in conditional statements to test whether or not
a particle is in one of these categories.

.. code-block:: python

    particles = [Particle('e-'), Particle('Fe-56'), Particle('alpha')]

    for particle in particles:
        if particle.element:
            print(f"{particle} corresponds to element {particle.element}")
        if particle.isotope:
            print(f"{particle} corresponds to isotope {particle.isotope}")
        if particle.ion:
            print(f"{particle} corresponds to ion {particle.ion}")


The `categories` attribute returns a `set` with the classification
categories corresponding to the particle.

>>> sorted(electron.categories)
['charged', 'electron', 'fermion', 'lepton', 'matter', 'stable']

Membership of a particle within a category may be checked using the
`~plasmapy.atomic.Particle.is_category` method.

>>> alpha.is_category('lepton')
False
>>> electron.is_category('fermion', 'lepton', 'charged')
True
>>> iron56.is_category(['element', 'isotope'])
True

The `require` keyword specifies categories that a particle must
belong to in order for `is_category` to return `True`.

>>> deuteron.is_category(require={'element', 'isotope', 'ion'})
True

The particle must belong to at least one of the categories specified
with the `any_of` keyword

The `any_of` keyword specifies categories of which the particle must
belong to at least one in order for `is_category` to return `True`.

>>> Fe56.is_category(any_of=['charged', 'uncharged'])
False

The particle

The `any_of` keyword specifies categories of which the particle must
belong to

Calling the `is_category` method with no arguments returns a set
containing all of the valid categories for any particle.

>>> sorted(proton.is_category())  # all valid categories
['actinide',
 'alkali metal',
 'alkaline earth metal',
 'antibaryon',
 'antilepton',
 'antimatter',
 'antineutrino',
 'baryon',
 'boson',
 'charged',
 'electron',
 'element',
 'fermion',
 'halogen',
 'ion',
 'isotope',
 'lanthanide',
 'lepton',
 'matter',
 'metal',
 'metalloid',
 'neutrino',
 'neutron',
 'noble gas',
 'nonmetal',
 'positron',
 'post-transition metal',
 'proton',
 'stable',
 'transition metal',
 'uncharged',
 'unstable']
