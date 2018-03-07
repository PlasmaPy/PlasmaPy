---------
Functions
---------

In addition to the `~plasmapy.atomic.Particle` class, the
`~plasmapy.atomic` subpackage has a functional interface.

Symbols and Names
-----------------

Several functions in `~plasmapy.atomic` provide string representations
of particles:

Various string representations of elements, isotopes, ions, and other
particles are available in the

The `~plasmapy.atomic.atomic_symbol`

>>> from plasmapy.atomic import *
>>> atomic_symbol('alpha')
'He'
>>> isotope_symbol('alpha')
'He-4'
>>> ion_symbol('alpha')
'He-4 2+'
>>> particle_symbol('alpha')
'He-4 2+
>>> element_name('alpha')
'helium'
>>> particle_symbol('electron')
'e-'

Particle Properties
-------------------

>>> atomic_number('iron')
26
>>> mass_number('T+')
3

>>> integer_charge('H-')
-1

>>> electric_charge('muon antineutrino')
<Quantity 0. C>

>>> standard_atomic_weight('Pb').to('u')
<Quantity 207.2 u>

>>> particle_mass('deuteron')
<Quantity 3.34358372e-27 kg>

Isotopes
--------

>>> isotopic_abundance('H-1')
0.999885
>>> isotopic_abundance('D')
0.000115

>>> known_isotopes('H')
['H-1', 'D', 'T', 'H-4', 'H-5', 'H-6', 'H-7']

>>> common_isotopes('Fe')
['Fe-56', 'Fe-54', 'Fe-57', 'Fe-58']

>>> stable_isotopes('Pb')
['Pb-204', 'Pb-206', 'Pb-207', 'Pb-208']

Stability
---------

>>> is_stable('e-')
True
>>> is_stable('T')
False

The `~plasmapy.atomic.half_life` function returns the particle's
half-life in seconds, if known.

>>> half_life('n')
<Quantity 881.5 s>

For stable particles (or particles not experimentally determined to be
unstable), `~plasmapy.atomic.half_life` returns infinity seconds.

>>> half_life('p+')
<Quantity inf s>

If the particle's half-life is not known to sufficient precision, then
`~plasmapy.atomic.half_life` returns a `str` with the estimated value
while issuing a `~plasmapy.utils.MissingAtomicDataWarning`.

Reduced Mass
------------

The reduced mass is

>>> reduced_mass('e-', 'p+')
<Quantity 9.10442514e-31 kg>

>>> reduced_mass('D+', 'T+')
<Quantity 2.00486597e-27 kg>
