.. _ionization-state-data-structures:

Ionization state data structures
********************************

The ionization state (or charge state) of a plasma refers to the
fraction of an element that is at each ionization level.  For example,
the ionization state of a pure helium plasma could be 5% He⁰⁺, 94% He¹⁺,
and 1% He²⁺.

The ionization state of a single element
========================================

We may use the `~plasmapy.particles.IonizationState` class
to represent the ionization state of a single element, such as for this
example.

>>> from plasmapy.particles import IonizationState
>>> ionization_state = IonizationState("He", [0.05, 0.94, 0.01])

The ionization state for helium may be accessed using the
``ionic_fractions`` attribute.  These ionic fractions correspond to the
``integer_charges`` attribute.

>>> ionization_state.ionic_fractions
array([0.05, 0.94, 0.01])
>>> ionization_state.integer_charges
array([0, 1, 2])

The ``Z_mean`` attribute returns the mean integer charge averaged
over all particles in that element.

>>> ionization_state.Z_mean
0.96

The ``Z_rms`` attribute returns the root mean square integer charge.

>>> ionization_state.Z_rms
0.9899...

The ``Z_most_abundant`` attribute returns a `list` of the most abundant
ion(s).  The `list` may contain more than one integer charge in case of
a tie.

>>> ionization_state.Z_most_abundant
[1]

The ``summarize`` method prints out the ionic fraction for the ions with
an abundance of at least 1%.

>>> ionization_state.summarize()
IonizationState instance for He with Z_mean = 0.96
----------------------------------------------------------------
He  0+: 0.050
He  1+: 0.940
He  2+: 0.010
----------------------------------------------------------------

The number density of the element may be specified through the
``n_elem`` keyword argument.

>>> import astropy.units as u
>>> ionization_state = IonizationState(
...     "He", [0.05, 0.94, 0.01], n_elem = 1e19 * u.m ** -3,
... )

The ``n_e`` attribute provides the electron number density as a
`~astropy.units.Quantity`.

>>> ionization_state.n_e
<Quantity 9.6e+18 1 / m3>

The ``number_densities`` attribute provides the number density of each
ion or neutral.

>>> ionization_state.number_densities
<Quantity [5.0e+17, 9.4e+18, 1.0e+17] 1 / m3>

Ionization states for multiple elements
=======================================

The `~plasmapy.particles.IonizationStateCollection` class may be used to
represent the ionization state for multiple elements. This can be used,
for example, to describe the various impurities in a fusion plasma or
the charge state distributions of different elements in the solar wind.

>>> from plasmapy.particles import IonizationStateCollection

The minimal input to `~plasmapy.particles.IonizationStateCollection` is a `list`
of the elements or isotopes to represent.  Integers in the `list` will
be treated as atomic numbers.

>>> states = IonizationStateCollection(["H", 2])

To set the ionic fractions for hydrogen, we may do item assignment.

>>> states["H"] = [0.9, 0.1]

We may use indexing to retrieve an `~plasmapy.particles.IonizationState`
instance for an element.

>>> states["H"]
<IonizationState instance for H>

The ionization states for all of the elements may be specified directly
as arguments to the class.

>>> states = IonizationStateCollection(
...     {"H": [0.01, 0.99], "He": [0.04, 0.95, 0.01]},
...     abundances={"H": 1, "He": 0.08},
...     n0 = 5e19 * u.m ** -3,
... )

The ionic fractions will be stored as a `dict`.

>>> states.ionic_fractions
{'H': array([0.01, 0.99]), 'He': array([0.04, 0.95, 0.01])}

The number density for each element is the product of the number
density scaling factor ``n0`` with that element's abundance.
The number density for each ion is the product of ``n0``, the
corresponding element's abundance, and the ionic fraction.

>>> states.n0
<Quantity 5.e+19 1 / m3>
>>> states.abundances
{'H': 1.0, 'He': 0.08}
>>> states.number_densities["H"]
<Quantity [5.00e+17, 4.95e+19] 1 / m3>

The
corresponding element's abundance, and the ionic fraction.

>>> states.n0
<Quantity 5.e+19 1 / m3>
>>> states.abundances
{'H': 1.0, 'He': 0.08}
>>> states.number_densities["H"]
<Quantity [5.00e+17, 4.95e+19] 1 / m3>

The
corresponding element's abundance, and the ionic fraction.

>>> states.n
<Quantity 5.e+19 1 / m3>
>>> states.abundances
{'H': 1.0, 'He': 0.08}
>>> states.number_densities["H"]
<Quantity [5.00e+17, 4.95e+19] 1 / m3>

The `~plasmapy.particles.IonizationStates.summarize` method may also be
used to get a summary of the ionization states.

>>> states.summarize()
----------------------------------------------------------------
H  1+: 0.990    n_i = 4.95e+19 m**-3
----------------------------------------------------------------
He  0+: 0.040    n_i = 1.60e+17 m**-3
He  1+: 0.950    n_i = 3.80e+18 m**-3
----------------------------------------------------------------
n_e = 5.34e+19 m**-3
T_e = 1.30e+04 K
----------------------------------------------------------------
