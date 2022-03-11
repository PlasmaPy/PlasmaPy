---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  name: python3
---

[angular frequency]: https://en.wikipedia.org/wiki/Angular_frequency
[Boltzmann constant]: https://en.wikipedia.org/wiki/Boltzmann_constant
[electron-volt]: https://en.wikipedia.org/wiki/Electronvolt
[equivalencies]: https://docs.astropy.org/en/stable/units/equivalencies.html
[frequency]: https://en.wikipedia.org/wiki/Frequency
[lose their units with some operations]: https://docs.astropy.org/en/stable/known_issues.html#quantities-lose-their-units-with-some-operations
[NumPy]: https://numpy.org
[SciPy]: https://scipy.org
[performance tips]: https://docs.astropy.org/en/stable/units/index.html#performance-tips
[physical type]: https://docs.astropy.org/en/stable/units/physical_types.html

# Using Astropy Units

```{code-cell} ipython3
:tags: [remove-cell]
%xmode minimal
```

In scientific computing, we often represent physical quantities as numbers.

```{code-cell} ipython3
distance_in_miles = 50
time_in_hours = 2
velocity_in_mph = distance_in_miles / time_in_hours
print(velocity_in_mph)
```

Representing a physical quantity as a number has risks. We might
unknowingly perform operations with different units, like
`time_in_seconds + time_in_hours`. We might even accidentally perform
operations with physically incompatible units, like `length + time`,
without catching our mistake. We can avoid these problems by using a
units package.

This notebook introduces {py:mod}`astropy.units` with an emphasis on the
functionality needed to work with {py:mod}`plasmapy.particles` and
{py:mod}`plasmapy.formulary`. We typically import this subpackage as `u`.

```{code-cell} ipython3
import astropy.units as u
```

## Contents

1. [Unit basics](#Unit-basics)
2. [Unit operations](#Unit-operations)
3. [Unit conversations](#Unit-conversions)
4. [Detaching units and values](#Detaching-units-and-values)
5. [Equivalencies](#Equivalencies)
6. [Physical constants](#Physical-constants)
7. [Units in PlasmaPy](#Units-in-PlasmaPy)
8. [Optimizing unit operations](#Optimizing-unit-operations)
9. [Physical Types](#Physical-types)

+++

## Unit basics

+++

We can create a physical quantity by multiplying or dividing a number or array with a unit.

```{code-cell} ipython3
distance = 60 * u.km
print(distance)
```

This operation creates a {py:class}`~astropy.units.Quantity`: a number, sequence, or array that has been assigned a physical unit.

```{code-cell} ipython3
type(distance)
```

We can also create an object by using the {py:class}`~astropy.units.Quantity` class itself.

```{code-cell} ipython3
time = u.Quantity(120, u.min)
```

We can create {py:class}`~astropy.units.Quantity` objects with compound
units.

```{code-cell} ipython3
88 * u.imperial.mile / u.hour
```

We can even create {py:class}`~astropy.units.Quantity` objects that are
explicitly dimensionless.

```{code-cell} ipython3
3 * u.dimensionless_unscaled
```

We can also create a {py:class}`~astropy.units.Quantity` based off of a
[NumPy] array or a list.

```{code-cell} ipython3
import numpy as np

np.array([2.5, 3.2, 1.1]) * u.kg
```

```{code-cell} ipython3
[2, 3, 4] * u.m / u.s
```

## Unit operations

Operations between {py:class}`~astropy.units.Quantity` objects handle
unit conversions automatically. We can add
{py:class}`~astropy.units.Quantity` objects together as long as their
units have the same physical type.

```{code-cell} ipython3
1 * u.m + 25 * u.cm
```

Units get handled automatically during operations like multiplication,
division, and exponentiation.

```{code-cell} ipython3
distance / time
```

```{code-cell} ipython3
distance ** 2
```

Attempting an operation between physically incompatible units gives us
an error, which we can use to find bugs in our code.

```{code-cell} ipython3
:tags: [raises-exception]

3 * u.m + 3 * u.s
```

{py:class}`~astropy.units.Quantity` objects behave very similarly to
[NumPy] arrays because {py:class}`~astropy.units.Quantity` is a subclass
of {py:class}`numpy.ndarray`.

```{code-cell} ipython3
balmer_series = [656.279, 486.135, 434.0472, 410.1734] * u.nm
Hα = balmer_series[0]
print(Hα)
```

```{code-cell} ipython3
np.max(balmer_series)
```

Most frequently encountered [NumPy] and [SciPy] functions can be used
with {py:class}`~astropy.units.Quantity` objects.  However,
{py:class}`~astropy.units.Quantity` objects [lose their units with some
operations].

+++

## Unit conversions

The {py:meth}`~astropy.units.Quantity.to` method allows us to convert a
{py:class}`~astropy.units.Quantity` to different units of the same
physical type. This method accepts strings that represent a unit
(including compound units) or a unit object.

```{code-cell} ipython3
velocity = 54 * u.km / u.hr
velocity.to("m/s")
```

```{code-cell} ipython3
velocity.to(u.m / u.s)
```

The {py:attr}`~astropy.units.Quantity.si` and
{py:attr}`~astropy.units.Quantity.cgs`
attributes convert the {py:class}`~astropy.units.Quantity` to SI or CGS
units, respectively.

```{code-cell} ipython3
velocity.si
```

```{code-cell} ipython3
velocity.cgs
```

## Detaching units and values

The {py:attr}`~astropy.units.Quantity.value` attribute of a
{py:class}`~astropy.units.Quantity` provides the number (as a [NumPy]
scalar) or [NumPy] array without the unit.

```{code-cell} ipython3
time.value
```

The {py:attr}`~astropy.units.Quantity.unit` attribute of a
{py:class}`~astropy.units.Quantity` provides the unit without the value.

```{code-cell} ipython3
time.unit
```

## Equivalencies

+++

Plasma scientists often use the [electron-volt] (eV) as a unit of
temperature. This is a shortcut for describing the thermal energy per
particle, or more accurately the temperature multiplied by the
[Boltzmann constant], $k_B$. Because an electron-volt is a unit of
energy rather than temperature, we cannot directly convert
electron-volts to kelvin.

```{code-cell} ipython3
:tags: [raises-exception]

u.eV.to("K")
```

To handle non-standard unit conversions, {py:mod}`astropy.units` allows
the use of [equivalencies]. The conversion from eV to K can be done by
using the {py:func}`~astropy.units.temperature_energy` equivalency.

```{code-cell} ipython3
(1 * u.eV).to("K", equivalencies=u.temperature_energy())
```

Radians are treated dimensionlessly when the
{py:func}`~astropy.units.dimensionless_angles` equivalency is in effect.
Note that this equivalency does not account for the multiplicative
factor of $2π$ that is used when converting between [frequency] and
[angular frequency].

```{code-cell} ipython3
(3.2 * u.rad / u.s).to("1 / s", equivalencies=u.dimensionless_angles())
```

## Physical constants

+++

We can use {py:mod}`astropy.constants` to access the most commonly
needed physical constants.

```{code-cell} ipython3
from astropy.constants import c, e, eps0, hbar, k_B, mu0

print(c)
```

A {py:class}`astropy.constants.Constant` behaves very similarly to a
{py:class}`~astropy.units.Quantity`. For example, we can use the
Boltzmann constant to mimic the behavior of
{py:func}`~astropy.units.temperature_energy`.

```{code-cell} ipython3
thermal_energy_per_particle = 0.6 * u.keV
temperature = thermal_energy_per_particle / k_B
print(temperature.to("MK"))
```

Electromagnetic constants often need the unit system to be specified.
Code within PlasmaPy uses SI units.

```{code-cell} ipython3
:tags: [raises-exception]
2 * e
```

```{code-cell} ipython3
2 * e.si
```

## Units in PlasmaPy

+++

Now we can show some uses of {py:mod}`astropy.units` in PlasmaPy,
starting with {py:mod}`plasmapy.particles`. Many of the attributes of
{py:class}`~plasmapy.particles.particle_class.Particle` and
{py:class}`~plasmapy.particles.particle_collections.ParticleList`
provide {py:class}`~astropy.units.Quantity` objects.

```{code-cell} ipython3
from plasmapy.particles import Particle, ParticleList

alpha = Particle("He-4 2+")
```

```{code-cell} ipython3
alpha.charge
```

```{code-cell} ipython3
ions = ParticleList(["O 1+", "O 2+", "O 3+"])
ions.mass
```

Similarly, {py:class}`~astropy.units.Quantity` objects are the expected
inputs and outputs of most functions in {py:mod}`plasmapy.formulary`. We
can use them to calculate some plasma parameters for a typical region of
the solar corona.

```{code-cell} ipython3
from plasmapy.formulary import Alfven_speed, Debye_length, gyrofrequency
```

```{code-cell} ipython3
B = 0.01 * u.T
n = 1e15 * u.m ** -3
proton = Particle("p+")
```

```{code-cell} ipython3
Alfven_speed(B=B, density=n, ion=proton).to("km /s")
```

```{code-cell} ipython3
gyrofrequency(B=B, particle="e-")
```

The `to_hz` keyword provides the frequency in hertz rather than radians
per second, and accounts for the factor of $2π$.

```{code-cell} ipython3
gyrofrequency(B=B, particle="e-", to_hz=True)
```

Formulary functions perform calculations based on SI units, but accept
input arguments in other units. Temperature can be given in units of
temperature (e.g., kelvin) or energy (e.g., electron-volts).

```{code-cell} ipython3
Debye_length(T_e=1e6 * u.K, n_e=1e9 * u.m ** -3)
```

```{code-cell} ipython3
Debye_length(T_e=86.17 * u.eV, n_e=1e3 * u.cm ** -3)
```

## Optimizing unit operations

Astropy's documentation includes [performance tips] for using
{py:mod}`astropy.units` in computationally intensive situations. For
example, putting compound units in parentheses reduces the need to make
multiple copies of the data.

```{code-cell} ipython3
volume = 0.62 * (u.barn * u.Mpc)
```

## Physical types

+++

A [physical type] corresponds to physical quantities with dimensionally
compatible units. Astropy has functionality that represents different
physical types. These physical type objects can be accessed using either
the {py:attr}`~astropy.units.UnitBase.physical_type` attribute of a unit
or {py:func}`astropy.units.get_physical_type`.

```{code-cell} ipython3
(u.m ** 2 / u.s).physical_type
```

```{code-cell} ipython3
u.get_physical_type("number density")
```

These physical type objects can be used for dimensional analysis.

```{code-cell} ipython3
energy_density = (u.J * u.m ** -3).physical_type
velocity = u.get_physical_type("velocity")
print(energy_density * velocity)
```
