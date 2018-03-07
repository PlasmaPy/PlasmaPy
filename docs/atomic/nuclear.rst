-----------------
Nuclear Reactions
-----------------

The binding energy of a nuclide may be accessed either as an
attribute of a `~plasmapy.atomic.Particle` instance, or by using the
`~plasmapy.atomic.nuclear_binding_energy` function.

>>> from plasmapy.atomic import Particle, nuclear_binding_energy
>>> D = Particle('deuterium')
>>> D.binding_energy
<Quantity 3.56414847e-13 J>
>>> nuclear_binding_energy('D').to('GeV')
<Quantity 2.22456652 MeV>

The energy released from a nuclear reaction may be found using the
`~plasmapy.atomic.nuclear_reaction_energy` function.  The input may be
a `str` representing the reaction.

>>> from plasmapy.atomic import nuclear_reaction_energy
>>> nuclear_reaction_energy('Be-8 -> 2*alpha')
<Quantity 1.47143078e-14 J>

The reaction may also be inputted using the `reactants` and `products`
keywords.

>>> nuclear_reaction_energy(reactants=['D', 'T'], products=['alpha', 'n'])
<Quantity 17.58932778 MeV>
