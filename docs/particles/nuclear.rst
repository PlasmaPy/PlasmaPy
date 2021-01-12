.. _particles-nuclear-reactions:

Nuclear Reactions
*****************

.. _particles-nuclear-binding-energy:

Binding energy
==============

The binding energy of a nuclide may be accessed either as an
attribute of a `~plasmapy.particles.Particle` object, or by using the
`~plasmapy.particles.nuclear_binding_energy` function.

>>> from plasmapy.particles import Particle, nuclear_binding_energy
>>> D = Particle('deuterium')
>>> D.binding_energy
<Quantity 3.56414847e-13 J>
>>> nuclear_binding_energy('D').to('GeV')
<Quantity 0.00222457 GeV>

.. _particles-nuclear-reaction-energy:

Nuclear reaction energy
=======================

The energy released from a nuclear reaction may be found using the
`~plasmapy.particles.nuclear_reaction_energy` function.  The input may be
a `str` representing the reaction.

>>> from plasmapy.particles import nuclear_reaction_energy
>>> nuclear_reaction_energy('Be-8 + alpha --> carbon-12')
<Quantity 1.18025735e-12 J>

The reaction may also be inputted using the ``reactants`` and
``products`` keywords.

>>> nuclear_reaction_energy(reactants=['D', 'T'], products=['alpha', 'n'])
<Quantity 2.81812097e-12 J>
