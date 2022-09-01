from __future__ import annotations

import astropy.units as u

from numbers import Integral
from typing import Iterable, Sequence, Union

from plasmapy.particles import CustomParticle, Particle
from plasmapy.particles.particle_collections import ParticleList

RealParticleLike = Union[str, Integral, Particle]
CustomParticleLike = Union[u.Quantity, CustomParticle]
SingleParticleLike = Union[RealParticleLike, CustomParticleLike]
ParticleListLike = Iterable[SingleParticleLike]
ParticleLike = Union[SingleParticleLike, ParticleListLike]

ParticleLike.__doc__ = r"""
An `object` is |particle-like| if it is a |Particle|, |CustomParticle|,
or |ParticleList|; or could be cast into one.

When used as a type hint annotation, |ParticleLike| indicates that an
argument should represent one or more physical particles. An `object`
representing a single particle could be a string, integer, |Quantity|
[|mass|], |Quantity| [|charge|], |Particle|, or |CustomParticle|.

For more information, please refer to |particle-like|.

An
`object` representing multiple particles


Objects that
can represent a single particle



Notes
-----
Real world particles are typically represented as instances of the
`~plasmapy.particles.particle_class.Particle` class in PlasmaPy.

>>> from plasmapy.particles import Particle
>>> Particle("proton")
Particle("p+")

All `~plasmapy.particles.particle_class.Particle` instances, and objects
that can be cast into `~plasmapy.particles.particle_class.Particle`
instances, are particle-like.

* **Elements**

    An element may also be represented by a string that contains the atomic
    symbol (case-sensitive) or the name of the element, or an integer
    representing the atomic number. The element iron can be represented as
    ``"Fe"``, ``"iron"``, ``"Iron"``, ``26``, or ``Particle("Fe")``.

* **Isotopes**

    An isotope may be represented by a string that contains an atomic symbol
    or element name, followed by a hyphen and the mass number (with no spaces
    in between). The isotope :sup:`56`\ Fe can be represented as
    ``"Fe-56"``, ``"iron-56"``, or ``Particle("Fe-56")``. :sup:`1`\ H can be
    represented by ``"protium"``, :sup:`2`\ H can be represented by ``"D"``
    or ``"deuterium"``, and :sup:`3`\ H can be represented by ``"T"`` or
    ``"tritium"``.

* **Ions**

    An ion or ionic level may be represented by a string that contains a
    representation of an element or isotope, followed by charge information.
    For example, ``"He 1+"``, ``"He+"``, ``"helium 1+"``, and ``"He II"``
    all represent singly ionized helium.

    Charge information is typically separated from the element or isotope by
    a space, and given as an integer paired with a plus or minus sign. The
    sign can either precede or follow the integer (e.g., ``"Fe 0+"`` or
    ``"Fe +0"``). The charge information can also be given as a series of
    plus signs or of minus signs that immediately follow the element or
    isotope (e.g., ``"Fe++"`` for Fe\ :sup:`2+`\ ).

    Ions can also be represented using Roman numeral notation, where the Roman
    numeral indicates the charge number plus one (e.g., ``"H I"`` represents
    H\ :sup:`0+` and ``"He-4 II"`` represents :sup:`4`\ He\ :sup:`1+`\ ).

    D\ :sup:`1+` can also be represented by ``"deuteron"``, T\ :sup:`1+` can
    be represented by ``"triton"``, and :sup:`4`\ He\ :sup:`2+` can be
    represented by ``"alpha"``.

* **Special particles**

    A special particle may be represented by a string that contains
    the name of the particle (case-insensitive) or a standard symbol for it
    (case-sensitive). A neutron can be represented as ``"n"`` or
    ``"neutron"``; a proton can be represented as ``"p+"``, ``"p"``, or
    ``"Proton"``; and an electron can be represented by ``"e-"``, ``"e"``,
    or ``"ELECTRON"``.

* **Custom particles**

    `~plasmapy.particles.particle_class.CustomParticle` instances are
    particle-like because particle properties are provided in physical units.

.. note::

    `~plasmapy.particles.particle_class.DimensionlessParticle`
    instances are *not* particle-like because, without normalization
    information, they do not uniquely identify a physical particle.

See Also
--------
~plasmapy.particles.particle_class.Particle
~plasmapy.particles.particle_class.CustomParticle
~plasmapy.particles.particle_collections.ParticleList
~plasmapy.particles.decorators.particle_input


Examples
--------
Using `ParticleLike` as a type hint annotation indicates that an
argument or variable should represent a physical particle.

>>> from plasmapy.particles import ParticleLike, Particle
>>> def is_electron(particle: ParticleLike):
...     return particle == Particle("e-")
"""
ParticleListLike = Union[ParticleList, Sequence[ParticleLike]]


ParticleListLike.__doc__ = r"""
An `object` is |particle-list-like| if it can be identified as a
|ParticleList| or cast into one.

When used as a type hint annotation, |ParticleListLike| indicates that
the corresponding argument should represent a sequence of physical
particles. Each item in a |ParticleListLike| object must be
|particle-like|.

Notes
-----
`~plasmapy.particles.particle_class.DimensionlessParticle` instances do
not uniquely represent a physical particle, and are thus not
|ParticleLike| and cannot be contained in a |ParticleListLike| object.

See Also
--------
~plasmapy.particles.particle_collections.ParticleList
~plasmapy.particles.particle_class.ParticleLike
~plasmapy.particles.decorators.particle_input

Examples
--------
Using |ParticleListLike| as a type hint annotation indicates that an
argument or variable should represent a sequence of |ParticleLike|
objects.

>>> from plasmapy.particles import ParticleList, ParticleListLike
>>> def contains_only_leptons(particles: ParticleListLike):
...     particle_list = ParticleList(particles)
...     return all(particle_list.is_category("lepton"))
>>> contains_only_leptons(["electron", "muon"])
True
"""
