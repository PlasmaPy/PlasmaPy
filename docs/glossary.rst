.. currentmodule:: plasmapy

********
Glossary
********

.. glossary::

   args
      An abbreviation for `positional arguments`_.

   -like
      Used to indicate an `object` of that type or that can instantiate
      that type.  For example, ``"He 2+"`` is :term:`particle-like`
      because it can be used to instantiate |Particle|.

   charge number
      The charge of a particle in units of elementary charge.

   keyword-only
      An argument that must be provided along with the name of that
      argument. If ``z`` is a keyword-only argument to ``f(z)``, then
      it can be provided as ``f(z=2)`` but not ``f(2)``.

   kwargs
      An abbreviation for `keyword arguments`_.

   particle-like
      An `object` is *particle-like* if it is a |Particle| or
      |CustomParticle|, or can be cast into one.

      An **element** may be represented by a string containing the
      atomic symbol (case-sensitive), the name of the element, or an
      integer representing the atomic number. The element iron can be
      represented as ``"Fe"``, ``"iron"``, or ``26``.

      An **isotope** may be represented by a string that contains an
      atomic symbol or element name, followed by a hyphen and the mass
      number (with no spaces in between). The isotope :sup:`56`\ Fe can
      be represented as ``"Fe-56"``, or ``"iron-56"``. :sup:`2`\ H can
      be represented by ``"D"`` or ``"deuterium"``, and :sup:`3`\ H can
      be represented by ``"T"`` or ``"tritium"``.

      An **ion** or **neutral atom** may be represented by a string that
      contains a representation of an element or isotope, followed by
      charge information which is typically an integer representing the
      charge number and a plus or minus sign to indicate the charge. For
      example, a deuteron may be represented as ``"D 1+"`` and
      :sup:`56`\ Fe\ :sup:`1+` may be represented as ``"Fe-56 1+"``.

      A **special particle** may be represented by a string that
      contains the name of the particle (case insensitive) or a
      standard symbol for it (case insensitive). A neutron can be
      represented as ``"n"`` or ``"neutron"``; a proton can be
      represented as ``"p+"``, ``"p"``, or ``"proton"``; and an electron
      can be represented by ``"e-"``, ``"e"``, or ``"electron"``.

      For more complete details, refer to |ParticleLike|.

   real number
      Any numeric type that represents a real number. This could include
      a  `float`, `int`, a dimensionless |Quantity|, or any of the
      `numpy.number` types. Note that if a PlasmaPy function expects a
      dimensional |Quantity| and a real number is provided, then the
      real number is often assumed to have the appropriate SI units.

.. _`keyword arguments`: https://docs.python.org/3/glossary.html#term-argument
.. _`positional arguments`: https://docs.python.org/3/glossary.html#term-argument
