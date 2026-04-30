.. currentmodule:: plasmapy

.. _glossary:

********
Glossary
********

.. glossary::
   :sorted:

   alias
   aliases
      An abbreviated version of a commonly used function.  For
      example, `~plasmapy.formulary.speeds.va_` is an alias for
      `~plasmapy.formulary.speeds.Alfven_speed`.  Aliases are
      named with a trailing underscore.

      For further details, please refer to the :ref:`contributor guide's
      section on aliases <aliases>`.

   args
      An abbreviation for `positional arguments`_.

   atom-like
      A |particle-like| `object` is **atom-like** if it is or could be
      cast into:

      * A |Particle| representing an element, isotope, or ionic level; or
      * A |ParticleList| including only elements, isotopes, or ionic levels.

      For example, ``"p+"``, ``"He-4"``, ``"deuterium"``, ``"O 0+"``,
      ``Particle("Fe-56 16+)"``, ``["He-4 1+", "He-4 2+"]``,  and
      integers representing atomic numbers are all atom-like.

      Examples of objects that are |particle-like| but *not* atom-like
      include ``"neutron"``, ``"e-"``, and ``["e-", "e+"]``.
      Additionally, ``["He-4", "e-"]`` is not atom-like because this
      `list` contains an item that is not atom-like.

      Please refer to the glossary entry for |particle-like| for a full
      description of valid representations of elements, isotopes, and
      ions.

   -like
      Used to indicate an `object` of that type or that can instantiate
      that type.  For example, ``"He 2+"`` is |particle-like| because it
      can be used to instantiate |Particle|.

   charge number
      The electrical charge of a particle in units of the elementary
      charge. The charge number of an ion or neutral particle is usually
      denoted as ``Z``.

   collective behavior
      In |plasma| science, collective behavior refers to particles being
      able to interact with each other from a distance.

   fit-function
   fit-functions
      Any instance of a subclass of
      `~plasmapy.analysis.fit_functions.AbstractFitFunction`.  Also see
      module `~plasmapy.analysis.fit_functions`.

   force-free
      In plasma physics, a **force-free magnetic field** is a magnetic
      field where the Lorentz force is zero, meaning the magnetic
      pressure greatly exceeds the plasma pressure, allowing
      non-magnetic forces to be neglected. This condition is often
      approximated in the Sun's corona.

      For more details, visit the `Wikipedia page <https://en.wikipedia.org/wiki/Force-free_magnetic_field>`_.

   gas
      A compressible fluid consisting of particles that interact with
      each other predominantly through short-range binary collisions.

      Gas and plasma are distinct states of matter with fundamentally
      different properties. The overwhelming majority of particles in a
      gas are neutral and thus do not exert long-range forces on each
      other. While a real gas may contain trace amounts of charged
      particles, their effects are negligible. In contrast, a plasma
      contains enough charged particles able to exert long-range forces
      on each other to result in |collective behavior|.

      .. caution::

         The term "gas" is frequently used to encompass fluids in either
         of the gas or plasma states of matter, in particular in the
         astronomical community.

   integration test
      An **integration test** verifies that multiple software
      components work together as intended.

      Compared to a :term:`unit test`, an integration test is typically
      harder to write, slower to run, more difficult to maintain, and
      less useful at pinpointing the specific cause of a problem.
      However, integration tests are able to find problems that unit
      tests cannot. In particular, integration tests are able to find
      problems at the interfaces between different modules. On average,
      each integration test covers more lines of code than each related
      :term:`unit test`. Because unit tests and integration tests
      complement each other, both are important constituents of a test
      suite.

   keyword-only
      An |argument| or |parameter| is *keyword-only* when the |argument|
      must be provided with the name of the corresponding |parameter|.

      If ``z`` is a keyword-only |parameter| to ``f(z)``, then the
      |argument| ``2`` can be provided as ``f(z=2)`` but not ``f(2)``.

   kwargs
      An abbreviation for `keyword arguments`_.

   lite-function
   lite-functions
      An optimized version of an existing PlasmaPy function intended
      for applications where computational efficiency is most important.
      While most `~plasmapy.formulary` functions accept |Quantity|
      objects created using `astropy.units`, lite-functions accept
      numbers and |array_like| inputs that are implicitly assumed to be
      in SI units. The name of a lite-function ends with ``_lite``. A
      lite-function can be accessed as the ``lite`` attribute of the
      corresponding regular function.

      .. caution::

         Unlike most `~plasmapy.formulary` functions, no validations are
         performed on the arguments provided to a lite-function for the
         sake of computational efficiency. When using lite-functions, it
         is vital to double-check your implementation!

      For further details, please refer to the :ref:`contributor guide's
      section on lite-functions <lite-functions>`.

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
      charge number and a plus or minus sign to indicate the electrical
      charge. For example, a deuteron may be represented as ``"D 1+"``
      and :sup:`56`\ Fe\ :sup:`1+` may be represented as ``"Fe-56 1+"``.

      A **special particle** may be represented by a string that
      contains the name of the particle (case insensitive) or a
      standard symbol for it (case insensitive). A neutron can be
      represented as ``"n"`` or ``"neutron"``; a proton can be
      represented as ``"p+"``, ``"p"``, or ``"proton"``; and an electron
      can be represented by ``"e-"``, ``"e"``, or ``"electron"``.

      |DimensionlessParticle| instances are not particle-like because,
      without normalization information, they do not uniquely identify a
      physical particle.

      For more complete details, refer to |ParticleLike|.

   particle-list-like
      An `object` is *particle-list-like* if it is a |ParticleList|, or
      can be cast into one.

      For more complete details, refer to |ParticleListLike|.

   plasma
      A (usually |quasineutral|) fluid that contains enough charged
      particles (usually ions and electrons) for the fluid to exhibit
      |collective behavior|. A plasma can also contain neutral atoms,
      neutral molecules, and (usually charged) dust particles.

      *Plasma* usually refers to a |quasineutral| fluid that contains
      ions and free electrons, often along with neutral atoms, neutral
      molecules, and/or charged dust grains.


, and that may also
      contain neutral atoms and molecules,

      While neutral gas particles interact via


      *Plasma* usually refers to a fluid contain

      A fluid that contains a large enough fraction of charged particles
      (usually ions and electrons) for fluid behavior to depend on the
      forces

      A fluid that contains a large enough fraction of charged particles
      to cause the dynamics of the fluid


to exhibit
      |collective behavior|.

      *Plasma* usually refers to

      A fluid that contains a significant fraction of charged particles
      so that fluid behavior depends on

A fluid that contains an appreciable fraction of charged particles


      A fluid that contains an appreciable fraction of charged particles
      and exhibits |collective behavior|.


      A fluid that contains an appreciable fraction of charged particles.

      A fluid that contains a large enough fraction of charged particles
      to exhibit |collective behavior|.

      *Plasma* most often refers to a |quasineutral| fluid that contains

      A fluid that contains a large enough fraction of charged particles
      for its behavior to be influenced by long-range forces

      A fluid that contains enough electrically charged particles for
      its behavior to be influenced by long-range forces exert
      long-range forces exerted by charged particles on each other.

      The term *plasma* is most often used to refer to a |quasineutral|
      fluid that contains an appreciable fraction of electrons and
      positively charged ions

      The term *plasma* typically refers to a plasma containing unbound
      electrons and positively charged

      The term "plasma" typically refers to a |quasineutral| fluid
      containing an appreciable fraction of ions and electrons.

      The term "plasma" typically refers to a |quasineutral| fluid
      containing charged particles and neutral particles that exhibits
      |collective behavior|. The charged particles typically include
      ions and electrons, and might include charged dust particulates.
      The neutral particles include atoms and molecules.

      The term "plasma" is also used when describing a non-neutral
      plasma, electron-positron plasma, or quark-gluon plasma.

      While plasma is often referred to as "gas" or "ionized gas" in
      astronomy literature and films like _Star Trek VI_, plasma is a
      state of matter that is distinct from the gas state of matter.


      .. A fluid that contains an appreciable fraction of charged particles
      and exhibits

      .. A fluid that contains enough charged particles for the fluid to
      exhibit |collective behavior|.

      .. A fluid that contains an appreciable fraction of particles that
      exert long-range forces on each other.

      .. A fluid that contains enough charged particles for the fluid to
      exhibit |collective behavior|.

      .. A |quasineutral| fluid that contains enough charged particles for
      the fluid to exhibit |collective behavior| :cite:p:`chen:2016`.

      .. A :term:`quasineutral` fluid containing charged and neutral

      .. The term *plasma* most frequently refers to a fluid that contains

      .. A fluid that contains enough particles that exert long-range
      forces fraction of particles that

      .. A fluid that exhibits |collective behavior| resulting from
      particles that exert long-range forces on each other.

      .. A fluid that contains enough particles that exert long-range
      forces on each other to result in |collective behavior|.

      .. A fluid whose behavior is influenced to an appreciable degree by
      particles exerting long-range forces on each other

      .. A fluid containing enough particles that exert long-range forces
      on each other to result in |collective behavior|.

      .. A fluid that contains enough particles that exert long-range
      forces on each other to influence fluid behavior.

      .. More specifically, a plasma is a fluid containing

      .. A fluid that contains enough charged particles

      .. A fluid containing an appreciable fraction of particles that exert
      long-range forces on

      .. A fluid whose behavior is influenced by particles exerting
      long-range forces on each other.

      .. containing enough charged particles exerting long-range
      forces on each other

      .. A fluid whose behavior is appreciably influenced by long-range forces exerted by particles exerting long-range forces on each other.

      .. A fluid that contains enough particles exerting long-range forces on each other

      .. A fluid that contains enough particles that exert long-range
      forces on each other to influence the fluid's macroscopic

      .. A fluid whose macroscopic behavior is influenced by charegd

      .. A fluid containing an appreciable fraction of particles that exert
      long-range forces on each other

      .. A fluid that contains appreciable fraction of particles exert
      long-range forces on each other

      .. whose macroscopic behavior is influenced by long-range
      forces exerted by particles

      .. A fluid that contains charged particles where

      .. A fluid whose macroscopic behavior is influenced by long-range
      forces exerted by particles on other particles.

      .. A fluid that contains enough particles exerting long-range forces
      on each other fo the fluid to exhibit |collective behavior|.

      .. Start with the broadest definition to include positron-electron
         plasmas, quark-gluon plasmas, strongly coupled plasmas, and

      .. A fluid containing enough particles exerting long-range forces on
      each other for the fluid to exhibit |collective behavior|.

      .. A fluid containing particles the exert long-range forces on each
      other resulting in |collective behavior|.

      .. A fluid that exhibits |collective behavior| due to particles
      exerting long-range forces on each other.

      .. A fluid that exhibits |collective behavior| resulting from
      particles exerting long-range forces on each other.

      .. A fluid that exhibits |collective behavior|.

      .. A fluid that contains enough particles able to exert a force on
      each other from a distance for the fluid to engage in
      |collective behavior|.

      .. A fluid in which an appreciable fraction of the particles exert
      long-range forces on each other.

      .. A fluid in which an appreciable fraction of the particles exert
      long-range forces on each other such that the fluid exhibits
      |collective behavior|.

      .. A fluid that contains an appreciable fraction of particles that
      can exert forces on each other from a distance.

      .. A fluid that contains an appreciable fraction of particles that
      exert forces on each other outside of binary collisions.

      .. A fluid that exhibits |collective behavior| due to particles
      that exerting forces on each other from a distance.

      .. Above are brainstorms for the first sentence(s). Below is
         content for the rest of the definition.

   quasineutral
      A fluid containing an appreciable fraction of charged particles is
      quasineutral when the macroscopic charge density is roughly zero.

   real number
      Any numeric type that represents a real number. This could include
      a `float`, `int`, a dimensionless |Quantity|, or any of the
      `numpy.number` types. Note that if a PlasmaPy function expects a
      dimensional |Quantity| and a real number is provided, then the
      real number is often assumed to have the appropriate SI units.

   temperature
      Most functions in PlasmaPy accept temperature, :math:`T`, as a
      `~astropy.units.Quantity` with units of temperature (e.g., kelvin)
      or energy (e.g., electron-volts). A value for energy that is
      provided will be divided by the Boltzmann constant, :math:`k_B`,
      to be converted into units of temperature.

   unit test
      A **unit test** verifies a single unit of behavior, does it
      quickly, and does it in isolation from other tests
      :cite:p:`khorikov:2020`.

      Unit tests are intended to provide fast feedback that help pinpoint
      the locations of errors. Unit tests often abide by the following
      pattern :cite:p:`osherove:2013`:

      1. *Arrange*: gather inputs and get the system to the state in which
         the test is expected to run.

      2. *Act*: make the system under test undertake the operation that is
         being tested.

      3. *Assert*: verify that the actual outcome of the *act* phase matches
         the expected outcome.

      In a unit test for a function, the *arrange* phase involves
      collecting or constructing the inputs for the function. The *act*
      phase occurs when the function is called with those inputs. The
      *assert* phase is when the value returned by the function is
      compared to the expected result.

.. _`keyword arguments`: https://docs.python.org/3/glossary.html#term-argument
.. _`positional arguments`: https://docs.python.org/3/glossary.html#term-argument

.. The definitions of "plasma", "gas", and "collective behavior" in this
   file are made available under either PlasmaPy's license or the
   CC0 1.0 Universal license, which is available on the World Wide Web
   at: https://creativecommons.org/publicdomain/zero/1.0
