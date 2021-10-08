.. currentmodule:: plasmapy

********
Glossary
********

.. glossary::
   :sorted:

   args
      An abbreviation for `positional arguments`_.

   -like
      Used to indicate an `object` of that type or that can instantiate
      that type.  For example, ``"He 2+"`` is :term:`particle-like`
      because it can be used to instantiate |Particle|.

   charge number
      The charge of a particle in units of elementary charge.

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

      |DimensionlessParticle| instances are not particle-like because,
      without normalization information, they do not uniquely identify a
      physical particle.

      For more complete details, refer to |ParticleLike|.

   real number
      Any numeric type that represents a real number. This could include
      a `float`, `int`, a dimensionless |Quantity|, or any of the
      `numpy.number` types. Note that if a PlasmaPy function expects a
      dimensional |Quantity| and a real number is provided, then the
      real number is often assumed to have the appropriate SI units.

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

.. |argument| replace:: :term:`argument`
.. |parameter| replace:: :term:`parameter`

.. _`keyword arguments`: https://docs.python.org/3/glossary.html#term-argument
.. _`positional arguments`: https://docs.python.org/3/glossary.html#term-argument
