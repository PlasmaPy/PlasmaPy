.. _particles-decorators:

Decorators
**********

.. _particles-decorators-particle-input:

Passing `~plasmapy.particles.Particle` objects to functions and methods
=======================================================================

When calculating plasma parameters, we very frequently need to access
the properties of the particles that make up that plasma. The
`~plasmapy.particles.decorators.particle_input` decorator allows functions
and methods to easily access properties of different particles.

The `~plasmapy.particles.particle_input` decorator takes valid
representations of particles given in arguments to functions and passes
through the corresponding `~plasmapy.particles.Particle` object.  The
arguments must be annotated with `~plasmapy.particles.Particle`
so that the decorator knows to create the `~plasmapy.particles.Particle`
object.  The decorated function can then access particle properties by using
`~plasmapy.particles.Particle` attributes.  This decorator will raise an
`~plasmapy.particles.exceptions.InvalidParticleError` if the input does
not correspond to a valid particle.

Here is an example of a decorated function.

.. code-block:: python

  from plasmapy.particles import Particle, particle_input

  @particle_input
  def particle_mass(particle: Particle):
      return particle.mass

This function can now accept either `~plasmapy.particles.Particle` objects
or valid representations of particles.

>>> particle_mass('p+')  # string input
<Quantity 1.67262192e-27 kg>
>>> proton = Particle("proton")
>>> particle_mass(proton)  # Particle object input
<Quantity 1.67262192e-27 kg>

If only one positional or keyword argument is annotated with
`~plasmapy.particles.Particle`, then the keywords ``mass_numb`` and ``Z``
may be used when the decorated function is called.

.. code-block:: python

  @particle_input
  def integer_charge(particle: Particle, Z: int = None, mass_numb: int = None) -> int:
      return particle.integer_charge

The above example includes optional type hint annotations for ``Z`` and
``mass_numb`` and the returned value.  The
`~plasmapy.particles.particle_input` decorator may be used in methods in
classes as well:

.. code-block:: python

  class ExampleClass:
      @particle_input
      def particle_symbol(self, particle: Particle) -> str:
          return particle.symbol

On occasion it is necessary for a function to accept only certain
categories of particles.  The `~plasmapy.particles.particle_input`
decorator enables several ways to allow this.

If an annotated keyword is named ``element``, ``isotope``, or ``ion``;
then `~plasmapy.particles.particle_input` will raise an
`~plasmapy.particles.exceptions.InvalidElementError`,
`~plasmapy.particles.exceptions.InvalidIsotopeError`, or
`~plasmapy.particles.exceptions.InvalidIonError` if the particle is not
associated with an element, isotope, or ion; respectively.

.. code-block:: python

  @particle_input
  def capitalized_element_name(element: Particle):
      return element.element_name

  @particle_input
  def number_of_neutrons(isotope: Particle):
      return isotope.mass_number - isotope.atomic_number

  @particle_input
  def number_of_bound_electrons(ion: Particle):
      return ion.atomic_number - ion.integer_charge

The keywords ``require``, ``any_of``, and ``exclude`` to the
decorator allow further customization of the particle categories
allowed as inputs.  These keywords are used as in
`~plasmapy.particles.Particle.is_category`.

.. code-block:: python

  @particle_input(require='charged')
  def sign_of_charge(charged_particle: Particle):
      """Require a charged particle."""
      return '+' if charged_particle.integer_charge > 0 else '-'

  @particle_input(any_of=['charged', 'uncharged'])
  def integer_charge(particle: Particle) -> int:
      """Accept only particles with charge information."""
      return particle.integer_charge

  @particle_input(exclude={'antineutrino', 'neutrino'})
  def particle_mass(particle: Particle):
      """
      Exclude neutrinos/antineutrinos because these particles have
      weakly constrained masses.
      """
      return particle.mass
