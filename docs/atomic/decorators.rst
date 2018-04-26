.. _atomic-decorators

Decorators
**********

.. _atomic-decorators

Passing `~plasmapy.atomic.Particle` Instances to Functions and Methods
======================================================================

When calculating plasma parameters, we very frequently need to access
the properties of the particles that make up that plasma. The
`~plasmapy.atomic.particle_input` decorator allows functions and
methods to easily the properties of a particle.

The `~plasmapy.atomic.particle_input` decorator takes valid
representations of particles given in arguments to functions and passes
through the corresponding instance of the `~plasmapy.atomic.Particle`
class.  The arguments must be annotated with `~plasmapy.atomic.Particle`
so that the decorator knows to create the `~plasmapy.atomic.Particle`
instance.  The function can then access particle properties by using
`~plasmapy.atomic.Particle` attributes.  This decorator will raise an
`~plasmapy.utils.InvalidParticleError` if the input does not correspond
to a valid particle.

.. code-block:: python

  from plasmapy.atomic import Particle, particle_input

  @particle_input
  def particle_mass(particle: Particle):
      return particle.mass

If only one positional or keyword argument is annotated with
`~plasmapy.atomic.Particle`, then the keywords `mass_numb` and `Z` may
be used when the decorated function is called.

.. code-block:: python

  @particle_input
  def integer_charge(particle: Particle, Z: int = None, mass_numb: int = None) -> int:
      return particle.integer_charge

The above example includes optional type hint annotations for `Z` and
`mass_numb` and the returned value.  The
`~plasmapy.atomic.particle_input` decorator may be used in methods in
classes as well:

.. code-block:: python

  class ExampleClass:
      @particle_input
      def particle_symbol(self, particle: Particle) -> str:
          return particle.particle

On occasion it is necessary for a function to accept only certain
categories of particles.  The `~plasmapy.atomic.particle_input`
decorator enables several ways to allow this.

If an annotated keyword is named `element`, `isotope`, or `ion`; then
`~plasmapy.atomic.particle_input` will raise an
`~plasmapy.utils.InvalidElementError`,
`~plasmapy.utils.InvalidIsotopeError`, or
`~plasmapy.utils.InvalidIonError` if the particle is not associated with
an element, isotope, or ion; respectively.

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

The keywords `require`, `any_of`, and `exclude` to the decorator allow
further customization of the particle categories allowed as inputs.
These keywords are used as in `~plasmapy.atomic.Particle.is_category`.

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
