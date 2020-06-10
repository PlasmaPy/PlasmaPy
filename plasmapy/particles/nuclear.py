"""Functions that are related to nuclear reactions."""

import re
from typing import List, Optional, Union

from astropy import units as u

from plasmapy.particles.exceptions import AtomicError, InvalidParticleError
from plasmapy.particles.particle_class import Particle
from plasmapy.particles.particle_input import particle_input

__all__ = ["nuclear_binding_energy", "nuclear_reaction_energy", "mass_energy"]


@particle_input(any_of={"isotope", "baryon"})
def nuclear_binding_energy(
    particle: Particle, mass_numb: Optional[int] = None
) -> u.Quantity:
    """
    Return the nuclear binding energy associated with an isotope.

    Parameters
    ----------
    particle: `str`, `int`, or `~plasmapy.particles.Particle`
        A Particle object, a string representing an element or isotope,
        or an integer representing the atomic number of an element.

    mass_numb: `int`, optional
        The mass number of an isotope, which is required if and only
        if the first argument can only be used to determine the
        element and not the isotope.

    Returns
    -------
    binding_energy: `~astropy.units.Quantity`
        The binding energy of the nucleus in units of joules.

    Raises
    ------
    `~plasmapy.utils.InvalidParticleError`
        If the inputs do not correspond to a valid particle.

    `~plasmapy.utils.AtomicError`
        If the inputs do not correspond to a valid isotope or nucleon.

    `TypeError`
        If the inputs are not of the correct types.

    See Also
    --------
    `~plasmapy.particles.nuclear_reaction_energy` : Return the change in
        binding energy during nuclear fusion or fission reactions.

    `~plasmapy.particles.mass_energy` : Return the mass energy of a
        nucleon or particle.

    Examples
    --------
    >>> from astropy import units as u
    >>> nuclear_binding_energy('Fe-56').to(u.MeV)
    <Quantity 492.25957 MeV>
    >>> nuclear_binding_energy(26, 56)
    <Quantity 7.8868678e-11 J>
    >>> nuclear_binding_energy('p')  # proton
    <Quantity 0. J>
    >>> from astropy import units as u
    >>> before = nuclear_binding_energy("D") + nuclear_binding_energy("T")
    >>> after = nuclear_binding_energy("alpha")
    >>> (after - before).to(u.MeV)  # released energy from D + T --> alpha + n
    <Quantity 17.589 MeV>

    """
    return particle.binding_energy.to(u.J)


@particle_input
def mass_energy(particle: Particle, mass_numb: Optional[int] = None) -> u.Quantity:
    """
    Return a particle's mass energy.  If the particle is an isotope or
    nuclide, return the nuclear mass energy only.

    Parameters
    ----------
    particle: `str`, `int`, or `~plasmapy.particles.Particle`
        A Particle object, a string representing an element or isotope,
        or an integer representing the atomic number of an element.

    mass_numb: `int` (optional)
        The mass number of an isotope, which is required if and only
        if the first argument can only be used to determine the
        element and not the isotope.

    Returns
    -------
    mass_energy: `~astropy.units.Quantity`
        The mass energy of the particle (or in the case of ) in units of joules.

    Raises
    ------
    `~plasmapy.utils.InvalidParticleError`
        If the inputs do not correspond to a valid particle.

    `~plasmapy.utils.AtomicError`
        If the inputs do not correspond to a valid isotope or nucleon.

    `TypeError`
        If the inputs are not of the correct types.

    >>> mass_energy('He-4')
    <Quantity 5.9719e-10 J>

    """
    return particle.mass_energy


def nuclear_reaction_energy(*args, **kwargs):
    """
    Return the released energy from a nuclear reaction.

    Parameters
    ----------
    reaction: `str` (optional, positional argument only)
        A string representing the reaction, like
        ``"D + T --> alpha + n"`` or ``"Be-8 --> 2 * He-4"``.

    reactants: `list`, `tuple`, or `str`, optional, keyword-only
        A `list` or `tuple` containing the reactants of a nuclear
        reaction (e.g., ``['D', 'T']``), or a string representing the
        sole reactant.

    products: `list`, `tuple`, or `str`, optional, keyword-only
        A list or tuple containing the products of a nuclear reaction
        (e.g., ``['alpha', 'n']``), or a string representing the sole
        product.

    Returns
    -------
    energy: `~astropy.units.Quantity`
        The difference between the mass energy of the reactants and
        the mass energy of the products in a nuclear reaction.  This
        quantity will be positive if the reaction is exothermic
        (releases energy) and negative if the reaction is endothermic
        (absorbs energy).

    Raises
    ------
    `AtomicError`:
        If the reaction is not valid, there is insufficient
        information to determine an isotope, the baryon number is
        not conserved, or the charge is not conserved.

    `TypeError`:
        If the positional input for the reaction is not a string, or
        reactants and/or products is not of an appropriate type.

    See Also
    --------
    `~plasmapy.particles.nuclear_binding_energy` : finds the binding energy
        of an isotope

    Notes
    -----
    This function requires either a string containing the nuclear
    reaction, or reactants and products as two keyword-only lists
    containing strings representing the isotopes and other particles
    participating in the reaction.

    Examples
    --------
    >>> from astropy import units as u

    >>> nuclear_reaction_energy("D + T --> alpha + n")
    <Quantity 2.8181e-12 J>

    >>> triple_alpha1 = '2*He-4 --> Be-8'
    >>> triple_alpha2 = 'Be-8 + alpha --> carbon-12'
    >>> energy_triplealpha1 = nuclear_reaction_energy(triple_alpha1)
    >>> energy_triplealpha2 = nuclear_reaction_energy(triple_alpha2)
    >>> print(energy_triplealpha1, energy_triplealpha2)
    -1.471430e-14 J 1.1802573e-12 J
    >>> energy_triplealpha2.to(u.MeV)
    <Quantity 7.3665870 MeV>

    >>> nuclear_reaction_energy(reactants=['n'], products=['p+', 'e-'])
    <Quantity 1.25343e-13 J>

    """

    # TODO: Allow for neutrinos, under the assumption that they have no mass.

    # TODO: Add check for lepton number conservation; however, we might wish
    # to have violation of lepton number issuing a warning since these are
    # often omitted from nuclear reactions when calculating the energy since
    # the mass is tiny.

    errmsg = f"Invalid nuclear reaction."

    def process_particles_list(
        unformatted_particles_list: List[Union[str, Particle]]
    ) -> List[Particle]:
        """
        Take an unformatted list of particles and puts each
        particle into standard form, while allowing an integer and
        asterisk immediately preceding a particle to act as a
        multiplier.  A string argument will be treated as a list
        containing that string as its sole item.
        """

        if isinstance(unformatted_particles_list, str):
            unformatted_particles_list = [unformatted_particles_list]

        if not isinstance(unformatted_particles_list, (list, tuple)):
            raise TypeError(
                "The input to process_particles_list should be a "
                "string, list, or tuple."
            )

        particles = []

        for original_item in unformatted_particles_list:

            try:
                item = original_item.strip()

                if item.count("*") == 1 and item[0].isdigit():
                    multiplier_str, item = item.split("*")
                    multiplier = int(multiplier_str)
                else:
                    multiplier = 1

                try:
                    particle = Particle(item)
                except (InvalidParticleError) as exc:
                    raise AtomicError(errmsg) from exc

                if particle.element and not particle.isotope:
                    raise AtomicError(errmsg)

                [particles.append(particle) for i in range(multiplier)]

            except Exception:
                raise AtomicError(
                    f"{original_item} is not a valid reactant or "
                    "product in a nuclear reaction."
                ) from None

        return particles

    def total_baryon_number(particles: List[Particle]) -> int:
        """
        Find the total number of baryons minus the number of
        antibaryons in a list of particles.
        """
        total_baryon_number = 0
        for particle in particles:
            total_baryon_number += particle.baryon_number
        return total_baryon_number

    def total_charge(particles: List[Particle]) -> int:
        """
        Find the total integer charge in a list of nuclides
        (excluding bound electrons) and other particles.
        """
        total_charge = 0
        for particle in particles:
            if particle.isotope:
                total_charge += particle.atomic_number
            elif not particle.element:
                total_charge += particle.integer_charge
        return total_charge

    def add_mass_energy(particles: List[Particle]) -> u.Quantity:
        """
        Find the total mass energy from a list of particles, while
        taking the masses of the fully ionized isotopes.
        """
        total_mass_energy = 0.0 * u.J
        for particle in particles:
            total_mass_energy += particle.mass_energy
        return total_mass_energy.to(u.J)

    input_err_msg = (
        "The inputs to nuclear_reaction_energy should be either "
        "a string representing a nuclear reaction (e.g., "
        "'D + T -> He-4 + n') or the keywords 'reactants' and "
        "'products' as lists with the nucleons or particles "
        "involved in the reaction (e.g., reactants=['D', 'T'] "
        "and products=['He-4', 'n']."
    )

    reaction_string_is_input = args and not kwargs and len(args) == 1

    reactants_products_are_inputs = kwargs and not args and len(kwargs) == 2

    if reaction_string_is_input == reactants_products_are_inputs:
        raise AtomicError(input_err_msg)

    if reaction_string_is_input:

        reaction = args[0]

        if not isinstance(reaction, str):
            raise TypeError(input_err_msg)
        elif "->" not in reaction:
            raise AtomicError(
                f"The reaction '{reaction}' is missing a '->'"
                " or '-->' between the reactants and products."
            )

        try:
            LHS_string, RHS_string = re.split("-+>", reaction)
            LHS_list = re.split(r" \+ ", LHS_string)
            RHS_list = re.split(r" \+ ", RHS_string)
            reactants = process_particles_list(LHS_list)
            products = process_particles_list(RHS_list)
        except Exception as ex:
            raise AtomicError(f"{reaction} is not a valid nuclear reaction.") from ex

    elif reactants_products_are_inputs:

        try:
            reactants = process_particles_list(kwargs["reactants"])
            products = process_particles_list(kwargs["products"])
        except TypeError as t:
            raise TypeError(input_err_msg) from t
        except Exception as e:
            raise AtomicError(errmsg) from e

    if total_baryon_number(reactants) != total_baryon_number(products):
        raise AtomicError(
            "The baryon number is not conserved for "
            f"reactants = {reactants} and products = {products}."
        )

    if total_charge(reactants) != total_charge(products):
        raise AtomicError(
            "Total charge is not conserved for reactants = "
            f"{reactants} and products = {products}."
        )

    released_energy = add_mass_energy(reactants) - add_mass_energy(products)

    return released_energy
