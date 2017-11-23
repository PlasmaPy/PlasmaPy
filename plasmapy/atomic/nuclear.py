"""Functions that are related to nuclear reactions."""

from astropy import units, constants
import re
from .atomic import (isotope_symbol, mass_number, isotope_mass, atomic_number,
                     charge_state, _is_neutron, _is_electron, _is_positron)


def nuclear_binding_energy(argument, mass_numb=None):
    r"""Returns the nuclear binding energy associated with an isotope.

    Parameters
    ----------
    argument: string or integer
        A string representing an element or isotope, or an integer
        representing the atomic number of an element.

    mass_numb: integer, optional
        The mass number of an isotope, which is required if and only
        if the first argument can only be used to determine the
        element and not the isotope.

    Returns
    -------
    binding_energy: Quantity
        The binding energy of the nucleus in units of Joules.

    See also
    --------
    nuclear_reaction_energy : Returns the change in binding energy
        during nuclear fusion or fission reactions.

    Examples
    --------
    >>> nuclear_binding_energy('Fe-56')
    <Quantity 7.674003137600576e-11 J>
    >>> nuclear_binding_energy(26, 56)
    <Quantity 7.674003137600576e-11 J>
    >>> nuclear_binding_energy('p')  # proton
    <Quantity 0.0 J>

    >>> from astropy import units as u
    >>> before = nuclear_binding_energy("D") + nuclear_binding_energy("T")
    >>> after = nuclear_binding_energy("alpha")
    >>> (after - before).to(u.MeV)  # released energy from D + T --> alpha + n
    <Quantity 17.58929687252852 MeV>

    """

    if _is_neutron(argument) and mass_numb is None or mass_numb == 1:
        return 0.0 * units.J

    isotopic_symbol = isotope_symbol(argument, mass_numb)

    isotopic_mass = isotope_mass(isotopic_symbol)
    number_of_protons = atomic_number(argument)

    if mass_numb is None:
        mass_numb = mass_number(argument)
    number_of_neutrons = mass_numb - number_of_protons

    if number_of_protons == 1 and number_of_neutrons == 0:
        binding_energy = 0.0 * units.J
    else:
        mass_of_nucleons = (number_of_protons * constants.m_p +
                            number_of_neutrons * constants.m_n)
        mass_defect = mass_of_nucleons - isotopic_mass
        binding_energy = mass_defect * constants.c**2

    return binding_energy.to(units.J)


def nuclear_reaction_energy(*args, **kwargs):
    r"""Returns the released energy from a nuclear reaction.

    Parameters
    ----------
    reaction: string, optional
        A string representing the reaction, like "D + T --> alpha + n"

    reactants: list, optional
        A list of strings representing the reactants of a nuclear
        reaction

    products: list, optional
        A list of strings representing the products of a nuclear
        reaction

    Returns
    -------
    energy: Quantity
        The change in nuclear binding energy, which will be positive
        if the reaction releases energy and negative if the reaction
        is energetically unfavorable.

    Raises
    ------
    ValueError:
        If the reaction is not valid, there is insufficient
        information to determine an isotope, the number of nucleons is
        not conserved, or the charge is not conserved.

    TypeError:
        If the positional input for reaction is not a string.

    See also
    --------
    nuclear_binding_energy : finds the binding energy of an isotope

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
    <Quantity 2.8181160225476198e-12 J>
    >>> triple_alpha1 = 'alpha + He-4 --> Be-8'
    >>> triple_alpha2 = 'Be-8 + alpha --> carbon-12'
    >>> energy_triplealpha1 = nuclear_reaction_energy(triple_alpha1)
    >>> energy_triplealpha2 = nuclear_reaction_energy(triple_alpha2)
    >>> print(energy_triplealpha1, energy_triplealpha2)
    -1.4714307834595232e-14 J 1.1802573526724632e-12 J
    >>> energy_triplealpha2.to(u.MeV)
    <Quantity 7.366587037595994 MeV>
    >>> nuclear_reaction_energy(reactants=['n'], products=['p', 'e-'])
    <Quantity 1.25343510874046e-13 J>

    """

    def _get_species(side_of_reaction):
        r"""Parse a side of a reaction to get a list of the isotopes
        and other participants in the reaction."""
        species = []
        split_list = re.split(' \+ ', side_of_reaction)
        for item in split_list:
            item = item.strip()
            try:
                if _is_electron(item):
                    species.append('e-')
                elif _is_positron(item):
                    species.append('e+')
                else:
                    isotope = isotope_symbol(item)
                    species.append(isotope)
            except Exception:
                try:
                    multiplier_string = ''
                    while item[0].isdigit():
                        multiplier_string += item[0]
                        item = item[1:]
                        isotope = isotope_symbol(item)
                    for i in range(0, int(multiplier_string)):
                        species.append(isotope)
                except Exception:
                    raise ValueError("'"+str(item)+"' is an invalid isotope")
        return species

    def _nucleon_number(species_list):
        r"""Finds the total number of nucleons in a list of species."""
        nucleon_number = 0
        for species in species_list:
            if not (_is_electron(species) or _is_positron(species)):
                nucleon_number += mass_number(species)
        return nucleon_number

    def _total_charge(species_list):
        r"""Finds the total integer charge in a list of species,
        excluding bound electrons."""
        total_charge = 0
        for species in species_list:
            if _is_electron(species):
                total_charge -= 1
            elif _is_positron(species):
                total_charge += 1
            elif not _is_neutron(species):
                total_charge += atomic_number(species)
        return total_charge

    def _mass_energy(species_list):
        r"""Finds the total mass energy from a list of species."""
        total_mass = 0.0*units.kg
        for species in species_list:
            if _is_electron(species) or _is_positron(species):
                total_mass += constants.m_e
            elif _is_neutron(species):
                total_mass += constants.m_n
            else:
                total_mass += isotope_mass(species)
        return (total_mass * constants.c**2).to(units.J)

    input_err_msg = ("The inputs to nuclear_reaction_energy must be either a "
                     "a string representing a nuclear reaction (e.g., "
                     "'D + T -> He-4 + n') or the keywords 'reactants' and "
                     "'products' as lists with the nucleons or particles "
                     "involved in the reaction.")

    if kwargs and not args and len(kwargs) == 2:  # keyword inputs

        try:
            reactants = list(kwargs['reactants'])
            products = list(kwargs['products'])
        except:
            raise ValueError(input_err_msg)

    elif args and not kwargs and len(args) == 1:  # reaction string input

        try:
            LHS, RHS = re.split('-+>', args[0])
        except TypeError:
            raise TypeError("The left and right hand sides of the reaction "
                            "should be separated by '-->'")

        try:
            reactants = _get_species(LHS)
            products = _get_species(RHS)
        except Exception:
            raise ValueError(input_err_msg)

    else:
        raise ValueError(input_err_msg)

    try:
        if _nucleon_number(reactants) != _nucleon_number(products):
            raise ValueError("The number of nucleons is not conserved")

        if _total_charge(reactants) != _total_charge(products):
            raise ValueError("Total charge is not conserved")

        released_energy = _mass_energy(reactants) - _mass_energy(products)

    except Exception:
        raise ValueError("Invalid reactant(s) and/or product(s)")

    return released_energy.to(units.J)
