"""Functions that are related to nuclear reactions."""

from astropy import units, constants
import re
from .atomic import (isotope_symbol, mass_number, isotope_mass, atomic_number,
                     _is_neutron)


def nuclear_binding_energy(argument, mass_numb=None):
    """Returns the nuclear binding energy associated with an isotope.

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


def nuclear_reaction_energy(reaction):
    """Returns the released energy from a nuclear fusion reaction.

    Parameters
    ----------
    reaction: string
        A string representing the reaction, like "D + T --> alpha + n"

    Returns
    -------
    energy: Quantity
        The change in nuclear binding energy, which will be positive
        if the reaction releases energy and negative if the reaction
        is energetically unfavorable.

    Raises
    ------
    ValueError:
        If the input is not a valid reaction, there is insufficient
        information to determine an isotope, or if the number of
        nucleons is not conserved during the reaction.

    TypeError:
        If the input is not a string.

    See also
    --------
    nuclear_binding_energy : finds the binding energy of an isotope

    Examples
    --------
    >>> from astropy import units as u
    >>> nuclear_reaction_energy("D + T --> alpha + n")
    <Quantity 17.58929687252852 MeV>
    >>> triple_alpha1 = 'alpha + He-4 --> Be-8'
    >>> triple_alpha2 = 'Be-8 + alpha --> carbon-12'
    >>> energy_triplealpha1 = nuclear_reaction_energy(triple_alpha1)
    >>> energy_triplealpha2 = nuclear_reaction_energy(triple_alpha2)
    >>> print(energy_triplealpha1, '\\n', energy_triplealpha2)
    -0.0918394866304908 MeV 
     7.3665870375939875 MeV
    >>> energy_triplealpha1.to(u.J)
    <Quantity -1.4714307834564652e-14 J>
    >>> energy_triplealpha2.cgs
    <Quantity 1.1802573526721418e-05 erg>
    >>> nuclear_reaction_energy('alpha + alpha --> 2alpha')
    <Quantity 0.0 MeV>

    """

    def _get_isotopes_list(side):
        """Parse a side of a reaction to get a list of the isotopes."""
        pre_list = re.split(' \+ ', side)
        isotopes_list = []
        for item in pre_list:
            item = item.strip()
            try:
                isotope = isotope_symbol(item)
                isotopes_list.append(isotope)
            except Exception:
                try:
                    multiplier_string = ''
                    while item[0].isdigit():
                        multiplier_string += item[0]
                        item = item[1:]
                    isotope = isotope_symbol(item)
                    for i in range(0, int(multiplier_string)):
                        isotopes_list.append(isotope)
                except Exception:
                    raise
        return isotopes_list

    def _mass_number_of_list(isotopes_list):
        """Find the total number of nucleons in a list of isotopes."""
        mass_numb = 0
        for isotope in isotopes_list:
            mass_numb += mass_number(isotope)
        return mass_numb

    def _add_binding_energies(isotopes_list):
        """Finds the total binding energy from a list of isotopes."""
        total_binding_energy = 0.0*units.MeV
        for isotope in isotopes_list:
            total_binding_energy += nuclear_binding_energy(isotope)
        return total_binding_energy

    if not isinstance(reaction, str):
        raise TypeError("The input of nuclear_reaction_energy must be a string"
                        " representing a reaction (e.g., 'D + T --> He + n')")

    try:
        LHS, RHS = re.split('-+>', reaction)
    except Exception:
        raise ValueError("The left and right hand sides of the reaction "
                         "should be separated by '-->'")

    reactants = _get_isotopes_list(LHS)
    products = _get_isotopes_list(RHS)

    mass_num_reactants = _mass_number_of_list(reactants)
    mass_num_products = _mass_number_of_list(products)

    if mass_num_reactants != mass_num_products:
        raise ValueError("Mass numbers on LHS and RHS do not match for "
                         "reaction " + reaction)

    binding_energy_before = _add_binding_energies(reactants)
    binding_energy_after = _add_binding_energies(products)
    energy = binding_energy_after - binding_energy_before

    return energy.to(units.MeV)
