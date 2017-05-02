from astropy import units as u, constants as const
from .elements import atomic_symbols_list, atomic_symbols_dict, Elements
from .isotopes import Isotopes

import warnings


def element_symbol(argument):
    """Returns the atomic symbol.

    Parameters
    ----------
    argument: integer or string
        An integer representing the atomic number, or a string containing a
        representation of an element or an isotope.

    Returns
    -------
    symbol: string
        The atomic symbol of the element or isotope.

    Raises
    ------
    ValueError:
        If the argument cannot be used to identify the element.

    See also
    --------
    isotope_symbol : returns isotope symbol instead of atomic symbol

    element_name : returns the name of an element

    Notes
    -----
    This function returns the symbol of the element rather than the symbol of
    an isotope.  For example, 'deuterium', 'T', or 'hydrogen-2' will yield 'H';
    'alpha' will yield 'He'; and 'iron-56' or 'Fe-56' will yield 'Fe'.

    This function is case insensitive when there is no ambiguity associated
    with case.  However, this function will return 'H' for hydrogen if the 
    argument is lower case 'p' for proton, and will return capital 'P' if the
    argument is 'P' for phosphorus.

    Examples
    --------
    >>> element_symbol('helium')
    'He'
    >>> element_symbol(42)
    'Mo'
    >>> element_symbol('D')
    'H'
    >>> element_symbol('C-13')
    'C'
    >>> element_symbol('alpha')
    'He'
    >>> element_symbol('79')
    'Au'

    """

#    if str(argument) == '0':
#        warnings.warn('Assuming atomic number of zero corresponds to neutron')
#        return 'n'

    if type(argument) == str:


        if '-' in argument:  # check for isotope notation
            dash_position = argument.find('-')
            argument = argument[:dash_position]
        elif argument in ['n', 'neutron', 0]:
            return 'n'
        elif argument == 'p' or argument.lower() == 'proton':
            return 'H'
        if atomic_symbols_dict.keys().__contains__(argument.lower()):
            symbol = atomic_symbols_dict[argument.lower()]
        elif atomic_symbols_list.__contains__(argument.capitalize()):
            symbol = argument.capitalize()
        elif argument in ['p', 'p+']:
            symbol = 'H'
        elif argument.lower() in ['d', 't', 'protium', 'deuterium', 'tritium']:
            symbol = 'H'
        elif 'alpha' == argument.lower():
            symbol = 'He'
        elif argument.isalnum():
            argument = int(argument)
        else:
            raise ValueError(argument+" is an invalid argument for "
                             "element_symbol")

    if type(argument) == int:
        if 0 <= argument <= 118:
            symbol = atomic_symbols_list[argument]
        else:
            raise ValueError(str(argument)+" is an invalid atomic number in "
                             "element_symbol")

    return symbol


def isotope_symbol(argument, mass_numb=None):
    """Returns the symbol representing an isotope.

    Parameters
    ----------
    argument: integer or string
        An integer representing the atomic number, or a string containing a
        representation of an element or an isotope.

    mass_numb: integer or string
        An integer or string representing the mass number of the isotope.

    Returns
    -------
    symbol: string
        The isotopic symbol. The result will generally be returned as something
        like 'He-4' or 'Au-197', but will return 'D' for deuterium and 'T' for
        tritium.

    Raises
    ------
    ValueError:
        Invalid argument.

    TypeError:
        Invalid type.

    Warning:
        Resulting isotope not in database, which means it is an uncommon or
        very rare isotope.

    See also
    --------
    element_symbol : returns atomic symbol instead of isotopic symbol

    Notes
    -----
    This function returns the symbol of the element rather than the symbol of
    an isotope.  For example, 'deuterium', 'T', or 'hydrogen-2' will yield 'H';
    'alpha' will yield 'He'; and 'iron-56' or 'Fe-56' will yield 'Fe'.

    Examples
    --------
    >>> isotope_symbol('He, 4')
    'He-4'
    >>> isotope_symbol(79, 197)
    'Au-197
    >>> isotope_symbol('hydrogen-2')
    'D'
    >>> isotope_symbol('carbon-13')
    'C-13'
    >>> isotope_symbol('alpha')
    'He-4'

    """

    if argument in ['n', 'neutron', 'Neutron', 0, '0']:
        return 'n'

    if type(mass_numb) == str:
        try:
            mass_numb = int(mass_numb)
        except:
            raise TypeError("The second argument in isotope_symbol is an "
                            "invalid mass number.")
    elif type(argument) == int and type(mass_numb) == int:
        if argument > mass_numb:
            raise ValueError("The first argument of isotope_symbol, which "
                             "represents the atomic number when it is an "
                             "integer, cannot exceed the second (keyword) "
                             "argument which represents the mass number.")

    atomic_symbol = element_symbol(argument)

    if mass_numb is None:
        if type(argument) == str:
            if '-' in argument:
                dash_position = argument.find('-')
                mass_numb = argument[dash_position+1:]
            elif argument.lower() in ['protium', 'proton', 'p', 'p+']:
                mass_numb = 1
            elif argument.lower() in ['d', 'deuterium', 'h-2', 'hydrogen-2']:
                mass_numb = '2'
            elif argument.lower() in ['t', 'tritium', 'h-3', 'hydrogen-3']:
                mass_numb = '3'
            elif argument.lower() == 'alpha':
                mass_numb = 4
        else:
            raise TypeError("Argument in isotope_symbol must be a string "
                            "representing an isotope if no mass number is " 
                            "provided.")

    isotope = atomic_symbol + '-' + str(mass_numb)

    if isotope == 'H-2':
        isotope = 'D'
    elif isotope == 'H-3':
        isotope = 'T'

    if atomic_number(isotope) > int(mass_numb):
        raise ValueError("The mass number (" + str(mass_numb) + ") "
                         "cannot exceed the atomic number (" + 
                         str(atomic_number(isotope)) + ") in isotope_symbol.")

#    if isotope not in Isotopes.keys():
#        warnings.warn("No data is available for isotope " + isotope + 
#                      ", so this isotope is either unknown or rare.")

    return isotope


def atomic_number(argument):
    """Returns the atomic number (the number of protons in an atom)
    from an atomic symbol or name.

    Parameters
    ----------
    argument: string
        A string represnting the symbol or name of an element or isotope.

    Returns
    -------
    atomic_number: integer
        An integer representing the atomic number of the element or isotope.

    See also
    --------
    mass_number : returns the mass number (the total number of protons and
        neutrons) of an isotope.

    Examples
    --------
    >>> atomic_number("H")
    1
    >>> atomic_number("tritium")
    1
    >>> atomic_number("alpha")
    2
    >>> atomic_number("oganesson")
    118

    """

    atomic_symbol = element_symbol(argument)
    atomic_numb = Elements[atomic_symbol]['atomic_number']
    return atomic_numb

def is_isotope_stable(isotope, mass_numb=None):
    """Returns true for stable isotopes and false otherwise."""

    stable_isotopes = ['H-1', 'D', 'He-3', 'He-4', 'Li-6', 'Li-7', 'Be-9', 
                       'B-10', 'B-11', 'C-12', 'C-13', 'N-14', 'N-15', 'O-16',
                       'O-17', 'O-18', 'F-19', 'Ne-20', 'Ne-21', 'Ne-22', 
                       'Na-23', 'Mg-24', 'Mg-25', 'Mg-26', 'Al-27', 'Si-28', 
                       'Si-29', 'Si-30', 'P-31', 'S-32', 'S-33', 'S-34', 
                       'S-36', 'Cl-35', 'Cl-37', 'Ar-36', 'Ar-38', 'Ar-40', 
                       'K-39', 'K-41', 'Ca-40', 'Ca-42', 'Ca-43', 'Ca-44', 
                       'Ca-46', 'Sc-45', 'Ti-46', 'Ti-47', 'Ti-48', 'Ti-49',
                       'Ti-50', 'V-51', 'Cr-50', 'Cr-52', 'Cr-53', 'Cr-54', 
                       'Mn-55', 'Fe-54', 'Fe-56', 'Fe-57', 'Fe-58', 'Co-59',
                       'Ni-58', 'Ni-60', 'Ni-61', 'Ni-62', 'Ni-64', 'Cu-63',
                       'Cu-65', 'Zn-64', 'Zn-66', 'Zn-67', 'Zn-68', 'Zn-70', 
                       'Ga-69', 'Ga-71', 'Ge-70', 'Ge-72', 'Ge-73', 'Ge-74',
                       'As-75', 'Se-74', 'Se-76', 'Se-77', 'Se-78', 'Se-80', 
                       'Br-79', 'Br-81', 'Kr-78', 'Kr-80', 'Kr-82', 'Kr-83',
                       'Kr-84', 'Kr-86', 'Rb-85', 'Sr-84', 'Sr-86', 'Sr-87', 
                       'Sr-88', 'Y-89', 'Zr-90', 'Zr-91', 'Zr-92', 'Zr-94', 
                       'Nb-93', 'Mo-92', 'Mo-94', 'Mo-95', 'Mo-96', 'Mo-97', 
                       'Mo-98', 'Ru-96', 'Ru-98', 'Ru-99', 'Ru-100', 'Ru-101', 
                       'Ru-102', 'Ru-104', 'Rh-103', 'Pd-102', 'Pd-104', 
                       'Pd-105', 'Pd-106', 'Pd-108', 'Pd-110', 'Ag-107', 
                       'Ag-109', 'Cd-106', 'Cd-108', 'Cd-110', 'Cd-111', 
                       'Cd-112', 'Cd-114', 'In-113', 'Sn-112', 'Sn-114', 
                       'Sn-115', 'Sn-116', 'Sn-117', 'Sn-118', 'Sn-119', 
                       'Sn-120', 'Sn-122', 'Sn-124', 'Sb-121', 'Sb-123', 
                       'Te-120', 'Te-122', 'Te-123', 'Te-124', 'Te-125', 
                       'Te-126', 'I-127', 'Xe-124', 'Xe-126', 'Xe-128', 
                       'Xe-129', 'Xe-130', 'Xe-131', 'Xe-132', 'Xe-134', 
                       'Cs-133', 'Ba-132', 'Ba-134', 'Ba-135', 'Ba-136', 
                       'Ba-137', 'Ba-138', 'La-139', 'Ce-136', 'Ce-138', 
                       'Ce-140', 'Ce-142', 'Pr-141', 'Nd-142', 'Nd-143', 
                       'Nd-145', 'Nd-146', 'Nd-148', 'Sm-144', 'Sm-149', 
                       'Sm-150', 'Sm-152', 'Sm-154', 'Eu-153', 'Gd-154', 
                       'Gd-155', 'Gd-156', 'Gd-157', 'Gd-158', 'Gd-160', 
                       'Tb-159', 'Dy-156', 'Dy-158', 'Dy-160', 'Dy-161', 
                       'Dy-162', 'Dy-163', 'Dy-164', 'Ho-165', 'Er-162', 
                       'Er-164', 'Er-166', 'Er-167', 'Er-168', 'Er-170', 
                       'Tm-169', 'Yb-168', 'Yb-170', 'Yb-171', 'Yb-172', 
                       'Yb-173', 'Yb-174', 'Yb-176', 'Lu-175', 'Hf-176', 
                       'Hf-177', 'Hf-178', 'Hf-179', 'Hf-180', 'Ta-180', 
                       'Ta-181', 'W-182', 'W-183', 'W-184', 'W-186', 'Re-185', 
                       'Os-184', 'Os-187', 'Os-188', 'Os-189', 'Os-190', 
                       'Os-192', 'Ir-191', 'Ir-193', 'Pt-192', 'Pt-194', 
                       'Pt-195', 'Pt-196', 'Pt-198', 'Au-197', 'Hg-196', 
                       'Hg-198', 'Hg-199', 'Hg-200', 'Hg-201', 'Hg-202', 
                       'Hg-204', 'Tl-203', 'Tl-205', 'Pb-204', 'Pb-206', 
                       'Pb-207', 'Pb-208']
    
    if isotope_symbol(isotope, mass_numb) in stable_isotopes:
        return True
    else:
        return False


def mass_number(isotope):
    """Get the mass number (the number of protons and neutrals) of an isotope.

    Parameters
    ----------
    isotope : string
        A string representing the symbol or name of an isotope.

    Returns
    -------
    mass_number : integer
       An integer representing the mass number of the isotope.

    Raises
    ------

    Warning:
        If the isotope is not in the list of common isotopes.

    See also
    --------
    atomic_number : returns the number of protons in an isotope or element

    Examples
    --------
    >>> mass_number("tritium")
    3

    """

    if isotope == 'n':
        return 1

    symbol = isotope_symbol(isotope)
    try:
        mass_numb = Isotopes[symbol]["mass_number"]
    except:
        if symbol == 'D':
            mass_numb = 2
        elif symbol == 'T':
            mass_numb = 3
        elif '-' in symbol:
            dash_position = symbol.find('-')
            mass_numb_string = symbol[dash_position+1:]
            try:
                mass_numb = int(mass_numb_string)
            except:
                raise ValueError("Cannot find mass number for isotope " +
                                 str(isotope))
    return mass_numb


def element_name(argument):
    """Returns the name of an element.

    Parameters
    ----------
    argument : string or integer
        String representing an element or isotope, or an integer representing
        the atomic number of an element.

    Returns
    -------
    name : string
        The name of the element.

    See also
    --------
    element_symbol : returns the atomic symbol

    isotope_symbol : returns the symbol of an isotope

    Examples
    --------
    >>> element_name("H")
    'hydrogen'
    >>> element_name("T")
    'hydrogen'
    >>> element_name("alpha")
    'helium'
    >>> element_name()
    >>> element_name(42)
    'molybdenum'
    >>> element_name("C-12")
    'carbon'

    """
    atomic_symbol = element_symbol(argument)
    name = Elements[atomic_symbol]["name"]
    return name


def standard_atomic_weight(argument):
    """Returns the standard (conventional) atomic weight of an element based on
    the relative abundances of isotopes in terrestrial environments.

    Parameters
    ----------
    argument: string or integer
        A string representing the symbol or name of an element or isotope, or
        an integer representing an element's atomic number.

    Returns
    -------
    atomic_weight: astropy.units.Quantity with units of astropy.units.u
        The standard atomic weight of an element based on values frmo

    Raises
    ------
    ValueError:
        If the argument cannot be used to identify an element, or if the
        argument represents an isotope.

    UserWarning:
        If the standard atomic weight is not provided for an element (e.g., if
        an element does not naturally occur).

    See also
    --------
    isotope_mass : returns the atomic mass of an isotope.

    ion_mass : returns the mass of an ion of an element or isotope, accounting
        for reduction of mass from the neutral state due to different numbers
        of electrons.

    Notes
    -----
    The relative abundances of different isotopes of an element sometimes vary
    naturally in different locations within the terrestrial environment.  The
    CIAAW provides ranges for these element, which include H, Li, B, C, N, O,
    Mg, Si, S, Cl, Br, Tl.  This function provides a single value from the
    CIAWW 2015 standard values when a single value is given, and the lower
    accuracy conventional value given by Meija et al. (2013,
    doi:10.1515/pac-2015-0305) for the elements where a range is given.

    Examples
    --------
    >>> from astropy import units
    >>> standard_atomic_weight("H")
    <Quantity 1.008 u>
    >>> standard_atomic_weight("H").to(units.kg)
    <Quantity 1.673823232368e-27 kg>  # accounts for small amount of deuterium
    >>> isotope_mass("H-1")
    <Quantity 1.6735326915759943e-27 kg>  # pure hydrogen-1
    >>> standard_atomic_weight(82)
    <Quantity 238.02891 u>
    >>> standard_atomic_weight("lead")
    <Quantity 238.02891 u>
    >>> standard_atomic_weight(118) is None
    UserWarning: No standard atomic weight is available for Og. Returning None.
    True

    """

    if type(argument) == int:
        try:
            atomic_symbol = element_symbol(argument)
        except:
            raise ValueError(str(argument)+' is an invalid atomic number')
    elif type(argument) == str:
        if argument == 'P':
            atomic_symbol = 'P'
        elif (not argument.isalpha() or argument.lower() in
              ['p', 'protium', 'd', 'deuterium', 't', 'tritium', 'alpha']):
            raise ValueError("Use isotope_mass to get atomic mass of isotopes")
        else:
            atomic_symbol = element_symbol(argument)
    atomic_weight = Elements[atomic_symbol]["atomic_mass"]

    if atomic_weight is None:
        raise ValueError("No standard atomic weight is available for " +
                         atomic_symbol)
    return atomic_weight


def isotope_mass(argument, mass_numb=None):
    """Return the mass of an isotope.

    Parameters
    ----------
    argument : string or integer
        A string representing an element or isotope or an integer representing
        an atomic number

    mass_numb : integer (optional)
        The mass number of the isotope.

    Returns
    -------
    isotope_mass : Quantity
        The atomic mass of a neutral atom of an isotope.

    Raises
    ------
    UserWarning:
        Redundant (but consistent) isotope information is provided.

    ValueError:
        Contradictory or insufficient isotope information is provided.

    See also
    --------
    standard_atomic_weight : returns atomic weight of an element based on
        terrestrial abundances of isotopes

    ion_mass : returns the mass of an ion of an element or isotope, accounting
        for loss of electrons

    Notes
    -----
    The masses of rare isotopes may be unavailable.

    Examples
    --------
    >>>from astropy import units as u
    >>>isotope_mass("H-1")
    <Quantity 1.00782503223 u>
    >>>isotope_mass("H-1").to(units.kg)
    <Quantity 1.6735326915759943e-27 kg>
    >>> isotope_mass("He", 4)
    <Quantity 4.00260325413 u>
    >>> isotope_mass(2, 4)

    """

    if argument == "alpha":
        raise ValueError("Use ion_mass for mass of an alpha particle")

    if type(argument) == str and type(mass_numb) is int:
        if '-' in argument or argument in ['D', 'T']:
            if mass_numb == mass_number(argument):
                warnings.warn("Redundant isotope information in isotope_mass")
            else:
                raise ValueError("Contradictory or insufficient information" +
                                 " relating to isotope and mass number\n" +
                                 "(first argument = " + str(argument) + "   " +
                                 "second argument = " + str(mass_numb) + ")")

    try:
        symbol = isotope_symbol(argument, mass_numb)
        atomic_mass = Isotopes[symbol]['atomic_mass']
    except:
        if symbol == 'n':
            return 1.00866491588*u.u
        raise ValueError("Isotope mass not available for " + str(argument))

    return atomic_mass


def ion_mass(argument, Z=1, mass_numb=None):
    """Returns the mass of an ion.

 by finding the atomic mass of a neutral
    isotope (when available) or the standard atomic mass of an element based
    on terrestrial abundances, and then accounting for the change in mass due
    ionization.

    Parameters
    ----------
    argument: string or integer
        A string representing an element or isotope, or an integer representing
        the atomic number of an element.

    Z: integer (optional)
        The ionization state of the ion (defaulting to a charge of Z = 1)

    mass_numb: integer (optional)
        The mass number of an isotope.

    Returns
    -------
    m_i: Quantity
        The mass of a single ion of the isotope or element with charge state Z.

    Raises
    ------
    ValueError
        If the ionization state exceeds the atomic number

    See also
    --------
    standard_atomic_mass : returns the conventional atomic mass of an element
        based on terrestrial values and assuming the atom is neutral.

    isotope_mass : returns the mass of an isotope (if available) assuming the
        atom is neutral.

    Notes
    -----
    This function in general finds the mass of an isotope (or the standard
    atomic weight based on terrestrial values if a unique isotope cannot be
    identified), and then substracts the mass of Z electrons.  If Z is not
    provided as an input, then this function assumes that the ion is singly
    ionized.

    Specific values are returns for protons, deuterons, tritons, alpha
    particles, and positrons.

    Calling ion_mass('H') does not return the mass of a proton but instea
    uses hydrogen's standard atomic weight based on terrestrial values.  To
    get the mass of a proton, use ion_mass('p').

    Examples
    --------
    >>> ion_mass('p')  # proton
    <Constant name='Proton mass' value=1.672621777e-27 uncertainty=7.4e-35
    unit='kg' reference='CODATA 2010'>
    >>> ion_mass('H') == ion_mass('p')
    False
    >>> ion_mass('P')  # phosphorus
    <Quantity 5.143222638917872e-26 kg>
    >>> ion_mass('He-4', 2)
    <Quantity 6.64465723e-27 kg>
    >>> ion_mass('T')
    <Quantity 3.343583719e-27 kg>
    >>> ion_mass(26, Z=1, mass_numb=56)
    <Quantity 9.288122788133088e-26 kg>
    >>> ion_mass('Fe-56')
    <Quantity 9.288122788133088e-26 kg>
    """

    if str(argument).lower() in ['e+', 'positron']:
        return const.m_e

    if atomic_number(argument) < Z:
        raise ValueError("The ionization state cannot exceed the atomic number")

    if argument == 'alpha' or element_symbol(argument) == 'He' and Z == 2:
        return 6.644657230e-27*u.kg
    elif argument in ['p', 'p+'] or str(argument).lower() in \
            ['proton', 'protium']:
        return const.m_p
    elif atomic_number(argument) == 1:
        if type(argument) == str and '-1' in str(argument):
            return const.m_p
        elif argument == 1 and mass_numb == 1:
            return const.m_p

    try:
        isotope_symb = isotope_symbol(argument, mass_numb)
        if isotope_symb == 'D' and Z == 1:
            return 3.343583719e-27 * u.kg
        elif isotope_symb == 'T' and Z == 1:
            return 5.007356665e-27 * u.kg
        atomic_mass = isotope_mass(isotope_symb)
    except:
        try:
            atomic_mass = standard_atomic_weight(argument)
        except:
            errormessage = "No isotope mass or standard atomic weight is " +\
                "available to get ion mass for " + str(argument)
            if type(mass_numb) is int:
                errormessage += " with mass number " + str(mass_numb)
            
            raise ValueError(errormessage)

    m_i = (atomic_mass - Z*const.m_e).to(u.kg)

    return m_i


def nuclear_binding_energy(argument, mass_numb=None):
    """Returns the nuclear binding energy associated with an isotope.

    Parameters
    ----------
    argument: string or integer
        A string representing an element or isotope, or an integer representing
        the atomic number of an element.

    mass_numb: integer
        The mass number of an isotope, which is required if and only if the
        first argument can only be used 

    Returns
    -------
    binding_energy: Quantity
        The binding energy of the nucleus.

    Raises
    ------
    ValueError
        If the isotope cannot be uniquely determined from the inputs.

    Examples
    --------
    >>> nuclear_binding_energy('Fe-56')
    <Quantity 7.674002331034521e-11 J>
    >>> nuclear_binding_energy(26, 56)
    <Quantity 7.674002331034521e-11 J>
    >>> nuclear_binding_energy('p')  # proton
    <Quantity 0.0 J>

    >>> from astropy import units as u
    >>> before = nuclear_binding_energy("D") + nuclear_binding_energy("T")
    >>> after = nuclear_binding_energy("alpha")
    >>> (after - before).to(u.MeV)  # released energy from D + T --> alpha + n
    <Quantity 17.589296207151556 MeV>  
    """

    if argument == 'n' and mass_numb is None or mass_numb == 1:
        return 0 * u.J

    isotopic_symbol = isotope_symbol(argument, mass_numb)
    isotopic_mass = isotope_mass(isotopic_symbol)
    number_of_protons = atomic_number(argument)

    if mass_numb is None:
        mass_numb = mass_number(argument)
    number_of_neutrons = mass_numb - number_of_protons
    
    if number_of_protons == 1 and number_of_neutrons == 0:
        binding_energy = 0 * u.J
    else:
        mass_of_nucleons = (number_of_protons * const.m_p + 
                            number_of_neutrons * const.m_n)
        mass_defect = mass_of_nucleons - isotopic_mass        
        binding_energy = mass_defect * const.c**2
    
    return binding_energy.to(u.J)

def energy_from_nuclear_reaction(reaction):
    """Returns the released energy from a nuclear fusion reaction.

    Parameters
    ----------
    reaction: string
        A string representing the reaction, like "D + T --> alpha + n"

    """
    
    import re, string

    def _get_isotopes_list(side):
        print(side)
        pre_list = re.split(' \+ ', side)
        print(pre_list)
        isotopes_list = []
        for item in pre_list:
            item = item.strip()
            try:
                if item == 'n':
                    symbol = 'n'
                else:
                    symbol = isotope_symbol(item)
                isotopes_list.append(symbol)
            except:
                try:
                    multiplier_string = ''
                    while item[0].isdigit():
                        multiplier_string += item[0]
                        item = item[1:]
                    symbol = isotope_symbol(item)
                    for i in range(0,int(multiplier_string)):
                        isotopes_list.append(symbol)
                except:
                    raise
        return isotopes_list

    def _mass_number_of_list(isotopes_list):
        mass_numb = 0
        for symbol in isotopes_list:
            mass_numb += mass_number(symbol)
        return mass_numb

    def _add_binding_energies(isotopes_list):
        total_binding_energy = 0.0*u.MeV
        for isotope in isotopes_list:
            total_binding_energy += nuclear_binding_energy(isotope)
        return total_binding_energy

    if type(reaction) != str:
        raise TypeError("The input of energy_from_nuclear_reaction "
                        "must be a string representing the reaction (e.g., "
                        "'D + T --> He + n')")

    try:
        LHS, RHS = re.split('-+>', reaction)
    except:
        raise ValueError("The left and right hand sides of the reaction "
                         "should be separated by '-->'")
    
    reactants = _get_isotopes_list(LHS)
    products = _get_isotopes_list(RHS)

    mass_num_reactants = _mass_number_of_list(reactants)
    mass_num_products = _mass_number_of_list(products)

    if mass_num_reactants != mass_num_products:
        raise ValueError("Mass numbers on LHS and RHS do not match for "
                         "reaction "+reaction)
 
    binding_energy_before = _add_binding_energies(reactants) 
    binding_energy_after = _add_binding_energies(products)
    energy = binding_energy_after - binding_energy_before
    return energy.to(u.MeV)
