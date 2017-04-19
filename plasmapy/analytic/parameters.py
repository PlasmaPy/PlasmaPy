import numpy as np
from astropy import units as u
from astropy.constants import m_p, m_e, c, mu0, k_B
from astropy.units import UnitsError


def ion_mass(ion): # temporary; this needs improvement
    IonMasses = {'H':m_p, 'D':2*m_p, 'He3':3*m_p, 'He':4*m_p}
    try:
        m_i = IonMasses[ion]
    except:
        raise KeyError("No mass available for ion = '"+ion+"'")
    return m_i


def atomic_number(element):
    

    """Calculate the Alfven speed of a plasma.
    
    Parameters
    ----------


    Returns
    -------


    Raises
    ------


    See also
    --------

    
    Notes
    -----



    References
    ----------
    
    
    Examples
    --------

    >>>
    
    """






def Alfven_speed(B, density, element="H"):
    """Calculate the Alfven speed of a two-component plasma.
    
    Parameters
    ----------
    B : Astropy Quantity with units of magnetic field strength
    Magnetic field strength.

    density : Astropy Quantity with units of 1 / m**3 or kg / m**3
    Either the ion number density or the plasma mass density.

    element : string
    Symbol of ion element such as 'H', 'He'
    
    Returns
    -------
    V_A : Astropy Quantity with units of velocity    
    Alfven velocity.

    Raises
    ------

    ValueError
        If the Alfven velocity exceeds the speed of light

    Notes
    -----
    The Alfven velocity :math:`V_A` is the typical propagation speed
    of magnetic disturbances in a plasma, and is given by:

    .. math:: 
    V_A = \frac{B}{\sqrt{\mu_0\rho}}

    where the mass density is :math:`\rho \approx n_i m_i`.
    
    Examples
    --------

    >>> from astropy import units as u
    >>> import plasmapy
    >>> B = 140 * u.G
    >>> n = 5e19 * u.m**-3
    >>> rho = n * plasmapy.constants.m_i
    >>> element = 'H'
    >>> V_A1 = Alfven_speed(B,n,element)
    >>> print(V_A1)
    43185.62686969444 m / s
    >>> V_A2 = Alfven_speed(B,rho,element)
    >>> print(V_A2)
    43185.62686969444 m / s
    >>> print(V_A2.cgs)
    4318562.686969444 cm / s
    >>> print(V_A2.to(u.km/u.s))
    43.18562686969444 km / s

    """

    try:
        if density.si.unit=="m**-3":
            m_i = ion_mass(element)
            rho = density*m_i
        elif density.si.unit=="kg / m3":
            rho = density
    except:
        raise UnitsError("Second argument must have units of a number density or mass density")

    V_A = np.abs(B)/np.sqrt(mu0*rho)

    if V_A >= c:
        raise ValueError("Alfven velocity is greater than speed of light")

    return V_A.si


@u.quantity_input
def electron_gyrofrequency(B: u.T):
    """Calculate the electron gyrofrequency in units of radians per second"""
    omega_ce = e*B/m_e
    return omega_ce.si


def ion_gyrofrequency(B, ion="H"):
    m_i = ion_mass['H']
    try:
        if B.si.unit=="T":
            omega_ci = e*B/m_i
    except:
        raise UnitsError
    return omega_ci.si


def electron_plasma_frequency(n_e):
    """Calculates the electron plasma frequency"""
    omega_pe = e*np.sqrt(n_e/(eps0*m_e))
    return omega_pe.si


def ion_plasma_frequency(n_i,Z=1,ion='H'):
    """Calculates the ion plasma frequency"""
    m_i = ion_mass(ion)
    omega_pi = Z*e*np.sqrt(n_i/(eps0*m_e))
    return omega_pi.si


def Debye_length(T_e,n_e):
    lambda_D = e*np.sqrt(eps*k_B*T_e/n_e)
    return lambda_D.si

def Debye_number(T_e,n_e):
    lambda_D = Debye_length(T_e,n_e)
    N_D = (4/3)*pi*lambda_D**3
    return N_D.si


def ion_inertial_length(n_i,Z=1,ion='H'):
    omega_pi = ion_plasma_frequency(n_i,Z=Z,ion=ion)
    d_i = c/omega_pi
    return d_i.si


def electron_inertial_length(n_i,ion):
    omega_pe = electron_plasma_frequency(n_e)
    d_e = c/omega_pe
    return d_e.si


#def electron_thermal_velocity(T_e):
    


#def ion_thermal_velocity(T_i,ion):

def ion_sound_velocity(T_i,ion):
    m_i = ion_mass(ion)
    V_S = np.sqrt(k_B*T_i/m_i)
    return V_S.si


#def transverse_Spitzer_resistivity():
#def parallel_Spitzer_resistivity():
def magnetic_pressure(B):
    p_B = B**2/(2*mu0)
    return p_B

#def plasma_pressure(n,T):

#def plasma_beta(n,T,etc)


#def Coulomb_logarithm():
#def Reynolds_number():
#def magnetic_Reynolds_number():
#def Lundquist_number():
#def magnetic_Prandtl_number():


