__all__ = ["cold_plasma_function_solution"]

import astropy.units as u
import numpy as np

import astropy.constants as const

from plasmapy.utils.decorators import validate_quantities
from plasmapy.formulary import parameters as pfp

from sympy import Symbol, solve


@validate_quantities(
    B={"can_be_negative": False},
    theta={"can_be_negative": False},
)

def cold_plasma_function_solution(
    *,
    B: u.T,
    k: u.rad/u.m,
    omega_p: u.rad/u.s,
    omega_e:u.rad/u.s,
    theta: u.rad,
):
 
    r"""
    Using the solution provided by Bellan 2012, calculate the numerical
    solution to the Stix, cold plasma method (:math:`\omega`) dispersion
    relation presented by Stix 1992.  This dispersion relation also
    assummes a uniform magnetic field :math:`\mathbf{B_o}`, theta is angle 
    between the magnetic field and the normal surface of the wave vector, 
    and k is the wave vector in direction of propagtion normal to the surface.
    For more information see the **Notes** section below.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to :math:`T`.
    k : `~astropy.units.Quantity`, single valued or 1-D array
        Wavenumber in units convertible to :math:`rad / m`.  Either single
        valued or 1-D array of length :math:`N`.
    omega_e : `~astropy.units.Quantity`
        Electron frequency in units convertible to :math:`s^{-1}`.
    omega_p : `~astropy.units.Quantity`
        Plasma frequency in units convertible to :math:`s^{-1}`.
    theta : `~astropy.units.Quantity`, single valued or 1-D array
        The angle of propagation of the wave with respect to the magnetic field,
        :math:`\cos^{-1}(k_z / k)`, in units must be convertible to :math:`deg`.
        Either single valued or 1-D array of size :math:`M`.
    

    Returns
    -------
    omega : array[float, `~astropy.units.Quantity`]
        An array of computed wave frequencies in units :math:`rad/s`.  

    Raises
    ------
    TypeError
        If applicable arguments are not instances of `~astropy.units.Quantity` or
        cannot be converted into one.

    TypeError
        If ``omega_e`` or ``omega_p``  are not of type `int` or `float`.

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the expected
        units.

    ValueError
        If any of ``B``, ``k`` or ``theta`` is negative.

    ValueError
        If ``k`` is negative or zero.

    ValueError
        If ``B`` is not single valued
        `astropy.units.Quantity` (i.e. an array).

    ValueError
        If ``k`` or ``theta`` are not single valued or a 1-D array.


    Notes
    -----

    The cold plasma waves equation presented by Strix 1992 [2]_ (equation 9
    of Bellan 2012 [1]_) is:

    .. math::
       #---

    where

    .. math::
       #---

    Following section 5 of Bellan 2012 [1]_ the exact roots of the above dispersion
    equation can be derived and expressed as one analytical solution (equation 38
    of Bellan 2012 [1]_):

    .. math::
        #-----

    where :math:`j = 0` represents the fast mode, :math:`j = 1` represents the
    Alfvén mode, and :math:`j = 2` represents the acoustic mode.  Additionally,

    .. math::
        #----

    References
    ----------
    .. [1] PM Bellan, Improved basis set for low frequency plasma waves, 2012,
       JGR, 117, A12219, doi: `10.1029/2012JA017856
       <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2012JA017856>`_.

    .. [2] TH Strix, 1992, Waves in Plasmas, Illustrated, 
        Springer Science & Business Media, 1992, New York 
       Part C, doi: `10.1088/0368-3281/5/2/304
       <https://doi.org/10.1088/0368-3281/5/2/304>`_

    Examples
    --------
    >>> from astropy import units as u
    >>> from plasmapy.dispersion import cold_plasma_function_solution
    >>> inputs = {
    ...     "k": 0.01 * u.rad / u.m,
    ...     "theta": 30 * u.deg,
    ...     "B": 8.3e-9 * u.T,
    ...     "n_i": 5e6 * u.m ** -3,
    ...     "T_e": 1.6e6 * u.K,
    ...     "T_i": 4.0e5 * u.K,
    ...     "ion": "p+",
    ... }
    >>> omegas = two_fluid_dispersion_solution(**inputs)
    >>> omegas
    {'fast_mode': <Quantity 1520.57... rad / s>,
     'alfven_mode': <Quantity 1261.75... rad / s>,
     'acoustic_mode': <Quantity 0.688152... rad / s>}

    >>> inputs = {
    ...     "k": [1e-7, 2e-7] * u.rad / u.m,
    ...     "theta": [10, 20] * u.deg,
    ...     "B": 8.3e-9 * u.T,
    ...     "n_i": 5e6 * u.m ** -3,
    ...     "T_e": 1.6e6 * u.K,
    ...     "T_i": 4.0e5 * u.K,
    ...     "ion": "He+",
    ... }
    >>> omegas = two_fluid_dispersion_solution(**inputs)
    >>> omegas['fast_mode']
    <Quantity [[0.00767..., 0.00779... ],
               [0.01534..., 0.01558...]] rad / s>
    """
    print('test')
    # validate arguments
    for arg_name in ("B"):
        val = locals()[arg_name].squeeze()
        if not (val.dim ==0)
            raise ValueError(
                f"Argument '{arg_name}' must a float or an integer "
                f"shape {val.shape}."
            )
        locals()[arg_name] = val


    #validate arguments
    omega_e = omega_e.squeeze()
    omega_p = omega_p.squeeze()
    if not (omega_e.ndim == 0 or omega_e.ndim == 1):
        raise ValueError(
            f"Argument 'omega_e' needs to be a single valued or 1D array astropy Quantity,"
            f" got array of shape {omega_e.shape}."
        )
    if not (omega_p.ndim == 0 or omega_p.ndim == 1):for 
        raise ValueError(
            f"Argument 'omega_p' needs to be a single valued or 1D array astropy Quantity,"
            f" got array of shape {omega_p.shape}."
        )
   

    
    #for arg_name in ("omega_e", "omega_p"):
    #    print(arg_name)
    #    if not isinstance(locals()[arg_name], (int, np.integer, float, np.floating)):
    #        raise TypeError(
    #            f"Expected int or float for argument '{arg_name}', but got "
    #            f"{type(locals()[arg_name])}."
    #        )

    # validate argument k
    k = k.squeeze()
    if not (k.ndim == 0 or k.ndim == 1):
        raise ValueError(
            f"Argument 'k' needs to be a single valued or 1D array astropy Quantity,"
            f" got array of shape {k.shape}."
        )
    if np.any(k <= 0):
        raise ValueError("Argument 'k' can not be a or have negative values.")

    # validate argument theta
    theta = theta.squeeze()
    theta = theta.to(u.radian)
    if not (theta.ndim == 0 or theta.ndim == 1):
        raise ValueError(
            f"Argument 'theta' needs to be a single valued or 1D array astropy "
            f"Quantity, got array of shape {theta.shape}."
        )
    
    #Set our variable to solve for
    w = Symbol('omega')*u.rad/u.s

    ck = k*const.c

    component_frequency = np.tile(0*u.rad/u.s, 2)
    
    print(B)
    print('file sync works')
    component_frequency[0] = pfp.gyrofrequency(B=B, particle='H+', signed=False)
    print(component_frequency)
    component_frequency[1] = pfp.gyrofrequency(B=B, particle='e-', signed=True)

    
    if k.ndim == 0:
        k_int = True
    elif k.ndim == 1:
        k_int = False
    else:
        raise ValueError(
            f"Error in input k can only be a float or an array of floats."
            f" k of shape {k.shape}."
        )
    lengths = []
    if omega_p.ndim == 0 and omega_e.ndim == 0:
        omega_int = True
        lengths.append(1)
        lengths.append(1)
    elif omega_p.ndim == 1 and omega_e.ndim == 1:
        omega_int = False
        lengths.append(len(omega_p))
        lengths.append(len(omega_e))
    else:
        raise ValueError(
            f"Error in omega input, input dimensions must be the same and input can only be a float or an array of floats."
            f"omega_p of shape {omega_p.shape}, omega_e of shape {omega_e.shape}."
        )
    
    sum_len = min(lengths)
    print(sum_len)
    plasma_proton = np.tile(0*u.rad/u.s, sum_len)
    plasma_electron = np.tile(0*u.rad/u.s, sum_len)

    #Inital conditions
    S_Sum = np.tile(0*u.rad/u.s, sum_len)
    S_Sum[0] = 1*u.rad/u.s
    D_Sum = np.tile(0*u.rad/u.s, sum_len)
    P_Sum = np.tile(0*u.rad/u.s, sum_len)
    P_Sum[0] = 1*u.rad/u.s

    if omega_int == False:  
        for i in range(sum_len):
            plasma_proton[i] = omega_p[i]
            plasma_electron[i] = omega_e[i]
    elif omega_int == True:
        plasma_proton[0] = omega_p
        plasma_electron[0] = omega_e
    else:
        raise ValueError(
            f"Error in omega input, omega_int has no been set correctly."
            f"omega_int of shape {omega_int.shape}."
        )
     

    #Calculate S, D and P

    prot_w = (component_frequency[0]**2)
    print(prot_w)
    elec_w = (component_frequency[1]**2)
    print(w**2, prot_w)
        
    for i in range(sum_len):
        S_Sum  =+  (plasma_proton[i]**2)/(w**2 - prot_w) #+ (plasma_electron[i]**2)/(omega**2 - component_frequency[1]**2)
        print('b')
    
    #for i in range(0,sum_len):
    #    P_Sum =+  (plasma_proton[i]**2)/(w**2) + (plasma_electron[i]**2)/(w**2)
    #
    #for i in range(0,sum_len):
    #    D_Sum =+ (plasma_proton[i]**2/w)*(component_frequency[0]/(w**2-component_frequency[0]**2)) + (plasma_electron[i]**2/w)*(component_frequency[1]/(w**2-component_frequency[1]**2))
    
    #Define our values for SPD
    S = S0 - S_Sum
    P = P0 - P_Sum
    D = D_Sum
    
    #Define R and L
    R = S + D
    L = S - D
    
    #Set up eqation to solve
    
    A = S*(np.sin(theta)**2) + P*(np.cos(theta)**2)
    B = R*L*(np.sin(theta)**2) +P*S*(1 + (np.cos(theta)**2))
    C = P*R*L

    #eq = A*((ck/w)**4) - B*((ck/w)**2) + C

    #sol = solve(eq, w)
    
    #sol_len = len(sol)
    
    #for i in range(0, sol_len):
    #    sol[i] = round(sol[i], 4)
    
    dummy_set = set(sol)
    omega_final = list(dummy_set)
    omega_final.sort()
    print(omega_final)

    return 
    
    #End#    


    return omega_final

