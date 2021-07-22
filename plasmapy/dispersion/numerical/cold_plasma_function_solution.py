__all__ = ["cold_plasma_function_solution"]

import astropy.units as u
import numpy as np

import astropy.constants as const

from plasmapy.utils.decorators import validate_quantities

import pandas as pd
from sympy import var, solve



@validate_quantities(
    B={"can_be_negative": False},
    theta={"can_be_negative": False},
)
def cold_plasma_function_solution(
    *,
    B: u.T,
    k: u.rad / u.m,
    omega_p: 1/u.s,
    omega_e: 1/u.s,
    theta: u.rad,
):
    r"""
    Using the solution provided by Bellan 2012, calculate the analytical
    solution to the two fluid, low-frequency (:math:`\omega/kc \ll 1`) dispersion
    relation presented by Stringer 1963.  This dispersion relation also
    assummes a uniform magnetic field :math:`\mathbf{B_o}`, no D.C. electric
    field :math:`\mathbf{E_o}=0`, and quasi-neutrality.  For more information
    see the **Notes** section below.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to :math:`T`.
    k : `~astropy.units.Quantity`, single valued or 1-D array
        Wavenumber in units convertible to :math:`rad / m`.  Either single
        valued or 1-D array of length :math:`N`.
    omega_e : `~astropy.units.Quantity`
        Ion number density in units convertible to :math:`m^{-3}`.
    omega_p : `~astropy.units.Quantity`
        The electron temperature in units of :math:`K` or :math:`eV`.
    T_i : `~astropy.units.Quantity`
        The ion temperature in units of :math:`K` or :math:`eV`.
    theta : `~astropy.units.Quantity`, single valued or 1-D array
        The angle of propagation of the wave with respect to the magnetic field,
        :math:`\cos^{-1}(k_z / k)`, in units must be convertible to :math:`deg`.
        Either single valued or 1-D array of size :math:`M`.
    

    Returns
    -------
    omega : Dict[str, `~astropy.units.Quantity`]
        A dictionary of computed wave frequencies in units :math:`rad/s`.  The
        dictionary contains three keys: ``'fast_mode'`` for the fast mode,
        ``'alfven_mode'`` for the Alfvén mode, and ``'acoustic_mode'`` for the
        ion-acoustic mode.  The value for each key will be a :math:`N x M` array.

    Raises
    ------
    TypeError
        If applicable arguments are not instances of `~astropy.units.Quantity` or
        cannot be converted into one.

    TypeError
        If ``ion`` is not of type or convertible to `~plasmapy.particles.Particle`.

    TypeError
        If ``gamma_e``, ``gamma_i``, or``z_mean`` are not of type `int` or `float`.

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the expected
        units.

    ValueError
        If any of ``B``, ``k``, ``n_i``, ``T_e``, or ``T_i`` is negative.

    ValueError
        If ``k`` is negative or zero.

    ValueError
        If ``ion`` is not of category ion or element.

    ValueError
        If ``B``, ``n_i``, ``T_e``, or ``T_I`` are not single valued
        `astropy.units.Quantity` (i.e. an array).

    ValueError
        If ``k`` or ``theta`` are not single valued or a 1-D array.

    Warns
    -----
    : `~plasmapy.utils.exceptions.PhysicsWarning`
        When the computed wave frequencies violate the low-frequency
        (:math:`\omega/kc \ll 1`) assumption of the dispersion relation.

    Notes
    -----

    The complete dispersion equation presented by Springer 1963 [2]_ (equation 1
    of Bellan 2012 [1]_) is:

    .. math::
        \left( \cos^2 \theta - Q \frac{\omega^2}{k^2 {v_A}^2} \right) &
        \left[
            \left( \cos^2 \theta - \frac{\omega^2}{k^2 {c_s}^2} \right)
            - Q \frac{\omega^2}{k^2 {v_A}^2} \left(
                1 - \frac{\omega^2}{k^2 {c_s}^2}
            \right)
        \right] \\
            &= \left(1 - \frac{\omega^2}{k^2 {c_s}^2} \right)
              \frac{\omega^2}{{\omega_{ci}}^2} \cos^2 \theta

    where

    .. math::
        Q &= 1 + k^2 c^2/{\omega_{pe}}^2 \\
        \cos \theta &= \frac{k_z}{k} \\
        \mathbf{B_o} &= B_{o} \mathbf{\hat{z}}

    :math:`\omega` is the wave frequency, :math:`k` is the wavenumber, :math:`v_A`
    is the Alfvén velocity, :math:`c_s` is the sound speed, :math:`\omega_{ci}` is
    the ion gyrofrequency, and :math:`\omega_{pe}` is the electron plasma frequency.
    This relation does additionally assume low-frequency waves
    :math:`\omega/kc \ll 1`, no D.C. electric field :math:`\mathbf{E_o}=0` and
    quasi-neutrality.

    Following section 5 of Bellan 2012 [1]_ the exact roots of the above dispersion
    equation can be derived and expressed as one analytical solution (equation 38
    of Bellan 2012 [1]_):

    .. math::
        \frac{\omega}{\omega_{ci}} = \sqrt{
            2 \Lambda \sqrt{-\frac{P}{3}} \cos\left(
                \frac{1}{3} \cos^{-1}\left(
                    \frac{3q}{2p} \sqrt{-\frac{3}{p}}
                \right)
                - \frac{2 \pi}{3}j
            \right)
            + \frac{\Lambda A}{3}
        }

    where :math:`j = 0` represents the fast mode, :math:`j = 1` represents the
    Alfvén mode, and :math:`j = 2` represents the acoustic mode.  Additionally,

    .. math::
        p &= \frac{3B-A^2}{3} \; , \; q = \frac{9AB-2A^3-27C}{27} \\
        A &= \frac{Q + Q^2 \beta + Q \alpha + \alpha \Lambda}{Q^2} \;
            , \; B = \alpha \frac{1 + 2 Q \beta + \Lambda \beta}{Q^2} \;
            , \; C = \frac{\alpha^2 \beta}{Q^2} \\
        \alpha &= \cos^2 \theta \;
            , \; \beta = \left( \frac{c_s}{v_A}\right)^2 \;
            , \; \Lambda = \left( \frac{k v_{A}}{\omega_{ci}}\right)^2

    References
    ----------
    .. [1] PM Bellan, Improved basis set for low frequency plasma waves, 2012,
       JGR, 117, A12219, doi: `10.1029/2012JA017856
       <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2012JA017856>`_.

    .. [2] TE Stringer, Low-frequency waves in an unbounded plasma, 1963, JNE,
       Part C, doi: `10.1088/0368-3281/5/2/304
       <https://doi.org/10.1088/0368-3281/5/2/304>`_

    Examples
    --------
    >>> from astropy import units as u
    >>> from plasmapy.dispersion import two_fluid_dispersion
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

    # validate arguments
    for arg_name in ("B"):
        val = locals()[arg_name].squeeze()
        if val.shape != ():
            raise ValueError(
                f"Argument '{arg_name}' must a single value and not an array of "
                f"shape {val.shape}."
            )
        locals()[arg_name] = val

    # validate arguments
    for arg_name in ("omega_e", "omega_p"):
        if not isinstance(locals()[arg_name], (int, np.integer, float, np.floating)):
            raise TypeError(
                f"Expected int or float for argument '{arg_name}', but got "
                f"{type(locals()[arg_name])}."
            )

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
            f"Quantity, got array of shape {k.shape}."
        )


    #Start#
    
    plasma_proton = np.zeros(1)  
    plasma_electron = np.zeros(1) 
    component_frequency = np.zeros(2)
    
    lenr = 1

    #Set our variable to solve for
    w = var('omega')

    #Inital conditions
    S_Sum = 0
    S0 = 1
    D_Sum = 0
    P_Sum = 0
    P0 = 1
    
    ck = k*const.c

    plasma_proton[0] = omega_p
    plasma_electron[0] = omega_e
    component_frequency[0] = (abs(const.e)*B)/((const.m_e)*const.c) #Proton
    component_frequency[1] = ((const.e)*B)/((const.m_e)*const.c) #Elelctron

    #Calculate S, D and P
    for i in range(0,lenr):
        S_Sum  =+  (plasma_proton[i]**2)/(w**2 - component_frequency[0]**2) + (plasma_electron[i]**2)/(w**2 - component_frequency[1]**2)
    
    for i in range(0,lenr):
        P_Sum =+  (plasma_proton[i]**2)/(w**2) + (plasma_electron[i]**2)/(w**2)
    
    for i in range(0,lenr):
        D_Sum =+ (plasma_proton[i]**2/w)*(component_frequency[0]/(w**2-component_frequency[0]**2)) + (plasma_electron[i]**2/w)*(component_frequency[1]/(w**2-component_frequency[1]**2))
    
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

    eq = A*((ck/w)**4) - B*((ck/w)**2) + C

    sol = solve(eq, w)
    
    sol_len = len(sol)
    
    for i in range(0, sol_len):
        sol[i] = round(sol[i], 4)
    
    dummy_set = set(sol)
    final_list = list(dummy_set)
    final_list.sort()
    print(final_list)

    
    #End#    


    return final_list

inputs = {
   "k": 0.01 * u.rad / u.m,
   "theta": 30 * u.deg,
   "B": 8.3e-9 * u.T,
   "omega_p": 1.6e6 * 1/u.s,
   "omega_e": 4.0e5 * 1/u.s,
}
