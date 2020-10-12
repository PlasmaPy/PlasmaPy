__all__ = ["two_fluid_dispersion_solution"]

import astropy.units as u
import numpy as np

from astropy.constants.si import c, e, k_B, m_e, m_p, mu_0

import plasmapy.formulary.parameters as pfp

from plasmapy.utils.decorators import validate_quantities


@validate_quantities(
    n={"can_be_negative": False},
    B={"can_be_negative": False},
    T_i={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    theta={"can_be_negative":True},
    m_i={"can_be_negative":False},
    m_e={"can_be_negative":False}
)
def two_fluid_dispersion_solution(n: u.m ** -3, B: u.T, T_i: u.K, T_e: u.K, k: u.m ** -1, theta: u.deg = 45 * u.deg, z=1, gamma_e=1, gamma_i=3) :

    r"""
    Computes the wave frequency in the low frequency regime, corresponding to the plasma dispersion relation.

    **Aliases:** `tfds_`

    Parameters
    ----------
    n : ~astropy.units.Quantity
        Number density.
    z : float or integer
        Average ionization number which defaults to 1 for protons.
    B : ~astropy.units.Quantity
        Magnetic field
    T_i : ~astropy.units.Quantity
        The ion temperature
    T_e : ~astropy.units.Quantity
        The electron temperature
    thetha : ~astropy.units.Quantity
        angle of propagation
    k : ~astropy.units.Quantity
        Wave vector
    m_e : ~astropy.units.Quantity
        Mass of the electron
    m_i : ~astropy.units.Quantity
        Mass of the ion
    gamma_e : float or int
        The adiabatic index for electrons, which defaults to 1.  This
        value assumes that the electrons are able to equalize their
        temperature rapidly enough that the electrons are effectively
        isothermal.
    gamma_i : float or int
        The adiabatic index for ions, which defaults to 3. This value
        assumes that ion motion has only one degree of freedom, namely
        along magnetic field lines.

    Examples
    --------
    #TODO: Add some examples here

    Returns
    -------
    omega : ~astropy.units.Quantiry
        Wave frequency, in 1/s
    v_ph : ~astropy.units.Quantity
        Phase velocity of the wave, in m/s
    
    Notes
    -----
    Gets the solution for wave dispersion relation based on equation 38 of
    Bellan2012 (doi:10.1029/2012JA017856)

    .. math::

        \frac{\omega}{\omega_{ci}} = \sqrt(2\Lambda \sqrt(-\frac{p}{3}) cos \left( \frac{1}{3} cos^{-1}\left( \frac{3q}{2p} \sqrt(-\frac{3}{p} \right) - \frac{2\pi}{3}j \right) + \frac{A\Lambda}{3})
    
    Where, 
        j = 0 ==> fast mode
        j = 1 ==> Alfven mode
        j = 2 ==> Acoustic mode

    References
    ----------
    - PM bellan, Improved basis set for low frequency plasma waves, 2012, JGR, 117, A12219, doi:10.1029/2012JA017856
    - TE Stringer, Low-frequency waves in an unbounded plasma, 1963, JNE, Part C, doi:10.1088/0368-3281/5/2/304

    """

    # Required derived parameters
    c_s = pfp.ion_sound_speed(T_e=T_e, T_i=T_i, n_e=z * n, gamma_e =
    gamma_e, gamma_i = gamma_i, ion='p+')
    v_A = pfp.Alfven_speed( B, n, ion)
    omega_ci = pfp.gyrofrequency(B=B, particle='p+', signed=False, Z=z)
    omega_pe = pfp.plasma_frequency(n=n, particle='e-', z_mean=z)

    alpha = (np.cos(theta)**2).value
    beta = (c_s**2/v_A**2).value
    Lambda = (k**2 * v_A**2/omega_ci**2).value

    Q = 1 + (k**2 * c**2/omega_pe**2).value

    A = (Q + Q**2 * beta + Q * alpha + alpha * Lambda)/Q**2
    B = alpha * (1 + 2 * Q * beta + Lambda * beta)/Q**2
    C = alpha**2 * beta/Q**2

    p = (3 * B - A**2)/3
    q = (9 * A * B - 2 * A**3 - 27 * C)/27


    keys = ['fast_mode', 'alfven_mode', 'acoustic_mode']
    omega    = dict.fromkeys(keys)


    for (j,key) in zip(range(3), keys):

        # The solution corresponding to equation 38
        omega[key] = omega_ci * np.sqrt(2 * Lambda * np.sqrt(-p/3)
        * np.cos(1/3 * np.arccos(3 * q/(2 * p) * np.sqrt(-3/p)) - 2 * np.pi/3
        * j) + Lambda * A/3)

    return omega

tfds_ = two_fluid_dispersion_solution
""" Alias to :func:`two_fluid_dispersion_solution`. """