__all__ = ["two_fluid_dispersion_solution"]

import astropy.units as u
import numpy as np

from astropy.constants.si import c, e, k_B, m_e, m_p, mu0

import plasmapy.formulary.parameters as pfp

from plasmapy.utils.decorators import validate_quantities


@validate_quantities(
    B={"can_be_negative": False},
    m_e={"can_be_negative": False},
    m_i={"can_be_negative": False},
    n={"can_be_negative": False},
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    T_i={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    theta={"can_be_negative": True},
)
def two_fluid_dispersion_solution(
    B: u.T,
    k: u.m ** -1,
    n: u.m ** -3,
    T_e: u.K,
    T_i: u.K,
    gamma_e=1,
    gamma_i=3,
    ion="p+",
    m_e: u.kg = m_e,
    m_i: u.kg = m_p,
    theta: u.deg = 45 * u.deg,
    z=1,
):
    r"""
    Return a dictionary dictionary of frequencies corresponding to the
    solutions of dispersion relation in low frequency regime.

    **Aliases:** `tfds_`

    Parameters
    ----------
    B : ~astropy.units.Quantity
        Magnetic field
    gamma_e : float or int
        The adiabatic index for electrons, which defaults to 1.  This
        value assumes that the electrons are able to equalize their
        temperature rapidly enough that the electrons are effectively
        isothermal.
    gamma_i : float or int
        The adiabatic index for ions, which defaults to 3. This value
        assumes that ion motion has only one degree of freedom, namely
        along magnetic field lines.
    ion : string, optional
        Representation of the ion species (e.g., `'p'` for protons,
        `'D+'` for deuterium, or 'He-4 +1' for singly ionized
        helium-4), which defaults to protons.  If no charge state
        information is provided, then the ions are assumed to be
        singly charged.
    k : ~astropy.units.Quantity
        Wave number
    m_e : ~astropy.units.Quantity
        Mass of negative ion which defaults to electron
    m_i : ~astropy.units.Quantity
        Mass of positive ion which defaults to proton
    n : ~astropy.units.Quantity
        Number density.
    T_e : ~astropy.units.Quantity
        The electron temperature
    T_i : ~astropy.units.Quantity
        The ion temperature
    theta : ~astropy.units.Quantity
        Angle of propagation which defaults to 45 degrees.
    z : float or integer
        Average ionization number which defaults to 1 for protons.

    Returns
    -------
    omega : ~astropy.units.Quantity
        A dictionary of Wave frequencies corresponding to three modes namely
        a) Ion-acoustic mode, b) Alfven mode and c) Fast mode, in 1/s units

    Notes
    -----
    Computes the solution for wave dispersion relation based on equation 38 of
    Bellan2012JGR (doi:10.1029/2012JA017856)

    .. math::

        \frac{\omega}{\omega_{ci}} = \sqrt(2\Lambda \sqrt(-\frac{p}{3}) cos
        \left( \frac{1}{3} cos^{-1}\left( \frac{3q}{2p} \sqrt(-\frac{3}{p}
        \right) - \frac{2\pi}{3}j \right) + \frac{A\Lambda}{3})

    Where,
        j = 0 ==> fast mode
        j = 1 ==> Alfven mode
        j = 2 ==> Acoustic mode

    The above equation is derived from the general wave equation in the low
    frequency regime, where both electrons and ions play significant role (in
    the high frequency regime ions do not play any significant role and the
    process is dominated by electron dynamics).

    The complete dispersion equation is thus written as (from equation (1) of
    Bellan2012JGR):

    .. math::
        \left( cos^2\theta - Q\frac{\omega^2}{k^2 {v_A}^2} \right) \left[
        \left( cos^2\theta - \frac{\omega^2}{k^2 {c_s}^2 \right) -
        Q\frac{\omega^2}{k^2 {v_A}^2 \left( 1 - \frac{\omega^2}{k^2 {c_s}^2
        \right) \right] = \left(1 - \frac{\omega^2}{k^2 {c_s}^2 \right)
        \frac{\omega^2}{{\omega_{ci}}^2}cos^2\theta

    Here,
    .. math::
        Q = 1 + k^2 c^2/{\omega_{pe}}^2
    \omega_{ci} is the proton gyrofrequency

    References
    ----------
    .. [Bellan2012JGR] PM bellan, Improved basis set for low frequency plasma
    waves, 2012, JGR, 117, A12219, doi:10.1029/2012JA017856
    .. [Stringer1963JNE] TE Stringer, Low-frequency waves in an unbounded
    plasma, 1963, JNE, Part C, doi:10.1088/0368-3281/5/2/304
    .. [Rogers2001PRL] Rogers, B. N.; Denton, R. E.; Drake, J. F. & Shay, M. A.
    Role of Dispersive Waves in Collisionless Magnetic Reconnection, prl, 2001,
    87, 195004

    Examples
    --------
    >>> from astropy import units as u
    >>> from plasmapy.dispersion import two_fluid_dispersion_solution as tfds
    >>> k = 0.01 * u.m ** -1
    >>> theta = 30 * u.deg
    >>> B = 8.3E-9 * u.T
    >>> T_e = 1.6e6 * u.K
    >>> T_i = 4.e5 * u.K
    >>> z = 1
    >>> omega = tfds(B, k, n, T_e, T_i, theta, z)
    {'fast_mode': <Quantity [[1520.57692532]] rad / s>,
    'alfven_mode': <Quantity [[1260.01546096]] rad / s>,
    'acoustic_mode': <Quantity [[0.68815159]] rad / s>}

    >>> k = np.linspace(10**-7, 10**-2, 1E4) * u.m ** -1
    >>> theta = np.linspace(5, 85, 100) * u.deg
    >>> n = 5 * u.cm ** -3
    >>> B = 8.3E-9 * u.T
    >>> T_e = 1.6e6 * u.K
    >>> T_i = 4.e5 * u.K
    >>> z = 1
    >>> c = 3.e8 * u.m/u.s
    >>> c_s = pfp.ion_sound_speed(T_e=T_e, T_i=T_i, n_e=z * n)
    >>> v_A = pfp.Alfven_speed( B, n, ion='p+')
    >>> omega_ci = pfp.gyrofrequency(B=B, particle='p+', signed=False, Z=z)
    >>> omega = two_fluid_dispersion_solution(n=n, B=B, T_e=T_e, T_i=T_i,
        theta=theta, z=z, k=k)
    >>> omega['fast_mode'][:,40]
        [0.016117629, 0.17733531, 0.33868854, … , 1520.3016, 1520.4535,
        1520.6055] rad/s

    """
    # Required derived parameters
    # Compute the ion sound speed using the function from
    # plasmapy.formulary.parameters
    c_s = pfp.ion_sound_speed(
        T_e=T_e, T_i=T_i, n_e=z * n, gamma_e=gamma_e, gamma_i=gamma_i, ion=ion
    )

    # Compute the ion Alfven speed using the function from
    # plasmapy.formulary.parameters
    v_A = pfp.Alfven_speed(B, n, ion=ion)

    # Compute the ion gyro frequency using the function from
    # plasmapy.formulary.parameters
    omega_ci = pfp.gyrofrequency(B=B, particle="p+", signed=False, Z=z)

    # Compute the electron plasma frequency using the function from
    # plasmapy.formulary.parameters
    omega_pe = pfp.plasma_frequency(n=n, particle="e-", z_mean=z)

    # Compute the dimensionless parameters corresponding to equation 32 of
    # Bellan2012JGR
    alpha = (np.cos(theta.to("rad")) ** 2).value
    beta = (c_s ** 2 / v_A ** 2).value

    # Create a meshgrid of direction of propagation and the wavenumber
    alphav, kv = np.meshgrid(alpha, k)

    Lambda = (kv ** 2 * v_A ** 2 / omega_ci ** 2).value

    # Compute the dimensionless parameters corresponding to equation 2 of
    # Bellan2012JGR
    Q = 1 + (kv ** 2 * c ** 2 / omega_pe ** 2).value

    # Compute the dimensionless parameters corresponding to equation 35 of
    # Bellan2012JGR
    A = (Q + Q ** 2 * beta + Q * alphav + alphav * Lambda) / Q ** 2
    B = alphav * (1 + 2 * Q * beta + Lambda * beta) / Q ** 2
    C = alphav ** 2 * beta / Q ** 2

    # Compute the dimensionless parameters corresponding to equation 36 of
    # Bellan2012JGR
    p = (3 * B - A ** 2) / 3
    q = (9 * A * B - 2 * A ** 3 - 27 * C) / 27

    # List out the three wave modes for which this function gives the
    # frequencies
    keys = ["fast_mode", "alfven_mode", "acoustic_mode"]

    # Create a dictionary with the wave mode names as its keys
    omega = dict.fromkeys(keys)

    # Compute the value of  omega for each key and for each value of wavenumber
    # and direction of propagation
    for (j, key) in zip(range(3), keys):
        # The solution corresponding to equation 38
        omega[key] = omega_ci * np.sqrt(
            2
            * Lambda
            * np.sqrt(-p / 3)
            * np.cos(
                1 / 3 * np.arccos(3 * q / (2 * p) * np.sqrt(-3 / p)) - 2 * np.pi / 3 * j
            )
            + Lambda * A / 3
        )

    return omega


tfds_ = two_fluid_dispersion_solution
""" Alias to :func:`two_fluid_dispersion_solution`. """
