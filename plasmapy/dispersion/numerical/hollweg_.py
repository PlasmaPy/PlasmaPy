"""
This module contains functionality for calculating various numerical
solutions to Hollweg's two fluid dispersion relation
"""

import astropy.units as u
import numpy as np
import warnings

from astropy.constants.si import c
from typing import Union

from plasmapy.formulary import parameters as pfp
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import ChargeError
from plasmapy.utils.decorators import validate_quantities
from plasmapy.utils.exceptions import PhysicsWarning


@validate_quantities(
    B={"can_be_negative": False},
    n_i={"can_be_negative": False},
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    T_i={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def hollweg(
    *,
    B: u.T,
    ion: Union[str, Particle],
    k: u.rad / u.m,
    n_i: u.m ** -3,
    T_e: u.K,
    T_i: u.K,
    theta: u.deg,
    gamma_e: Union[float, int] = 1,
    gamma_i: Union[float, int] = 3,
    z_mean: Union[float, int] = None,
):

    r"""
    Using the equation provided by Bellan 2012, this function
    calculates the numerical solution to the two fluid dispersion relation
    presented by Hollweg 1999. This dispersion relation assumes
    :math:`\omega/\omega_{\rm ci} \ll 1`, a uniform magnetic field
    :math: `\mathbf{B_0}`, and quasi-neutrality.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to :math:`T`.
    ion : `str` or `~plasmapy.particles.particle_class.Particle`
        Representation of the ion species (e.g., ``'p'`` for protons, ``'D+'``
        for deuterium, ``'He-4 +1'`` for singly ionized helium-4, etc.). If no
        charge state information is provided, then the ions are assumed to be
        singly ionized.
    k : `~astropy.units.Quantity`, single valued or 1-D array
        Wavenumber in units convertible to :math:`rad / m`.  Either single
        valued or 1-D array of length :math:`N`.
    n_i : `~astropy.units.Quantity`
        Ion number density in units convertible to :math:`m^{-3}`.
    T_e : `~astropy.units.Quantity`
        The electron temperature in units of :math:`K` or :math:`eV`.
    T_i : `~astropy.units.Quantity`
        The ion temperature in units of :math:`K` or :math:`eV`.
    theta : `~astropy.units.Quantity`, single valued or 1-D array
        The angle of propagation of the wave with respect to the magnetic field,
        :math:`\cos^{-1}(k_z / k)`, in units must be convertible to :math:`deg`.
        Either single valued or 1-D array of size :math:`M`.
    gamma_e : `float` or `int`, optional
        The adiabatic index for electrons, which defaults to 1.  This
        value assumes that the electrons are able to equalize their
        temperature rapidly enough that the electrons are effectively
        isothermal.
    gamma_i : `float` or `int`, optional
        The adiabatic index for ions, which defaults to 3. This value
        assumes that ion motion has only one degree of freedom, namely
        along magnetic field lines.
    z_mean : `float` or int, optional
        The average ionization state (arithmetic mean) of the ``ion`` composing
        the plasma.  Will override any charge state defined by argument ``ion``.

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
        When the computed wave frequencies violate the
        :math:`\omega/\omega_{\rm ci} \ll 1` assumption of the dispersion relation.

    Notes
    -----

    The equation presented in Hollweg 1999 [2]_ (equation 3 in Bellan 2012
    [1]_) is:

    .. math::
        \left( \frac{\omega^2}{k_{\rm z}^2 v_{\rm A}^2} - 1 \right) &
        \left[
            \omega^2 \left( \omega^2 - k^2 v_{\rm A}^2 \right)
            - \beta k^2 v_{\rm A}^2 \left(
                \omega^2 - k_{\rm z}^2 v_{\rm A}^2
            \right)
        \right] \\
        &= \omega^2 \left(\omega^2 - k^2 v_{\rm A}^2 \right) k_{\rm x}^2
        \left(
            \frac{c_{\rm s}^2}{\omega_{\rm ci}^2}
            - \frac{c^2}{\omega_{\rm pe}^2} \frac{\omega^2}{k_{\rm z}^2v_{\rm A}^2}
        \right)

    where

    .. math::
        k_{\rm x} = \mathbf{k} \cdot \hat{x}

    References
    ----------
    .. [1] PM Bellan, Improved basis set for low frequency plasma waves, 2012,
       JGR, 117, A12219, doi: `10.1029/2012JA017856
       <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2012JA017856>`_.

    .. [2] JV Hollweg, Kinetic Alfven wave revisited, 1999, JGR, 104(A7),
       14811–14819, doi: `10.1029/1998JA900132
       <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/1998JA900132>`_

    Examples
    --------
    >>> from astropy import units as u
    >>> from plasmapy.dispersion.numerical import hollweg_
    >>> inputs = {
    ...    "k": np.logspace(-7, -2, 2) * u.rad / u.m,
    ...    "theta": 88 * u.deg,
    ...    "n_i": 5 * u.cm ** -3,
    ...    "B": 2.2e-8 * u.T,
    ...    "T_e": 1.6e6 * u.K,
    ...    "T_i": 4.0e5 * u.K,
    ...    "ion": Particle("p+"),
    ... }
    >>> omegas = hollweg(**inputs)
    >>> omegas
    {'fast_mode': <Quantity [2.62911663e-02, 2.27876968e+03] rad / s>,
     'alfven_mode': <Quantity [7.48765909e-04, 2.13800404e+03] rad / s>,
     'acoustic_mode': <Quantity [0.00043295, 0.07358991] rad / s>}
    """

    # validate argument ion
    if not isinstance(ion, Particle):
        try:
            ion = Particle(ion)
        except TypeError:
            raise TypeError(
                f"For argument 'ion' expected type {Particle} but got {type(ion)}."
            )
    if not (ion.is_ion or ion.is_category("element")):
        raise ValueError("The particle passed for 'ion' must be an ion or element.")

    # validate z_mean
    if z_mean is None:
        try:
            z_mean = abs(ion.integer_charge)
        except ChargeError:
            z_mean = 1
    else:
        if not isinstance(z_mean, (int, np.integer, float, np.floating)):
            raise TypeError(
                f"Expected int or float for argument 'z_mean', but got {type(z_mean)}."
            )
        z_mean = abs(z_mean)

    # validate arguments
    for arg_name in ("B", "n_i", "T_e", "T_i"):
        val = locals()[arg_name].squeeze()
        if val.shape != ():
            raise ValueError(
                f"Argument '{arg_name}' must a single value and not an array of "
                f"shape {val.shape}."
            )
        locals()[arg_name] = val

    # validate arguments
    for arg_name in ("gamma_e", "gamma_i"):
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

    # Calc needed plasma parameters
    n_e = z_mean * n_i
    c_s = pfp.ion_sound_speed(
        T_e=T_e,
        T_i=T_i,
        ion=ion,
        n_e=n_e,
        gamma_e=gamma_e,
        gamma_i=gamma_i,
        z_mean=z_mean,
    )
    v_A = pfp.Alfven_speed(B, n_i, ion=ion, z_mean=z_mean)
    omega_ci = pfp.gyrofrequency(B=B, particle=ion, signed=False, Z=z_mean)
    omega_pe = pfp.plasma_frequency(n=n_e, particle="e-")

    # Parameters kx and kz
    kz = np.cos(theta.value) * k
    kx = np.sqrt(k ** 2 - kz ** 2)

    # Bellan2012JGR beta param equation 3
    beta = (c_s / v_A) ** 2

    # Parameters D, F, sigma, and alpha to simplify equation 3
    D = (c_s / omega_ci) ** 2
    F = (c / omega_pe) ** 2
    sigma = (kz * v_A) ** 2
    alpha = (k * v_A) ** 2

    # Polynomial coefficients: c3*x^3 + c2*x^2 + c1*x + c0 = 0
    c3 = F * kx ** 2 + 1
    c2 = -sigma * ((alpha / sigma) * (1 + beta + F * kx ** 2) + D * kx ** 2 + 1)
    c1 = sigma * alpha * (1 + 2 * beta + D * kx ** 2)
    c0 = -beta * alpha * sigma ** 2

    omega = {}
    fast_mode = []
    alfven_mode = []
    acoustic_mode = []

    # If a single k value is given
    if np.isscalar(k.value) is True:
        w = np.emath.sqrt(np.roots([c3.value, c2.value, c1.value, c0.value]))
        fast_mode = np.max(w)
        alfven_mode = np.median(w)
        acoustic_mode = np.min(w)

    # If mutliple k values are given
    else:
        # a0*x^3 + a1*x^2 + a2*x^3 + a3 = 0
        for (a0, a1, a2, a3) in zip(c3, c2, c1, c0):

            w = np.emath.sqrt(np.roots([a0.value, a1.value, a2.value, a3.value]))
            fast_mode.append(np.max(w))
            alfven_mode.append(np.median(w))
            acoustic_mode.append(np.min(w))

    omega["fast_mode"] = fast_mode * u.rad / u.s
    omega["alfven_mode"] = alfven_mode * u.rad / u.s
    omega["acoustic_mode"] = acoustic_mode * u.rad / u.s

    # check the low-frequency limit

    m1 = np.max(omega["fast_mode"])
    m2 = np.max(omega["alfven_mode"])
    m3 = np.max(omega["acoustic_mode"])

    w_max = max(m1, m2, m3)
    w_wci_max = w_max / omega_ci

    # dispersion relation is only valid in the regime w << w_ci
    if w_max / omega_ci > 0.1:
        warnings.warn(
            f"This solver is valid in the regime w/w_ci << 1. "
            f"A w value of {w_max:.2f} and a w/w_ci value of {w_wci_max:.2f} "
            f"were calculated which may affect the validity of the solution.",
            PhysicsWarning,
        )

    return omega
