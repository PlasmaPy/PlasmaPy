"""
This module contains functionality for calculating various numerical
solutions to Hollweg's two fluid dispersion relation
"""
__all__ = ["hollweg"]

import astropy.units as u
import numpy as np
import warnings

from astropy.constants.si import c
from typing import Union

from plasmapy.formulary.frequencies import gyrofrequency, plasma_frequency
from plasmapy.formulary.speeds import Alfven_speed, ion_sound_speed
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import ChargeError
from plasmapy.utils.decorators import validate_quantities
from plasmapy.utils.exceptions import PhysicsWarning

c_si_unitless = c.value


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
    n_i: u.m**-3,
    T_e: u.K,
    T_i: u.K,
    theta: u.rad,
    gamma_e: Union[float, int] = 1,
    gamma_i: Union[float, int] = 3,
    z_mean: Union[float, int] = None,
):
    r"""
    Calculate the two fluid dispersion relation presented by
    :cite:t:`hollweg:1999`, and discussed by :cite:t:`bellan:2012`.
    This is a numberical solver of equation 3 in :cite:t:`bellan:2012`.
    See the **Notes** section below for additional details.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to T.
    ion : `str` or `~plasmapy.particles.particle_class.Particle`
        Representation of the ion species (e.g., ``'p'`` for protons,
        ``'D+'`` for deuterium, ``'He-4 +1'`` for singly ionized
        helium-4, etc.). If no charge state information is provided,
        then the ions are assumed to be singly ionized.
    k : `~astropy.units.Quantity`, single valued or 1-D array
        Wavenumber in units convertible to rad/m.  Either single
        valued or 1-D array of length :math:`N`.
    n_i : `~astropy.units.Quantity`
        Ion number density in units convertible to m\ :sup:`-3`.
    T_e : `~astropy.units.Quantity`
        The electron temperature in units of K or eV.
    T_i : `~astropy.units.Quantity`
        The ion temperature in units of K or eV.
    theta : `~astropy.units.Quantity`, single valued or 1-D array
        The angle of propagation of the wave with respect to the
        magnetic field, :math:`\cos^{-1}(k_z / k)`, in units convertible
        to radians.  Either single valued or 1-D array of size
        :math:`M`.
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
        The average ionization state (arithmetic mean) of the ``ion``
        composing the plasma.  Will override any charge state defined
        by argument ``ion``.

    Returns
    -------
    omega : Dict[str, `~astropy.units.Quantity`]
        A dictionary of computed wave frequencies in units rad/s.  The
        dictionary contains three keys: ``'fast_mode'`` for the fast
        mode, ``'alfven_mode'`` for the Alfvén mode, and
        ``'acoustic_mode'`` for the ion-acoustic mode.  The value for
        each key will be a :math:`N x M` array.

    Raises
    ------
    TypeError
        If applicable arguments are not instances of
        `~astropy.units.Quantity` or cannot be converted into one.

    TypeError
        If ``ion`` is not of type or convertible to
        `~plasmapy.particles.particle_class.Particle`.

    TypeError
        If ``gamma_e``, ``gamma_i``, or ``z_mean`` are not of type `int`
        or `float`.

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the
        expected units.

    ValueError
        If any of ``B``, ``k``, ``n_i``, ``T_e``, or ``T_i`` is
        negative.

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
        When :math:`\omega / \omega_{\rm ci} > 0.1`, violation of the
        low-frequency assumption.

    : `~plasmapy.utils.exceptions.PhysicsWarning`
        When :math:`c_{\rm s} / v_{\rm A} > 0.1`, violation of low-β.

    : `~plasmapy.utils.exceptions.PhysicsWarning`
        When :math:`|θ - π/2| > 0.1`, violation of quasi-perpendicular
        propagation.

    Notes
    -----

    The dispersion relation presented in :cite:t:`hollweg:1999`
    (equation 3 in :cite:t:`bellan:2012`) is:

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
        \mathbf{B_o} &= B_{o} \mathbf{\hat{z}} \\
        \cos \theta &= \frac{k_z}{k} \\
        \mathbf{k} &= k_{\rm x} \hat{x} + k_{\rm z} \hat{z}

    :math:`\omega` is the wave frequency, :math:`k` is the wavenumber,
    :math:`v_{\rm A}` is the Alfvén velocity, :math:`c_{\rm s}` is the
    sound speed, :math:`\omega_{\rm ci}` is the ion gyrofrequency, and
    :math:`\omega_{\rm pe}` is the electron plasma frequency. In the
    derivation of this relation Hollweg assumed low-frequency waves
    :math:`\omega / \omega_{\rm ci} \ll 1`, no D.C. electric field
    :math:`\mathbf{E_o}=0`, and quasi-neutrality.

    :cite:t:`hollweg:1999` asserts this expression is valid for
    arbitrary :math:`c_{\rm s} / v_{\rm A}` (β) and
    :math:`k_{\rm z} / k` (θ).  Contrarily, :cite:t:`bellan:2012`
    states in §1.7 that due to the inconsistent retention of the
    :math:`\omega / \omega_{\rm ci} \ll 1` terms the expression can
    only be valid if both :math:`c_{\rm s} \ll v_{\rm A}` (low-β) and
    the wave propgation is nearly perpendicular to the magnetic field.

    This routine solves for ω for given :math:`k` values by numerically
    solving for the roots of the above expression.

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
    {'fast_mode': <Quantity [2.62911663e-02+0.j, 2.27876968e+03+0.j] rad / s>,
     'alfven_mode': <Quantity [7.48765909e-04+0.j, 2.13800404e+03+0.j] rad / s>,
     'acoustic_mode': <Quantity [0.00043295+0.j, 0.07358991+0.j] rad / s>}
    """

    # validate argument ion
    if not isinstance(ion, Particle):
        try:
            ion = Particle(ion)
        except TypeError as ex:
            raise TypeError(
                f"For argument 'ion' expected type {Particle} but got {type(ion)}."
            ) from ex
    if not ion.is_ion and not ion.is_category("element"):
        raise ValueError("The particle passed for 'ion' must be an ion or element.")

    # validate z_mean
    if z_mean is None:
        try:
            z_mean = abs(ion.charge_number)
        except ChargeError:
            z_mean = 1
    elif isinstance(z_mean, (int, np.integer, float, np.floating)):
        z_mean = abs(z_mean)
    else:
        raise TypeError(
            f"Expected int or float for argument 'z_mean', but got {type(z_mean)}."
        )

    # validate arguments
    for arg_name in ("B", "n_i", "T_e", "T_i"):
        val = locals()[arg_name].squeeze()
        if val.shape != ():
            raise ValueError(
                f"Argument '{arg_name}' must be single valued and not an array of "
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
    if k.ndim not in [0, 1]:
        raise ValueError(
            f"Argument 'k' needs to be single valued or a 1D array astropy Quantity,"
            f" got array of shape {k.shape}."
        )
    if np.any(k <= 0):
        raise ValueError("Argument 'k' can not be a or have negative values.")

    # validate argument theta
    theta = theta.squeeze()
    if theta.ndim not in (0, 1):
        raise ValueError(
            f"Argument 'theta' needs to be a single valued or 1D array astropy "
            f"Quantity, got array of shape {theta.shape}."
        )

    # Single k value case
    if np.isscalar(k.value):
        k = np.array([k.value]) * u.rad / u.m

    # Calc needed plasma parameters
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PhysicsWarning)
        n_e = z_mean * n_i
        c_s = ion_sound_speed(
            T_e=T_e,
            T_i=T_i,
            ion=ion,
            n_e=n_e,
            gamma_e=gamma_e,
            gamma_i=gamma_i,
            z_mean=z_mean,
        ).value
        v_A = Alfven_speed(B, n_i, ion=ion, z_mean=z_mean).value
        omega_ci = gyrofrequency(B=B, particle=ion, signed=False, Z=z_mean).value
        omega_pe = plasma_frequency(n=n_e, particle="e-").value

    cs_vA = c_s / v_A
    thetav, kv = np.meshgrid(theta.value, k.value)

    # Parameters kx and kz
    kz = np.cos(thetav) * kv
    kx = np.sin(thetav) * kv

    # Define helpful parameters
    beta = (c_s / v_A) ** 2
    alpha_A = (kv * v_A) ** 2
    alpha_s = (kv * c_s) ** 2  # == alpha_A * beta
    sigma = (kz * v_A) ** 2
    D = (c_s / omega_ci) ** 2
    F = (c_si_unitless / omega_pe) ** 2

    # Polynomial coefficients: c3*x^3 + c2*x^2 + c1*x + c0 = 0
    c3 = F * kx**2 + 1
    c2 = -alpha_A * (1 + beta + F * kx**2) - sigma * (1 + D * kx**2)
    c1 = sigma * alpha_A * (1 + 2 * beta + D * kx**2)
    c0 = -alpha_s * sigma**2

    # Find roots to polynomial
    coefficients = np.array([c3, c2, c1, c0], ndmin=3)
    nroots = coefficients.shape[0] - 1  # 3
    nks = coefficients.shape[1]
    nthetas = coefficients.shape[2]
    roots = np.empty((nroots, nks, nthetas), dtype=np.complex128)
    for ii in range(nks):
        for jj in range(nthetas):
            roots[:, ii, jj] = np.roots(coefficients[:, ii, jj])

    roots = np.sqrt(roots)
    roots = np.sort(roots, axis=0)

    # Warn about NOT low-beta
    if c_s / v_A > 0.1:
        warnings.warn(
            f"This solver is valid in the low-beta regime, "
            f"c_s/v_A << 1 according to Bellan, 2012, Sec. 1.7 "
            f"(see documentation for DOI). A c_s/v_A value of {cs_vA:.2f} "
            f"was calculated which may affect the validity of the solution.",
            PhysicsWarning,
        )

    # Warn about theta not nearly perpendicular
    theta_diff_max = np.amax(np.abs(thetav - np.pi / 2))
    if theta_diff_max > 0.1:
        warnings.warn(
            f"This solver is valid in the regime where propagation is "
            f"nearly perpendicular to B according to Bellan, 2012, Sec. 1.7 "
            f"(see documentation for DOI). A |theta - pi/2| value of "
            f"{theta_diff_max:.2f} was calculated which may affect the "
            f"validity of the solution.",
            PhysicsWarning,
        )

    # dispersion relation is only valid in the regime w << w_ci
    w_max = np.max(roots)
    w_wci_max = w_max / omega_ci
    if w_wci_max > 0.1:
        warnings.warn(
            f"This solver is valid in the regime w/w_ci << 1.  A w "
            f"value of {w_max:.2f} and a w/w_ci value of "
            f"{w_wci_max:.2f} were calculated which may affect the "
            f"validity of the solution.",
            PhysicsWarning,
        )

    return {
        "fast_mode": roots[2, :].squeeze() * u.rad / u.s,
        "alfven_mode": roots[1, :].squeeze() * u.rad / u.s,
        "acoustic_mode": roots[0, :].squeeze() * u.rad / u.s,
    }
