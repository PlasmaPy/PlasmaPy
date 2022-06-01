"""
This module contains functionality for calculating various numerical
solutions to the kinetic alfven dispersion relation
"""
__all__ = ["kinetic_alfven"]

import astropy.units as u
import numpy as np
import warnings

from astropy.constants.si import c
from typing import Union

from plasmapy.formulary import frequencies as pfp
from plasmapy.formulary import speeds as speed
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
def kinetic_alfven(
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
    Using the equation provided in :cite:t:`bellan:2012`, this function
    calculates the numerical solution to the two fluid dispersion
    relation presented by :cite:t:`hirose:2004`.
    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to :math:`T`.
    ion : `str` or `~plasmapy.particles.particle_class.Particle`
        Representation of the ion species (e.g., ``'p'`` for protons,
        ``'D+'`` for Deuterium, ``'He-4 +1'`` for singly ionized
        Helium-4, etc.). If no charge state information is provided,
        then the ions are assumed to be singly ionized.
    k : `~astropy.units.Quantity`, single valued or 1-D array
        Wavenumber in units convertible to :math:`rad / m`.  Either
        single valued or 1-D array of length :math:`N`.
    n_i : `~astropy.units.Quantity`
        Ion number density in units convertible to :math:`m^{-3}`.
    T_e : `~astropy.units.Quantity`
        The electron temperature in units of :math:`K` or :math:`eV`.
    T_i : `~astropy.units.Quantity`
        The ion temperature in units of :math:`K` or :math:`eV`.
    theta : `~astropy.units.Quantity`, single valued or 1-D array
        The angle of propagation of the wave with respect to the
        magnetic field, :math:`\cos^{-1}(k_z / k)`, in units must be
        convertible to :math:`rad`. Either single valued or 1-D array
        of size :math:`M`.
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
        A dictionary of computed wave frequencies in units
        :math:`rad/s`.  The dictionary contains three keys:
        ``'fast_mode'`` for the fast mode, ``'alfven_mode'`` for the
        Alfv\'{e}n mode, and ``'acoustic_mode'`` for the ion-acoustic
        mode.  The value for each key will be a :math:`N x M` array.
    Raises
    ------
    TypeError
        If applicable arguments are not instances of
        `~astropy.units.Quantity` or cannot be converted into one.
    TypeError
        If ``ion`` is not of type or convertible to
        `~plasmapy.particles.Particle`.
    TypeError
        If ``gamma_e``, ``gamma_i``, or``z_mean`` are not of type
        `int` or `float`.
    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the
        expected units.
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
    Notes
    -----
    Solves equation 5 in Bellan2012JGR (2x2 matrix method
    argued in Hasegawa and Uberoi 1982, Morales and Maggs 1997,
    and Lysak and Lotko 1996)
    ..math::
        \omega^2 = k_{\rm z}^2 v_{\rm A}^2 \left(1 + \frac{k_{\rm x}^2 &
        c_{\rm s}^2}{\omega_{\rm ci}^2} \right)
    Examples
    --------
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from plasmapy.particles import Particle
    >>> from plasmapy.dispersion.numerical import kinetic_alfven_
    >>> inputs = {
    ...    "k": np.logspace(-7,-2,2) * u.rad / u.m,
    ...    "theta": 30 * u.deg,
    ...    "B": 8.3e-9 * u.T,
    ...    "n_i": 5 * u.m ** -3,
    ...    "T_e": 1.6e6 * u.K,
    ...    "T_i": 4.0e5 * u.K,
    ...    "ion": Particle("p+"),
    ...}
    >>> omegas = kinetic_alfven(**inputs)
    [7.01005647e+00 6.70197761e+08] rad / s
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
            z_mean = abs(ion.charge_number)
        except ChargeError:
            z_mean = 1
    else:
        if not isinstance(z_mean, (int, np.integer, float, np.floating)):
            raise TypeError(
                "Expected int or float for argument 'z_mean', "
                f"instead got {type(z_mean)}."
            )
        z_mean = abs(z_mean)

    # validate arguments
    for arg_name in ("B", "n_i", "T_e", "T_i"):
        val = locals()[arg_name].squeeze()
        if val.shape != ():
            raise ValueError(
                f"Argument '{arg_name}' must a single value and "
                f"not an array of shape {val.shape}."
            )
        locals()[arg_name] = val

    # validate arguments
    for arg_name in ("gamma_e", "gamma_i"):
        if not isinstance(locals()[arg_name], (int, np.integer, float, np.floating)):
            raise TypeError(
                f"Expected int or float for argument '{arg_name}', "
                f"instead got {type(locals()[arg_name])}."
            )

    # validate argument k
    k = k.squeeze()
    if not (k.ndim == 0 or k.ndim == 1):
        raise ValueError(
            "Argument 'k' needs to be a single valued or 1D array "
            f"astropy Quantity, instead got array of shape {k.shape}."
        )
    if np.any(k <= 0):
        raise ValueError("Argument 'k' can not be a or have negative values.")

    # validate argument theta
    theta = theta.squeeze()
    theta = theta.to(u.radian)
    if not (theta.ndim == 0 or theta.ndim == 1):
        raise ValueError(
            "Argument 'theta' needs to be a single valued or 1D array "
            f"astropy Quantity, instead got array of shape {theta.shape}."
        )

    n_e = z_mean * n_i
    c_s = speed.ion_sound_speed(
        T_e=T_e,
        T_i=T_i,
        ion=ion,
        n_e=n_e,
        gamma_e=gamma_e,
        gamma_i=gamma_i,
        z_mean=z_mean,
    )
    v_A = speed.Alfven_speed(B, n_i, ion=ion, z_mean=z_mean)
    omega_ci = pfp.gyrofrequency(B=B, particle=ion, signed=False, Z=z_mean)

    # parameters kz
    kz = np.cos(theta.value) * k
    kx = np.sqrt(k ** 2 - kz ** 2)

    # parameters sigma, D, and F to simplify equation 3
    A = (kz * v_A) ** 2
    F = ((kx * c_s) / omega_ci) ** 2

    omega = np.sqrt(A * (1 + F))

    # thermal speeds for electrons and ions in plasma
    v_Te = speed.thermal_speed(T=T_e, particle="e-")
    v_Ti = speed.thermal_speed(T=T_i, particle=ion)

    # maximum value of omega
    w_max = np.max(omega)

    # maximum and minimum values for w/kz
    omega_kz = omega / kz

    omega_kz_max = np.max(omega_kz)
    omega_kz_min = np.min(omega_kz)

    # dispersion relation is only valid in v_Te >> w/kz >> v_Ti

    # maximum value for w/kz test
    if omega_kz_max / v_Te > 0.1 or v_Ti / omega_kz_max > 0.1:
        warnings.warn(
            "This calculation produced one or more invalid w/kz "
            "value(s), which violates the regime in which the "
            "dispersion relation is valid (v_Te >> w/kz >> v_Ti)",
            PhysicsWarning,
        )

    # minimum value for w/kz test
    elif omega_kz_min / v_Te > 0.1 or v_Ti / omega_kz_min > 0.1:
        warnings.warn(
            "This calculation produced one or more invalid w/kz "
            "value(s) which violates the regime in which the "
            "dispersion relation is valid (v_Te >> w/kz >> v_Ti)",
            PhysicsWarning,
        )

    # dispersion relation is only valid in the regime w << w_ci
    if w_max / omega_ci > 0.1:
        warnings.warn(
            "The calculation produced a high-frequency wave, "
            "which violates the low frequency assumption (w << w_ci)",
            PhysicsWarning,
        )

    return omega


inputs = {
    "k": np.logspace(-7, -2, 2) * u.rad / u.m,
    "theta": 30 * u.deg,
    "B": 8.3e-9 * u.T,
    "n_i": 5 * u.m ** -3,
    "T_e": 1.6e6 * u.K,
    "T_i": 4.0e5 * u.K,
    "ion": Particle("p+"),
}

omegas = kinetic_alfven(**inputs)
print(omegas)
