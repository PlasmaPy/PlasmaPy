"""
This module contains functionality for calculating various analytical
solutions to the two fluid dispersion relation.
"""
__all__ = ["two_fluid"]

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


@validate_quantities(
    B={"can_be_negative": False},
    n_i={"can_be_negative": False},
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    T_i={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def two_fluid(
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
    Using the solution provided by :cite:t:`bellan:2012`, calculate the
    analytical solution to the two fluid, low-frequency
    (:math:`\omega/kc \ll 1`) dispersion relation presented by
    :cite:t:`stringer:1963`.  This dispersion relation also assumes a
    uniform magnetic field :math:`\mathbf{B_o}`, no D.C. electric field
    :math:`\mathbf{E_o}=0`, and quasi-neutrality.  For more information
    see the **Notes** section below.

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
        Wavenumber in units convertible to rad/m`.  Either single
        valued or 1-D array of length :math:`N`.
    n_i : `~astropy.units.Quantity`
        Ion number density in units convertible to m\ :sup:`-3`.
    T_e : `~astropy.units.Quantity`
        The electron temperature in units of K or eV.
    T_i : `~astropy.units.Quantity`
        The ion temperature in units of K or eV.
    theta : `~astropy.units.Quantity`, single valued or 1-D array
        The angle of propagation of the wave with respect to the
        magnetic field, :math:`\cos^{-1}(k_z / k)`, in units must be
        convertible to radians. Either single valued or 1-D array of
        size :math:`M`.
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
        each key will be a :math:`N × M` array.

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

    The complete dispersion equation presented by :cite:t:`stringer:1963`
    (equation 1 of :cite:t:`bellan:2012`) is:

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

    :math:`\omega` is the wave frequency, :math:`k` is the wavenumber,
    :math:`v_A` is the Alfvén velocity, :math:`c_s` is the sound speed,
    :math:`\omega_{ci}` is the ion gyrofrequency, and
    :math:`\omega_{pe}` is the electron plasma frequency. This relation
    does additionally assume low-frequency waves
    :math:`\omega/kc \ll 1`, no D.C. electric field
    :math:`\mathbf{E_o}=0` and quasi-neutrality.

    Following section 5 of :cite:t:`bellan:2012` the exact roots of the
    above dispersion equation can be derived and expressed as one
    analytical solution (equation 38 of :cite:t:`bellan:2012`):

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

    where :math:`j = 0` represents the fast mode, :math:`j = 1`
    represents the Alfvén mode, and :math:`j = 2` represents the
    acoustic mode.  Additionally,

    .. math::
        p &= \frac{3B-A^2}{3} \; , \; q = \frac{9AB-2A^3-27C}{27} \\
        A &= \frac{Q + Q^2 \beta + Q \alpha + \alpha \Lambda}{Q^2} \;
            , \; B = \alpha \frac{1 + 2 Q \beta + \Lambda \beta}{Q^2} \;
            , \; C = \frac{\alpha^2 \beta}{Q^2} \\
        \alpha &= \cos^2 \theta \;
            , \; \beta = \left( \frac{c_s}{v_A}\right)^2 \;
            , \; \Lambda = \left( \frac{k v_{A}}{\omega_{ci}}\right)^2

    Examples
    --------
    >>> from astropy import units as u
    >>> from plasmapy.dispersion.analytical import two_fluid
    >>> inputs = {
    ...     "k": 0.01 * u.rad / u.m,
    ...     "theta": 30 * u.deg,
    ...     "B": 8.3e-9 * u.T,
    ...     "n_i": 5e6 * u.m ** -3,
    ...     "T_e": 1.6e6 * u.K,
    ...     "T_i": 4.0e5 * u.K,
    ...     "ion": "p+",
    ... }
    >>> omegas = two_fluid(**inputs)
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
    >>> omegas = two_fluid(**inputs)
    >>> omegas['fast_mode']
    <Quantity [[0.00767..., 0.00779... ],
               [0.01534..., 0.01558...]] rad / s>
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
    if k.ndim not in [0, 1]:
        raise ValueError(
            f"Argument 'k' needs to be a single valued or 1D array astropy Quantity,"
            f" got array of shape {k.shape}."
        )
    if np.any(k <= 0):
        raise ValueError("Argument 'k' can not be a or have negative values.")

    # validate argument theta
    theta = theta.squeeze()
    if theta.ndim not in [0, 1]:
        raise ValueError(
            f"Argument 'theta' needs to be a single valued or 1D array astropy "
            f"Quantity, got array of shape {k.shape}."
        )

    # Calc needed plasma parameters
    n_e = z_mean * n_i
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PhysicsWarning)
        c_s = ion_sound_speed(
            T_e=T_e,
            T_i=T_i,
            ion=ion,
            n_e=n_e,
            gamma_e=gamma_e,
            gamma_i=gamma_i,
            z_mean=z_mean,
        )
    v_A = Alfven_speed(B, n_i, ion=ion, z_mean=z_mean)
    omega_ci = gyrofrequency(B=B, particle=ion, signed=False, Z=z_mean)
    omega_pe = plasma_frequency(n=n_e, particle="e-")

    # Bellan2012JGR params equation 32
    alpha = np.cos(theta.value) ** 2
    beta = (c_s / v_A).to(u.dimensionless_unscaled).value ** 2
    alphav, kv = np.meshgrid(alpha, k.value)  # create grid
    Lambda = (kv * v_A.value / omega_ci.value) ** 2

    # Bellan2012JGR params equation 2
    Q = 1 + (kv * c.value / omega_pe.value) ** 2

    # Bellan2012JGR params equation 35
    A = ((1 + alphav) / Q) + beta + (alphav * Lambda / Q**2)
    B = alphav * (1 + 2 * Q * beta + Lambda * beta) / Q**2
    C = beta * (alphav / Q) ** 2

    # Bellan2012JGR params equation 36
    p = (3 * B - A**2) / 3
    q = (9 * A * B - 2 * A**3 - 27 * C) / 27

    # Bellan2012JGR params equation 38
    R = 2 * Lambda * np.emath.sqrt(-p / 3)
    S = 3 * q / (2 * p) * np.emath.sqrt(-3 / p)
    T = Lambda * A / 3
    omega = {}
    for ind, key in enumerate(("fast_mode", "alfven_mode", "acoustic_mode")):
        # The solution corresponding to equation 38
        w = omega_ci * np.emath.sqrt(
            R * np.cos(1 / 3 * np.emath.arccos(S) - 2 * np.pi / 3 * ind) + T
        )
        omega[key] = w.squeeze()

        # check for violation of dispersion relation assumptions
        # (i.e. low-frequency, w/kc << 0.1)
        wkc_max = np.max(w.value / (kv * c.value))
        if wkc_max > 0.1:
            warnings.warn(
                f"The {key} calculation produced a high-frequency wave (w/kc == "
                f"{wkc_max:.3f}), which violates the low-frequency (w/kc << 1) "
                f"assumption of the dispersion relation.",
                PhysicsWarning,
            )

    return omega
