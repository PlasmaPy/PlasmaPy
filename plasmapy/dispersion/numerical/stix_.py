"""
This module contains functionality for calculating the numerical
solutions to the Stix cold plasma function.
"""

__all__ = ["stix"]

import astropy.units as u
import numpy as np

from astropy.constants.si import c
from sympy import Symbol
from sympy.solvers import solve

from plasmapy.formulary.frequencies import gyrofrequency
from plasmapy.particles import Particle
from plasmapy.utils.decorators import validate_quantities

c_si_unitless = c.value


@validate_quantities(B={"can_be_negative": False}, k={"can_be_negative": False})
def stix(
    B: u.T,
    k: u.rad / u.m,
    ions: Particle,
    omega_ions: u.rad / u.s,
    theta: u.rad,
):
    r"""
    Calculate the cold plasma function solution by using :cite:t:`bellan:2012`,
    this uses the numerical method to find (:math:`\omega`) dispersion
    relation provided by :cite:t:`stringer:1963`. This dispersion relation also
    assumes a uniform magnetic field :math:`\mathbf{B_0}`, theta is the
    angle between the magnetic and the normal surface of the wave
    vector. For more information see the **Notes** section below.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        Value of the magnitude of the magnetic field in units convertible
        to :math:`T`.
    k : single value or 1 D array astropy `~astropy.units.Quantity`
        Value of the wavenumber in units convertible to radians / m.
        ions: single particle value or 1 D array of particles, ion(s) composing
        the plasma as expressed by chemical symbols.
        omega_ions: single value or 1 D array astropy `~astropy.units.Quantity`
        Frequency value for the associated ion in units convertible to
        radians / s.
    theta: single value or 1 D array astropy `~astropy.units.Quantity`
        Value of theta with respect to the magnetic field,
        :math:`\cos^{-1}(k_z/k)`, must be in units convertible to
        radians.

    Returns
    -------
    omegas : Dict[`str`, `~astropy.units.Quantity`]
        Presents the wavenumber used to find the value(s) of the cold
        plasma frequencies, (omega), dispersion solver(s) and then the
        solutions themselves.

    Raises
    ------
    `TypeError`
        If the argument is of an invalid type.
        `~astropy.units.UnitsError`
        If the argument is a `~astropy.units.Quantity` but is not
        dimensionless.

    `ValueError`
        If the number of frequencies for each ion isn't the same.

    `NoConvergence`
        If a solution cannot be found and the convergence failed to root.

    Notes
    -----
    The cold plasma function is defined by :cite:t:`stringer:1963`, this is equation  8
    of :cite:t:`bellan:2012` presented here:

    .. math::
        (S\sin^{2}(\theta) + P\cos^{2}(\theta))(ck/\omega)^{4} - [RL\sin^{2}() +
        PS(1 + \cos^{2}(theta))](ck/\omega)^{2} + PRL = 0

    where,

    .. math::
        \mathbf{n} = \frac{c \mathbf{k}}{\omega}

    .. math::
        S = 1 - \sum \frac{\omega^{2}_{p\sigma}}{\omega^{2} - \omega^{2}_{c\sigma}}

    .. math::
        P = 1 - \sum \frac{\omega^{2}_{p\sigma}}{\omega^{2} }

    .. math::
        D = \sum \frac{\omega_{c\sigma}}{\omega} \frac{\omega^{2}_{p\sigma}}{\omega^{2} - \omega_{c\sigma}^{2} }

    Following on section 1.6 of :cite:t:`bellan:2012` expresses following derived quantities
    as follows.

    .. math::
        R = S + D \hspace{1cm} L = S - D

    The equation is valid for all :math:`\omega` and :math:`k`
    providing that :math:`\frac{\omega}{k_{z}} >> \nu_{Te}` with
    :math:`\nu_{Ti}` and :math:`k_{x}r_{Le,i} << 1`.  The prediction of
    :math:`k \to 0` occurs when P, R or L cut off and predicts
    :math:`k \to \inf` for perpendicular propagation during wave
    resonance :math:`S \to 0`.
    Example
    -------
    >>> from astropy import units as u
    >>> from plasmapy.dispersion.numerical.stix_ import stix
    >>> inputs = {
    ...     "B": 8.3e-9 * u.T,
    ...     "k": 0.001 * u.rad / u.m,
    ...     "ions": [Particle("H+"), Particle("e-")],
    ...     "omega_ions": [4.0e5,2.0e5] * u.rad / u.s,
    ...     "theta": 30 * u.deg,
    >>> }
    >>> w = stix(**inputs)
    >>> print(w[0.001])

    """

    # validate B argument
    if B.ndim != 0:
        raise ValueError(
            f"Argument 'B' must be a scalar, got array of shape {B.shape}."
        )

    # validate ion frequency(ies) and find dimension
    if omega_ions.ndim == 0 and len(ions) == 0:
        omega_int = True
        lengths = 1
    elif omega_ions.ndim == 1 and len(ions) == len(omega_ions):
        omega_int = False
        lengths = min(len(omega_ions), len(ions))
    else:
        raise ValueError(
            f"Argument 'omega_ions' and 'ions' need to be the same length,"
            f"got value of shape {len(omega_ions)} and {len(ions)}."
        )

    # validate ion argument
    for i in range(lengths):
        if type(ions[i]) is not Particle:
            raise TypeError(
                f"Argument 'ions[i]' need to be particles of particle type"
                f"got value of type {type(ions[i])}."
            )

    # validate k argument and dimension
    k = k.squeeze()
    if not (k.ndim == 0 or k.ndim == 1):
        raise ValueError(
            f"Argument 'k' needs to be a single value or a 1D array astropy Quantity,"
            f"got a value of shape {k.shape}."
        )

    # validate ion frequencies
    omega_ions = omega_ions.squeeze()
    if not (omega_ions.ndim == 0 or omega_ions.ndim == 1):
        raise TypeError(
            f"Argument 'omega_e' needs to be a single value or a single valued "
            f"1D array astropy Quantity, got value of shape {omega_ions.shape}."
        )

    # validate theta value
    theta = theta.squeeze()
    theta = theta.to(u.radian)
    if not (theta.ndim == 0):
        raise TypeError(
            f"Argument 'theta' needs to be a single value astropy Quantity,"
            f"got value of shape {theta.shape}."
        )

    # validate k arguement and find the dimension
    k_dim = k.ndim
    if k_dim == 0:
        ck = np.zeros(1)
        val = k * c_si_unitless
        ck[0] = val.value
    elif k_dim == 1:
        ck = np.zeros(len(k))
        for i in range(len(k)):
            val = k[i] * c_si_unitless
            ck[i] = val.value
    else:
        raise TypeError(
            f"Argument 'k' needs to be a single value or 1D array astropy Quantity,"
            f"got value of shape {k.shape}."
        )

    # Generate frequencies from the given ions
    sum_len = lengths

    plasma_freq = np.zeros(sum_len)

    component_frequency = np.zeros(sum_len)
    for i in range(sum_len):
        component_frequency[i] = gyrofrequency(B=B, particle=ions[i], signed=True).value

    if omega_int is False:
        for i in range(sum_len):
            plasma_freq[i] = float(omega_ions[i].value)
    elif omega_int is True:
        plasma_freq[0] = float(omega_ions.value)
    else:
        raise TypeError(
            f"Argument 'omega_ions', quantity type could not be determined,"
            f"got value of shape {omega_ions.shape}."
        )

    # Stix method implemented
    w = Symbol("w")

    S = 1
    P = 1
    D = 0

    omegas = {}

    for i in range(sum_len):
        S += (plasma_freq[i] ** 2) / (w ** 2 + component_frequency[i] ** 2)
        P += (plasma_freq[i] / w) ** 2
        D += ((plasma_freq[i] ** 2) / (w ** 2 + component_frequency[i] ** 2)) * (
            component_frequency[i] / w
        )

    R = S + D
    L = S - D

    A = S * (np.sin(theta.value) ** 2) + P * (np.cos(theta.value) ** 2)
    B = R * L * (np.sin(theta.value) ** 2) + P * S * (1 + np.cos(theta.value) ** 2)
    C = P * R * L

    # solve the stix equation for single k value or an array
    for i in range(len(ck)):
        eq = A * ((ck[i] / w) ** 4) - B * ((ck[i] / w) ** 2) + C

        sol = solve(eq, w, warn=True)

        sol_omega = []

        for j in range(len(sol)):
            val = complex(sol[j]) * u.rad / u.s
            sol_omega.append(val)

        omegas[i] = sol_omega
        val = ck[i] / c_si_unitless
        omegas[val] = omegas.pop(i)

    return omegas
