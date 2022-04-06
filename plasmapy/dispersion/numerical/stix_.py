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

from plasmapy.formulary.frequencies import gyrofrequency, plasma_frequency
from plasmapy.particles import Particle, ParticleList
from plasmapy.utils.decorators import validate_quantities

c_si_unitless = c.value


@validate_quantities(
    B={"can_be_negative": False},
    n_i={"can_be_negative": False},
    k={"can_be_negative": False, "equivalencies": u.spectral()},
)
def stix(
    B: u.T,
    k: u.rad / u.m,
    ions: Particle,
    n_i: u.m ** -3,
    theta: u.rad,
):
    r"""
    Calculate the cold plasma function solution by using
    :cite:t:`bellan:2012`, this uses the numerical method to find
    (:math:`\omega`) dispersion relation provided by
    :cite:t:`stringer:1963`. This dispersion relation also assumes
    uniform magnetic field :math:`\mathbf{B_0}`, theta is the angle
    between the magnetic and the normal surface of the wave vector.
    For more information see the **Notes** section below.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        Value of the magnitude of the magnetic field in units convertible
        to :math:`T`.
    k : single value or 1 D array astropy `~astropy.units.Quantity`
        Value of the wavenumber in units convertible to radians / m.

    ions: a single or `list` of :term:`particle-like` object(s)
        epresentation of the ion species (e.g., ``"p"`` for protons,
        ``"D+"`` for deuterium, ``["H+", "He+"]`` for hydrogen and
        helium, etc.).  The charge state for each species must be
        specified.

    n_i: `~astropy.units.Quantity`, single valud or 1-D array
        Ion number density in units convertible to m\ :sup:`-3`.  Must
        be single valued or equal length to ``ions``.

    theta: `~astropy.units.Quantity`, single valued or 1-D array
        The angle of propagation of the wave with respect to the
        magnetic field, :math:`\cos^{-1}(k_z / k)`, in units convertible
        to radians.  Either single valued or 1-D array of size
        :math:`M`.

    Returns
    -------
    omegas : Dict[`str`, `~astropy.units.Quantity`]
        A dictionary of computed wave frequencies in units rad/s.  The
        dictionary contains three keys: ``"fast_mode"`` for the fast
        mode, ``"alfven_mode"`` for the Alfvén mode, and
        ``"acoustic_mode"`` for the ion-acoustic mode.  The value for
        each key will be a :math:`N x M` array.

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
        If a solution cannot be found and the convergence failed to
        root.

    Notes
    -----
    The cold plasma function is defined by :cite:t:`stringer:1963`,
    this is equation 8 of :cite:t:`bellan:2012` and is presented below.
    It is assumed that the zero-order quantities are uniform in space
    and static in time; while the first-order quantities are assumed to
    vary as :math:`e^{\left [ i (\textbf{k}\cdot\textbf{r} - \omega t)
    \right ]}` :cite:t:`stix:1992`.

    .. math::
        (S\sin^{2}(\theta) + P\cos^{2}(\theta))(ck/\omega)^{4}
            - [
                RL\sin^{2}(\theta) + PS(1 + \cos^{2}(\theta))
            ](ck/\omega)^{2} + PRL = 0

    where,

    .. math::
        \mathbf{n} = \frac{c \mathbf{k}}{\omega}

    .. math::
        S = 1 - \sum \frac{\omega^{2}_{p\sigma}}{\omega^{2} - \omega^{2}_{c\sigma}}

    .. math::
        P = 1 - \sum_{\sigma} \frac{\omega^{2}_{p\sigma}}{\omega^{2}}

    .. math::
        D = \sum_{\sigma}
            \frac{\omega_{c\sigma}}{\omega}
            \frac{\omega^{2}_{p\sigma}}{\omega^{2} - \omega_{c\sigma}^{2}}

    The Cold plasma assumption, Following on section 1.6 of
    :cite:t:`bellan:2012` expresses following derived quantities as
    follows.

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
    >>> from plasmapy.particles import Particle
    >>> from plasmapy.dispersion.numerical.stix_ import stix
    >>> inputs = {
    ...     "B": 8.3e-9 * u.T,
    ...     "k": 0.001 * u.rad / u.m,
    ...     "ions": [Particle("H+"), Particle("e-")],
    ...     "n_i": [4.0e5,2.0e5] * u.m**-3,
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

    # validate ions argument
    if not isinstance(ions, (list, tuple)):
        ions = [ions]
    ions = ParticleList(ions)

    if not all(failed := [ion.is_ion and ion.charge_number > 0 for ion in ions]):
        raise ValueError(f"The particle passed for 'ions' must only be positive. {failed}")

    # validate n_i argument
    if not (n_i.ndim == 0 or n_i.ndim == 1):
        raise ValueError(
            f"Argument 'n_i' must be a float or an array of floats,"
            f" instead got shape of {n_i.shape}"
        )

    # find dimension and validate
    if n_i.ndim == 0 and len(ions) == 1:
        omega_int = True
        lengths = 1
    elif n_i.ndim == 1 and len(ions) >= 1:
        omega_int = False
        lengths = min(len(n_i), len(ions))
    else:
        raise ValueError(
            f"Argument 'n_i' and 'ions' need to be the same length, "
            f"got value of shape {len(n_i)} and {len(ions)}."
        )

    species = ions + [Particle("e-")]
    arg_ = 0
    for i in range(lengths):
        arg_ = arg_ + n_i[i].val
    arg_ = arg_ * u.m ** -3
    n_i.append(arg_)

    # validate k argument and dimension
    k = k.squeeze()
    if not (k.ndim == 0 or k.ndim == 1):
        raise ValueError(
            f"Argument 'k' needs to be a single value or a 1D array astropy Quantity,"
            f"got a value of shape {k.shape}."
        )

    # validate ion densities
    n_i = n_i.squeeze()
    if not (n_i.ndim == 0 or n_i.ndim == 1):
        raise TypeError(
            f"Argument 'omega_e' needs to be a single value or a single valued "
            f"1D array astropy Quantity, got value of shape {n_i.shape}."
        )

    for density in n_i:
        arg_ = plasma_frequency(density)
        density = arg_

    # validate theta value
    theta = theta.squeeze()
    theta = theta.to(u.radian)
    if not (theta.ndim == 0):
        raise TypeError(
            f"Argument 'theta' needs to be a single value astropy Quantity,"
            f"got value of shape {theta.shape}."
        )

    # validate k argument and find the dimension
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
            plasma_freq[i] = float(n_i[i].value)
    elif omega_int is True:
        plasma_freq[0] = float(n_i.value)
    else:
        raise TypeError(
            f"Argument 'n_i', quantity type could not be determined,"
            f"got value of shape {n_i.shape}."
        )

    # Stix method implemented
    w = Symbol("w")

    S = 1
    P = 1
    D = 0

    omegas = {}

    for i in range(sum_len):
        S += 0#(plasma_freq[i] ** 2) / (w ** 2 + component_frequency[i] ** 2)
        P += 0#(plasma_freq[i] / w) ** 2
        D += 0#((plasma_freq[i] ** 2) / (w ** 2 + component_frequency[i] ** 2)) * (
            #component_frequency[i] / w
        #)

    R = S + D
    L = S - D

    A = S * (np.sin(theta.value) ** 2) + P * (np.cos(theta.value) ** 2)
    B = R * L * (np.sin(theta.value) ** 2) + P * S * (1 + np.cos(theta.value) ** 2)
    C = P * R * L

    print("Note: Solution computation time may vary.")

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
