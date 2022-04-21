"""
This module contains functionality for calculating the numerical
solutions to the Stix cold plasma function.
"""

__all__ = ["stix"]

import astropy.units as u
import numpy as np

from astropy.constants.si import c

from plasmapy.formulary.frequencies import gyrofrequency, plasma_frequency
from plasmapy.particles import Particle, ParticleList
from plasmapy.utils.decorators import validate_quantities

c_si_unitless = c.value


@validate_quantities(
    B={"can_be_negative": False},
    n_i={"can_be_negative": False},
    w={"can_be_negative": False},
)
def stix(
    B: u.T,
    w: u.rad / u.s,
    ions: Particle,
    n_i: u.m ** -3,
    theta: u.rad,
):
    r"""
    Calculate the cold plasma function solution by using
    :cite:t:`bellan:2012`, this uses the numerical method to find the
    wave number(s), (:math:`k`), for the dispersion relation provided by
    :cite:t:`stringer:1963`. This dispersion relation also assumes
    uniform magnetic field :math:`\mathbf{B_0}`, theta is the angle
    between the magnetic and the normal surface of the wave vector.
    For more information see the **Notes** section below.

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field magnitude in units convertible to T.

    w : `~astropy.units.Quantity`, single value omega or 1-D array in
        units convertible to rad/s.

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
        A dictionary of computed wave numbers in units rad/m.  The
        dictionary contains keys for each wave number, this will return
        an array  of value :math:`\theta x 4`.

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
        S = 1 - \sum \frac{\omega^{2}_{p\sigma}}{\omega^{2} -
            \omega^{2}_{c\sigma}}

    .. math::
        P = 1 - \sum_{\sigma} \frac{\omega^{2}_{p\sigma}}{\omega^{2}}

    .. math::
        D = \sum_{\sigma}
            \frac{\omega_{c\sigma}}{\omega}
            \frac{\omega^{2}_{p\sigma}}{\omega^{2} -
            \omega_{c\sigma}^{2}}

    The Cold plasma assumption, Following on section 1.6 of
    :cite:t:`bellan:2012` expresses following derived quantities as
    follows.

    .. math::
        R = S + D \hspace{1cm} L = S - D

    The equation is valid for all :math:`\omega` and :math:`k`
    providing that :math:`\frac{\omega}{k_{z}} >> \nu_{Te}` with
    :math:`\nu_{Ti}` and :math:`k_{x}r_{Le,i} << 1`.  The prediction of
    :math:`k \to 0` occurs when P, R or L cut off and predicts
    :math:`k \to \infty` for perpendicular propagation during wave
    resonance :math:`S \to 0`.

    Example
    -------
    >>> from astropy import units as u
    >>> from plasmapy.particles import Particle
    >>> from plasmapy.dispersion.numerical.stix_ import stix
    >>> inputs = {
    ...     "B": 8.3e-9 * u.T,
    ...     "w": 0.001 * u.rad / u.s,
    ...     "ions": [Particle("H+"), Particle("He+")],
    ...     "n_i": [4.0e5,2.0e5] * u.m**-3,
    ...     "theta": 30 * u.deg,
    >>> }
    >>> w = stix(**inputs)
    >>> print(w)

    """

    # Validate ions argument
    if not isinstance(ions, (list, tuple)):
        ions = [ions]
    ions = ParticleList(ions)

    if not all(failed := [ion.is_ion and ion.charge_number > 0 for ion in ions]):
        raise ValueError(
            "Particle(s) passed to 'ions' must be a positively charged"
            " ion. The following particle(s) is(are) not allowed "
            f"{[ion for ion, fail in zip(ions, failed) if not fail]}"
        )

    # Validate n_i argument
    if n_i.ndim not in (0, 1):
        raise ValueError(
            "Argument 'n_i' must be a single valued or a 1D array of "
            f"size 1 or {len(ions)}, instead got shape of {n_i.shape}"
        )
    elif n_i.ndim == 1 and n_i.size != len(ions):
        raise ValueError(
            "Argument 'n_i' and 'ions' need to be the same length, got"
            f" value of shape {len(ions)} and {len(n_i.shape)}."
        )

    n_i = n_i.value
    if n_i.ndim == 0:
        n_i = np.array([n_i] * len(ions))
    elif n_i.size == 1:
        n_i = np.repeat(n_i, len(ions))

    species = ions + [Particle("e-")]
    densities = np.zeros(n_i.size + 1)
    densities[:-1] = n_i
    densities[-1] = np.sum(n_i * ions.charge_number)

    # Validate B argument
    B = B.squeeze()
    if B.ndim != 0:
        raise ValueError(
            "Argument 'B' must be single valued and not an array of"
            f" shape  {B.shape}."
        )

    # Validate w argument and dimension
    w = w.squeeze()
    if not (w.ndim == 0 or w.ndim == 1):
        raise ValueError(
            f"Argument 'w' needs to be a single value or a 1D array "
            f" astropy Quantity, got a value of shape {w.shape}."
        )
    if np.any(w <= 0):
        raise ValueError(f"Argument 'w' can not have a negative value.")
    if np.isscalar(w.value):
        w = np.array([w.value])

    # Validate theta value
    theta = theta.squeeze()
    theta = theta.to(u.radian)
    if theta.ndim not in (0, 1):
        raise TypeError(
            f"Argument 'theta' needs to be a single value or 1D array "
            f" astropy Quantity, got array of shape {theta.shape}."
        )
    if np.isscalar(theta.value):
        theta = np.array([theta.value])

    # Generate mesh grid of w x theta
    w, theta = np.meshgrid(w, theta, indexing="ij")

    # Generate the plasma parameters needed
    wps = []
    wcs = []

    for par, dens in zip(species, densities.tolist()):
        wps.append(plasma_frequency(n=dens * u.m ** -3, particle=par).value)
        wcs.append(gyrofrequency(B=B, particle=par, signed=False).value)
    wps = np.array(wps)
    wcs = np.array(wcs)

    # Stix method implemented
    S = 1
    P = 1
    D = 0

    for wc, wp in zip(wcs, wps):
        S -= (wp ** 2) / (w ** 2 - wc ** 2)
        P -= (wp / w) ** 2
        D += ((wp ** 2) / (w ** 2 + wc ** 2)) * (wc / w)

    R = S + D
    L = S - D

    # Generate coefficients to solve
    A = (c_si_unitless / w) ** 4 * (S * (np.sin(theta) ** 2) + P * (np.cos(theta) ** 2))
    B = -((c_si_unitless / w) ** 2) * (
        R * L * (np.sin(theta) ** 2) + P * S * (1 + (np.cos(theta) ** 2))
    )
    C = P * R * L

    # Solve for k values
    k = np.empty(4, dtype=np.complex128)

    k[0] = np.emath.sort((-B + np.emath.sqrt(B ** 2 - 4 * A * C)) / (2 * A))
    k[1] = -k[0]
    k[2] = np.emath.sort((-B - np.emath.sqrt(B ** 2 - 4 * A * C)) / (2 * A))
    k[3] = -k[2]

    return k
