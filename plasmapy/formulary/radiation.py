"""
Functions for calculating quantities associated with electromagnetic
radiation.
"""

__all__ = [
    "thermal_bremsstrahlung",
]

import astropy.constants as const
import astropy.units as u
import numpy as np

from scipy.special import exp1

from plasmapy.formulary.frequencies import plasma_frequency
from plasmapy.particles import particle_input, ParticleLike
from plasmapy.utils.decorators import validate_quantities
from plasmapy.utils.exceptions import PhysicsError


@validate_quantities(
    frequencies={"can_be_negative": False},
    n_e={"can_be_negative": False},
    n_i={"can_be_negative": False},
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
@particle_input
def thermal_bremsstrahlung(
    frequencies: u.Hz,
    n_e: u.m**-3,
    T_e: u.K,
    n_i: u.m**-3 = None,
    ion_species: ParticleLike = "H+",
    kmax: u.m = None,
) -> np.ndarray:
    r"""
    Calculate the bremsstrahlung emission spectrum for a Maxwellian plasma
    in the Rayleigh-Jeans limit :math:`ℏ ω ≪ k_B T_e`

    .. math::
       \frac{dP}{dω} = \frac{8 \sqrt{2}}{3\sqrt{π}}
       \bigg ( \frac{e^2}{4 π ε_0} \bigg )^3
       \bigg ( m_e c^2 \bigg )^{-\frac{3}{2}}
       \bigg ( 1 - \frac{ω_{pe}^2}{ω^2} \bigg )^\frac{1}{2}
       \frac{Z_i^2 n_i n_e}{\sqrt(k_B T_e)}
       E_1(y)

    where :math:`E_1` is the exponential integral

    .. math::
        E_1 (y) = - \int_{-y}^∞ \frac{e^{-t}}{t}dt

    and :math:`y` is the dimensionless argument

    .. math::
        y = \frac{1}{2} \frac{ω^2 m_e}{k_{max}^2 k_B T_e}

    where :math:`k_{max}` is a maximum wavenumber approximated here as
    :math:`k_{max} = 1/λ_B` where  :math:`λ_B` is the electron
    de Broglie wavelength.

    Parameters
    ----------
    frequencies : `~astropy.units.Quantity`
        Array of frequencies over which the bremsstrahlung spectrum will be
        calculated (convertible to Hz).

    n_e : `~astropy.units.Quantity`
        Electron number density in the plasma (convertible to m\ :sup:`-3`\ ).

    T_e : `~astropy.units.Quantity`
        Temperature of the electrons (in K or convertible to eV).

    n_i : `~astropy.units.Quantity`, optional
        Ion number density in the plasma (convertible to m\ :sup:`-3`\ ). Defaults
        to the quasi-neutral condition :math:`n_i = n_e / Z`\ .

    ion_species : `str` or `~plasmapy.particles.particle_class.Particle`, optional
        An instance of `~plasmapy.particles.particle_class.Particle`, or a string
        convertible to `~plasmapy.particles.particle_class.Particle`.

    kmax :  `~astropy.units.Quantity`
        Cutoff wavenumber (convertible to radians per meter). Defaults
        to the inverse of the electron de Broglie wavelength.

    Returns
    -------
    spectrum : `~astropy.units.Quantity`
        Computed bremsstrahlung spectrum over the frequencies provided.

    Notes
    -----
    For details, see :cite:t:`bekefi:1966`\ .
    """

    # Default n_i is n_e/Z:
    if n_i is None:
        n_i = n_e / ion_species.charge_number

    # Default value of kmax is the electrom thermal de Broglie wavelength
    if kmax is None:
        kmax = (np.sqrt(const.m_e.si * const.k_B.si * T_e) / const.hbar.si).to(1 / u.m)

    # Convert frequencies to angular frequencies
    ω = (frequencies * 2 * np.pi * u.rad).to(u.rad / u.s)

    # Calculate the electron plasma frequency
    ω_pe = plasma_frequency(n=n_e, particle="e-")

    # Check that all ω < wpe (this formula is only valid in this limit)
    if np.min(ω) < ω_pe:
        raise PhysicsError(
            "Lowest frequency must be larger than the electron "
            f"plasma frequency {ω_pe:.1e}, but min(ω) = {np.min(ω):.1e}"
        )

    # Check that the parameters given fall within the Rayleigh-Jeans limit
    # hω << kT_e
    rj_const = (
        np.max(ω) * const.hbar.si / (2 * np.pi * u.rad * const.k_B.si * T_e)
    ).to(u.dimensionless_unscaled)
    if rj_const.value > 0.1:

        raise PhysicsError(
            "Rayleigh-Jeans limit not satisfied: "
            f"ℏω/kT_e = {rj_const.value:.2e} > 0.1. "
            "Try lower ω or higher T_e."
        )

    # Calculate the bremsstrahlung power spectral density in several steps
    c1 = (
        (8 / 3)
        * np.sqrt(2 / np.pi)
        * (const.e.si**2 / (4 * np.pi * const.eps0.si)) ** 3
        * 1
        / (const.m_e.si * const.c.si**2) ** 1.5
    )

    Zi = ion_species.charge_number
    c2 = (
        np.sqrt(1 - ω_pe**2 / ω**2)
        * Zi**2
        * n_i
        * n_e
        / np.sqrt(const.k_B.si * T_e)
    )

    # Dimensionless argument for exponential integral
    arg = 0.5 * ω**2 * const.m_e.si / (kmax**2 * const.k_B.si * T_e) / u.rad**2
    # Remove units, get ndarray of values
    arg = (arg.to(u.dimensionless_unscaled)).value

    return c1 * c2 * exp1(arg)
