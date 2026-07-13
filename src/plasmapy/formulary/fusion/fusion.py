"""Fusion reaction cross-sections, reactivities, and power densities.

This module provides functions for calculating nuclear fusion reaction
rates in thermal plasmas using parameterized cross-sections and
reactivities from the literature.

References
----------
H.-S. Bosch and G. M. Hale, "Improved formulas for fusion cross-sections
and thermal reactivities," Nuclear Fusion, vol. 32, no. 4, pp. 611-631, 1992.

W. M. Nevins and R. Swain, "The thermonuclear fusion rate coefficient for
p-11B reactions," Nuclear Fusion, vol. 40, no. 4, pp. 865-872, 2000.
"""

__all__ = [
    "fusion_cross_section",
    "reactivity",
    "fusion_reaction_rate",
    "fusion_power_density",
]

import math

import astropy.units as u
import numpy as np

from plasmapy.formulary.fusion.parameters import (
    _ReactionParameters,
    _lookup_reaction,
)
from plasmapy.utils.decorators import validate_quantities

_DT_PEAK_REACTIVITY_KEV = 70.0
_DT_PEAK_REACTIVITY_VALUE = 1.1e-22


def _fusion_cross_section_bosch_hale(
    E_keV: float, params: _ReactionParameters
) -> float:
    r"""Bosch-Hale center-of-mass cross-section in m^2.

    Parameters
    ----------
    E_keV : float
        Center-of-mass energy in keV.

    params : _ReactionParameters
        Reaction parameters.

    Returns
    -------
    float
        Cross-section in m^2.
    """
    E = E_keV
    use_high = (
        params.high_energy_cutoff_keV > 0
        and E >= params.high_energy_cutoff_keV
    )
    max_E = params.high_energy_max_keV if use_high else params.max_E_keV
    if E < params.min_E_keV or E > max_E:
        return 0.0
    if use_high:
        A1, A2, A3, A4, A5 = (
            params.high_energy_A1,
            params.high_energy_A2,
            params.high_energy_A3,
            params.high_energy_A4,
            params.high_energy_A5,
        )
        B1, B2, B3, B4 = (
            params.high_energy_B1,
            params.high_energy_B2,
            params.high_energy_B3,
            params.high_energy_B4,
        )
    else:
        A1, A2, A3, A4, A5 = (
            params.A1,
            params.A2,
            params.A3,
            params.A4,
            params.A5,
        )
        B1, B2, B3, B4 = (
            params.B1,
            params.B2,
            params.B3,
            params.B4,
        )

    S = (
        A1
        + E * (A2 + E * (A3 + E * (A4 + E * A5)))
    ) / (1 + E * (B1 + E * (B2 + E * (B3 + E * B4))))

    return (1.0 / E) * S * math.exp(-params.BG / math.sqrt(E)) * 1e-31


def _cross_section_pb11(E_keV: float, params: _ReactionParameters) -> float:
    """Nevins-Swain p-B11 center-of-mass cross-section in m^2.

    Parameters
    ----------
    E_keV : float
        Center-of-mass energy in keV.

    params : _ReactionParameters
        Reaction parameters for p-B11.

    Returns
    -------
    float
        Cross-section in m^2.
    """
    if E_keV < params.min_E_keV or E_keV > params.max_E_keV:
        return 0.0

    C0, C1, C2 = 197.0, 0.240, 2.31e-4
    AL, EL, delEL = 1.82e4, 148.0, 2.35
    D0, D1, D2, D5 = 330.0, 66.1, -20.3, -1.58
    A0, A1, A2, A3 = 2.57e6, 5.67e5, 1.34e5, 5.68e5
    E0, E1, E2, E3 = 581.3, 1083.0, 2405.0, 3344.0
    delE0, delE1, delE2, delE3 = 85.7, 234.0, 138.0, 309.0
    B = 4.38
    EG = 22.589

    if E_keV <= 400:
        S = (
            C0
            + C1 * E_keV
            + C2 * (E_keV**2)
            + AL / ((E_keV - EL) ** 2 + delEL**2)
        )
    elif 400 < E_keV < 642:
        x = (E_keV - 400) / 100
        S = D0 + D1 * x + D2 * x**2 + D5 * x**5
    else:
        S = (
            A0 / ((E_keV - E0) ** 2 + delE0**2)
            + A1 / ((E_keV - E1) ** 2 + delE1**2)
            + A2 / ((E_keV - E2) ** 2 + delE2**2)
            + A3 / ((E_keV - E3) ** 2 + delE3**2)
            + B
        )

    return (S * 1000 / E_keV) * math.exp(-math.sqrt(EG * 1000 / E_keV)) * 1e-28


def _reactivity_bosch_hale(T_keV: float, params: _ReactionParameters) -> float:
    r"""Bosch-Hale parameterized reactivity in m^3/s.

    Parameters
    ----------
    T_keV : float
        Ion temperature in keV.

    params : _ReactionParameters
        Reaction parameters.

    Returns
    -------
    float
        Reactivity in m^3/s.
    """
    if T_keV < params.min_T_keV or T_keV > params.max_T_keV:
        return 0.0

    C1, C2, C3, C4, C5, C6, C7 = (
        params.C1,
        params.C2,
        params.C3,
        params.C4,
        params.C5,
        params.C6,
        params.C7,
    )

    theta = T_keV / (
        1.0
        - (
            T_keV
            * (C2 + T_keV * (C4 + T_keV * C6))
            / (1.0 + T_keV * (C3 + T_keV * (C5 + T_keV * C7)))
        )
    )
    xi = (params.BG**2 / (4.0 * theta)) ** (1.0 / 3.0)
    r_cm3_s = (
        C1
        * theta
        * math.sqrt(xi / (params.m1c2_keV * params.m2c2_keV * T_keV**3))
        * math.exp(-3.0 * xi)
    )

    return r_cm3_s * 1e-6


def _reactivity_nevins_swain_pb11(T_keV: float, params: _ReactionParameters) -> float:
    """Nevins-Swain parameterized p-B11 reactivity in m^3/s.

    Parameters
    ----------
    T_keV : float
        Ion temperature in keV.

    params : _ReactionParameters
        Reaction parameters for p-B11.

    Returns
    -------
    float
        Reactivity in m^3/s.
    """
    if T_keV < params.min_T_keV or T_keV > params.max_T_keV:
        return 0.0

    EG = 22.589
    P1, P2, P3, P4, P5, P6, P7 = (
        params.C1,
        params.C2,
        params.C3,
        params.C4,
        params.C5,
        params.C6,
        params.C7,
    )
    Mrc2 = params.m2c2_keV

    R = 5.41e-21 * T_keV ** (-1.5) * math.exp(-148.0 / T_keV)

    theta = T_keV / (
        1.0
        - (
            T_keV
            * (P2 + T_keV * (P4 + T_keV * P6))
            / (1.0 + T_keV * (P3 + T_keV * (P5 + T_keV * P7)))
        )
    )
    xi = (EG * 1000 / (4.0 * theta)) ** (1.0 / 3.0)
    NR_high_T = (
        P1
        * theta
        * math.sqrt(xi / (Mrc2 * T_keV**3))
        * math.exp(-3.0 * xi)
    )

    return NR_high_T + R


@validate_quantities(
    E={"units": u.keV, "can_be_negative": False, "can_be_zero": False},
)
def fusion_cross_section(
    reaction: str | tuple[str, str],
    E: u.Quantity[u.keV],
) -> u.Quantity[u.m**2]:
    r"""Return the center-of-mass fusion cross-section.

    The cross-section :math:`\sigma(E)` represents the effective area
    for a fusion reaction to occur when two reactants collide with a
    given center-of-mass energy.

    Parameters
    ----------
    reaction : `str` or `tuple` of `str`
        The fusion reaction. Can be specified as:
        - A reaction string: ``"D + T --> He-4 + n"``
        - A short name: ``"T(d,n)4He"``
        - A tuple of reactants: ``("D", "T")``

        For D-D reactions, the full reaction string or short name
        must be used to distinguish the two branches.

    E : `~astropy.units.Quantity`
        Center-of-mass energy in units convertible to keV.

    Returns
    -------
    sigma : `~astropy.units.Quantity`
        Fusion cross-section in units of m^2.

    Raises
    ------
    `ValueError`
        If the reaction is not recognized or the energy is out of
        the valid range for the parameterization.

    Notes
    -----
    For D + T, D + He-3, and D + D reactions, this function uses the
    Bosch-Hale parameterization :cite:p:`bosch-hale:1992`.

    For p + B-11, this function uses the Nevins-Swain
    parameterization :cite:p:`nevins-swain:2000`.

    Examples
    --------
    >>> import astropy.units as u
    >>> fusion_cross_section("D + T --> He-4 + n", 10 * u.keV)
    <Quantity 2.702...e-30 m2>
    """
    params = _lookup_reaction(reaction)
    E_keV = E.to(u.keV).value

    if params.name == "11B(p,4He)4He4He":
        sigma_m2 = _cross_section_pb11(E_keV, params)
    elif params.BG > 0:
        sigma_m2 = _fusion_cross_section_bosch_hale(E_keV, params)
    else:
        raise ValueError(
            f"Cannot compute cross-section for {params.name}. "
            f"No cross-section parameterization available."
        )

    if sigma_m2 == 0.0:
        E_min = params.min_E_keV
        E_max = (
            params.high_energy_max_keV
            if params.high_energy_cutoff_keV
            else params.max_E_keV
        )
        raise ValueError(
            f"Energy {E_keV:.1f} keV is outside the valid range "
            f"[{E_min}, {E_max}] keV for {params.name}."
        )

    return sigma_m2 * u.m**2


@validate_quantities(
    T={"units": u.keV, "can_be_negative": False},
)
def reactivity(
    reaction: str | tuple[str, str],
    T: u.Quantity[u.keV],
) -> u.Quantity[u.m**3 / u.s]:
    r"""Return the Maxwellian-averaged fusion reactivity.

    The reactivity :math:`\langle \sigma v \rangle` is the cross-section
    averaged over a thermal (Maxwellian) velocity distribution at
    ion temperature :math:`T`.

    Parameters
    ----------
    reaction : `str` or `tuple` of `str`
        The fusion reaction. See
        `~plasmapy.formulary.fusion.fusion.fusion_cross_section` for details.

    T : `~astropy.units.Quantity`
        Ion temperature in units convertible to keV.

    Returns
    -------
    sigma_v : `~astropy.units.Quantity`
        Fusion reactivity in units of m^3 / s.

    Raises
    ------
    `ValueError`
        If the reaction is not recognized or the temperature is
        outside the valid range.

    Notes
    -----
    For D + T, D + He-3, and D + D reactions, this function uses the
    Bosch-Hale parameterized reactivity :cite:p:`bosch-hale:1992`.

    For p + B-11, this function uses the Nevins-Swain parameterized
    reactivity :cite:p:`nevins-swain:2000`.

    The D-T reactivity peaks at approximately 70 keV with a value of
    about :math:`1.1 \times 10^{-22}` m\ :sup:`3`\ /s.

    Examples
    --------
    >>> import astropy.units as u
    >>> reactivity("D + T --> He-4 + n", 10 * u.keV)
    <Quantity 5.249...e-26 m3 / s>
    """
    params = _lookup_reaction(reaction)
    T_keV = T.to(u.keV).value

    if params.name == "11B(p,4He)4He4He":
        sigma_v = _reactivity_nevins_swain_pb11(T_keV, params)
    elif params.BG > 0:
        sigma_v = _reactivity_bosch_hale(T_keV, params)
    else:
        raise ValueError(
            f"Cannot compute reactivity for {params.name}."
        )

    if sigma_v == 0.0:
        raise ValueError(
            f"Temperature {T_keV:.3f} keV is outside the valid range "
            f"[{params.min_T_keV}, {params.max_T_keV}] keV "
            f"for {params.name}."
        )

    return sigma_v * u.m**3 / u.s


@validate_quantities(
    T={"units": u.keV, "can_be_negative": False},
    n1={"units": u.m**-3, "can_be_negative": False},
    n2={"units": u.m**-3, "can_be_negative": False},
)
def fusion_reaction_rate(
    reaction: str | tuple[str, str],
    T: u.Quantity[u.keV],
    n1: u.Quantity[u.m**-3],
    n2: u.Quantity[u.m**-3],
) -> u.Quantity[u.m**-3 / u.s]:
    r"""Return the volumetric fusion reaction rate.

    The reaction rate :math:`R = n_1 n_2 \langle \sigma v \rangle`
    gives the number of fusion reactions per unit volume per second.

    Parameters
    ----------
    reaction : `str` or `tuple` of `str`
        The fusion reaction. See
        `~plasmapy.formulary.fusion.fusion.fusion_cross_section` for details.

    T : `~astropy.units.Quantity`
        Ion temperature in units convertible to keV.

    n1 : `~astropy.units.Quantity`
        Number density of the first reactant in units convertible to
        m\ :sup:`-3`.

    n2 : `~astropy.units.Quantity`
        Number density of the second reactant in units convertible to
        m\ :sup:`-3`.

    Returns
    -------
    rate : `~astropy.units.Quantity`
        Fusion reaction rate in units of m\ :sup:`-3`\ s\ :sup:`-1`.

    Notes
    -----
    For identical particles (D + D), the rate includes a factor of
    :math:`1/(1 + \delta_{12})` to avoid double-counting.

    Examples
    --------
    >>> import astropy.units as u
    >>> fusion_reaction_rate(
    ...     "D + T --> He-4 + n",
    ...     10 * u.keV,
    ...     1e20 * u.m**-3,
    ...     1e20 * u.m**-3,
    ... )
    <Quantity 5.249...e+14 m-3 / s>
    """
    params = _lookup_reaction(reaction)
    sigma_v = reactivity(reaction, T)
    n1_val = n1.to(u.m**-3)
    n2_val = n2.to(u.m**-3)
    factor = 1.0 if not params.identical_particles else 2.0
    return n1_val * n2_val * sigma_v / factor


@validate_quantities(
    T={"units": u.keV, "can_be_negative": False},
    n1={"units": u.m**-3, "can_be_negative": False},
    n2={"units": u.m**-3, "can_be_negative": False},
)
def fusion_power_density(
    reaction: str | tuple[str, str],
    T: u.Quantity[u.keV],
    n1: u.Quantity[u.m**-3],
    n2: u.Quantity[u.m**-3],
) -> u.Quantity[u.W / u.m**3]:
    r"""Return the volumetric fusion power density.

    The power density :math:`P = R \times Q` is the product of the
    reaction rate and the energy released per reaction.

    Parameters
    ----------
    reaction : `str` or `tuple` of `str`
        The fusion reaction. See
        `~plasmapy.formulary.fusion.fusion.fusion_cross_section` for details.

    T : `~astropy.units.Quantity`
        Ion temperature in units convertible to keV.

    n1 : `~astropy.units.Quantity`
        Number density of the first reactant in units convertible to
        m\ :sup:`-3`.

    n2 : `~astropy.units.Quantity`
        Number density of the second reactant in units convertible to
        m\ :sup:`-3`.

    Returns
    -------
    power : `~astropy.units.Quantity`
        Fusion power density in units of W / m\ :sup:`3`.

    Examples
    --------
    >>> import astropy.units as u
    >>> fusion_power_density(
    ...     "D + T --> He-4 + n",
    ...     10 * u.keV,
    ...     1e20 * u.m**-3,
    ...     1e20 * u.m**-3,
    ... )
    <Quantity 1.480...e+03 W / m3>
    """
    params = _lookup_reaction(reaction)
    rate = fusion_reaction_rate(reaction, T, n1, n2)
    Q_energy = params.Q_keV * u.keV
    return (rate * Q_energy).to(u.W / u.m**3)
