"""Cross sections and Maxwellian reactivities for nuclear fusion reactions."""

import json
from importlib.resources import as_file, files
from pathlib import Path

import astropy.constants as const
import astropy.units as u
import h5py
import numpy as np
from scipy.interpolate import CubicSpline

from plasmapy.utils.decorators import validate_quantities

"""
Opening Bosch and Hale Tables IV and Json Files
"""

DATA_DIR = files("plasmapy.utils.data").joinpath("config.json")

with as_file(DATA_DIR) as physical_path:
    with Path.open(physical_path / "bosch_hale_table_iv.json") as f:
        xs_coeffs = json.load(f)["reactions"]

    with Path.open(physical_path / "bosch_hale_table_v.json") as f:
        table_v = json.load(f)

    with Path.open(physical_path / "bosch_hale_table_vii.json") as f:
        rxty_coeffs = json.load(f)["reactions"]

    with Path.open(physical_path / "bosch_hale_table_viii.json") as f:
        table_viii = json.load(f)


@validate_quantities
def cross_section(energy: u.Quantity[u.keV], reaction: str, source: str) -> float:
    r"""
    Calculate the fusion cross section :math:`σ(E)` for a two-body fusion
    reaction.

    Parameters
    ----------
    energy : `~astropy.units.Quantity`
        The center-of-mass kinetic energy of the reactants, in units
        convertible to keV.

    reaction : `str`
        The fusion reaction to evaluate. All nine reactions below are
        available when ``source`` is ``"ENDF"``\ ; only the first four
        are available when ``source`` is ``"BH"``\ :

        - ``"D(t,n)A"`` — :math:`D + T → n + α`
        - ``"3He(d,p)A"`` — :math:`^3He + D → p + α`
        - ``"D(d,p)T"`` — :math:`D + D → p + T`
        - ``"D(d,n)3He"`` — :math:`D + D → n + ^3He`
        - ``"3He(3He,2p)A"`` — :math:`^3He + ^3He → 2p + α`
        - ``"3He(t,n+p)A"`` — :math:`^3He + T → n + p + α`
        - ``"3He(t,d)A"`` — :math:`^3He + T → D + α`
        - ``"T(t,2n)A"`` — :math:`T + T → 2n + α`
        - ``"11B(p,a)2A"`` — :math:`^{11}B + p → 3α`

    source : {``"BH"``, ``"ENDF"``}
        The data source used to evaluate the cross section:

        - ``"BH"`` — the Padé parametrization of :cite:t:`bosch:1992`,
          valid for the four reactions in their Table IV over the energy
          range given for each reaction (at most :math:`0` - :math:`4900`
          keV).
        - ``"ENDF"`` — a log-log cubic spline interpolation of the
          tabulated ENDF/B evaluated cross sections, available for all
          nine reactions.

    Returns
    -------
    sigma : `~astropy.units.Quantity`
        The fusion cross section, in units of millibarn when ``source``
        is ``"BH"`` and m\ :sup:`2` when ``source`` is ``"ENDF"``\ .

    Raises
    ------
    `ValueError`
        If ``reaction`` is not available for the requested ``source``\ ,
        if ``source`` is neither ``"BH"`` nor ``"ENDF"``\ , or if
        ``energy`` falls outside the Bosch-Hale validity range for
        ``reaction``\ .

    `~astropy.units.UnitTypeError`
        If ``energy`` does not have units convertible to keV.

    See Also
    --------
    reactivity

    Notes
    -----
    The Bosch-Hale parametrization writes the cross section in terms of
    the astrophysical :math:`S`\ -function and the Gamow penetrability
    factor as

    .. math::

        σ(E) = \frac{S(E)}{E \exp\left(B_G / \sqrt{E}\right)},

    where :math:`B_G` is the Gamow constant and :math:`S(E)` is a Padé
    approximant of the form

    .. math::

        S(E) = \frac{A_1 + E\left(A_2 + E\left(A_3 + E\left(A_4 +
        E A_5\right)\right)\right)}{1 + E\left(B_1 + E\left(B_2 +
        E\left(B_3 + E B_4\right)\right)\right)}.

    The ``"ENDF"`` interpolant is constructed in log-log space and is not
    extrapolated; energies outside the tabulated range return zero.

    Examples
    --------
    >>> import astropy.units as u
    >>> cross_section(100 * u.keV, "D(t,n)A", "BH")  # doctest: +SKIP
    <Quantity 3427.245 mbarn>
    >>> cross_section(100 * u.keV, "D(t,n)A", "ENDF")  # doctest: +SKIP
    <Quantity [4.96e-28] m2>
    """
    if source == "ENDF":
        if reaction not in ENDF_rxns:
            raise ValueError("Reaction in supported ENDF reactions")

        return _ENDF_cross_section(energy, reaction)

    if source == "BH":
        if reaction not in xs_coeffs:
            raise ValueError("does not have available Bosch and Hale coefficients")
        if not _in_BH_rxn_energy_range(energy, reaction):
            raise ValueError(
                f"{energy!r} is not in Bosch and Hale Cross Sectional energy range of 0 to 4900 keV"
            )

        return _BH_cross_section(energy, reaction)

    raise ValueError(f"Unknown source {source!r}; expected 'BH' or 'ENDF'.")


@validate_quantities
def reactivity(ion_temp: u.Quantity[u.keV], reaction: str, source: str) -> float:
    r"""
    Calculate the Maxwellian-averaged fusion reactivity :math:`⟨σv⟩(T)`
    for a two-body fusion reaction.

    Parameters
    ----------
    ion_temp : `~astropy.units.Quantity`
        The ion temperature, in units convertible to keV.

    reaction : `str`
        The fusion reaction to evaluate. All nine reactions listed in
        `cross_section` are available when ``source`` is ``"ENDF"`` ;
        only ``"D(t,n)A"`` , ``"3He(d,p)A"`` , ``"D(d,p)T"`` , and
        ``"D(d,n)3He"`` are available when ``source`` is ``"BH"`` .

    source : {``"BH"``, ``"ENDF"``}
        The data source used to evaluate the reactivity:

        - ``"BH"`` — the closed-form fit of :cite:t:`bosch:1992` (their
          Eqs. 12-14), valid for the four reactions in their Table VII
          over the temperature range given for each reaction (at most
          :math:`0` - :math:`190` keV).
        - ``"ENDF"`` — direct numerical integration of the ENDF
          :math:`σ(E)` interpolant over a Maxwellian distribution.

    Returns
    -------
    sv : `~astropy.units.Quantity`
        The Maxwellian-averaged reactivity, in units of cm\ :sup:`3`
        s\ :sup:`-1` when ``source`` is ``"BH"`` and m\ :sup:`3`
        s\ :sup:`-1` when ``source`` is ``"ENDF"``\ .

    Raises
    ------
    `ValueError`
        If ``reaction`` is not available for the requested ``source``\ ,
        if ``source`` is neither ``"BH"`` nor ``"ENDF"``\ , or if
        ``ion_temp`` falls outside the Bosch-Hale validity range for
        ``reaction``\ .

    `~astropy.units.UnitTypeError`
        If ``ion_temp`` does not have units convertible to keV.

    See Also
    --------
    cross_section

    Notes
    -----
    For a Maxwellian distribution of relative velocities, the reactivity
    is the average of :math:`σ(E) v` over the distribution,

    .. math::

        ⟨σv⟩(T) = \sqrt{\frac{8}{π μ}} \frac{1}{(k_B T)^{3/2}}
        \int_0^{∞} σ(E)\, E \exp\left(\frac{-E}{k_B T}\right) dE,

    where :math:`μ` is the reduced mass of the two reactants. The
    ``"ENDF"`` backend evaluates this integral numerically on a
    logarithmic energy grid.

    The ``"BH"`` backend instead evaluates the closed-form fit

    .. math::

        ⟨σv⟩(T) = C_1 θ \sqrt{\frac{ξ}{m_r c^2 T^3}} \exp(-3ξ),
        \qquad ξ = \left(\frac{B_G^2}{4θ}\right)^{1/3},

    with :math:`θ(T)` given by a Padé approximant in :math:`T`.

    The two backends should agree closely for the four overlapping
    reactions; a side-by-side comparison is the standard
    self-consistency check.

    Examples
    --------
    >>> import astropy.units as u
    >>> reactivity(10 * u.keV, "D(t,n)A", "BH")  # doctest: +SKIP
    <Quantity 1.13616547e-16 cm3 / s>
    >>> reactivity(10 * u.keV, "D(t,n)A", "ENDF")  # doctest: +SKIP
    <Quantity [1.13e-22] m3 / s>
    """
    if source == "ENDF":
        if reaction not in ENDF_rxns:
            raise ValueError("Reaction in supported ENDF reactions")

        return _ENDF_reactivity(ion_temp, reaction)

    if source == "BH":
        if reaction not in rxty_coeffs:
            raise ValueError("does not have available Bosch and Hale coefficients")
        if not _in_BH_rxn_ion_temp_range(ion_temp, reaction):
            raise ValueError(
                f"{ion_temp!r} is not in Bosch and Hale Reactivity energy range of 0 to 190 keV"
            )

        return _BH_reactivity(ion_temp, reaction)

    raise ValueError(f"Unknown source {source!r}; expected 'BH' or 'ENDF'.")


def _in_BH_rxn_energy_range(E, reaction):
    rxn = _xs_coeff(reaction)
    in_range = (rxn["E_min_keV"] * u.keV <= E) & (rxn["E_max_keV"] * u.keV >= E)
    return bool(np.all(in_range))


def _in_BH_rxn_ion_temp_range(T, reaction):
    rxn = _rxty_coeff(reaction)
    in_range = (rxn["T_min_keV"] * u.keV <= T) & (rxn["T_max_keV"] * u.keV >= T)
    return bool(np.all(in_range))


def _xs_coeff(r):
    return xs_coeffs[r]


def _rxty_coeff(r):
    return rxty_coeffs[r]


"""
Cross Section Function and its Helper Functions
"""


def _pade_polynomial(rxn, E):
    r"""
    Evaluate the Bosch-Hale Padé approximant for the S-function.
    """
    S_vals = rxn["A1"] + E * (
        rxn["A2"] + E * (rxn["A3"] + E * (rxn["A4"] + E * rxn["A5"]))
    )
    S_vals /= 1 + E * (rxn["B1"] + E * (rxn["B2"] + E * (rxn["B3"] + E * rxn["B4"])))
    return S_vals


def _parametrization_formula(S_Vals, energy, rxn):  # B&H Eq (8)
    r"""
    Combine the S-function with the Gamow penetrability to get a cross-section.
    """
    return S_Vals / (energy * np.exp(rxn["B_G"] / np.sqrt(energy)))


@u.quantity_input
def _BH_cross_section(energy: u.Quantity, reaction: str) -> float:
    r"""
    Compute the fusion cross-section using the Bosch-Hale Padé parametrization.
    """
    rxn = _xs_coeff(reaction)
    E_keV = energy.to(u.keV).value
    S = _pade_polynomial(rxn, E_keV)
    sigma = _parametrization_formula(S, E_keV, rxn)
    return sigma * u.mbarn


"""
Reactivity Functions and its helpers
"""


def _rxty_polynomial(T, r):
    r"""
    Evaluate the :math:`\theta(T)` Padé approximant for the Bosch-Hale
    reactivity.
    """
    theta = T * (r["C2"] + T * (r["C4"] + T * r["C6"]))
    theta /= 1 + T * (r["C3"] + T * (r["C5"] + T * r["C7"]))
    theta = 1 - theta
    return T / theta


def _find_xi(Theta, r):
    r"""
    Compute the :math:`\xi` factor for the Bosch-Hale reactivity formula.
    """
    xi = (r["B_G"]) ** 2
    xi /= 4 * Theta
    return (xi) ** (1 / 3)


def _find_reactivity(T, r, theta, xi):
    r"""
    Assemble the Bosch-Hale reactivity from its precomputed pieces.
    """
    rcty = np.sqrt(xi / (r["m_r_c2"] * T**3))
    rcty *= r["C1"] * theta * (np.exp(-3 * xi))
    return rcty


@u.quantity_input
def _BH_reactivity(ion_temp: u.Quantity, reaction: str) -> float:
    r"""
    Compute the Maxwellian fusion reactivity using the Bosch-Hale
    parametrization.
    """
    rxn = _rxty_coeff(reaction)
    T_keV = ion_temp.to(u.keV).value
    Theta = _rxty_polynomial(T_keV, rxn)
    xi = _find_xi(Theta, rxn)
    sv = _find_reactivity(T_keV, rxn, Theta, xi)
    return sv * (u.cm**3 / u.s)


amu = const.u.si.value  # atomic mass unit [kg]
e = const.e.si.value  # elementary charge [C]

h5_files = {
    "D(t,n)A": "T(D,n)A.h5",
    "3He(d,p)A": "He-3(D,p)A.h5",
    "D(d,p)T": "D(D,p)T.h5",
    "D(d,n)3He": "D(D,n)He-3.h5",
    "3He(3He,2p)A": "He-3(He-3,2p)A.h5",
    "3He(t,n+p)A": "He-3(T,n+p)A.h5",
    "3He(t,d)A": "He-3(T,D)A.h5",
    "T(t,2n)A": "T(T,2n)A.h5",
    "11B(p,a)2A": "B-11(p,He-4)2He-4.h5",
}

ENDF_rxns = list(h5_files)


def _load_h5(path):

    with h5py.File(path, "r") as f:
        sigma = f["SIG"][:]  # units = m^2
        E = f["energy"][:]  # units = eV, CoM reference frame
    return E, sigma


def _load_reaction(rxn_key):
    E, sigma = _load_h5(DATA_DIR / h5_files[rxn_key])
    return E * 1e-3, sigma


def _build_Xsec_interpolation(rxn_key):
    r"""
    Build a log-log cubic spline interpolant of the ENDF cross-section data.
    """
    E_kev, sigma = _load_reaction(rxn_key)

    E_unique, idx = np.unique(
        E_kev, return_index=True
    )  # raise an exception non unique values
    sigma_unique = sigma[idx]

    mask = sigma_unique > 0
    E_pos = E_unique[mask]
    sigma_pos = sigma_unique[mask]

    return CubicSpline(np.log(E_pos), np.log(sigma_pos), extrapolate=False)


@u.quantity_input(energy=u.keV)
def _ENDF_cross_section(energy, rxn_key):
    r"""
    Compute the fusion cross-section by interpolating tabulated ENDF data.
    """
    E_keV = np.atleast_1d(energy.to(u.keV).value)
    cs = _build_Xsec_interpolation(rxn_key)
    sigma = np.exp(cs(np.log(E_keV)))
    sigma = np.nan_to_num(sigma, nan=0.0)
    return sigma * u.m**2


_Egrid_keV = np.logspace(0, 5, 1000)  # internal, plain floats
Egrid = _Egrid_keV * u.keV  # for the API/plotting
masses = {
    "D": 2.014,
    "T": 3.016,
    "3He": 3.016,
    "11B": 11.009305167,
    "p": const.m_p.to_value(u.u),
}


def _find_mu(rxn):
    r"""
    Compute the reduced mass of the two reactants in a reaction key.
    """
    target, rest = rxn.split("(", 1)
    beam = rest.split(",", 1)[0]
    alias = {"d": "D", "t": "T", "p": "p", "3He": "3He"}
    beam = alias.get(beam, beam)
    m1 = masses[target] * amu
    m2 = masses[beam] * amu
    return m1 * m2 / (m1 + m2)


@u.quantity_input(T=u.keV)
def _ENDF_reactivity(T, rxn_key):
    r"""
    Compute the Maxwellian fusion reactivity by integrating ENDF:
    math:`\sigma(E)`.
    """
    T_keV = np.atleast_1d(T.to(u.keV).value)
    mu = _find_mu(rxn_key)

    cs = _build_Xsec_interpolation(rxn_key)
    sigma = np.nan_to_num(np.exp(cs(np.log(_Egrid_keV))), nan=0.0)

    E_col = _Egrid_keV[:, None]
    T_row = T_keV[None, :]
    integrand = sigma[:, None] * E_col * np.exp(-E_col / T_row)
    Integration = np.trapezoid(integrand, _Egrid_keV, axis=0)

    fac = 4.0 / np.sqrt(2.0 * np.pi * mu)
    fac = fac * (1000.0 * e) ** 2 / (1000.0 * e * T_keV) ** 1.5

    return (fac * Integration) * (u.m**3 / u.s)
