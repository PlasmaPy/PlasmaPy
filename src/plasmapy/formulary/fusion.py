from plasmapy.utils import data 
import astropy.units as u
import numpy as np
from scipy.interpolate import CubicSpline
import json
from pathlib import Path
import h5py

"""
Opening Bosch and Hale Tables IV and Json Files
"""

DATA_DIR = data.__file__

with open(DATA_DIR / "bosch_hale_table_iv.json") as f:
    xs_coeffs = json.load(f)["reactions"]

with open(DATA_DIR / "bosch_hale_table_v.json") as f:
    table_v = json.load(f)

with open(DATA_DIR / "bosch_hale_table_vii.json") as f:
    rxty_coeffs = json.load(f)["reactions"]

with open(DATA_DIR / "bosch_hale_table_viii.json") as f:
    table_viii = json.load(f)

E_tab = np.array([row[0] for row in table_v["data"]]) * u.keV
sigma_tab = {
    name: np.array([row[i+1] for row in table_v["data"]]) * u.mbarn
    for i, name in enumerate(table_v["columns"][1:])
}

T_tab = np.array([row[0] for row in table_viii["data"]]) *u.keV
sv_tab = {
    name: np.array([row[i+1] for row in table_viii["data"]]) * (u.cm**3/ u.s)
    for i, name in enumerate(table_viii["columns"][1:])
}

available_reactions = [
    'D(t,n)A', 
    '3He(d,p)A', 
    'D(d,p)T', 
    'D(d,n)3He', 
    '3He(3He,2p)A', 
    '3He(t,n+p)A', 
    '3He(t,d)A', 
    'T(t,2n)A', 
    '11B(p,a)2A'
]

def cross_section(energy:u.Quantity, reaction:str, source:str):
    r"""
    Fusion cross-section :math:`\sigma(E)` from either the Bosch-Hale fit or ENDF data.

    Parameters
    ----------
    energy : `~astropy.units.Quantity`
        Center-of-mass kinetic energy. Must have units of energy (e.g. keV).
    reaction : str
        Reaction key. Must be one of the entries in ``available_reactions``:
        ``"D-T"``, ``"D-D_a"``, ``"D-D_b"``, ``"D-3He"``, ``"p-11B"``,
        ``"T-T"``, ``"T-3He_a"``, ``"T-3He_b"``, ``"3He-3He"``.
    type : {"BH", "ENDF"}
        Data source:

        - ``"BH"``   : Bosch-Hale (1992) Padé parametrization, valid for
          the four reactions in Table IV up to ~ 4900 keV.
        - ``"ENDF"`` : log-log cubic spline interpolation of tabulated
          ENDF/B data, available for all nine reactions above.

    Returns
    -------
    sigma : `~astropy.units.Quantity`
        Cross-section. Units are millibarn (``BH``) or m\ :sup:`2` (``ENDF``).

    Raises
    ------
    ValueError
        If ``reaction`` is not in ``available_reactions``, if ``type`` is
        not ``"BH"`` or ``"ENDF"``, or if ``energy`` falls outside the
        Bosch-Hale validity range (0 to 4900 keV).

    Examples
    --------
    >>> import astropy.units as u
    >>> cross_section(100 * u.keV, "D-T", "BH")    # doctest: +SKIP
    <Quantity 3427.245 mbarn>
    >>> cross_section(100 * u.keV, "D-T", "ENDF")  # doctest: +SKIP
    <Quantity [4.96e-28] m2>
    """
    if reaction not in available_reactions:
        raise ValueError(f"{reaction!r} is not in the available reactions allowed")
    if source not in ("BH", "ENDF"):
        raise ValueError(f"Unknown source {source!r}; expected 'BH' or 'ENDF'.")

    if source == "ENDF":
        return _ENDF_cross_section(energy, reaction)

    # BH
    if source == "BH":
        if reaction not in xs_coeffs:
            raise ValueError("does not have available Bosch and Hale coefficients")
        if not ((0 * u.keV < energy) & (energy < 4900 * u.keV)).all():
            raise ValueError(f"{energy!r} is not in Bosch and Hale Cross Sectional energy range of 0 to 4900 keV")
        return _BH_cross_section(energy, reaction)

def reactivity(ion_temp:u.Quantity, reaction:str, source:str):
    r"""
    Maxwellian reactivity :math:`\langle \sigma v \rangle (T)` from either Bosch-Hale or ENDF.

    Parameters
    ----------
    ion_temp : `~astropy.units.Quantity`
        Ion temperature. Must have units of energy (e.g. keV).
    reaction : str
        Reaction key. Must be one of the entries in ``available_reactions``.
    type : {"BH", "ENDF"}
        Data source:

        - ``"BH"``   : Closed-form Bosch-Hale (1992) reactivity fit
          (Eqs. 12-14), valid for the four reactions in Table VII up
          to ~ 190 keV.
        - ``"ENDF"`` : Direct numerical integration of the ENDF
          :math:`\sigma(E)` interpolant over a Maxwellian.

    Returns
    -------
    sv : `~astropy.units.Quantity`
        Maxwellian-averaged reactivity in cm\ :sup:`3`/s.

    Raises
    ------
    ValueError
        If ``reaction`` is not in ``available_reactions``, if ``type`` is
        not ``"BH"`` or ``"ENDF"``, or if ``ion_temp`` falls outside the
        Bosch-Hale validity range (0 to 190 keV).

    Notes
    -----
    The two backends should agree closely for the four overlapping
    reactions (D-T, D-\ :sup:`3`\ He, D(d,p)T, D(d,n)\ :sup:`3`\ He);
    a side-by-side comparison plot is the standard self-consistency check.

    Examples
    --------
    >>> import astropy.units as u
    >>> reactivity(10 * u.keV, "D-T", "BH")    # doctest: +SKIP
    <Quantity 1.13616547e-16 cm3 / s>
    >>> reactivity(10 * u.keV, "D-T", "ENDF")  # doctest: +SKIP
    <Quantity [1.13e-16] cm3 / s>
    """

    if reaction not in available_reactions:
        raise ValueError(f"{reaction} is not in the available reactions allowed")
    if source not in ("BH", "ENDF"):
        raise ValueError(f"Unknown source {source!r}; expected 'BH' or 'ENDF'.")
    
    if source == "ENDF":
        return _ENDF_reactivity(ion_temp, reaction)
    
    if not ((0 * u.keV < ion_temp) & (ion_temp < 190 * u.keV)).all():
        raise ValueError(f"{ion_temp!r} is not in Bosch and Hale Reactivity energy range of 0 to 190 keV")
    return _BH_reactivity(ion_temp, reaction)


mbarn = u.def_unit("mbarn", 1e-3 * u.barn)

"""
Cross Section Function and its Helper Functions
"""

def _get_xs_co(r):
    return xs_coeffs[r]

def _get_rxty(r):
    return rxty_coeffs[r]

def _pade_polynomial(rxn, e):
    r"""
    Evaluate the Bosch-Hale Padé approximant for the S-function.
    """
    S_vals = rxn["A1"] + e*(rxn["A2"] + e*(rxn["A3"] + e*(rxn["A4"] + e*rxn["A5"])))
    S_vals /= 1 + e*(rxn["B1"] + e*(rxn["B2"] + e*(rxn["B3"] + e*rxn["B4"])))
    return S_vals

def _parametrization_formula(S_Vals, energy, rxn): #B&H Eq (8)
    r"""
    Combine the S-function with the Gamow penetrability to get a cross-section.
    """
    sigma = S_Vals/(energy*np.exp(rxn["B_G"]/np.sqrt(energy)))
    return sigma

@u.quantity_input
def _BH_cross_section(energy:u.Quantity, reaction:str):
    r"""
    Compute the fusion cross-section using the Bosch-Hale Padé parametrization.
    """
    rxn = _get_xs_co(reaction)
    E_keV = energy.to(u.keV).value
    S = _pade_polynomial(rxn, E_keV)
    sigma = _parametrization_formula(S, E_keV, rxn)
    return sigma * u.mbarn

"""
Reactivity Functions and its helpers
"""

def _get_polynomial(T, r):
    r"""
    Evaluate the :math:`\theta(T)` Padé approximant for the Bosch-Hale reactivity.
    """
    theta = T*(r["C2"] + T*(r["C4"] + T*r["C6"]))
    theta /= 1 + T*(r["C3"] + T*(r["C5"] + T*r["C7"]))
    theta = 1 - theta
    theta = T/theta
    return theta

def _get_xi(Theta, r):
    r"""
    Compute the :math:`\xi` factor for the Bosch-Hale reactivity formula.
    """
    xi = (r["B_G"])**2
    xi /= 4*Theta
    xi = (xi)**(1/3)
    return xi

def _get_reactivity(T, r, theta, xi):
    r"""
    Assemble the Bosch-Hale reactivity from its precomputed pieces.
    """
    rcty = np.sqrt(xi/((r["m_r_c2"]*T**3)))
    rcty *= r["C1"]*theta*(np.exp(-3*xi))
    return rcty

@u.quantity_input
def _BH_reactivity(ion_temp:u.Quantity, reaction:str):
    r"""
    Compute the Maxwellian fusion reactivity using the Bosch-Hale parametrization.
    """
    rxn = _get_rxty(reaction)
    T_keV = ion_temp.to(u.keV).value
    Theta = _get_polynomial(T_keV, rxn)
    xi = _get_xi(Theta, rxn)
    sv = _get_reactivity(T_keV, rxn, Theta, xi)
    return sv * (u.cm**3 / u.s)

amu = 1.66053906660e-27           # atomic mass unit [kg]
e   = 1.602176634e-19             # elementary charge [C]

h5_files = {
    'D(t,n)A': 'T(D,n)A.h5',
    '3He(d,p)A': 'He-3(D,p)A.h5',
    'D(d,p)T': 'D(D,p)T.h5',
    'D(d,n)3He': 'D(D,n)He-3.h5',
    '3He(3He,2p)A': 'He-3(He-3,2p)A.h5',
    '3He(t,n+p)A': 'He-3(T,n+p)A.h5',
    '3He(t,d)A': 'He-3(T,D)A.h5',
    'T(t,2n)A': 'T(T,2n)A.h5',
    '11B(p,a)2A': 'B-11(p,He-4)2He-4.h5'
}

new_rxns = list(h5_files)

def _load_h5(path):
     
    with h5py.File(path, "r") as f:
                sigma = f["SIG"][:] # units = m^2
                E = f["energy"][:] # units = eV, CoM reference frame
    return E, sigma

def _load_reaction(rxn_key):
    E, sigma = _load_h5(DATA_DIR / h5_files[rxn_key])
    return E * 1e-3, sigma


def _build_Xsec_interpolation(rxn_key):
    r"""
    Build a log-log cubic spline interpolant of the ENDF cross-section data.
    """
    E_kev, sigma = _load_reaction(rxn_key)

    E_unique, idx = np.unique(E_kev, return_index=True) #raise an exception non unique values
    sigma_unique = sigma[idx]

    mask = sigma_unique > 0
    E_pos = E_unique[mask]
    sigma_pos = sigma_unique[mask]

    spline = CubicSpline(np.log(E_pos), np.log(sigma_pos), extrapolate=False)
    return spline

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

_Egrid_keV = np.logspace(0, 5, 1000)        # internal, plain floats
Egrid = _Egrid_keV * u.keV                  # for the API/plotting
masses = {'D': 2.014, 'T': 3.016, '3He': 3.016, '11B': 11.009305167,
          'p': 1.007276466620409}

def _find_mu(rxn):
    r"""
    Compute the reduced mass of the two reactants in a reaction key.
    """
    target, rest = rxn.split('(', 1)
    beam = rest.split(',', 1)[0]
    alias = {'d': 'D', 't': 'T', 'p': 'p', '3He': '3He'}
    beam = alias.get(beam, beam)
    m1 = masses[target] * amu
    m2 = masses[beam] * amu
    return m1 * m2 / (m1 + m2)

@u.quantity_input(T=u.keV)
def _ENDF_reactivity(T, rxn_key):
    r"""
    Compute the Maxwellian fusion reactivity by integrating ENDF :math:`\sigma(E)`.
    """
    T_keV = np.atleast_1d(T.to(u.keV).value)
    mu = _find_mu(rxn_key)

    cs = _build_Xsec_interpolation(rxn_key)
    sigma = np.nan_to_num(np.exp(cs(np.log(_Egrid_keV))), nan=0.0)

    E_col = _Egrid_keV[:, None]
    T_row = T_keV[None, :]
    integrand = sigma[:, None] * E_col * np.exp(-E_col / T_row)
    I = np.trapezoid(integrand, _Egrid_keV, axis=0)

    fac = 4.0 / np.sqrt(2.0 * np.pi * mu)
    fac = fac * (1000.0 * e)**2 / (1000.0 * e * T_keV)**1.5

    return (fac * I) * (u.m**3 / u.s) 