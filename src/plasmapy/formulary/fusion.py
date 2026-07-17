"""Cross sections and Maxwellian reactivities for nuclear fusion reactions."""

__all__ = ["fusion_cross_section", "fusion_reactivity"]

import json
from importlib.resources import as_file, files
from pathlib import Path

import astropy.units as u
import numpy as np

from plasmapy.utils.decorators import validate_quantities

"""
Opening Bosch and Hale Tables IV and Json Files
"""

DATA_DIR = files("plasmapy.utils.data")

with as_file(DATA_DIR) as physical_path:
    with Path.open(physical_path / "bosch_hale_ENDF_xs_table.json") as f:
        xs_coeffs = json.load(f)["reactions"]

    with Path.open(physical_path / "bosch_hale_table_v.json") as f:
        table_v = json.load(f)

    with Path.open(physical_path / "bosch_hale_ENDF_rxty_table.json") as f:
        rxty_coeffs = json.load(f)["reactions"]

    with Path.open(physical_path / "bosch_hale_table_viii.json") as f:
        table_viii = json.load(f)


@validate_quantities
def fusion_cross_section(
    energy: u.Quantity[u.keV], reaction: str
) -> u.Quantity[u.m**2]:
    r"""
    Calculate the fusion cross section :math:`σ(E)` for a two reactant fusion
    reaction.

    The cross section is evaluated from the Bosch-Hale Padé
    parametrization :cite:p:`bosch:1992`. For the reactions tabulated by
    Bosch and Hale the published coefficients are used; for the remaining
    reactions the same functional form is fit (Bosch and Hale eq 8 and 9)
    to ENDF/B evaluated cross sections (see Notes).

    Parameters
    ----------
    energy : `~astropy.units.Quantity`
        The center-of-mass kinetic energy of the reactants, in units
        convertible to keV.

    reaction : `str`
        The fusion reaction to evaluate, given as one of the following
        keys:

        - ``"D(t,n)A"`` — :math:`D + T → n + α` (Bosch and Hale fit)
                        — :energy range: 0.5keV to 550keV
        - ``"3He(d,p)A"`` — :math:`^3He + D → p + α` (Bosch and Hale fit)
                          — :energy range: 0.3keV to 900keV
        - ``"D(d,p)T"`` — :math:`D + D → p + T` (Bosch and Hale fit)
                        — :energy range: 0.5keV to 5000keV
        - ``"D(d,n)3He"`` — :math:`D + D → n + ^3He` (Bosch and Hale fit)
                          — :energy range: 0.5keV to 4900keV
        - ``"3He(3He,2p)A"`` — :math:`^3He + ^3He → 2p + α` (ENDF data fit)
                             — :energy range: 1.0keV to 10000keV
        - ``"3He(t,n+p)A"`` — :math:`^3He + T → n + p + α` (ENDF data fit)
                            — :energy range: 1.0keV to 10000keV
        - ``"3He(t,d)A"`` — :math:`^3He + T → D + α` (ENDF data fit)
                          — :energy range: 1.0keV to 10000keV
        - ``"T(t,2n)A"`` — :math:`T + T → 2n + α` (ENDF data fit)
                         — :energy range: 0.5keV to 9000keV
        - ``"11B(p,a)2A"`` — :math:`^{11}B + p → 3α` (ENDF data fit)
                           — :energy range: 200keV to 5000keV

        The first four use the coefficients published in Table IV of
        :cite:t:`bosch:1992`\ ; the remainder use coefficients obtained
        by fitting the identical Padé form to ENDF/B data, as described
        in the Notes. A reaction is available only if its coefficients
        are present in the loaded table.

    Returns
    -------
    sigma : `~astropy.units.Quantity`
        The fusion cross section, in units of millibarn.

    Raises
    ------
    `ValueError`
        If ``reaction`` has no available coefficients, or if ``energy``
        falls outside the validity range for ``reaction``.

    `~astropy.units.UnitTypeError`
        If ``energy`` does not have units convertible to keV.

    See Also
    --------
    fusion_reactivity

    Notes
    -----
    The Bosch-Hale parametrization writes the cross section in terms of
    the astrophysical :math:`S`\ -function and the Gamow penetrability
    factor as

    .. math::

        σ(E) = \frac{S(E)}{E \exp\left(B_G / \sqrt{E}\right)},

    where :math:`B_G` is the :wikipedia:`Gamow factor` for the reaction and
    :math:`S(E)` is the slowly varying astrophysical :math:`S`\ -function,
    represented by the Padé approximant

    .. math::

        S(E) = \frac{A_1 + E\left(A_2 + E\left(A_3 + E\left(A_4 +
        E A_5\right)\right)\right)}{1 + E\left(B_1 + E\left(B_2 +
        E\left(B_3 + E B_4\right)\right)\right)}.

    For the four reactions in Table IV of :cite:t:`bosch:1992` the
    coefficients :math:`A_i` and :math:`B_j` are taken directly from the
    published fit. For the remaining reactions the coefficients are
    obtained by fitting this same functional form to evaluated cross
    sections from the ENDF/B library, as demonstrated below:

    Starting from the tabulated ENDF cross section :math:`σ_{ENDF}(E)`,
    the strong exponential energy dependence is first removed by
    inverting the relation above to recover the :math:`S`\ -function,

    .. math::

        S(E) = σ_{ENDF}(E)\, E \exp\left(B_G / \sqrt{E}\right),

    using the analytically known :math:`B_G` for that reaction. Because
    :math:`S(E)` is smooth and slowly varying — unlike :math:`σ(E)`,
    which spans many orders of magnitude near threshold — it is far
    better conditioned for a rational-function fit. The Padé coefficients
    are then found by nonlinear least squares
    (`scipy.optimize.curve_fit`) applied to the :math:`S(E)` samples,

    .. math::

        \min_{A_i,\, B_j} \sum_k
        \left[S_{\mathrm{model}}(E_k; A_i, B_j) - S(E_k)\right]^2 .

    A few numerical steps keep this fit well behaved. Energies at or
    below the reaction threshold, where :math:`σ_{ENDF}(E) = 0`, are
    masked out before the :math:`S`\ -function is formed, since the
    inversion is undefined there. The nonlinear solver is seeded from a
    linearized polynomial pre-fit — cross-multiplying the Padé form,

    .. math::

        S(E)\left[1 + E\left(B_1 + \cdots\right)\right]
        = A_1 + E\left(A_2 + \cdots\right),

    gives a problem that is linear in the :math:`A_i` and :math:`B_j`,
    whose least-squares solution provides a good starting point near the
    true minimum. Finally, the numerator and denominator coefficients
    differ by many orders of magnitude, so they are scaled before fitting
    to keep the Jacobian columns comparable and avoid an
    ill-conditioned solve. The resulting :math:`A_i` and :math:`B_j` are
    stored in the same format as the published coefficients and evaluated
    through exactly the same code path, so ``fusion_cross_section`` behaves
    identically for fitted and published reactions.

    Examples
    --------
    >>> import astropy.units as u
    >>> fusion_cross_section(100 * u.keV, "D(t,n)A")  # doctest: +SKIP
    <Quantity 3427.245 mbarn>
    """
    if reaction not in xs_coeffs:
        raise ValueError(
            f"{reaction!r} is not one of the available reactions: "
            f"{', '.join(xs_coeffs)}"
        )
    if not _in_BH_rxn_energy_range(energy, reaction):
        rxn = _xs_coeff(reaction)
        raise ValueError(
            f"{energy!r} is not in the {reaction!r} energy range "
            f"of {rxn['E_min_keV']} to {rxn['E_max_keV']} keV"
        )

    return _BH_cross_section(energy, reaction)


@validate_quantities
def fusion_reactivity(
    ion_temp: u.Quantity[u.keV], reaction: str
) -> u.Quantity[u.m**3 / u.s]:
    r"""
    Calculate the Maxwellian-averaged fusion reactivity :math:`⟨σv⟩(T)`
    for a two-body fusion reaction.

    The reactivity is evaluated from the closed-form fit of
    :cite:t:`bosch:1992` (their Eqs. 12-14), which reproduces the
    Maxwellian average of the cross section without numerical
    integration. For the reactions tabulated by Bosch and Hale the
    published coefficients are used; for the remaining reactions the same
    functional form is fit to reactivities derived from the ENDF/B cross
    sections (see Notes).

    Parameters
    ----------
    ion_temp : `~astropy.units.Quantity`
        The ion temperature, in units convertible to keV.

    reaction : `str`
        The fusion reaction to evaluate, given as one of the following
        keys:

        - ``"D(t,n)A"`` — :math:`D + T → n + α` (Bosch and Hale fit)
                        — :ion temp range: 0.2keV to 100keV
        - ``"3He(d,p)A"`` — :math:`^3He + D → p + α` (Bosch and Hale fit)
                          — :ion temp range: 0.5keV to 190keV
        - ``"D(d,p)T"`` — :math:`D + D → p + T` (Bosch and Hale fit)
                        — :ion temp range: 0.2keV to 100keV
        - ``"D(d,n)3He"`` — :math:`D + D → n + ^3He` (Bosch and Hale fit)
                          — :ion temp range: 0.2keV to 100keV
        - ``"3He(3He,2p)A"`` — :math:`^3He + ^3He → 2p + α` (ENDF data fit)
                             — :ion temp range: 8.27keV to 100keV
        - ``"3He(t,n+p)A"`` — :math:`^3He + T → n + p + α` (ENDF data fit)
                            — :ion temp range: 1.0keV to 100keV
        - ``"3He(t,d)A"`` — :math:`^3He + T → D + α` (ENDF data fit)
                          — :ion temp range: 1.0keV to 100keV
        - ``"T(t,2n)A"`` — :math:`T + T → 2n + α` (ENDF data fit)
                         — :ion temp range: 1.0keV to 100keV
        - ``"11B(p,a)2A"`` — :math:`^{11}B + p → 3α` (ENDF data fit)
                           — :ion temp range: 50keV to 500keV

        The first four use the reactivity coefficients published in
        Table VII of :cite:t:`bosch:1992`\ ; the remainder use
        coefficients obtained by fitting the identical closed form to
        reactivities derived from ENDF/B data, as described in the Notes.
        A reaction is available only if its coefficients are present in
        the loaded table.

    Returns
    -------
    sv : `~astropy.units.Quantity`
        The Maxwellian-averaged reactivity, in units of cm\ :sup:`3`
        s\ :sup:`-1`\ .

    Raises
    ------
    `ValueError`
        If ``reaction`` has no available reactivity coefficients, or if
        ``ion_temp`` falls outside the validity range for ``reaction``
        (at most :math:`0` - :math:`190` keV).

    `~astropy.units.UnitTypeError`
        If ``ion_temp`` does not have units convertible to keV.

    See Also
    --------
    fusion_cross_section

    Notes
    -----
    For a Maxwellian distribution of relative velocities, the reactivity
    is the average of :math:`σ(E) v` over the distribution,

    .. math::

        ⟨σv⟩(T) = \sqrt{\frac{8}{π μ}} \frac{1}{(k_B T)^{3/2}}
        \int_0^{∞} σ(E)\, E \exp\left(\frac{-E}{k_B T}\right) dE,

    where :math:`μ` is the reduced mass of the two reactants. Rather than
    evaluating this integral, this function uses the closed-form fit of
    :cite:t:`bosch:1992`,

    .. math::

        ⟨σv⟩(T) = C_1 θ \sqrt{\frac{ξ}{m_r c^2 T^3}} \exp(-3ξ),
        \qquad ξ = \left(\frac{B_G^2}{4θ}\right)^{1/3},

    where the temperature-dependent factor :math:`θ(T)` is itself a Padé
    approximant in :math:`T`,

    .. math::

        θ(T) = \frac{T}{\,1 - \dfrac{T\left(C_2 + T\left(C_4 +
        T C_6\right)\right)}{1 + T\left(C_3 + T\left(C_5 +
        T C_7\right)\right)}\,}.

    This fit reproduces the Maxwellian integral of :math:`σ(E)` to well
    within a percent over its stated temperature range, so it can be
    compared directly against a numerical integration of the cross
    section as a self-consistency check.

    For the four reactions in Table VII of :cite:t:`bosch:1992` the
    coefficients :math:`C_1, \ldots, C_7` are taken directly from the
    published fit. For the remaining reactions they are obtained the same
    way as the cross-section coefficients (see `fusion_cross_section`). Reference
    reactivities are first generated by numerically integrating the fitted
    cross section over a Maxwellian (the integral above) on a
    temperature grid, and the closed-form fit is then matched to those
    values by nonlinear least squares (`scipy.optimize.curve_fit`),

    .. math::

        \min_{C_1, \ldots, C_7} \sum_k
        \left[⟨σv⟩_{\mathrm{model}}(T_k; C_1, \ldots, C_7)
        - ⟨σv⟩(T_k)\right]^2 ,

    with the reduced-mass energy :math:`m_r c^2` and the :wikipedia:`Gamow
    factor` :math:`B_G` held at their known physical values, and with the same
    seeding and coefficient-scaling safeguards used for the cross-section
    fit. The resulting :math:`C_i` are stored in the same format as the
    published coefficients and evaluated through exactly the same code
    path, so ``fusion_reactivity`` behaves identically for fitted and published
    reactions.

    Examples
    --------
    >>> import astropy.units as u
    >>> fusion_reactivity(10 * u.keV, "D(t,n)A")  # doctest: +SKIP
    <Quantity 1.13616547e-16 cm3 / s>
    """
    if reaction not in rxty_coeffs:
        raise ValueError(
            f"{reaction!r} is not one of the available reactions: "
            f"{', '.join(rxty_coeffs)}"
        )
    if not _in_BH_rxn_ion_temp_range(ion_temp, reaction):
        rxn = _rxty_coeff(reaction)
        raise ValueError(
            f"{ion_temp!r} is not in the {reaction!r} ion temp range "
            f"of {rxn['T_min_keV']} to {rxn['T_max_keV']} keV"
        )
    return _BH_reactivity(ion_temp, reaction)


def _in_BH_rxn_energy_range(E, reaction):
    r"""Return whether ``E`` lies within the cross-section validity range."""
    rxn = _xs_coeff(reaction)
    in_range = (rxn["E_min_keV"] * u.keV <= E) & (rxn["E_max_keV"] * u.keV >= E)
    return bool(np.all(in_range))


def _in_BH_rxn_ion_temp_range(T, reaction):
    r"""Return whether ``T`` lies within the reactivity validity range."""
    rxn = _rxty_coeff(reaction)
    in_range = (rxn["T_min_keV"] * u.keV <= T) & (rxn["T_max_keV"] * u.keV >= T)
    return bool(np.all(in_range))


def _xs_coeff(r):
    r"""Return the cross-section coefficient block for reaction ``r``."""
    return xs_coeffs[r]


def _rxty_coeff(r):
    r"""Return the reactivity coefficient block for reaction ``r``."""
    return rxty_coeffs[r]


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
