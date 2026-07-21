"""Cross sections and Maxwellian reactivities for nuclear fusion reactions."""

__all__ = ["fusion_cross_section", "fusion_reactivity"]

# The overall approach for evaluating the fusion cross
# sections was informed in part by this scipython.com notebook:
# https://scipython.com/blog/nuclear-fusion-cross-sections/

import json
import warnings
from functools import cache
from importlib.resources import as_file, files
from pathlib import Path

import astropy.units as u
import numpy as np

from plasmapy.utils.decorators import validate_quantities

"""
Opening Bosch and Hale Tables IV and Json Files
"""

_DATA_DIR = files("plasmapy.utils.data")


@cache
def _load_reactions(name):
    with as_file(_DATA_DIR) as physical_path, Path.open(physical_path / name) as f:
        return json.load(f)["reactions"]


@validate_quantities
def fusion_cross_section(
    energy: u.Quantity[u.keV], reaction: str, out_of_range: str = "raise"
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
        The fusion reaction to evaluate. Must be one of the reaction keys
        listed in the table under Notes (for example, ``"D(t,n)A"``).

    Returns
    -------
    sigma : `~astropy.units.Quantity`
        The fusion cross section, in SI units of square meters

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
    The available reactions, the source of their coefficients, and the
    energy range over which each fit is valid are:

    .. list-table::
       :header-rows: 1
       :widths: 18 34 22 26

       * - Key
         - Reaction
         - Coefficient source
         - Valid energy range
       * - ``"D(t,n)A"``
         - :math:`D + T \rightarrow n + \alpha`
         - Bosch-Hale (Table IV)
         - 0.5 – 550 keV
       * - ``"3He(d,p)A"``
         - :math:`{}^{3}\mathrm{He} + D \rightarrow p + \alpha`
         - Bosch-Hale (Table IV)
         - 0.3 – 900 keV
       * - ``"D(d,p)T"``
         - :math:`D + D \rightarrow p + T`
         - Bosch-Hale (Table IV)
         - 0.5 – 5000 keV
       * - ``"D(d,n)3He"``
         - :math:`D + D \rightarrow n + {}^{3}\mathrm{He}`
         - Bosch-Hale (Table IV)
         - 0.5 – 4900 keV
       * - ``"3He(3He,2p)A"``
         - :math:`{}^{3}\mathrm{He} + {}^{3}\mathrm{He} \rightarrow 2p + \alpha`
         - ENDF/B fit
         - 1.0 – 10000 keV
       * - ``"3He(t,n+p)A"``
         - :math:`{}^{3}\mathrm{He} + T \rightarrow n + p + \alpha`
         - ENDF/B fit
         - 1.0 – 10000 keV
       * - ``"3He(t,d)A"``
         - :math:`{}^{3}\mathrm{He} + T \rightarrow D + \alpha`
         - ENDF/B fit
         - 1.0 – 10000 keV
       * - ``"T(t,2n)A"``
         - :math:`T + T \rightarrow 2n + \alpha`
         - ENDF/B fit
         - 0.5 – 9000 keV
       * - ``"11B(p,a)2A"``
         - :math:`{}^{11}\mathrm{B} + p \rightarrow 3\alpha`
         - ENDF/B fit
         - 200 – 5000 keV

    The four Bosch-Hale reactions use the coefficients published in
    Table IV of :cite:t:`bosch:1992`; the ENDF/B rows use coefficients
    obtained by fitting the identical Padé form to ENDF/B data, as
    described below. A reaction is available only if its coefficients are
    present in the loaded table.

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
    xs_coeff = _load_reactions("xs_pade_polynomial_coefficients.json")
    if reaction not in xs_coeff:
        raise ValueError(
            f"{reaction!r} is not one of the available reactions: {', '.join(xs_coeff)}"
        )
    if out_of_range not in ("raise", "nan"):
        raise ValueError(f"out_of_range must be 'raise' or 'nan', got {out_of_range!r}")
    rxn = xs_coeff[reaction]

    E_keV = np.asarray(energy.to(u.keV).value, dtype=float)
    scalar = E_keV.ndim == 0
    E_arr = np.atleast_1d(E_keV)
    in_range = (E_arr >= rxn["E_min_keV"]) & (E_arr <= rxn["E_max_keV"])

    if out_of_range == "raise" and not in_range.all():
        raise ValueError(
            f"{energy!r} is not in the {reaction!r} energy range "
            f"of {rxn['E_min_keV']} to {rxn['E_max_keV']} keV"
        )

    sigma = np.full(E_arr.shape, np.nan)
    E_in = E_arr[in_range]
    sigma[in_range] = _xs_pade_polynomial(rxn, E_in) / (
        E_in * np.exp(rxn["B_G"] / np.sqrt(E_in))
    )

    if out_of_range == "nan" and E_arr.size and np.all(np.isnan(sigma)):
        warnings.warn(
            f"all input energies are outside the {reaction!r} range "
            f"({rxn['E_min_keV']} to {rxn['E_max_keV']} keV); returning all NaN",
            stacklevel=2,
        )

    result = sigma[0] if scalar else sigma
    return result * u.mbarn


@validate_quantities
def fusion_reactivity(
    ion_temp: u.Quantity[u.keV], reaction: str, out_of_range: str = "raise"
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
        The fusion reaction to evaluate. Must be one of the reaction keys
        listed in the table under Notes (for example, ``"D(t,n)A"``).

    Returns
    -------
    sv : `~astropy.units.Quantity`
        The Maxwellian-averaged reactivity, in SI units of
        m\ :sup:`3`\  s\ :sup:`-1`\ .

    Raises
    ------
    `ValueError`
        If ``reaction`` has no available reactivity coefficients, or if
        ``ion_temp`` falls outside the validity range for ``reaction``.

    `~astropy.units.UnitTypeError`
        If ``ion_temp`` does not have units convertible to keV.

    See Also
    --------
    fusion_cross_section

    Notes
    -----
    The available reactions, the source of their coefficients, and the
    ion-temperature range over which each fit is valid are:

    .. list-table::
       :header-rows: 1
       :widths: 18 34 22 26

       * - Key
         - Reaction
         - Coefficient source
         - Valid ion-temperature range
       * - ``"D(t,n)A"``
         - :math:`D + T \rightarrow n + \alpha`
         - Bosch-Hale (Table VII)
         - 0.2 – 100 keV
       * - ``"3He(d,p)A"``
         - :math:`{}^{3}\mathrm{He} + D \rightarrow p + \alpha`
         - Bosch-Hale (Table VII)
         - 0.5 – 190 keV
       * - ``"D(d,p)T"``
         - :math:`D + D \rightarrow p + T`
         - Bosch-Hale (Table VII)
         - 0.2 – 100 keV
       * - ``"D(d,n)3He"``
         - :math:`D + D \rightarrow n + {}^{3}\mathrm{He}`
         - Bosch-Hale (Table VII)
         - 0.2 – 100 keV
       * - ``"3He(3He,2p)A"``
         - :math:`{}^{3}\mathrm{He} + {}^{3}\mathrm{He} \rightarrow 2p + \alpha`
         - ENDF/B fit
         - 8.27 – 100 keV
       * - ``"3He(t,n+p)A"``
         - :math:`{}^{3}\mathrm{He} + T \rightarrow n + p + \alpha`
         - ENDF/B fit
         - 1.0 – 100 keV
       * - ``"3He(t,d)A"``
         - :math:`{}^{3}\mathrm{He} + T \rightarrow D + \alpha`
         - ENDF/B fit
         - 1.0 – 100 keV
       * - ``"T(t,2n)A"``
         - :math:`T + T \rightarrow 2n + \alpha`
         - ENDF/B fit
         - 1.0 – 100 keV
       * - ``"11B(p,a)2A"``
         - :math:`{}^{11}\mathrm{B} + p \rightarrow 3\alpha`
         - ENDF/B fit
         - 50 – 500 keV

    The four Bosch-Hale reactions use the reactivity coefficients
    published in Table VII of :cite:t:`bosch:1992`; the ENDF/B rows use
    coefficients obtained by fitting the identical closed form to
    reactivities derived from ENDF/B data, as described below. A reaction
    is available only if its coefficients are present in the loaded table.

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
    rxty_coeff = _load_reactions("rxty_pade_polynomial_coefficients.json")
    if reaction not in rxty_coeff:
        raise ValueError(
            f"{reaction!r} is not one of the available reactions: "
            f"{', '.join(rxty_coeff)}"
        )
    if out_of_range not in ("raise", "nan"):
        raise ValueError(f"out_of_range must be 'raise' or 'nan', got {out_of_range!r}")
    rxn = rxty_coeff[reaction]

    T_keV = np.asarray(ion_temp.to(u.keV).value, dtype=float)
    scalar = T_keV.ndim == 0
    T_arr = np.atleast_1d(T_keV)
    in_range = (T_arr >= rxn["T_min_keV"]) & (T_arr <= rxn["T_max_keV"])

    if out_of_range == "raise" and not in_range.all():
        raise ValueError(
            f"{ion_temp!r} is not in the {reaction!r} ion temp range "
            f"of {rxn['T_min_keV']} to {rxn['T_max_keV']} keV"
        )

    sv = np.full(T_arr.shape, np.nan)
    T_in = T_arr[in_range]
    Theta = _rxty_pade_polynomial(rxn, T_in)
    xi = ((rxn["B_G"] ** 2) / (4 * Theta)) ** (1 / 3)
    sv[in_range] = (
        np.sqrt(xi / (rxn["m_r_c2"] * T_in**3)) * rxn["C1"] * Theta * np.exp(-3 * xi)
    )

    if out_of_range == "nan" and T_arr.size and np.all(np.isnan(sv)):
        warnings.warn(
            f"all input temperatures are outside the {reaction!r} range "
            f"({rxn['T_min_keV']} to {rxn['T_max_keV']} keV); returning all NaN",
            stacklevel=2,
        )

    result = sv[0] if scalar else sv
    return result * (u.cm**3 / u.s)


def _xs_pade_polynomial(rxn, E):
    r"""
    Evaluate the Bosch-Hale Padé approximant for the S-function.
    """
    S_vals = rxn["A1"] + E * (
        rxn["A2"] + E * (rxn["A3"] + E * (rxn["A4"] + E * rxn["A5"]))
    )
    S_vals /= 1 + E * (rxn["B1"] + E * (rxn["B2"] + E * (rxn["B3"] + E * rxn["B4"])))
    return S_vals


def _rxty_pade_polynomial(rxn, T):
    r"""
    Evaluate the :math:`\theta(T)` Padé approximant for the Bosch-Hale
    reactivity.
    """
    theta = T * (rxn["C2"] + T * (rxn["C4"] + T * rxn["C6"]))
    theta /= 1 + T * (rxn["C3"] + T * (rxn["C5"] + T * rxn["C7"]))
    theta = 1 - theta
    return T / theta
