"""Tests for functionality contained in `plasmapy.formulary.fusion`."""

import json
from functools import cache
from importlib.resources import as_file, files
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary import fusion
from plasmapy.particles import Particle

#: Reactions with Bosch-Hale Padé coefficients (Tables IV and VII).
_BH_REACTIONS = ("D(t,n)A", "3He(d,p)A", "D(d,p)T", "D(d,n)3He")

_CM3_PER_S = u.cm**3 / u.s

#: The published tables are quoted to four significant figures.
_TABLE_RTOL = 1e-3

_DATA_DIR = files("plasmapy.utils.data")


@cache
def _load_table(name):
    with as_file(_DATA_DIR) as physical_path, Path.open(physical_path / name) as f:
        return json.load(f)


_TABLE_V = _load_table("bosch_hale_table_v.json")
_TABLE_VIII = _load_table("bosch_hale_table_viii.json")
_XS_COEFF = fusion._load_reactions("xs_pade_polynomial_coefficients.json")
_RXTY_COEFF = fusion._load_reactions("rxty_pade_polynomial_coefficients.json")

#: All reactions with tabulated ENDF cross-sections.
_XS_REACTIONS = list(_XS_COEFF)
_RXTY_REACTIONS = list(_RXTY_COEFF)
_BH_RXTY_REACTIONS = [r for r in _BH_REACTIONS if r in _RXTY_COEFF]


def _table_params(table, x_unit, coeff):
    """
    Flatten a ``{"columns": [...], "data": [[x, y1, y2, ...]]}`` table into
    ``(x, reaction, expected)`` parameter triples.

    ``None`` entries and reactions absent from ``coeff`` are emitted as
    skipped parameters rather than dropped, so an unverified or unavailable
    cell stays visible in the pytest report instead of quietly vanishing
    from the run — or blowing up with a ``KeyError`` in the raw helper.
    """
    reactions = table["columns"][1:]
    params = []
    for row in table["data"]:
        x, values = row[0], row[1:]
        for reaction, expected in zip(reactions, values, strict=True):
            case = (x * x_unit, reaction, expected)
            if reaction not in coeff:
                mark = pytest.mark.skip(
                    reason=f"{reaction} has no tabulated coefficients"
                )
            elif expected is None:
                mark = pytest.mark.skip(
                    reason="unverified table entry; see known_issues"
                )
            else:
                params.append(pytest.param(*case))
                continue
            params.append(pytest.param(*case, marks=mark))
    return params


_TABLE_V_PARAMS = _table_params(_TABLE_V, u.keV, _XS_COEFF)
_TABLE_VIII_PARAMS = _table_params(_TABLE_VIII, u.keV, _RXTY_COEFF)


class TestBoschHaleCrossSection:
    """Padé parametrization of sigma(E), Bosch & Hale Eq. (8) and Table IV."""

    @pytest.mark.parametrize(("energy", "reaction", "expected_mbarn"), _TABLE_V_PARAMS)
    def test_reproduces_table_v(self, energy, reaction, expected_mbarn):
        """Every cell of Table V is reproduced to its quoted precision."""
        sigma = fusion.fusion_cross_section(energy, reaction)
        assert_quantity_allclose(
            sigma,
            expected_mbarn * u.mbarn,
            rtol=_TABLE_RTOL,
            err_msg=f"Table V mismatch for {reaction} at {energy}",
        )

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_energy_unit_independence(self, reaction):
        """A CM energy of 300 keV is a CM energy of 0.3 MeV."""
        in_keV = fusion.fusion_cross_section(300 * u.keV, reaction)
        in_MeV = fusion.fusion_cross_section(0.3 * u.MeV, reaction)
        assert_quantity_allclose(in_keV, in_MeV, rtol=1e-12)

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_scalar_input_gives_scalar_output(self, reaction):
        sigma = fusion.fusion_cross_section(300 * u.keV, reaction)
        assert sigma.isscalar

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_array_input_preserves_shape(self, reaction):
        energy = np.array([300, 350, 400]) * u.keV
        sigma = fusion.fusion_cross_section(energy, reaction)
        assert sigma.shape == energy.shape

    def test_dt_resonance_position(self):
        """
        The D(t,n)α cross-section peaks at the 3/2+ resonance of 5-He, near a
        CM energy of 64 keV.  This is the sharpest physical feature in the fit
        and catches sign or ordering errors in the Padé coefficients.
        """
        energy = np.linspace(40, 90, 2001) * u.keV
        sigma = fusion.fusion_cross_section(energy, "D(t,n)A")
        peak = energy[np.argmax(sigma)]
        assert_quantity_allclose(peak, 64.7 * u.keV, atol=1.0 * u.keV)

    def test_pade_denominator_is_unity_when_b_coefficients_vanish(self):
        """
        D(d,p)T has B1 = ... = B4 = 0, so S(E) collapses to a bare quartic.
        Comparing against the polynomial written out by hand checks the Horner
        nesting in ``_pade_polynomial``.
        """
        rxn = _XS_COEFF["D(d,p)T"]
        E = np.array([1.0, 10.0, 100.0, 1000.0])
        expected = (
            rxn["A1"]
            + rxn["A2"] * E
            + rxn["A3"] * E**2
            + rxn["A4"] * E**3
            + rxn["A5"] * E**4
        )
        assert np.allclose(fusion._xs_pade_polynomial(rxn, E), expected, rtol=1e-12)

    def test_gamow_factor_suppresses_low_energy(self):
        """Sigma ∝ exp(-B_G/√E)/E, so sigma must fall by orders of magnitude as E → 0."""
        low = fusion.fusion_cross_section(1 * u.keV, "D(t,n)A")
        high = fusion.fusion_cross_section(10 * u.keV, "D(t,n)A")
        assert low < 1e-3 * high


class TestBoschHaleReactivity:
    """Closed-form Maxwellian <sv>(T), Bosch & Hale Eqs. (12)-(14) corrected table VIII values."""

    @pytest.mark.parametrize(
        ("ion_temp", "reaction", "expected_cm3_per_s"), _TABLE_VIII_PARAMS
    )
    def test_reproduces_table_viii(self, ion_temp, reaction, expected_cm3_per_s):
        """Every verified cell of Table VIII is reproduced to its precision."""
        sv = fusion.fusion_reactivity(ion_temp, reaction)
        assert_quantity_allclose(
            sv,
            expected_cm3_per_s * u.cm**3 / u.s,
            rtol=_TABLE_RTOL,
            err_msg=f"Table VIII mismatch for {reaction} at {ion_temp}",
        )

    @pytest.mark.parametrize("reaction", _BH_RXTY_REACTIONS)
    def test_returns_volumetric_rate(self, reaction):
        sv = fusion.fusion_reactivity(60 * u.keV, reaction)
        assert sv.unit.is_equivalent(_CM3_PER_S)

    @pytest.mark.parametrize("reaction", _BH_RXTY_REACTIONS)
    def test_temperature_unit_independence(self, reaction):
        in_keV = fusion.fusion_reactivity(60 * u.keV, reaction)
        in_MeV = fusion.fusion_reactivity(0.06 * u.MeV, reaction)
        assert_quantity_allclose(in_keV, in_MeV, rtol=1e-12)

    @pytest.mark.parametrize("reaction", _BH_RXTY_REACTIONS)
    def test_scalar_input_gives_scalar_output(self, reaction):
        assert fusion.fusion_reactivity(60 * u.keV, reaction).isscalar

    @pytest.mark.parametrize("reaction", _BH_RXTY_REACTIONS)
    def test_array_input_preserves_shape(self, reaction):
        ion_temp = np.array([60, 70, 80]) * u.keV
        assert fusion.fusion_reactivity(ion_temp, reaction).shape == ion_temp.shape

    @pytest.mark.parametrize("reaction", _BH_RXTY_REACTIONS)
    def test_monotonic_below_the_peak(self, reaction):
        """
        <σv> rises monotonically with T until the sigma(E) resonance is passed.
        D(t,n)α turns over near 64 keV, so the check stops at 50 keV.
        """
        ion_temp = np.logspace(np.log10(0.5), np.log10(50), 200) * u.keV
        sv = fusion.fusion_reactivity(ion_temp, reaction).value
        assert np.all(np.diff(sv) > 0), f"{reaction} <sigmav> is not increasing in T"

    def test_dt_dominates_dd_at_ignition_temperatures(self):
        """D-T is the easy reaction: ~200x the D(d,p)T rate at 10 keV."""
        dt = fusion.fusion_reactivity(10 * u.keV, "D(t,n)A")
        dd = fusion.fusion_reactivity(10 * u.keV, "D(d,p)T")
        assert dt / dd > 100

    def test_dd_branches_are_nearly_equal(self):
        """
        The two D-D branches share an entrance channel and differ by only a few
        percent in <σv> at fusion-relevant temperatures.
        """
        ion_temp = np.array([2.0, 10.0, 50.0]) * u.keV
        p_branch = fusion.fusion_reactivity(ion_temp, "D(d,p)T")
        n_branch = fusion.fusion_reactivity(ion_temp, "D(d,n)3He")
        assert_quantity_allclose(p_branch, n_branch, rtol=0.15)


class TestEnergyRangePredicate:
    """The per-reaction energy window is enforced by masking, never by raising."""

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_midpoint_is_in_range(self, reaction):
        c = _XS_COEFF[reaction]
        mid = 0.5 * (c["E_min_keV"] + c["E_max_keV"]) * u.keV
        sigma = fusion.fusion_cross_section(mid, reaction)
        assert np.isfinite(sigma.value)
        assert sigma > 0

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_endpoints_are_inclusive(self, reaction):
        """The check is <= / >=, so both bounds must count as in-range."""
        c = _XS_COEFF[reaction]
        assert fusion.fusion_cross_section(c["E_min_keV"] * u.keV, reaction) > 0
        assert fusion.fusion_cross_section(c["E_max_keV"] * u.keV, reaction) > 0

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_out_of_range_is_nan(self, reaction):
        """Below E_min and above E_max both give NaN rather than an exception."""
        c = _XS_COEFF[reaction]
        for energy in (0.5 * c["E_min_keV"], 2.0 * c["E_max_keV"]):
            with pytest.warns(UserWarning, match="all NaN"):
                sigma = fusion.fusion_cross_section(energy * u.keV, reaction)
            assert np.isnan(sigma.value)

    def test_array_all_in_range_is_true(self):
        c = _XS_COEFF["D(t,n)A"]
        E = np.linspace(c["E_min_keV"], c["E_max_keV"], 5) * u.keV
        sigma = fusion.fusion_cross_section(E, "D(t,n)A")
        assert sigma.shape == E.shape
        assert np.all(sigma.value > 0)

    def test_array_one_bad_element_is_masked_elementwise(self):
        """One out-of-range element is masked without disturbing its neighbours."""
        c = _XS_COEFF["D(t,n)A"]
        E = np.array([c["E_min_keV"], 2.0 * c["E_max_keV"]]) * u.keV
        sigma = fusion.fusion_cross_section(E, "D(t,n)A")
        assert np.isfinite(sigma.value[0])
        assert np.isnan(sigma.value[1])


class TestAvailableReactions:
    """The public listings mirror the loaded coefficient tables."""

    def test_cross_section_listing_matches_coefficients(self):
        assert fusion.available_cross_section_reactions() == list(_XS_COEFF)

    def test_reactivity_listing_matches_coefficients(self):
        assert fusion.available_reactivity_reactions() == list(_RXTY_COEFF)

    @pytest.mark.parametrize(
        "lister",
        [
            fusion.available_cross_section_reactions,
            fusion.available_reactivity_reactions,
        ],
    )
    def test_returns_a_fresh_list(self, lister):
        """Mutating the result must not corrupt the cached coefficient dict."""
        lister().append("Fe(p,gamma)Co")
        assert "Fe(p,gamma)Co" not in lister()

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_every_listed_xs_reaction_evaluates(self, reaction):
        """Anything advertised must actually work in public function."""
        assert reaction in fusion.available_cross_section_reactions()
        c = _XS_COEFF[reaction]
        E = 0.5 * (c["E_min_keV"] + c["E_max_keV"]) * u.keV
        assert fusion.fusion_cross_section(E, reaction) > 0

    @pytest.mark.parametrize("reaction", _RXTY_REACTIONS)
    def test_every_listed_rxty_reaction_evaluates(self, reaction):
        """Anything advertised must actually work in public function."""
        assert reaction in fusion.available_reactivity_reactions()
        r = _RXTY_COEFF[reaction]
        T = 0.5 * (r["T_min_keV"] + r["T_max_keV"]) * u.keV
        assert fusion.fusion_reactivity(T, reaction) > 0


class TestCrossSectionInputValidation:
    """Public-API contract for ``fusion_cross_section``: unknown reactions,
    unit coercion, wrong units, and the SI (m^2) output unit."""

    def test_unknown_reaction_raises(self):
        with pytest.raises(ValueError):
            fusion.fusion_cross_section(100 * u.keV, "Fe(p,gamma)Co")

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_returns_square_meters(self, reaction):
        sigma = fusion.fusion_cross_section(300 * u.keV, reaction)
        assert sigma.unit == u.m**2

    def test_bare_number_warns_and_assumes_keV(self):
        with pytest.warns(u.UnitsWarning):
            sigma = fusion.fusion_cross_section(300, "D(t,n)A")  # in-window value
        assert_quantity_allclose(
            sigma, fusion.fusion_cross_section(300 * u.keV, "D(t,n)A")
        )

    def test_wrong_units_raise(self):
        with pytest.raises(u.UnitsError):
            fusion.fusion_cross_section(100 * u.s, "D(t,n)A")


class TestReactivityInputValidation:
    """Public-API contract for ``fusion_reactivity``: unknown reactions,
    unit coercion, wrong units, and the SI (m^3/s) output unit."""

    def test_unknown_reaction_raises(self):
        with pytest.raises(ValueError):
            fusion.fusion_reactivity(10 * u.keV, "Fe(p,gamma)Co")

    @pytest.mark.parametrize("reaction", _RXTY_REACTIONS)
    def test_output_is_si_volumetric_rate(self, reaction):
        sigma = fusion.fusion_reactivity(60 * u.keV, reaction)
        assert sigma.unit == u.m**3 / u.s

    def test_validity_range_is_per_reaction(self):
        """
        0.3 keV is inside the D(t,n)α window but below the 3He(d,p)α floor, so
        the same temperature evaluates for one reaction and masks for the other.
        """
        assert fusion.fusion_reactivity(0.3 * u.keV, "D(t,n)A").value > 0
        with pytest.warns(UserWarning, match="all NaN"):
            sv = fusion.fusion_reactivity(0.3 * u.keV, "3He(d,p)A")
        assert np.isnan(sv.value)

    def test_bare_number_warns_and_assumes_keV(self):
        with pytest.warns(u.UnitsWarning):
            sigma = fusion.fusion_reactivity(60, "D(t,n)A")  # in-window value
        assert_quantity_allclose(sigma, fusion.fusion_reactivity(60 * u.keV, "D(t,n)A"))

    def test_wrong_units_raise(self):
        with pytest.raises(u.UnitsError):
            fusion.fusion_reactivity(10 * u.s, "D(t,n)A")


class TestCrossSectionOutOfRange:
    """Out-of-range energies are masked to NaN; there is no opt-out."""

    def test_out_of_range_never_raises(self):
        """No policy argument exists: an out-of-range energy is not an error."""
        c = _XS_COEFF["D(t,n)A"]
        E = np.array([c["E_min_keV"], 2.0 * c["E_max_keV"]]) * u.keV
        sigma = fusion.fusion_cross_section(E, "D(t,n)A")
        assert np.isnan(sigma.value[1])

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_masks_out_of_range(self, reaction):
        """Out-of-range slots become NaN; in-range ones are still evaluated."""
        c = _XS_COEFF[reaction]
        lo, hi = c["E_min_keV"], c["E_max_keV"]
        E = np.array([0.5 * lo, 0.5 * (lo + hi), 2.0 * hi]) * u.keV
        sigma = fusion.fusion_cross_section(E, reaction)
        assert sigma.shape == (3,)
        assert sigma.unit == u.m**2
        v = sigma.value
        assert np.isnan(v[0])
        assert np.isnan(v[2])
        assert np.isfinite(v[1])
        assert v[1] > 0

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_endpoints_are_in_range(self, reaction):
        """Bounds are inclusive: E_min and E_max evaluate rather than mask."""
        c = _XS_COEFF[reaction]
        E = np.array([c["E_min_keV"], c["E_max_keV"]]) * u.keV
        sigma = fusion.fusion_cross_section(E, reaction)
        assert np.all(np.isfinite(sigma.value))

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_scalar_out_of_range_returns_nan(self, reaction):
        """A scalar out-of-range value returns a NaN scalar instead of raising."""
        c = _XS_COEFF[reaction]
        with pytest.warns(UserWarning, match="all NaN"):
            sigma = fusion.fusion_cross_section(2.0 * c["E_max_keV"] * u.keV, reaction)
        assert sigma.isscalar
        assert np.isnan(sigma.value)
        assert sigma.unit == u.m**2

    def test_all_out_of_range_warns(self):
        c = _XS_COEFF["D(t,n)A"]
        E = np.array([2.0, 4.0]) * c["E_max_keV"] * u.keV
        with pytest.warns(UserWarning, match="all NaN"):
            sigma = fusion.fusion_cross_section(E, "D(t,n)A")
        assert np.all(np.isnan(sigma.value))

    def test_preserves_2d_shape(self):
        """Masking is elementwise and keeps a non-1D input's shape."""
        c = _XS_COEFF["D(t,n)A"]
        lo, hi = c["E_min_keV"], c["E_max_keV"]
        mid = 0.5 * (lo + hi)
        grid = np.array([[0.5 * lo, mid], [mid, 2.0 * hi]]) * u.keV
        sigma = fusion.fusion_cross_section(grid, "D(t,n)A")
        assert sigma.shape == (2, 2)
        v = sigma.value
        assert np.isnan(v[0, 0])
        assert np.isnan(v[1, 1])
        assert np.isfinite(v[0, 1])
        assert np.isfinite(v[1, 0])


class TestReactivityOutOfRange:
    """Out-of-range ion temperatures are masked to NaN; there is no opt-out."""

    def test_out_of_range_never_raises(self):
        """No policy argument exists: an out-of-range temperature is not an error."""
        c = _RXTY_COEFF["D(t,n)A"]
        T = np.array([c["T_min_keV"], 2.0 * c["T_max_keV"]]) * u.keV
        sv = fusion.fusion_reactivity(T, "D(t,n)A")
        assert np.isnan(sv.value[1])

    @pytest.mark.parametrize("reaction", _RXTY_REACTIONS)
    def test_masks_out_of_range(self, reaction):
        c = _RXTY_COEFF[reaction]
        lo, hi = c["T_min_keV"], c["T_max_keV"]
        T = np.array([0.5 * lo, 0.5 * (lo + hi), 2.0 * hi]) * u.keV
        sv = fusion.fusion_reactivity(T, reaction)
        assert sv.shape == (3,)
        assert sv.unit == u.m**3 / u.s
        v = sv.value
        assert np.isnan(v[0])
        assert np.isnan(v[2])
        assert np.isfinite(v[1])
        assert v[1] > 0

    @pytest.mark.parametrize("reaction", _RXTY_REACTIONS)
    def test_endpoints_are_in_range(self, reaction):
        c = _RXTY_COEFF[reaction]
        T = np.array([c["T_min_keV"], c["T_max_keV"]]) * u.keV
        sv = fusion.fusion_reactivity(T, reaction)
        assert np.all(np.isfinite(sv.value))

    @pytest.mark.parametrize("reaction", _RXTY_REACTIONS)
    def test_scalar_out_of_range_returns_nan(self, reaction):
        c = _RXTY_COEFF[reaction]
        with pytest.warns(UserWarning, match="all NaN"):
            sv = fusion.fusion_reactivity(2.0 * c["T_max_keV"] * u.keV, reaction)
        assert sv.isscalar
        assert np.isnan(sv.value)
        assert sv.unit == u.m**3 / u.s

    def test_all_out_of_range_warns(self):
        c = _RXTY_COEFF["D(t,n)A"]
        T = np.array([2.0, 4.0]) * c["T_max_keV"] * u.keV
        with pytest.warns(UserWarning, match="all NaN"):
            sv = fusion.fusion_reactivity(T, "D(t,n)A")
        assert np.all(np.isnan(sv.value))

    def test_preserves_2d_shape(self):
        c = _RXTY_COEFF["D(t,n)A"]
        lo, hi = c["T_min_keV"], c["T_max_keV"]
        mid = 0.5 * (lo + hi)
        grid = np.array([[0.5 * lo, mid], [mid, 2.0 * hi]]) * u.keV
        sv = fusion.fusion_reactivity(grid, "D(t,n)A")
        assert sv.shape == (2, 2)
        v = sv.value
        assert np.isnan(v[0, 0])
        assert np.isnan(v[1, 1])
        assert np.isfinite(v[0, 1])
        assert np.isfinite(v[1, 0])


def _cm_over_lab(reaction):
    targ, _, rest = reaction.lower().partition("(")
    target = Particle(fusion._NUCLIDES[targ])
    projectile = Particle(fusion._NUCLIDES[rest.partition(",")[0]])
    return (target.mass / (target.mass + projectile.mass)).value


#: Reactions with identical projectile and target, so E_cm = E_lab / 2 exactly.
_SYMMETRIC_REACTIONS = ("D(d,p)T", "D(d,n)3He", "3He(3He,2p)A", "T(t,2n)A")


class TestCrossSectionReferenceFrame:
    """``reference_frame="lab"`` maps a lab projectile energy to the CM energy."""

    def test_cm_is_the_default(self):
        """Omitting reference_frame equals asking for the CM frame."""
        default = fusion.fusion_cross_section(300 * u.keV, "D(t,n)A")
        cm = fusion.fusion_cross_section(300 * u.keV, "D(t,n)A", reference_frame="CM")
        assert_quantity_allclose(default, cm, rtol=1e-12)

    def test_unknown_frame_raises(self):
        with pytest.raises(ValueError, match="reference_frame"):
            fusion.fusion_cross_section(300 * u.keV, "D(t,n)A", reference_frame="beam")

    @pytest.mark.parametrize("reaction", _XS_REACTIONS)
    def test_lab_matches_hand_converted_cm(self, reaction):
        """sigma(E_lab, lab) == sigma(E_lab * m_X/(m_X+m_a), CM)."""
        ratio = _cm_over_lab(reaction)
        c = _XS_COEFF[reaction]
        # a CM energy safely inside the window, and the lab energy mapping to it
        E_cm = c["E_min_keV"] + 0.25 * (c["E_max_keV"] - c["E_min_keV"])
        E_lab = E_cm / ratio
        lab = fusion.fusion_cross_section(
            E_lab * u.keV, reaction, reference_frame="lab"
        )
        cm = fusion.fusion_cross_section(E_lab * ratio * u.keV, reaction)
        assert_quantity_allclose(lab, cm, rtol=1e-10)

    @pytest.mark.parametrize("reaction", _SYMMETRIC_REACTIONS)
    def test_symmetric_reactions_halve_the_energy(self, reaction):
        """Equal masses give E_cm = E_lab / 2 exactly."""
        c = _XS_COEFF[reaction]
        E_cm = 0.5 * (c["E_min_keV"] + c["E_max_keV"]) * u.keV
        lab = fusion.fusion_cross_section(2 * E_cm, reaction, reference_frame="lab")
        cm = fusion.fusion_cross_section(E_cm, reaction)
        assert_quantity_allclose(lab, cm, rtol=1e-10)

    def test_lab_lands_at_lower_cm_energy(self):
        """m_X/(m_X+m_a) < 1, so a lab energy evaluates at a lower CM energy."""
        ratio = _cm_over_lab("D(t,n)A")
        assert 0 < ratio < 1
        E = 120 * u.keV
        lab = fusion.fusion_cross_section(E, "D(t,n)A", reference_frame="lab")
        cm_shifted = fusion.fusion_cross_section(E * ratio, "D(t,n)A")
        assert_quantity_allclose(lab, cm_shifted, rtol=1e-10)

    def test_lab_differs_from_cm_when_masses_differ(self):
        """For unequal masses the frame choice actually changes the result."""
        E = 100 * u.keV
        lab = fusion.fusion_cross_section(E, "3He(d,p)A", reference_frame="lab")
        cm = fusion.fusion_cross_section(E, "3He(d,p)A")
        assert abs(lab - cm) > 1e-6 * cm

    def test_lab_preserves_array_shape(self):
        ratio = _cm_over_lab("D(t,n)A")
        E_cm = np.array([10.0, 60.0, 120.0])
        E_lab = (E_cm / ratio) * u.keV
        lab = fusion.fusion_cross_section(E_lab, "D(t,n)A", reference_frame="lab")
        assert lab.shape == E_lab.shape
        assert_quantity_allclose(
            lab, fusion.fusion_cross_section(E_cm * u.keV, "D(t,n)A"), rtol=1e-10
        )

    def test_unmapped_species_raises_value_error(self, monkeypatch):
        """
        A reaction present in the coefficient file but absent from the nuclide
        map must fail the validator, not fall through to a bare ``KeyError``
        on the dict index.
        """
        patched = {**_XS_COEFF, "Li(y,n)A": _XS_COEFF["D(t,n)A"]}
        monkeypatch.setattr(fusion, "_load_reactions", lambda _name: patched)
        with pytest.raises(ValueError, match=r"No nuclide mapping for 'li'"):
            fusion.fusion_cross_section(300 * u.keV, "Li(y,n)A", reference_frame="lab")
