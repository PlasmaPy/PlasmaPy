import json
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary import fusion

#: Reactions with Bosch-Hale Padé coefficients (Tables IV and VII).
BH_REACTIONS = ("D(t,n)A", "3He(d,p)A", "D(d,p)T", "D(d,n)3He")

#: All reactions with tabulated ENDF cross-sections.
ENDF_REACTIONS = (
    "D(t,n)A",
    "3He(d,p)A",
    "D(d,p)T",
    "D(d,n)3He",
    "3He(3He,2p)A",
    "3He(t,n+p)A",
    "3He(t,d)A",
    "T(t,2n)A",
    "11B(p,a)2A",
)

CM3_PER_S = u.cm**3 / u.s

#: The published tables are quoted to four significant figures.
TABLE_RTOL = 1e-3


def _load_table(name):
    with Path.open(fusion.DATA_DIR / name) as f:
        return json.load(f)


TABLE_V = _load_table("bosch_hale_table_v.json")
TABLE_VIII = _load_table("bosch_hale_table_viii.json")


def _table_params(table, x_unit):
    """
    Flatten a ``{"columns": [...], "data": [[x, y1, y2, ...]]}`` table into
    ``(x, reaction, expected)`` parameter triples.

    ``None`` entries are emitted as skipped parameters rather than dropped, so
    that an unverified cell stays visible in the pytest report instead of
    quietly vanishing from the run.
    """
    reactions = table["columns"][1:]
    params = []
    for row in table["data"]:
        x, values = row[0], row[1:]
        for reaction, expected in zip(reactions, values, strict=True):
            case = (x * x_unit, reaction, expected)
            if expected is None:
                mark = pytest.mark.skip(
                    reason="unverified table entry; see known_issues"
                )
                params.append(pytest.param(*case, marks=mark))
            else:
                params.append(pytest.param(*case))
    return params


TABLE_V_PARAMS = _table_params(TABLE_V, u.keV)
TABLE_VIII_PARAMS = _table_params(TABLE_VIII, u.keV)


class TestBoschHaleCrossSection:
    """Padé parametrization of sigma(E), Bosch & Hale Eq. (8) and Table IV."""

    @pytest.mark.parametrize(("energy", "reaction", "expected_mbarn"), TABLE_V_PARAMS)
    def test_reproduces_table_v(self, energy, reaction, expected_mbarn):
        """Every cell of Table V is reproduced to its quoted precision."""
        sigma = fusion._BH_cross_section(energy, reaction)
        assert_quantity_allclose(
            sigma,
            expected_mbarn * u.mbarn,
            rtol=TABLE_RTOL,
            err_msg=f"Table V mismatch for {reaction} at {energy}",
        )

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_returns_millibarns(self, reaction):
        sigma = fusion._BH_cross_section(300 * u.keV, reaction)
        assert sigma.unit.physical_type == "area"

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_energy_unit_independence(self, reaction):
        """A CM energy of 300 keV is a CM energy of 0.3 MeV."""
        in_keV = fusion._BH_cross_section(300 * u.keV, reaction)
        in_MeV = fusion._BH_cross_section(0.3 * u.MeV, reaction)
        assert_quantity_allclose(in_keV, in_MeV, rtol=1e-12)

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_scalar_input_gives_scalar_output(self, reaction):
        sigma = fusion._BH_cross_section(300 * u.keV, reaction)
        assert sigma.isscalar

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_array_input_preserves_shape(self, reaction):
        energy = np.array([300, 350, 400]) * u.keV
        sigma = fusion._BH_cross_section(energy, reaction)
        assert sigma.shape == energy.shape

    def test_dt_resonance_position(self):
        """
        The D(t,n)α cross-section peaks at the 3/2+ resonance of 5-He, near a
        CM energy of 64 keV.  This is the sharpest physical feature in the fit
        and catches sign or ordering errors in the Padé coefficients.
        """
        energy = np.linspace(40, 90, 2001) * u.keV
        sigma = fusion._BH_cross_section(energy, "D(t,n)A")
        peak = energy[np.argmax(sigma)]
        assert_quantity_allclose(peak, 64.7 * u.keV, atol=1.0 * u.keV)

    def test_pade_denominator_is_unity_when_b_coefficients_vanish(self):
        """
        D(d,p)T has B1 = ... = B4 = 0, so S(E) collapses to a bare quartic.
        Comparing against the polynomial written out by hand checks the Horner
        nesting in ``_pade_polynomial``.
        """
        rxn = fusion._xs_coeff("D(d,p)T")
        E = np.array([1.0, 10.0, 100.0, 1000.0])
        expected = (
            rxn["A1"]
            + rxn["A2"] * E
            + rxn["A3"] * E**2
            + rxn["A4"] * E**3
            + rxn["A5"] * E**4
        )
        assert np.allclose(fusion._pade_polynomial(rxn, E), expected, rtol=1e-12)

    def test_gamow_factor_suppresses_low_energy(self):
        """Sigma ∝ exp(-B_G/√E)/E, so sigma must fall by orders of magnitude as E → 0."""
        low = fusion._BH_cross_section(1 * u.keV, "D(t,n)A")
        high = fusion._BH_cross_section(10 * u.keV, "D(t,n)A")
        assert low < 1e-3 * high


class TestBoschHaleReactivity:
    """Closed-form Maxwellian <sv>(T), Bosch & Hale Eqs. (12)-(14)."""

    @pytest.mark.parametrize(
        ("ion_temp", "reaction", "expected_cm3_per_s"), TABLE_VIII_PARAMS
    )
    def test_reproduces_table_viii(self, ion_temp, reaction, expected_cm3_per_s):
        """Every verified cell of Table VIII is reproduced to its precision."""
        sv = fusion._BH_reactivity(ion_temp, reaction)
        assert_quantity_allclose(
            sv,
            expected_cm3_per_s * u.cm**3 / u.s,
            rtol=TABLE_RTOL,
            err_msg=f"Table VIII mismatch for {reaction} at {ion_temp}",
        )

    @pytest.mark.parametrize("reaction", BH_REACTIONS)
    def test_returns_volumetric_rate(self, reaction):
        sv = fusion._BH_reactivity(60 * u.keV, reaction)
        assert sv.unit.is_equivalent(CM3_PER_S)

    @pytest.mark.parametrize("reaction", BH_REACTIONS)
    def test_temperature_unit_independence(self, reaction):
        in_keV = fusion._BH_reactivity(60 * u.keV, reaction)
        in_MeV = fusion._BH_reactivity(0.06 * u.MeV, reaction)
        assert_quantity_allclose(in_keV, in_MeV, rtol=1e-12)

    @pytest.mark.parametrize("reaction", BH_REACTIONS)
    def test_scalar_input_gives_scalar_output(self, reaction):
        assert fusion._BH_reactivity(60 * u.keV, reaction).isscalar

    @pytest.mark.parametrize("reaction", BH_REACTIONS)
    def test_array_input_preserves_shape(self, reaction):
        ion_temp = np.array([60, 70, 80]) * u.keV
        assert fusion._BH_reactivity(ion_temp, reaction).shape == ion_temp.shape

    @pytest.mark.parametrize("reaction", BH_REACTIONS)
    def test_monotonic_below_the_peak(self, reaction):
        """
        <σv> rises monotonically with T until the sigma(E) resonance is passed.
        D(t,n)α turns over near 64 keV, so the check stops at 50 keV.
        """
        ion_temp = np.logspace(np.log10(0.5), np.log10(50), 200) * u.keV
        sv = fusion._BH_reactivity(ion_temp, reaction).value
        assert np.all(np.diff(sv) > 0), f"{reaction} <sigmav> is not increasing in T"

    def test_dt_dominates_dd_at_ignition_temperatures(self):
        """D-T is the easy reaction: ~200x the D(d,p)T rate at 10 keV."""
        dt = fusion._BH_reactivity(10 * u.keV, "D(t,n)A")
        dd = fusion._BH_reactivity(10 * u.keV, "D(d,p)T")
        assert dt / dd > 100

    def test_dd_branches_are_nearly_equal(self):
        """
        The two D-D branches share an entrance channel and differ by only a few
        percent in <σv> at fusion-relevant temperatures.
        """
        ion_temp = np.array([2.0, 10.0, 50.0]) * u.keV
        p_branch = fusion._BH_reactivity(ion_temp, "D(d,p)T")
        n_branch = fusion._BH_reactivity(ion_temp, "D(d,n)3He")
        assert_quantity_allclose(p_branch, n_branch, rtol=0.15)


class TestEnergyRangePredicate:
    """``_in_BH_rxn_energy_range`` returns a bool; it never raises."""

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_midpoint_is_in_range(self, reaction):
        c = fusion._xs_coeff(reaction)
        mid = 0.5 * (c["E_min_keV"] + c["E_max_keV"]) * u.keV
        assert fusion._in_BH_rxn_energy_range(mid, reaction) is True

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_endpoints_are_inclusive(self, reaction):
        """The check is <= / >=, so both bounds must count as in-range."""
        c = fusion._xs_coeff(reaction)
        assert fusion._in_BH_rxn_energy_range(c["E_min_keV"] * u.keV, reaction) is True
        assert fusion._in_BH_rxn_energy_range(c["E_max_keV"] * u.keV, reaction) is True

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_out_of_range_is_false(self, reaction):
        c = fusion._xs_coeff(reaction)
        assert (
            fusion._in_BH_rxn_energy_range(0.5 * c["E_min_keV"] * u.keV, reaction)
            is False
        )
        assert (
            fusion._in_BH_rxn_energy_range(2.0 * c["E_max_keV"] * u.keV, reaction)
            is False
        )

    def test_array_all_in_range_is_true(self):
        c = fusion._xs_coeff("D(t,n)A")
        E = np.linspace(c["E_min_keV"], c["E_max_keV"], 5) * u.keV
        assert fusion._in_BH_rxn_energy_range(E, "D(t,n)A") is True

    def test_array_one_bad_element_is_false(self):
        """All-or-nothing: one out-of-range element flips the whole array False."""
        c = fusion._xs_coeff("D(t,n)A")
        E = np.array([c["E_min_keV"], 2.0 * c["E_max_keV"]]) * u.keV
        assert fusion._in_BH_rxn_energy_range(E, "D(t,n)A") is False


class TestCrossSectionDispatch:
    """``cross_section`` routes to a backend, validates inputs, or raises."""

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_matches_backend(self, reaction):
        assert_quantity_allclose(
            fusion.cross_section(300 * u.keV, reaction),
            fusion._BH_cross_section(300 * u.keV, reaction),
        )

    def test_unknown_reaction_raises(self):
        with pytest.raises(ValueError):
            fusion.cross_section(100 * u.keV, "Fe(p,gamma)Co")

    @pytest.mark.parametrize(
        ("reaction", "energy"),
        [
            ("D(t,n)A", 0.1 * u.keV),  # below E_min = 0.5 keV
            ("D(t,n)A", 600 * u.keV),  # above E_max = 550 keV
            ("3He(d,p)A", 0.2 * u.keV),  # below E_min = 0.3 keV
            ("3He(d,p)A", 1000 * u.keV),  # above E_max = 900 keV
            ("D(d,p)T", 6000 * u.keV),  # above E_max = 5000 keV
            ("D(d,n)3He", 5000 * u.keV),  # above E_max = 4900 keV
            ("3He(3He,2p)A", 11000 * u.keV),  # above E_max = 10000 keV
            ("3He(3He,2p)A", 0.1 * u.keV),  # below E_min = 1 keV
            ("3He(t,n+p)A", 11000 * u.keV),  # above E_max = 10000 keV
            ("3He(t,n+p)A", 0.1 * u.keV),  # below E_min = 1 keV
            ("3He(t,d)A", 11000 * u.keV),  # above E_max = 10000 keV
            ("3He(t,d)A", 0.1 * u.keV),  # below E_min = 0.5 keV
            ("T(t,2n)A", 10000 * u.keV),  # above E_max = 9000 keV
            ("T(t,2n)A", 0.1 * u.keV),  # below E_min = 0.5 keV
            ("11B(p,a)2A", 6000 * u.keV),  # above E_max = 5000 keV
            ("11B(p,a)2A", 100 * u.keV),  # below E_min = 200 keV
        ],
    )
    def test_energy_outside_validity_range_raises(self, reaction, energy):
        with pytest.raises(ValueError, match="energy range"):
            fusion.cross_section(energy, reaction)

    def test_partially_out_of_range_array_raises(self):
        """The range check is all or nothing: one bad element rejects the array."""
        energy = np.array([10.0, 100.0, 10_000.0]) * u.keV
        with pytest.raises(ValueError, match="energy range"):
            fusion.cross_section(energy, "D(t,n)A")

    def test_bare_number_warns_and_assumes_keV(self):
        with pytest.warns(u.UnitsWarning):
            sigma = fusion.cross_section(300, "D(t,n)A")  # in-window value
        assert_quantity_allclose(sigma, fusion.cross_section(300 * u.keV, "D(t,n)A"))

    def test_wrong_units_raise(self):
        with pytest.raises(u.UnitsError):
            fusion.cross_section(100 * u.s, "D(t,n)A")

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_public_cross_section_is_finite_and_positive(self, reaction):
        sigma = fusion.cross_section(300 * u.keV, reaction)
        assert np.isfinite(sigma.value)
        assert sigma.value > 0


class TestReactivityDispatch:
    """``reactivity`` routes to a backend, validates inputs, or raises."""

    @pytest.mark.parametrize("reaction", ENDF_REACTIONS)
    def test_bh_source_matches_backend(self, reaction):
        assert_quantity_allclose(
            fusion.reactivity(60 * u.keV, reaction),
            fusion._BH_reactivity(60 * u.keV, reaction),
        )

    def test_unknown_reaction_raises(self):
        with pytest.raises(ValueError):
            fusion.reactivity(10 * u.keV, "Fe(p,gamma)Co")

    @pytest.mark.parametrize(
        ("reaction", "ion_temp"),
        [
            ("D(t,n)A", 0.1 * u.keV),  # below T_min = 0.2 keV
            ("D(t,n)A", 150 * u.keV),  # above T_max = 100 keV
            ("3He(d,p)A", 0.3 * u.keV),  # below T_min = 0.5 keV
            ("3He(d,p)A", 200 * u.keV),  # above T_max = 190 keV
            ("D(d,p)T", 101 * u.keV),  # above T_max = 100 keV
            ("D(d,n)3He", 0.1 * u.keV),  # below T_min = 0.2 keV
            ("3He(3He,2p)A", 110 * u.keV),  # above T_max = 100 keV
            ("3He(3He,2p)A", 0.1 * u.keV),  # below T_min = 8.27 keV
            ("3He(t,n+p)A", 110 * u.keV),  # above T_max = 100 keV
            ("3He(t,n+p)A", 0.1 * u.keV),  # below T_min = 1 keV
            ("3He(t,d)A", 110 * u.keV),  # above T_max = 100 keV
            ("3He(t,d)A", 0.1 * u.keV),  # below T_min = 1 keV
            ("T(t,2n)A", 110 * u.keV),  # above T_max = 100 keV
            ("T(t,2n)A", 0.1 * u.keV),  # below T_min = 1 keV
            ("11B(p,a)2A", 510 * u.keV),  # above T_max = 500 keV
            ("11B(p,a)2A", 40 * u.keV),  # below T_min = 50 keV
        ],
    )
    def test_temperature_outside_bh_validity_range_raises(self, reaction, ion_temp):
        with pytest.raises(ValueError, match="energy range"):
            fusion.reactivity(ion_temp, reaction)

    def test_validity_range_is_per_reaction(self):
        """
        0.3 keV is inside the D(t,n)α window but below the 3He(d,p)α floor.
        This is why the Table VIII regression tests call ``_BH_reactivity``.
        """
        assert fusion.reactivity(0.3 * u.keV, "D(t,n)A").value > 0
        with pytest.raises(ValueError):
            fusion.reactivity(0.3 * u.keV, "3He(d,p)A")

    def test_bare_number_warns_and_assumes_keV(self):
        with pytest.warns(u.UnitsWarning):
            sigma = fusion.reactivity(60, "D(t,n)A")  # in-window value
        assert_quantity_allclose(sigma, fusion.reactivity(60 * u.keV, "D(t,n)A"))

    def test_wrong_units_raise(self):
        with pytest.raises(u.UnitsError):
            fusion.reactivity(10 * u.s, "D(t,n)A")
