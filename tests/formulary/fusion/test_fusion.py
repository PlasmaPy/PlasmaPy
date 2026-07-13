"""Tests for `plasmapy.formulary.fusion.fusion`."""

import astropy.units as u
import numpy as np
import pytest
from astropy.tests.helper import assert_quantity_allclose

from plasmapy.formulary.fusion.fusion import (
    fusion_cross_section,
    fusion_power_density,
    fusion_reaction_rate,
    reactivity,
)
from plasmapy.formulary.fusion.parameters import (
    DDn,
    DDp,
    DHe3,
    DT,
    _lookup_reaction,
)

# Bosch-Hale reference values computed from the parameterization
# Cross-section reference values in barns (1 barn = 1e-28 m^2)
DT_XSEC_10_KEV = 2.702072e-2  # barns at 10 keV
DT_XSEC_50_KEV = 4.218617  # barns at 50 keV
DT_XSEC_100_KEV = 3.427245  # barns at 100 keV
DDP_XSEC_50_KEV = 1.557252e-2  # barns at 50 keV
DDN_XSEC_50_KEV = 1.648744e-2  # barns at 50 keV
DHE3_XSEC_50_KEV = 8.687707e-3  # barns at 50 keV

# Reactivity reference values in m^3/s
DT_REACT_10_KEV = 5.249405e-26
DT_REACT_20_KEV = 2.000675e-25
DT_REACT_50_KEV = 3.996121e-25
DT_REACT_70_KEV = 4.129763e-25
DT_REACT_100_KEV = 3.903058e-25
DDP_REACT_10_KEV = 2.984959e-28
DDN_REACT_50_KEV = 5.849748e-27
DHE3_REACT_50_KEV = 2.566188e-26

# NRL Formulary Table 1 reference: D-T cross-section at 10 keV is ~0.03 barns
# (NRL Formulary 2019, page 45, Table 1)
# Bosch-Hale gives 0.027 barns at 10 keV


def test_dt_cross_section_10_keV() -> None:
    r"""D-T cross-section at 10 keV should be ~0.027 barns."""
    sigma = fusion_cross_section("D + T --> He-4 + n", 10 * u.keV)
    assert sigma.unit.is_equivalent(u.barn)
    assert np.isclose(sigma.to(u.barn).value, DT_XSEC_10_KEV, rtol=1e-3)


def test_dt_cross_section_50_keV() -> None:
    r"""D-T cross-section at 50 keV should be ~4.2 barns."""
    sigma = fusion_cross_section("D + T --> He-4 + n", 50 * u.keV)
    assert np.isclose(sigma.to(u.barn).value, DT_XSEC_50_KEV, rtol=1e-3)


def test_dt_cross_section_100_keV() -> None:
    r"""D-T cross-section at 100 keV should be ~3.4 barns."""
    sigma = fusion_cross_section("D + T --> He-4 + n", 100 * u.keV)
    assert np.isclose(sigma.to(u.barn).value, DT_XSEC_100_KEV, rtol=1e-3)


def test_dt_reactivity_10_keV() -> None:
    r"""D-T reactivity at 10 keV."""
    sv = reactivity("D + T --> He-4 + n", 10 * u.keV)
    assert sv.unit.is_equivalent(u.m**3 / u.s)
    assert np.isclose(sv.value, DT_REACT_10_KEV, rtol=1e-3)


def test_dt_reactivity_20_keV() -> None:
    r"""D-T reactivity at 20 keV should be ~2e-25 m^3/s."""
    sv = reactivity("D + T --> He-4 + n", 20 * u.keV)
    assert np.isclose(sv.value, DT_REACT_20_KEV, rtol=1e-3)


def test_dt_reactivity_peaks_around_70_keV() -> None:
    r"""D-T reactivity should peak at approximately 70 keV."""
    sv = reactivity("D + T --> He-4 + n", 70 * u.keV)
    assert np.isclose(sv.value, DT_REACT_70_KEV, rtol=1e-3)


def test_dt_reactivity_100_keV() -> None:
    r"""D-T reactivity at 100 keV should be ~3.9e-25 m^3/s."""
    sv = reactivity("D + T --> He-4 + n", 100 * u.keV)
    assert np.isclose(sv.value, DT_REACT_100_KEV, rtol=1e-3)


def test_dd_cross_section_branches() -> None:
    r"""D-D cross-sections for both branches at 50 keV."""
    sigma_ddp = fusion_cross_section("D(d,p)T", 50 * u.keV)
    sigma_ddn = fusion_cross_section("D(d,n)3He", 50 * u.keV)
    assert sigma_ddp.unit.is_equivalent(u.barn)
    assert sigma_ddn.unit.is_equivalent(u.barn)
    assert np.isclose(sigma_ddp.to(u.barn).value, DDP_XSEC_50_KEV, rtol=1e-3)
    assert np.isclose(sigma_ddn.to(u.barn).value, DDN_XSEC_50_KEV, rtol=1e-3)


def test_dd_reactivity_string() -> None:
    r"""D-D reactivity using full reaction string."""
    sv = reactivity("D(d,p)T", 50 * u.keV)
    assert sv.unit.is_equivalent(u.m**3 / u.s)


def test_dd_reactivity_invalid_tuple() -> None:
    r"""D-D with just reactants should raise an error."""
    with pytest.raises(ValueError, match="two reaction branches"):
        reactivity(("D", "D"), 50 * u.keV)


def test_dhe3_cross_section() -> None:
    r"""D-He3 cross-section at 50 keV."""
    sigma = fusion_cross_section("D + He-3 --> He-4 + p", 50 * u.keV)
    assert sigma.unit.is_equivalent(u.barn)
    assert np.isclose(sigma.to(u.barn).value, DHE3_XSEC_50_KEV, rtol=1e-3)


def test_dhe3_reactivity() -> None:
    r"""D-He3 reactivity at 50 keV."""
    sv = reactivity("D + He-3 --> He-4 + p", 50 * u.keV)
    assert sv.unit.is_equivalent(u.m**3 / u.s)
    assert np.isclose(sv.value, DHE3_REACT_50_KEV, rtol=1e-3)


def test_reactivity_short_name() -> None:
    r"""Reactivity with short reaction names."""
    sv1 = reactivity("T(d,n)4He", 20 * u.keV)
    sv2 = reactivity("D + T --> He-4 + n", 20 * u.keV)
    assert_quantity_allclose(sv1, sv2)


def test_alpha_alias() -> None:
    r"""Reaction with 'alpha' naming."""
    sv1 = reactivity("D + T --> alpha + n", 20 * u.keV)
    sv2 = reactivity("D + T --> He-4 + n", 20 * u.keV)
    assert_quantity_allclose(sv1, sv2)


def test_reactants_tuple_dt() -> None:
    r"""Tuple of reactants should work for unambiguous reactions."""
    sv = reactivity(("D", "T"), 20 * u.keV)
    assert np.isfinite(sv.value)


def test_reaction_rate_dt() -> None:
    r"""D-T reaction rate at typical fusion conditions."""
    n = 1e20 * u.m**-3
    rate = fusion_reaction_rate("D + T --> He-4 + n", 20 * u.keV, n, n)
    assert rate.unit.is_equivalent(u.m**-3 / u.s)
    assert np.isfinite(rate.value)
    assert rate.value > 0


def test_reaction_rate_dd_identical_particles() -> None:
    r"""D-D rate should include the 1/2 factor for identical particles."""
    n = 1e20 * u.m**-3
    rate = fusion_reaction_rate("D(d,p)T", 50 * u.keV, n, n)
    # Compare with non-identical case: should be roughly half of n^2*sigma_v
    sv = reactivity("D(d,p)T", 50 * u.keV)
    expected = 0.5 * n * n * sv
    assert_quantity_allclose(rate, expected, rtol=1e-6)


def test_power_density_dt() -> None:
    r"""D-T power density."""
    n = 1e20 * u.m**-3
    P = fusion_power_density("D + T --> He-4 + n", 20 * u.keV, n, n)
    assert P.unit.is_equivalent(u.W / u.m**3)
    assert np.isfinite(P.value)
    assert P.value > 0


def test_cross_section_out_of_range() -> None:
    r"""Energy outside the valid range should raise."""
    with pytest.raises(ValueError, match="outside the valid range"):
        fusion_cross_section("D + T --> He-4 + n", 1e6 * u.keV)


def test_reactivity_out_of_range() -> None:
    r"""Temperature outside the valid range should raise."""
    with pytest.raises(ValueError, match="outside the valid range"):
        reactivity("D + T --> He-4 + n", 1e6 * u.keV)


def test_invalid_reaction() -> None:
    r"""Invalid reaction string should raise."""
    with pytest.raises(ValueError, match="Unknown reaction"):
        fusion_cross_section("C + O --> Si", 50 * u.keV)


def test_pb11_reactivity() -> None:
    r"""p-B11 reactivity should be finite and positive."""
    sv = reactivity("p + B-11 --> He-4 + He-4 + He-4", 300 * u.keV)
    assert sv.unit.is_equivalent(u.m**3 / u.s)
    assert np.isfinite(sv.value)
    assert sv.value > 0


def test_pb11_reactivity_short_name() -> None:
    r"""p-B11 reactivity via short name."""
    sv = reactivity("11B(p,4He)4He4He", 300 * u.keV)
    assert np.isfinite(sv.value)


def test_3_alpha_alias() -> None:
    r"""p-B11 reaction via '3 alpha' alias."""
    sv = reactivity("p + B-11 --> 3 alpha", 300 * u.keV)
    assert np.isfinite(sv.value)


def test_cross_section_unit() -> None:
    r"""Cross-section unit should be m^2."""
    sigma = fusion_cross_section("D + T --> He-4 + n", 50 * u.keV)
    assert sigma.unit == u.m**2


def test_reactivity_unit() -> None:
    r"""Reactivity unit should be m^3/s."""
    sv = reactivity("D + T --> He-4 + n", 20 * u.keV)
    assert sv.unit == u.m**3 / u.s


def test_reaction_rate_unit() -> None:
    r"""Reaction rate unit should be m^{-3} s^{-1}."""
    n = 1e20 * u.m**-3
    rate = fusion_reaction_rate("D + T --> He-4 + n", 20 * u.keV, n, n)
    assert rate.unit == u.m**-3 / u.s


def test_power_density_unit() -> None:
    r"""Power density unit should be W/m^3."""
    n = 1e20 * u.m**-3
    P = fusion_power_density("D + T --> He-4 + n", 20 * u.keV, n, n)
    assert P.unit.is_equivalent(u.W / u.m**3)


def test_lookup_dt() -> None:
    r"""Lookup D-T reaction via short name."""
    params = _lookup_reaction("T(d,n)4He")
    assert params.name == "T(d,n)4He"
    assert params.reactant1 == "D"
    assert params.reactant2 == "T"


def test_lookup_ddp() -> None:
    r"""Lookup D-D proton branch."""
    params = _lookup_reaction("D(d,p)T")
    assert params.name == "D(d,p)T"


def test_lookup_ddn() -> None:
    r"""Lookup D-D neutron branch."""
    params = _lookup_reaction("D(d,n)3He")
    assert params.name == "D(d,n)3He"


def test_lookup_reaction_string_dt() -> None:
    r"""Lookup via reaction string."""
    params = _lookup_reaction("D + T --> He-4 + n")
    assert params.name == "T(d,n)4He"


def test_lookup_reaction_string_dd_ambiguous() -> None:
    r"""Lookup D-D via reaction string."""
    params = _lookup_reaction("D + D --> T + p")
    assert params.name == "D(d,p)T"


def test_lookup_reactants_tuple() -> None:
    r"""Lookup via reactant tuple for unambiguous reactions."""
    params = _lookup_reaction(("D", "T"))
    assert params.name == "T(d,n)4He"
    params = _lookup_reaction(("He-3", "D"))
    assert params.name == "3He(d,p)4He"


def test_lookup_reactants_dd_raises() -> None:
    r"""Lookup D-D via tuple should raise."""
    with pytest.raises(ValueError, match="two reaction branches"):
        _lookup_reaction(("D", "D"))
