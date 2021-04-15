import astropy.units as u
import datetime
import hypothesis
import numpy as np
import pytest

from astropy.tests.helper import assert_quantity_allclose
from hypothesis import example, given, settings
from hypothesis import strategies as st

from plasmapy.particles import IonizationStateCollection
from plasmapy.transport.flows import FlowCalculator

all_species = IonizationStateCollection(
    {"H": [0, 1], "C": [0, 1 / 1.1, 0.1 / 1.1, 0, 0, 0, 0],},
    n0=1e20 * u.m ** -3,
    abundances={"H": 1, "C": 0.11},
    T_e=10 * u.eV,
)
hydrogen = all_species["H"]
carbon_states = all_species["C"]

density_gradient = {
    "H 1+": 1e18 * u.m ** -3 / u.m,
    "C 1+": 1e18 * u.m ** -3 / u.m,
}
temperature_gradient = {
    "H 1+": -10 * u.K / u.m,
    "C 1+": -10 * u.K / u.m,
}


@pytest.fixture(scope="module")
def fc(flux_surface):
    fc = FlowCalculator(
        all_species, flux_surface, density_gradient, temperature_gradient
    )
    return fc


def test_get_flows(fc, num_regression):
    for ion, r in fc.flows.items():
        assert np.isfinite(r).all(), ion
    num_regression.check({key: value.si.value for key, value in fc.flows.items()})


def test_fluxes_partial(fc, num_regression):
    d_partial = {}
    for key, fluxes in [
        ("BP", fc._fluxes_BP),
        ("PS", fc._fluxes_PS),
        ("CL", fc._fluxes_CL),
    ]:
        for ion, (Γ, q) in fluxes.items():
            assert np.isfinite(Γ).all(), ion
            assert np.isfinite(q).all(), ion
            d_partial[f"Γ_{key}_{ion}"] = Γ.si.value
            d_partial[f"q_{key}_{ion}"] = q.si.value
    num_regression.check(d_partial)


@pytest.mark.xfail(reason="units are off in fluxes")
def test_diffusion_coefficient(fc, num_regression):
    d = {}
    for ion, D in fc.diffusion_coefficient.items():
        assert np.isfinite(D).all(), ion
        d[ion] = D.si.value
        assert D.unit.physical_type == "diffusion coefficient"
    num_regression.check(d)


@pytest.mark.xfail(reason="units are off in fluxes")
def test_thermal_coefficient(fc, num_regression):
    d = {}
    for ion, χ in fc.thermal_conductivity.items():
        assert np.isfinite(χ).all(), ion
        d[ion] = χ.si.value
        assert χ.unit.physical_type == "heat conductivity"
    num_regression.check(d)


@pytest.mark.xfail(reason="Units are off; need a tesla in the denominator")
def test_bootstrap_current(fc, num_regression):
    Ib = fc.bootstrap_current
    assert np.isfinite(Ib).all(), ion
    d = {"bootstrap current": Ib.si.value}
    num_regression.check(d)
    assert Ib.unit.physical_type == "electric current density"


@pytest.mark.xfail(reason="units are off")
def test_fluxes(fc, num_regression):
    d = {}
    for ion, (Γ, q) in fc.fluxes.items():
        assert np.isfinite(Γ).all(), ion
        assert np.isfinite(q).all(), ion
        d[f"Γ_{ion}"] = Γ.si.value
        d[f"q_{ion}"] = q.si.value
    num_regression.check(d)
