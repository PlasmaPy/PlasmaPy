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


def test_fluxes(fc, num_regression):
    d = {}
    for ion, (Γ, q) in fc.fluxes.items():
        assert np.isfinite(Γ).all(), ion
        assert np.isfinite(q).all(), ion
        d[f"Γ_{ion}"] = Γ.si.value
        d[f"q_{ion}"] = q.si.value
    num_regression.check(d)
