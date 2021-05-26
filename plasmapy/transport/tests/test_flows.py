import astropy
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

# if "profile" not in globals():
#     def profile(func):
#         return func


all_species = IonizationStateCollection(
    {
        "H": [0, 1],
        "C": [0, 1 / 1.1, 0.1 / 1.1, 0, 0, 0, 0],
    },
    n0=1e20 * u.m ** -3,
    abundances={"H": 1, "C": 0.11},
    T_e=10 * u.eV,
)
hydrogen = all_species["H"]
carbon_states = all_species["C"]

density_gradient = {
    "H": np.ones(2) * 1e18 * u.m ** -3 / u.m,
    "C": np.ones(7) * 1e18 * u.m ** -3 / u.m,
}
temperature_gradient = {
    "H": np.ones(2) * -10 * u.K / u.m,
    "C": np.ones(7) * -10 * u.K / u.m,
}


# profile
@pytest.fixture(scope="module")
def fc(flux_surface):
    fc = FlowCalculator(
        all_species,
        flux_surface,
        density_gradient,
        temperature_gradient,
        mu_N=1000,
    )
    return fc


# profile
def test_get_flows(fc, num_regression):
    for ion, r in fc._charge_state_flows.items():
        if "0" in ion:
            continue
        assert np.isfinite(r).all(), ion
    num_regression.check(
        {key: value.si.value for key, value in fc._charge_state_flows.items()}
    )


@pytest.mark.parametrize(
    "key",
    [
        "BP",
        "CL",
        "PS",
    ],
)
def test_fluxes_partial(fc, key, num_regression):
    fluxes = getattr(fc, f"_fluxes_{key}")
    d_partial = {}
    for ion, (Γ, q) in fluxes.items():
        if "0" in ion:
            continue
        assert np.isfinite(Γ).all(), ion
        assert np.isfinite(q).all(), ion
        d_partial[f"Γ_{key}_{ion}"] = Γ.si.value
        d_partial[f"q_{key}_{ion}"] = q.si.value
    num_regression.check(d_partial)


def test_diffusion_coefficient(fc, num_regression):
    d = {}
    for ion, D in fc.diffusion_coefficient.items():
        assert np.isfinite(D).all(), ion
        d[ion] = D.si.value
        assert D.unit.si == (u.m ** 2 / u.s)
    num_regression.check(d)


def test_thermal_coefficient(fc, num_regression):
    d = {}
    for ion, χ in fc.thermal_conductivity.items():
        if "0" in ion:
            continue
        assert np.isfinite(χ).all(), ion
        d[ion] = χ.si.value
        assert χ.si.unit == u.Unit("W / (K m)")
    num_regression.check(d)


@pytest.mark.xfail(
    raises=u.UnitsError,
    reason="need a tesla in the denominator",
)
def test_bootstrap_current(fc, num_regression):
    Ib = fc.bootstrap_current
    assert np.isfinite(Ib), ion
    current_density_unit = u.MA / u.m ** 2

    # if this crashes, we have replaced the current issue that Ib is actually <B * I_b> with another one
    (Ib.unit / u.T).to(current_density_unit)

    assert_quantity_allclose(Ib, -0.007 * current_density_unit)


def test_fluxes(fc, num_regression):
    d = {}
    for ion, (Γ, q) in fc.fluxes.items():
        if "0" in ion:
            continue
        assert np.isfinite(Γ).all(), ion
        assert np.isfinite(q).all(), ion
        d[f"Γ_{ion}"] = Γ.si.value
        d[f"q_{ion}"] = q.si.value
    num_regression.check(d)


def test_particle_velocities_heat_fluxes(fc, num_regression):
    d = {}
    for ion, particle_velocities in fc.local_flow_velocities.items():
        if "0" in ion:
            continue
        d[f"u_{ion}"] = particle_velocities.ravel().value
        d[f"q_{ion}"] = fc.local_heat_flux_components[ion].ravel().value
    num_regression.check(d)
