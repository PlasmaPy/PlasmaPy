import astropy
import astropy.units as u
import datetime
import hypothesis
import numpy as np
import pytest

from astropy.tests.helper import assert_quantity_allclose
from hypothesis import example, given, settings
from hypothesis import strategies as st

from plasmapy.particles import IonizationStateCollection, Particle
from plasmapy.transport.flows import FlowCalculator
from plasmapy.transport.Houlberg1997 import ExtendedParticleList
import warnings

@pytest.fixture(scope="module",
    params = [
        ExtendedParticleList(
            [
                Particle("C 1+"),
                Particle("p+"),
            ],
            10 * u.eV,
            u.Quantity([1e20 / 1.11,  1e20], u.m ** -3),
            dn=u.Quantity([1e18, 1e18]) * u.m ** -4,
            dT=u.Quantity([-10, -10]) * u.K / u.m,
        ),
        ExtendedParticleList(
            [
                Particle("C 1+"),
                Particle("C 2+"),
                Particle("p+"),
            ],
            10 * u.eV,
            u.Quantity([1e20 / 1.11, 0.1e20 / 1.11, 1e20], u.m ** -3),
            dn=u.Quantity([1e18, 1e18, 1e18]) * u.m ** -4,
            dT=u.Quantity([-10, -10, -10]) * u.K / u.m,
        ),
        ExtendedParticleList(
            [
                Particle("C 1+"),
                Particle("C 2+"),
                Particle("C 3+"),
                Particle("p+"),
            ],
            10 * u.eV,
            u.Quantity([1e20 / 1.11, 0.1e20 / 1.11, 0.01e20 / 1.11, 1e20], u.m ** -3),
            dn=u.Quantity([1e18, 1e18, 1e18, 1e18]) * u.m ** -4,
            dT=u.Quantity([-10, -10, -10, -10]) * u.K / u.m,
        ),
    ],
    ids = ["2particles", "3particles", "4particles"],
)
def all_species(request):
    return request.param

@pytest.fixture(scope="module")
def fc(all_species, flux_surface):

    fc = FlowCalculator(
        all_species,
        flux_surface,
        mu_N=1000,
    )
    return fc


def test_get_flows(fc, num_regression):
    assert fc._charge_state_flows.to(u.V / u.m)
    num_regression.check(
        {"_charge_state_flows": fc._charge_state_flows.si.value.ravel()}
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
    num_regression.check(
        dict(
            particle=fluxes.particle_flux.si.value,
            heat=fluxes.heat_flux.si.value,
        )
    )


def test_diffusion_coefficient(fc, num_regression):
    D = fc.diffusion_coefficient
    assert D.to(u.m ** 2 / u.s)
    num_regression.check({"diffusion_coefficient": D.si.value})


def test_thermal_coefficient(fc, num_regression):
    Chi = fc.thermal_conductivity
    assert Chi.to(u.Unit("W / (K m)"))
    num_regression.check({"thermal_conductivity": Chi.si.value})


def test_bootstrap_current(fc, num_regression):
    Ib = fc.bootstrap_current
    assert np.isfinite(Ib), ion
    current_density_unit = u.MA / u.m ** 2

    # if this crashes, we have replaced the current issue that Ib is actually <B * I_b> with another one
    (Ib.unit / u.T).to(current_density_unit)

    num_regression.check({"bootstrap_current": Ib.si.value})


def test_fluxes(fc, num_regression):
    Γ, q = fc.fluxes
    d = dict(
        Γ=Γ.si.value,
        q=q.si.value,
    )
    num_regression.check(d)


def test_particle_velocities_heat_fluxes(fc, num_regression):
    d = {}
    d[f"u"] = fc.local_flow_velocities.ravel().value
    d[f"q"] = fc.local_heat_flux_components.ravel().value
    num_regression.check(d)


# @pytest.fixture(scope="module")
# def fc_with_electrons(flux_surface):
#     all_species = IonizationStateCollection(
#         {
#             "H": [0, 1],
#             "C": [0, 1 / 1.1, 0.1 / 1.1, 0, 0, 0, 0],
#             "e-": [1],
#         },
#         n0=1e20 * u.m ** -3,
#         abundances={"H": 1, "C": 0.11},
#         T_e=10 * u.eV,
#     )
#     hydrogen = all_species["H"]
#     carbon_states = all_species["C"]

#     density_gradient = {
#         "H": np.ones(2) * 1e18 * u.m ** -3 / u.m,
#         "C": np.ones(7) * 1e18 * u.m ** -3 / u.m,
#     }
#     temperature_gradient = {
#         "H": np.ones(2) * -10 * u.K / u.m,
#         "C": np.ones(7) * -10 * u.K / u.m,
#     }

#     fc = FlowCalculator(
#         all_species,
#         flux_surface,
#         density_gradient,
#         temperature_gradient,
#         mu_N=1000,
#     )
#     return fc


@pytest.mark.slow
def test_integration():
    import astropy.units as u
    import numpy as np
    import pandas as pd
    import xarray

    data_df = pd.read_csv("/home/dominik/home-import/IFPILM/Magisterka/HoulbergNSTX.csv", index_col=0).iloc[
        ::10, :
    ]

    T_i = T_e = T_C6 = data_df["ToverT0"] * 0.5 * 1000
    dT_i = dT_e = dT_C6 = data_df["ToverT0_derivative"] * 0.5 * 1000
    n_e, dn_e = data_df["n_e"], data_df["n_e_derivative"]
    n_i, dn_i = data_df["n_D"], data_df["n_D_derivative"]
    n_C6, dn_C6 = data_df["n_C"], data_df["n_C_derivative"]

    rho = data_df.x.values

    ## Multiple flux surfaces - radial grid

    from plasmaboundaries import NSTX_single_null

    NSTX_single_null

    NSTX_Bt0 = 0.3 * u.T
    NSTX_R0 = 0.8 * u.m
    NSTX_a0 = 0.64 * u.m
    NSTX_I = 1 * u.MA
    from plasmapy.plasma.symbolicequilibrium import SymbolicEquilibrium

    params = {"aspect_ratio": 1.25, "A": -0.05, "elongation": 2, "triangularity": 0.3}
    # TODO this is still not taken in

    eq = SymbolicEquilibrium(
        **NSTX_single_null,
        B0=NSTX_Bt0.si.value,  # TODO handle quantity input
        config="single-null",
    )
    rminmaxstep = (
        0.1,
        1.9,
        0.001,
    )  # these definitely, unfortunately, need to be moved into SymbolicEquilibrium
    zminmaxstep = (-2, 2, 0.001)

    surfaces = list(eq.get_multiple_flux_surfaces(
                rho_values=rho, rminmaxstep=rminmaxstep, zminmaxstep=zminmaxstep
            ))

    ## `FlowCalculator`

    rho_to_surface = {key: value for key, value in surfaces}

    import xarray

    dataset_H1 = xarray.Dataset(
        {
            "T": ("rho", T_i),
            "gradT": ("rho", dT_i),
            "n": ("rho", n_i),
            "gradn": ("rho", dn_i),
        },
        coords={"rho": rho, "particle": "H 1+"},
        attrs={
            "T unit": u.eV,
            "n unit": u.m ** -3,
            "gradT unit": u.eV / u.m,
            "gradn unit": u.m ** -3 / u.m,
        },
    )

    dataset_C6 = xarray.Dataset(
        {
            "T": ("rho", T_C6),
            "gradT": ("rho", dT_C6),
            "n": ("rho", n_C6),
            "gradn": ("rho", dn_C6),
        },
        coords={"rho": rho, "particle": "C 6+"},
        attrs={
            "T unit": u.eV,
            "n unit": u.m ** -3,
            "gradT unit": u.eV / u.m,
            "gradn unit": u.m ** -3 / u.m,
        },
    )
    dataset_e = xarray.Dataset(
        {
            "T": ("rho", T_e),
            "gradT": ("rho", dT_e),
            "n": ("rho", n_e),
            "gradn": ("rho", dn_e),
        },
        coords={"rho": rho, "particle": "e-"},
        attrs={
            "T unit": u.eV,
            "n unit": u.m ** -3,
            "gradT unit": u.eV / u.m,
            "gradn unit": u.m ** -3 / u.m,
        },
    )

    dataset = xarray.concat([dataset_H1, dataset_C6, dataset_e], dim="particle")
    # dataset["rho"] = ("rho", rho)
    dataset["charges"] = ("particle", [1, 6, -1])
    dataset["charge_density"] = "rho", (dataset.charges * dataset.n).sum(dim="particle")
    dataset

    dataset["psi"] = eq.rho_to_psi(rho)
    dataset

    fcs = []

    for i, (ψ, surface) in enumerate(surfaces):
        if surface is None:
            print(f"Skipping {i}-th surface at {ψ}")
            continue
        try:
            with pytest.warns(UserWarning, match="was not sorted"):
                fcs.append(
                    FlowCalculator.from_xarray_surface(dataset.isel(rho=i), surface)
                )
        except ImportError as e:
            display(e)


    datasets = [fc.to_dataset() for fc in fcs]

    results = xarray.concat(datasets, dim="rho")
    scaling = (fcs[0].bootstrap_current.unit / NSTX_Bt0).to(u.MA / u.m ** 2)
    results = results.assign(
        bootstrap_current_normalized=results.bootstrap_current * scaling
    )

    results.diffusion_coefficient.sel(particle="C 6+")

    import pandas as pd

    df = pd.read_csv("/home/dominik/home-import/IFPILM/Magisterka/NSTXplot1.csv")

    df.plot.line(x="x")
    results.bootstrap_current_normalized.plot.line(
        x="rho", label="My bootstrap current"
    )
