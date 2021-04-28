import plasmaboundaries

from plasmapy.plasma.symbolicequilibrium import SymbolicEquilibrium
from plasmapy.transport.flows import FlowCalculator

params = plasmaboundaries.ITER.copy()
equilibrium = SymbolicEquilibrium(**params, B0=5.2, config="single-null")
flux_surface = equilibrium.get_flux_surface(psi_value=-0.01)
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

fc = FlowCalculator(
    all_species,
    flux_surface,
    density_gradient,
    temperature_gradient,
    mu_N=1000,
)
fc._charge_state_flows
