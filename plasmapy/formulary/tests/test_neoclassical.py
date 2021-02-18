import astropy.units as u
import pytest

from plasmapy.formulary.neoclassical import (
    charge_weighting_factor,
    effective_momentum_relaxation_rate,
    M_matrix,
    N_matrix,
)
from plasmapy.particles import IonizationState, Particle, proton

hydrogen = IonizationState("p+", n_elem=1e20 * u.m ** -3, T_e=10 * u.eV)
pure_carbon = IonizationState("C+", n_elem=1e19 * u.m ** -3, T_e=10 * u.eV)
carbon_states = IonizationState(
    "C", [0, 1 / 1.1, 0.1 / 1.1, 0, 0, 0, 0], n_elem=1.1e19 * u.m ** -3, T_e=10 * u.eV
)


@pytest.mark.parametrize(
    ["function", "shape"], [(N_matrix, (3, 3)), (M_matrix, (3, 3)),],
)
def test_matrix_between_elements(function, shape, num_regression):
    data = function(hydrogen, pure_carbon)
    assert data.shape == shape, data.shape
    num_regression.check({function.__name__: data.flatten()})


@pytest.mark.parametrize("function", [effective_momentum_relaxation_rate,])
def test_number_between_ionization_states(function, num_regression):
    data = function(hydrogen, carbon_states)
    num_regression.check({function.__name__: data.si.value})


def test_weighted_ionization_factor(num_regression):
    data = charge_weighting_factor(1, carbon_states)
    num_regression.check({"xi": data.si.value})
