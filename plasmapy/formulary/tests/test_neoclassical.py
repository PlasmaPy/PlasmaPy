import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.neoclassical import (
    charge_weighting_factor,
    effective_momentum_relaxation_rate,
    L_friction_coefficient,
    M_matrix,
    M_script,
    N_matrix,
    N_script,
    pitch_angle_diffusion_rate,
)
from plasmapy.particles import IonizationState, Particle, proton

hydrogen = IonizationState("p+", n_elem=1e20 * u.m ** -3, T_e=10 * u.eV)
pure_carbon = IonizationState("C+", n_elem=1e19 * u.m ** -3, T_e=10 * u.eV)
carbon_states = IonizationState(
    "C", [0, 1 / 1.1, 0.1 / 1.1, 0, 0, 0, 0], n_elem=1.1e19 * u.m ** -3, T_e=10 * u.eV
)
all_species = [hydrogen, carbon_states]


@pytest.mark.parametrize(
    ["function", "shape"],
    [(N_matrix, (3, 3)), (M_matrix, (3, 3)), (N_script, (3, 3)),],
)
def test_matrix_between_elements(function, shape, num_regression):
    data = function(hydrogen, pure_carbon)
    try:
        data = data.si.value
    except AttributeError:
        pass  # we're already a numpy array
    assert data.shape == shape, data.shape
    num_regression.check({function.__name__: data.flatten()})


@pytest.mark.parametrize(
    ["function", "shape"], [(M_script, (3, 3)),],
)
def test_matrix_between_element_and_all_species(function, shape, num_regression):
    data = function(hydrogen, all_species)
    assert data.shape == shape, data.shape
    num_regression.check({function.__name__: data.si.value.flatten()})


@pytest.mark.parametrize("function", [effective_momentum_relaxation_rate,])
def test_number_between_ionization_states(function, num_regression):
    data = function(hydrogen, carbon_states)
    num_regression.check({function.__name__: data.si.value})


def test_weighted_ionization_factor(num_regression):
    data = charge_weighting_factor(1, carbon_states)
    num_regression.check({"xi": data.si.value})
    assert data == 1 - charge_weighting_factor(2, carbon_states)


def test_L_friction_coefficient(num_regression):
    data = L_friction_coefficient(hydrogen, 1, carbon_states, 1, all_species)
    num_regression.check({"L": data.si.value.flatten()})


def test_pitch_angle_diffusion_rate(num_regression):
    x = np.logspace(-6, 6, 5000)
    ν_D_ai = pitch_angle_diffusion_rate(x, 1, carbon_states, all_species)
    num_regression.check(
        {"x": x, "ν_D_ai": ν_D_ai.si.value}, tolerances={"ν_D_ai": 1e-4}
    )
