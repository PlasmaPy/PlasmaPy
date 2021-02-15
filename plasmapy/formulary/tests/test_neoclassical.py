import astropy.units as u
import pytest

from plasmapy.formulary.neoclassical import (
    collision_frequency_ai_bj,
    effective_momentum_relaxation_rate,
    M_matrix,
    N_matrix,
)
from plasmapy.particles import Particle, proton, Species

hydrogen = Species(proton, 1e20 * u.m ** -3, temperature=10 * u.eV)
carbon = Species(Particle("C+"), 1e19 * u.m ** -3, temperature=10 * u.eV)
carbon2 = Species(Particle("C++"), 1e18 * u.m ** -3, temperature=10 * u.eV)
hydrogen_states = [hydrogen]
carbon_states = [carbon, carbon2]


@pytest.mark.parametrize(
    ["function", "shape"], [(N_matrix, (3, 3)), (M_matrix, (3, 3)),],
)
def test_matrix_between_elements(function, shape, num_regression):
    data = function(hydrogen, carbon)
    assert data.shape == shape, data.shape
    num_regression.check({function.__name__: data.flatten()})


@pytest.mark.parametrize(
    "function", [collision_frequency_ai_bj,],
)
def test_number(function, num_regression):
    data = function(hydrogen, carbon)
    num_regression.check({function.__name__: data.si.value})


@pytest.mark.parametrize("function", [effective_momentum_relaxation_rate,])
def test_number_between_ionization_states(function, num_regression):
    data = function(hydrogen_states, carbon_states)
    num_regression.check({function.__name__: data.si.value})
