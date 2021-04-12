import astropy.units as u
import datetime
import hypothesis
import numpy as np
import pytest

from astropy.tests.helper import assert_quantity_allclose
from hypothesis import example, given, settings
from hypothesis import strategies as st

from plasmapy.formulary.neoclassical import (
    effective_momentum_relaxation_rate,
    get_flows,
    K_B_ai,
    K_ps_ai,
    M_matrix,
    M_script,
    N_matrix,
    N_script,
    pitch_angle_diffusion_rate,
    ν_T_ai,
)
from plasmapy.particles import IonizationStateCollection

all_species = IonizationStateCollection(
    {
        "H": [0, 1],
        #      "D": [0, 1],   raises ParticleError, why?
        "C": [0, 1 / 1.1, 0.1 / 1.1, 0, 0, 0, 0],
    },
    n0=1e20 * u.m ** -3,
    abundances={"H": 1, "C": 0.11},
    T_e=10 * u.eV,
)
hydrogen = all_species["H"]
carbon_states = all_species["C"]


@pytest.mark.parametrize(
    ["function", "shape"],
    [(N_matrix, (3, 3)), (M_matrix, (3, 3)), (N_script, (3, 3)),],
)
def test_matrix_between_elements(function, shape, num_regression):
    data = function(hydrogen, carbon_states)
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


def test_pitch_angle_diffusion_rate_and_banana_vsicosity(num_regression, flux_surface):
    x = np.logspace(-3, 6, 5000)
    ν_D_ai = pitch_angle_diffusion_rate(x, 1, carbon_states, all_species)
    k = K_B_ai(x, 1, carbon_states, all_species, flux_surface)
    num_regression.check(
        {"x": x, "ν_D_ai": ν_D_ai.si.value, "K_B_ai": k.si.value},
        tolerances={"ν_D_ai": {"rtol": 1e-4}, "K_B_ai": {"rtol": 1e-4}},
    )


@given(
    x=st.floats(
        min_value=1e-3,
        max_value=1e3,
        exclude_min=True,
        allow_nan=False,
        allow_infinity=False,
    )
)
@example(x=1)
@example(x=0.1)
@example(x=0.01)
@settings(deadline=datetime.timedelta(milliseconds=1000))
def test_ν_T_ai(x):
    result = ν_T_ai(x, 1, hydrogen, all_species)
    assert result > 0
    assert np.isfinite(result)


@given(
    x=st.floats(
        min_value=1e-3,
        max_value=1e3,
        exclude_min=True,
        allow_nan=False,
        allow_infinity=False,
    )
)
@example(x=270.08574852208653)
@example(x=684.765468434412)
@settings(deadline=datetime.timedelta(milliseconds=1000))
def test_K_ps_ai(x, flux_surface):
    result = K_ps_ai(x, 1, hydrogen, all_species, flux_surface)
    assert result > 0
    assert np.isfinite(result)
    second_result = K_ps_ai(x, 1, hydrogen, all_species, flux_surface)
    assert_quantity_allclose(result, second_result)


def test_get_flows(flux_surface):
    result = get_flows(
        all_species,
        flux_surface,
        density_gradient={
            "H 1+": 1e18 * u.m ** -3 / u.m,
            "C 1+": 1e18 * u.m ** -3 / u.m,
        },
        temperature_gradient={"H 1+": -10 * u.K / u.m, "C 1+": -10 * u.K / u.m,},
    )
    for r in result.values():
        assert np.isfinite(r).all()
