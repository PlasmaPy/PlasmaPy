import astropy.units as u
import datetime
import hypothesis
import numpy as np
import pytest

from astropy.tests.helper import assert_quantity_allclose
from hypothesis import example, given, settings
from hypothesis import strategies as st

from plasmapy.particles import IonizationStateCollection, Particle
from plasmapy.transport.Houlberg1997 import ExtendedParticleList

all_species = ExtendedParticleList(
    [Particle("C 1+"), Particle("C 2+"), Particle("C 3+"), Particle("p+"), ],
    10 * u.eV,
    u.Quantity([1e20/1.11, 0.1e20/1.11, 0.01e20/1.11, 1e20], u.m**-3),
)
x = np.logspace(-3, 6, 50)

@pytest.mark.parametrize(
    ["function", "shape"],
    [
        ("N_matrix", (3, 3, 4, 4)),
        ("M_matrix", (3, 3, 4, 4)),
        ("N_script", (3, 3, 4, 4)),
        ("M_script", (3, 3, 4)),
        ("effective_momentum_relaxation_rate", (4,)),
        ("ξ", (4,)),
    ],
)
def test_matrix_between_elements(function, shape, num_regression):
    data = getattr(all_species, function)
    try:
        data = data.si.value
    except AttributeError:
        pass  # we're already a numpy array
    assert data.shape == shape, data.shape
    num_regression.check({function: data.flatten()})


@pytest.mark.parametrize(
    ["function", "shape", "rtol", "args"],
    [
        ("pitch_angle_diffusion_rate", (50, 4), 1e-6, []),
        ("ν_T_ai", (50, 4), 1e-6, []),
        ("K_B_ai", (50, 4), 1e-6, [0.6234941403639689]),
    ],
)
def test_function_of_relative_velocity(num_regression, function, shape, rtol, args):
    if args:
        data = getattr(all_species, function)(x, *args)
    else:
        data = getattr(all_species, function)(x)

    try:
        data = data.si.value
    except AttributeError:
        pass  # we're already a numpy array
    assert data.shape == shape, data.shape
    num_regression.check({function: data.flatten()})

@pytest.mark.parametrize(
    ["function", "shape", "rtol", "args"],
    [
        ("K_ps_ai", (50, 4), 1e-6, []),
    ],
)
def test_function_on_flux_surface(num_regression, flux_surface, function, shape, rtol, args):
    if args:
        data = getattr(all_species, function)(x, flux_surface, *args)
    else:
        data = getattr(all_species, function)(x, flux_surface)

    try:
        data = data.si.value
    except AttributeError:
        pass  # we're already a numpy array
    assert data.shape == shape, data.shape
    num_regression.check({function: data.flatten()})


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
    x = np.array([x])
    result = ν_T_ai(x, hydrogen, all_species)[1]
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
    x = np.array([x])
    result = K_ps_ai(x, hydrogen, all_species, flux_surface)
    assert (result[1:] > 0).all()
    assert np.isfinite(result[1:]).all()
    second_result = K_ps_ai(x, hydrogen, all_species, flux_surface)
    assert_quantity_allclose(result, second_result)


def test_mu(flux_surface, num_regression):
    μ = mu_hat(hydrogen, all_species, flux_surface)
    num_regression.check({"mu": μ.si.value.flatten()})
