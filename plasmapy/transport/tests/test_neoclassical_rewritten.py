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
    [
        Particle("C 1+"),
        Particle("C 2+"),
        Particle("C 3+"),
        Particle("p+"),
    ],
    10 * u.eV,
    u.Quantity([1e20 / 1.11, 0.1e20 / 1.11, 0.01e20 / 1.11, 1e20], u.m ** -3),
)
x = np.logspace(-3, 6, 50)


@pytest.mark.parametrize(
    ["function", "shape"],
    [
        ("N_matrix", (3, 3, 2, 2)),
        ("M_matrix", (3, 3, 2, 2)),
        ("N_script", (3, 3, 2, 2)),
        ("M_script", (3, 3, 2)),
        ("effective_momentum_relaxation_rate", (2, 2)),
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
        ("K", (50, 4), 1e-6, []),
    ],
)
def test_function_on_flux_surface(
    num_regression, flux_surface, function, shape, rtol, args
):
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


def test_mu_hat(num_regression, flux_surface):
    data = all_species.mu_hat(flux_surface).si
    assert data.unit == u.Unit("kg / (m3 s)")
    assert data.shape == (3, 3, 4)
    num_regression.check({"mu_hat": data.value.flatten()})


def test_split_sum(num_regression):
    N = len(all_species)
    M = 9
    assert (all_species.split_index == np.array([0, 3])).all()

    arr = np.arange(N * M).reshape(N, M)
    split = all_species.split_isotopes(arr, axis=0)
    assert len(split) == 2

    assert all_species.num_isotopes == 2

    splitsum = all_species.compress(arr, axis=0)
    assert splitsum.shape == (2, M)


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
    result = all_species.ν_T_ai(x)
    assert np.all(result > 0)
    assert np.isfinite(result).all()


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
    result = all_species.K_ps_ai(x, flux_surface)
    # assert (result > 0).all()
    assert np.isfinite(result).all()
    second_result = all_species.K_ps_ai(x, flux_surface)
    assert_quantity_allclose(result, second_result)
