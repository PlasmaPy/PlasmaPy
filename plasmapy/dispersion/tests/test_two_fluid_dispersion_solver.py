"""Tests for the two fluid dispersion solver"""

import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.two_fluid_dispersion_solver import (
two_fluid_dispersion_solution,
)

k = 0.01 * u.m ** -1
theta = 30 * u.deg
B = 8.3E-9 * u.T
n = 5 * u.cm ** -3
T_e = 1.6e6 * u.K
T_i = 4.e5 * u.K
z = 1

B_neg = -1 * u.T
n_neg = -5 * u.cm ** -3
T_e_neg = -1.6e6 * u.K
T_i_neg = -4.e5 * u.K
z_neg = -1

k_arr = np.linspace(10**-7, 10**-2, 1E4) * u.m ** -1
theta_arr = np.linspace(5, 85, 100) * u.deg
n = 5 * u.cm ** -3
B = 8.3E-9 * u.T
T_e = 1.6e6 * u.K
T_i = 4.e5 * u.K
z = 1
c = 3.e8 * u.m/u.s

def test_two_fluid_dispersion_solution():
    r"""Test the two fluid analytical dispersion solution"""

    sol = two_fluid_dispersion_solution(B=B, k=k, n=n, T_e=T_e, T_i=T_i,
    theta=theta, z=z)

    sol_theta_0 = two_fluid_dispersion_solution(B=B, k=k, n=n, T_e=T_e, T_i=T_i,
    theta=0 * u.deg, z=z)

    sol_theta_90 = two_fluid_dispersion_solution(B=B, k=k, n=n, T_e=T_e, T_i=T_i,
    theta=90 * u.deg, z=z)

    sol_B_0 = two_fluid_dispersion_solution(B=0 * u.T, k=k, n=n, T_e=T_e, T_i=T_i,
    theta=theta, z=z)

    sol_k_0 = two_fluid_dispersion_solution(B=B, k=0 * u.m ** -1, n=n, T_e=T_e,
    T_i=T_i, theta=theta, z=z)


    assert np.isclose(sol['fast_mode'].value[0,0], 1520.57692, rtol=1e-6)
    assert np.isclose(sol['alfven_mode'].value[0,0], 1260.01546, rtol=1e-6)
    assert np.isclose(sol['acoustic_mode'].value[0,0], 0.68815, rtol=1e-6)
    assert np.isnan(sol_theta_0['alfven_mode'].value[0,0], 1455.23097, 
    rtol=1e-6)
    assert np.isnan(sol_theta_90['acoustic_mode'].value[0,0], np.nan)

    assert np.isnan(sol_B_0['fast_mode'].value[0,0], np.nan)
    assert np.isnan(sol_B_0['alfven_mode'].value[0,0], np.nan)
    assert np.isnan(sol_B_0['acoustic_mode'].value[0,0], np.nan)

    assert np.isclose(sol_k_0['fast_mode'].value[0,0], 0, rtol=1e-6)
    assert np.isclose(sol_k_0['alfven_mode'].value[0,0], 0, rtol=1e-6)
    assert np.isclose(sol_k_0['acoustic_mode'].value[0,0], 0, rtol=1e-6)

    with pytest.raises(ValueError):
        two_fluid_dispersion_solution(B=B_neg, k=k, n=n, T_e=T_e, T_i=T_i,
        theta=theta, z=z)

with pytest.raises(ValueError):
        two_fluid_dispersion_solution(B=B, k=k, n=n_neg, T_e=T_e, T_i=T_i,
        theta=theta, z=z)

with pytest.raises(ValueError):
        two_fluid_dispersion_solution(B=B, k=k, n=n, T_e=T_e_neg, T_i=T_i,
        theta=theta, z=z)

with pytest.raises(ValueError):
        two_fluid_dispersion_solution(B=B, k=k, n=n, T_e=T_e, T_i=T_i_neg,
        theta=theta, z=z)

with pytest.raises(ValueError):
        two_fluid_dispersion_solution(B=B, k=k, n=n, T_e=T_e, T_i=T_i,
        theta=theta, z=z_neg)

# Cases where one or more of the inputs is an array
sol_k_arr = two_fluid_dispersion_solution(B=B, k=k_arr, n=n, T_e=T_e, T_i=T_i,
theta=theta, z=z)
sol_theta_arr = two_fluid_dispersion_solution(B=B, k=k, n=n, T_e=T_e, T_i=T_i,
theta=theta_arr, z=z)
sol_k_arr0 = two_fluid_dispersion_solution(B=B, k=k_arr[0], n=n, T_e=T_e, T_i=T_i,
theta=theta, z=z)
sol_theta_arr0 = two_fluid_dispersion_solution(B=B, k=k, n=n, T_e=T_e, T_i=T_i,
theta=theta_arr[0], z=z)
sol_k_arr1 = two_fluid_dispersion_solution(B=B, k=k_arr[-1], n=n, T_e=T_e, T_i=T_i,
theta=theta, z=z)
sol_theta_arr1 = two_fluid_dispersion_solution(B=B, k=k, n=n, T_e=T_e, T_i=T_i,
theta=theta_arr[-1], z=z)

assert np.isclose(sol_k_arr['fast_mode'].value[0,0],
sol_k_arr0['fast_mode'].value[0,0])
assert np.isclose(sol_theta_arr['fast_mode'].value[0,0],
sol_theta_arr0['fast_mode'].value[0,0])

assert np.isclose(sol_k_arr['alfven_mode'].value[0,0],
sol_k_arr0['alfven_mode'].value[0,0])
assert np.isclose(sol_theta_arr['alfven_mode'].value[0,0],
sol_theta_arr0['alfven_mode'].value[0,0])

assert np.isclose(sol_k_arr['acoustic_mode'].value[0,0],
sol_k_arr0['acoustic_mode'].value[0,0])
assert np.isclose(sol_theta_arr['acoustic_mode'].value[0,0],
sol_theta_arr0['acoustic_mode'].value[0,0])

assert np.isclose(sol_k_arr['fast_mode'].value[-1,0],
sol_k_arr1['fast_mode'].value[0,0])
assert np.isclose(sol_theta_arr['fast_mode'].value[0,-1],
sol_theta_arr1['fast_mode'].value[0,0])

assert np.isclose(sol_k_arr['alfven_mode'].value[-1,0],
sol_k_arr1['alfven_mode'].value[0,0])
assert np.isclose(sol_theta_arr['alfven_mode'].value[0,-1],
sol_theta_arr1['alfven_mode'].value[0,0])

assert np.isclose(sol_k_arr['acoustic_mode'].value[-1,0],
sol_k_arr1['acoustic_mode'].value[0,0])
assert np.isclose(sol_theta_arr['acoustic_mode'].value[0,-1],
sol_theta_arr1['acoustic_mode'].value[0,0])
