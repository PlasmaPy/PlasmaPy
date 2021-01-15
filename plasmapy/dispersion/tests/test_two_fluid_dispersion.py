"""Tests for the two fluid dispersion"""

import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.two_fluid_dispersion import (
    tfds_,
    two_fluid_dispersion_solution,
)
from plasmapy.particles import Particle
from plasmapy.formulary import parameters as pfp
from plasmapy.utils.exceptions import PhysicsWarning

k = 0.0001 * u.rad/u.m
B = 8.3e-9 * u.T
n_i = 5.0e6 * u.m ** -3
theta = np.pi / 4 * u.rad
T_e = 1.6e6 * u.K
T_i = 4.0e5 * u.K
z_mean = 1
ion = 'p+'

B_neg = -1 * u.T
n_neg = -5.0e6 * u.m ** -3
T_e_neg = -1.6e6 * u.K
T_i_neg = -4.0e5 * u.K

k_arr = np.linspace(10 ** -7, 10 ** -2, 10000) * u.rad/u.m
theta_arr = np.linspace(5, 85, 100) * u.deg
c = 3.0e8 * u.m / u.s


@pytest.mark.skip
def test_two_fluid_dispersion():
    r"""Test the two fluid analytical dispersion solution"""

    sol = two_fluid_dispersion_solution(
        B=B, ion=ion, k=k, n_i=n_i, T_e=T_e, T_i=T_i, theta=theta, z_mean=z_mean
    )

    sol_theta_0 = two_fluid_dispersion_solution(
        B=B, ion=ion, k=k, n_i=n_i, T_e=T_e, T_i=T_i, theta=0 * u.deg, z_mean=z_mean
    )

    sol_theta_90 = two_fluid_dispersion_solution(
        B=B, ion=ion, k=k, n_i=n_i, T_e=T_e, T_i=T_i, theta=90 * u.deg, z_mean=z_mean
    )

    sol_k_0 = two_fluid_dispersion_solution(
        B=B, ion=ion, k=0 * u.rad/u.m, n_i=n_i, T_e=T_e, T_i=T_i, theta=theta, z_mean=z_mean
    )

    assert np.isclose(sol["fast_mode"].value, 56.00705004, rtol=1e-6)
    assert np.isclose(sol["alfven_mode"].value, 15.13224648, rtol=1e-6)
    assert np.isclose(sol["acoustic_mode"].value, 0.55619463, rtol=1e-6)
    assert np.isclose(sol_theta_0["alfven_mode"].value, 15.20273637, rtol=1e-6)
    assert np.isclose(
        np.real(sol_theta_90["acoustic_mode"].value), 0.0, rtol=1e-6
    )

    assert np.isclose(sol_k_0["fast_mode"].value, 0, rtol=1e-6)
    assert np.isclose(sol_k_0["alfven_mode"].value, 0, rtol=1e-6)
    assert np.isclose(sol_k_0["acoustic_mode"].value, 0, rtol=1e-6)

    with pytest.raises(ValueError):
        two_fluid_dispersion_solution(
            B=B_neg, ion=ion, k=k, n_i=n_i, T_e=T_e, T_i=T_i, theta=theta, z_mean=z_mean
        )

    with pytest.raises(ValueError):
        two_fluid_dispersion_solution(
            B=B, ion=ion, k=k, n_i=n_neg, T_e=T_e, T_i=T_i, theta=theta, z_mean=z_mean
        )

    with pytest.raises(ValueError):
        two_fluid_dispersion_solution(
            B=B, ion=ion, k=k, n_i=n_i, T_e=T_e_neg, T_i=T_i, theta=theta, z_mean=z_mean
        )

    with pytest.raises(ValueError):
        two_fluid_dispersion_solution(
            B=B, ion=ion, k=k, n_i=n_i, T_e=T_e, T_i=T_i_neg, theta=theta, z_mean=z_mean
        )

    # Cases where one or more of the inputs is an array
    sol_k_arr = two_fluid_dispersion_solution(
        B=B, ion=ion, k=k_arr, n_i=n_i, T_e=T_e, T_i=T_i, theta=theta, z_mean=z_mean
    )
    sol_theta_arr = two_fluid_dispersion_solution(
        B=B, ion=ion, k=k, n_i=n_i, T_e=T_e, T_i=T_i, theta=theta_arr, z_mean=z_mean
    )
    sol_k_arr0 = two_fluid_dispersion_solution(
        B=B, ion=ion, k=k_arr[0], n_i=n_i, T_e=T_e, T_i=T_i, theta=theta, z_mean=z_mean
    )
    sol_theta_arr0 = two_fluid_dispersion_solution(
        B=B, ion=ion, k=k, n_i=n_i, T_e=T_e, T_i=T_i, theta=theta_arr[0], z_mean=z_mean
    )
    sol_k_arr1 = two_fluid_dispersion_solution(
        B=B, ion=ion, k=k_arr[-1], n_i=n_i, T_e=T_e, T_i=T_i, theta=theta, z_mean=z_mean
    )
    sol_theta_arr1 = two_fluid_dispersion_solution(
        B=B, ion=ion, k=k, n_i=n_i, T_e=T_e, T_i=T_i, theta=theta_arr[-1], z_mean=z_mean
    )


    assert np.isclose(
        sol_k_arr["fast_mode"].value[0], sol_k_arr0["fast_mode"].value
    )
    assert np.isclose(
        sol_theta_arr["fast_mode"].value[0], sol_theta_arr0["fast_mode"].value
    )

    assert np.isclose(
        sol_k_arr["alfven_mode"].value[0], sol_k_arr0["alfven_mode"].value
    )
    assert np.isclose(
        sol_theta_arr["alfven_mode"].value[0],
        sol_theta_arr0["alfven_mode"].value,
    )

    assert np.isclose(
        sol_k_arr["acoustic_mode"].value[0], sol_k_arr0["acoustic_mode"].value
    )
    assert np.isclose(
        sol_theta_arr["acoustic_mode"].value[0],
        sol_theta_arr0["acoustic_mode"].value,
    )

    assert np.isclose(
        sol_k_arr["fast_mode"].value[-1], sol_k_arr1["fast_mode"].value
    )

    assert np.isclose(
        sol_theta_arr["fast_mode"].value[-1], sol_theta_arr1["fast_mode"].value
    )

    assert np.isclose(
        sol_k_arr["alfven_mode"].value[-1], sol_k_arr1["alfven_mode"].value
    )
    assert np.isclose(
        sol_theta_arr["alfven_mode"][-1].value,
        sol_theta_arr1["alfven_mode"].value,
    )

    assert np.isclose(
        sol_k_arr["acoustic_mode"].value[-1], sol_k_arr1["acoustic_mode"].value
    )
    assert np.isclose(
        sol_theta_arr["acoustic_mode"].value[-1],
        sol_theta_arr1["acoustic_mode"].value,
    )


class TestTwoFluidDispersionSolution:
    _kwargs_single_valued = {
        "B": 8.3e-9 * u.T,
        "ion": "p+",
        "k": 0.0001 * u.rad/u.m,
        "n_i": 5.0e6 * u.m ** -3,
        "T_e": 1.6e6 * u.K,
        "T_i": 4.0e5 * u.K,
        "theta": 45 * u.deg,
    }
    _kwargs_bellan2012 = {
        "B": 400e-4 * u.T,
        "ion": Particle("He+"),
        "n_i": 6.358e19 * u.m ** -3,
        "T_e": 20 * u.eV,
        "T_i": 10 * u.eV,
        "k": (2 * np.pi * u.rad) / (0.56547 * u.m),
    }

    def test_alias(self):
        """Test the associated alias is defined correctly."""
        assert tfds_ is two_fluid_dispersion_solution

    @pytest.mark.parametrize(
        "kwargs, _error",
        [
            ({**_kwargs_single_valued, "B": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "B": [8e-9, 8.5e-9] * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": -1 * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": 5 * u.m}, u.UnitTypeError),
            ({**_kwargs_single_valued, "ion": {"not": "a particle"}}, TypeError),
            ({**_kwargs_single_valued, "ion": "e-"}, ValueError),
            ({**_kwargs_single_valued, "ion": "He", "z_mean": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "k": np.ones((3, 2)) * u.rad / u.m}, ValueError),
            ({**_kwargs_single_valued, "k": 0 * u.rad / u.m}, ValueError),
            ({**_kwargs_single_valued, "k": -1.0 * u.rad / u.m}, ValueError),
            ({**_kwargs_single_valued, "k": 5 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "n_i": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "n_i": [5e6, 6e6] * u.m ** -3}, ValueError),
            ({**_kwargs_single_valued, "n_i": -5e6 * u.m ** -3}, ValueError),
            ({**_kwargs_single_valued, "n_i": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "T_e": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "T_e": [1.4e6, 1.7e6] * u.K}, ValueError),
            ({**_kwargs_single_valued, "T_e": -10 * u.eV}, ValueError),
            ({**_kwargs_single_valued, "T_e": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "T_i": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "T_i": [4e5, 5e5] * u.K}, ValueError),
            ({**_kwargs_single_valued, "T_i": -1 * u.eV}, ValueError),
            ({**_kwargs_single_valued, "T_i": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "theta": np.ones((3, 2)) * u.deg}, ValueError),
            ({**_kwargs_single_valued, "theta": 5 * u.eV}, u.UnitTypeError),
            ({**_kwargs_single_valued, "gamma_e": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "gamma_i": "wrong type"}, TypeError),
        ],
    )
    def test_raises(self, kwargs, _error):
        """Test scenarios that raise an `Exception`."""
        with pytest.raises(_error):
            two_fluid_dispersion_solution(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, _warning",
        [
            (
                {
                    "B": 8.3e-7 * u.T,
                    "ion": "p+",
                    "k": 0.0001 * u.rad/u.m,
                    "n_i": 3.0e6 * u.m ** -3,
                    "T_e": 1.6e6 * u.K,
                    "T_i": 4.0e5 * u.K,
                    "theta": 5 * u.deg,
                },
                PhysicsWarning,
            ),
        ],
    )
    def test_warns(self, kwargs, _warning):
        """Test scenarios the issue a `Warning`."""
        with pytest.warns(_warning):
            two_fluid_dispersion_solution(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            (
                {**_kwargs_bellan2012, "theta": 0 * u.deg},
                {
                    "fast_mode": 1.8631944,
                    "alfven_mode": 0.5366538,
                    "acoustic_mode": 0.4000832,
                },
            ),
            (
                {**_kwargs_bellan2012, "theta": 90 * u.deg},
                {
                    "fast_mode": 1.4000284,
                    "alfven_mode": 0.0,
                    "acoustic_mode": 0.0,
                },
            ),
        ],
    )
    def test_on_bellan2012_vals(self, kwargs, expected):
        """
        Test calculated values based on Figure 1 of Bellan 2012
        (DOI: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2012JA017856).
        """
        # theta and k values need to be single valued for this test to function
        # correctly

        cs = pfp.cs_(kwargs["T_e"], kwargs["T_i"], kwargs["ion"])
        va = pfp.va_(kwargs["B"], kwargs["n_i"], ion=kwargs["ion"])
        wci = pfp.wc_(kwargs["B"], kwargs["ion"])

        beta = (cs / va).value ** 2
        if not np.isclose(beta, 0.4, atol=1e-4):
            pytest.fail(
                f"The Bellan 2012 paper requires a 'beta' value of 0.4 and the test "
                f"parameters yielded {beta:.6f}."
            )

        Lambda = (kwargs["k"] * va / wci).value ** 2
        if not np.isclose(Lambda, 0.4, atol=1e-4):
            pytest.fail(
                f"The Bellan 2012 paper requires a 'Lambda' value of 0.4 and the test "
                f"parameters yielded {Lambda:.6f}."
            )

        ws = two_fluid_dispersion_solution(**kwargs)
        for mode, val in ws.items():
            norm = (np.absolute(val) / (kwargs["k"] * va)).value ** 2
            assert np.isclose(norm, expected[mode])

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            (
                    {
                        **_kwargs_bellan2012,
                        "ion": Particle("He"),
                        "z_mean": 2.0,
                        "theta": 0 * u.deg,
                    },
                    {
                        **_kwargs_bellan2012,
                        "ion": Particle("He +2"),
                        "theta": 0 * u.deg,
                    },
            ),
            (
                    {
                        **_kwargs_bellan2012,
                        "ion": Particle("He"),
                        "theta": 0 * u.deg,
                    },
                    {
                        **_kwargs_bellan2012,
                        "ion": Particle("He+"),
                        "theta": 0 * u.deg,
                    },
            ),
        ],
    )
    def test_z_mean_override(self, kwargs, expected):
        """Test overriding behavior of kw 'z_mean'."""
        ws = two_fluid_dispersion_solution(**kwargs)
        ws_expected = two_fluid_dispersion_solution(**expected)

        for mode in ws:
            assert np.isclose(ws[mode], ws_expected[mode], atol=0, rtol=1.7e-4)
