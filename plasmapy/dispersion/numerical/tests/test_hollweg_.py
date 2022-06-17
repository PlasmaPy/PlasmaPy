"""Tests for the hollweg dispersion solution."""

import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.numerical.hollweg_ import hollweg
from plasmapy.formulary import speeds
from plasmapy.particles import Particle
from plasmapy.utils.exceptions import PhysicsWarning


class TestHollweg:
    _kwargs_single_valued = {
        "k": 0.01 * u.rad / u.m,
        "theta": 88 * u.deg,
        "n_i": 5 * u.cm**-3,
        "B": 2.2e-8 * u.T,
        "T_e": 1.6e6 * u.K,
        "T_i": 4.0e5 * u.K,
        "ion": Particle("p+"),
    }

    _kwargs_hollweg1999 = {
        "theta": 90 * u.deg,
        "n_i": 5 * u.cm**-3,
        "T_e": 1.6e6 * u.K,
        "T_i": 4.0e5 * u.K,
        "ion": Particle("p+"),
    }

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
            ({**_kwargs_single_valued, "n_i": [5e6, 6e6] * u.m**-3}, ValueError),
            ({**_kwargs_single_valued, "n_i": -5e6 * u.m**-3}, ValueError),
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
            hollweg(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, _warning",
        [
            # w/w_ci << 1 PhysicsWarning
            (
                {
                    "k": 0.01 * u.rad / u.m,
                    "theta": 88 * u.deg,
                    "n_i": 0.05 * u.cm**-3,
                    "B": 2.2e-8 * u.T,
                    "T_e": 1.6e6 * u.K,
                    "T_i": 4.0e5 * u.K,
                    "ion": Particle("p+"),
                },
                PhysicsWarning,
            ),
            # c_s/v_A << 1 PhysicsWarning
            (
                {
                    "k": 10e-8 * u.rad / u.m,
                    "theta": 88 * u.deg,
                    "n_i": 5 * u.cm**-3,
                    "B": 6.98e-8 * u.T,
                    "T_e": 1.6e6 * u.K,
                    "T_i": 4.0e5 * u.K,
                    "ion": Particle("p+"),
                },
                PhysicsWarning,
            ),
            # theta nearly perpendicular PhysicsWarning
            (
                {
                    "k": 10e-8 * u.rad / u.m,
                    "theta": 84 * u.deg,
                    "n_i": 1 * u.cm**-3,
                    "B": 6.98e-8 * u.T,
                    "T_e": 1.6e6 * u.K,
                    "T_i": 4.0e5 * u.K,
                    "ion": Particle("p+"),
                },
                PhysicsWarning,
            ),
        ],
    )
    def test_warning(self, kwargs, _warning):
        """Test scenarios that raise a `Warning`."""
        with pytest.warns(_warning):
            hollweg(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            # k is an array, theta is single valued
            (
                {
                    **_kwargs_single_valued,
                    "k": np.logspace(-7, -2, 2) * u.rad / u.m,
                },
                {
                    "fast_mode": [2.62911663e-02 + 0.0j, 2.27876968e03 + 0.0j],
                    "alfven_mode": [7.48765909e-04 + 0.0j, 2.13800404e03 + 0.0j],
                    "acoustic_mode": [0.00043295 + 0.0j, 0.07358991 + 0.0j],
                },
            ),
            # theta is an array, k is single valued
            (
                {**_kwargs_single_valued, "theta": [87, 88] * u.deg},
                {
                    "fast_mode": [3406.43522969 + 0.0j, 2278.76967883 + 0.0j],
                    "alfven_mode": [2144.81200575 + 0.0j, 2138.00403666 + 0.0j],
                    "acoustic_mode": [0.11044097 + 0.0j, 0.07358991 + 0.0j],
                },
            ),
            # k and theta are an array
            (
                {
                    **_kwargs_single_valued,
                    "k": np.logspace(-7, -2, 2),
                    "theta": [86, 87, 88] * u.deg,
                },
                {
                    "fast_mode": [
                        [
                            2.62804756e-02 + 0.0j,
                            2.62867114e-02 + 0.0j,
                            2.62911663e-02 + 0.0j,
                        ],
                        [
                            4.53954617e03 + 0.0j,
                            3.40643523e03 + 0.0j,
                            2.27876968e03 + 0.0j,
                        ],
                    ],
                    "alfven_mode": [
                        [
                            1.49661942e-03 + 0.0j,
                            1.12286371e-03 + 0.0j,
                            7.48765909e-04 + 0.0j,
                        ],
                        [
                            2.14516382e03 + 0.0j,
                            2.14481201e03 + 0.0j,
                            2.13800404e03 + 0.0j,
                        ],
                    ],
                    "acoustic_mode": [
                        [
                            0.00086572 + 0.0j,
                            0.00064937 + 0.0j,
                            0.00043295 + 0.0j,
                        ],
                        [
                            0.14735951 + 0.0j,
                            0.11044097 + 0.0j,
                            0.07358991 + 0.0j,
                        ],
                    ],
                },
            ),
        ],
    )
    def test_handle_k_theta_arrays(self, kwargs, expected):
        """Test scenarios involving k and theta arrays."""
        ws = hollweg(**kwargs)
        for mode, val in ws.items():
            assert np.allclose(val.value, expected[mode])

    @pytest.mark.parametrize(
        "kwargs, expected, desired_beta",
        [
            (  # beta = 1/20 for kx*L = 0
                {**_kwargs_hollweg1999, "k": 1e-14 * u.rad / u.m, "B": 6.971e-8 * u.T},
                1 + 0j,
                1 / 20,
            ),
            (  # beta = 1/20 for kx*L = 1
                {
                    **_kwargs_hollweg1999,
                    "k": 0.0000439223874624874 * u.rad / u.m,
                    "B": 6.971e-8 * u.T,
                },
                1.4018 + 0j,
                1 / 20,
            ),
            (  # beta = 1/2 for kx*L = 0
                {**_kwargs_hollweg1999, "k": 1e-14 * u.rad / u.m, "B": 2.205e-8 * u.T},
                1 + 0j,
                0.5,
            ),
            (  # beta = 1/2 for kx*L = 1
                {
                    **_kwargs_hollweg1999,
                    "k": 0.000013893109303101 * u.rad / u.m,
                    "B": 2.205e-8 * u.T,
                },
                1.3536 + 0j,
                0.5,
            ),
            (  # beta = 2 for kx*L = 0
                {
                    **_kwargs_hollweg1999,
                    "k": 1e-14 * u.rad / u.m,
                    "B": 1.10232e-8 * u.T,
                },
                1 + 0j,
                2,
            ),
            (  # beta = 2 for kx*L = 1
                {
                    **_kwargs_hollweg1999,
                    "k": 0.00000691190063354451 * u.rad / u.m,
                    "B": 1.10232e-8 * u.T,
                },
                1.2607 + 0j,
                2,
            ),
            (  # beta = 1/2000 for kx*L = 0
                {
                    **_kwargs_hollweg1999,
                    "k": 1e-14 * u.rad / u.m,
                    "B": 6.97178e-7 * u.T,
                },
                1 + 0j,
                1 / 2000,
            ),
            (  # beta = 1/2000 for kx*L = 1
                {
                    **_kwargs_hollweg1999,
                    "k": 0.000439273010336778 * u.rad / u.m,
                    "B": 6.97178e-7 * u.T,
                },
                0.98750 + 0j,
                1 / 2000,
            ),
        ],
    )
    def test_hollweg1999_vals(self, kwargs, expected, desired_beta):
        """
        Test calculated values based on Figure 2 of Hollweg1999
        (DOI: https://doi.org/10.1029/1998JA900132) using eqn 3 of
        Bellan 2012
        (DOI: https://doi.org/10.1029/2012JA017856).

        The WebPlotDigitizer software was used to determine the test
        parameters for k, B, and expected omega from Figure 2 of
        Hollweg1999.

        - GitHub: https://github.com/ankitrohatgi/WebPlotDigitizer
        - Web Version: https://automeris.io/WebPlotDigitizer/
        """
        # k values need to be single valued for this test to function correctly

        cs = speeds.cs_(kwargs["T_e"], kwargs["T_i"], kwargs["ion"]).value
        va = speeds.va_(kwargs["B"], kwargs["n_i"], ion=kwargs["ion"]).value

        beta = (cs / va) ** 2
        if not np.isclose(beta, desired_beta, atol=2e-4):
            pytest.fail(
                f"The Holweg 1999 paper requires a 'beta' value of {desired_beta:0.5f} "
                f"and the test parameters yielded {beta:.6f}."
            )

        kz = (np.cos(kwargs["theta"]) * kwargs["k"]).value

        w_alfven = (hollweg(**kwargs)["alfven_mode"]).value
        big_omega = np.abs(w_alfven / (kz * va))

        assert np.allclose(big_omega, expected, atol=1e-2)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            (
                {
                    **_kwargs_single_valued,
                    "ion": Particle("He"),
                    "z_mean": 2.0,
                },
                {**_kwargs_single_valued, "ion": Particle("He +2")},
            ),
            #
            # z_mean defaults to 1
            (
                {**_kwargs_single_valued, "ion": Particle("He")},
                {**_kwargs_single_valued, "ion": Particle("He+")},
            ),
        ],
    )
    def test_z_mean_override(self, kwargs, expected):
        """Test overriding behavior of kw 'z_mean'."""
        ws = hollweg(**kwargs)
        ws_expected = hollweg(**expected)

        for mode in ws:
            assert np.isclose(ws[mode], ws_expected[mode], atol=1e-5, rtol=1.7e-4)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({**_kwargs_single_valued}, {"shape": ()}),
            (
                {
                    **_kwargs_single_valued,
                    "k": [1, 2, 3] * u.rad / u.m,
                },
                {"shape": (3,)},
            ),
            (
                {
                    **_kwargs_single_valued,
                    "k": [1, 2, 3] * u.rad / u.m,
                    "theta": [50, 77] * u.deg,
                },
                {"shape": (3, 2)},
            ),
            (
                {
                    **_kwargs_single_valued,
                    "theta": [50, 77] * u.deg,
                },
                {"shape": (2,)},
            ),
        ],
    )
    def test_return_structure(self, kwargs, expected):
        """Test the structure of the returned values."""
        ws = hollweg(**kwargs)

        assert isinstance(ws, dict)
        assert len({"acoustic_mode", "alfven_mode", "fast_mode"} - set(ws.keys())) == 0

        for mode, val in ws.items():
            assert isinstance(val, u.Quantity)
            assert val.unit == u.rad / u.s
            assert val.shape == expected["shape"]
