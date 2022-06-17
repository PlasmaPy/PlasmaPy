"""Tests for the two fluid dispersion solution."""

import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.analytical.two_fluid_ import two_fluid
from plasmapy.formulary.frequencies import wc_
from plasmapy.formulary.speeds import cs_, va_
from plasmapy.particles import Particle
from plasmapy.utils.exceptions import PhysicsWarning


class TestTwoFluid:
    _kwargs_single_valued = {
        "B": 8.3e-9 * u.T,
        "ion": "p+",
        "k": 0.0001 * u.rad / u.m,
        "n_i": 5.0e6 * u.m**-3,
        "T_e": 1.6e6 * u.K,
        "T_i": 4.0e5 * u.K,
        "theta": 45 * u.deg,
    }
    _kwargs_bellan2012 = {
        "B": 400e-4 * u.T,
        "ion": Particle("He+"),
        "n_i": 6.358e19 * u.m**-3,
        "T_e": 20 * u.eV,
        "T_i": 10 * u.eV,
        "k": (2 * np.pi * u.rad) / (0.56547 * u.m),
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
            two_fluid(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, _warning",
        [
            # violates the low-frequency assumption (w/kc << 1)
            (
                {
                    "B": 8.3e-7 * u.T,
                    "ion": "p+",
                    "k": 0.0001 * u.rad / u.m,
                    "n_i": 3.0e6 * u.m**-3,
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
            two_fluid(**kwargs)

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
                {"fast_mode": 1.4000284, "alfven_mode": 0.0, "acoustic_mode": 0.0},
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

        cs = cs_(kwargs["T_e"], kwargs["T_i"], kwargs["ion"])
        va = va_(kwargs["B"], kwargs["n_i"], ion=kwargs["ion"])
        wci = wc_(kwargs["B"], kwargs["ion"])

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

        ws = two_fluid(**kwargs)
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
                {**_kwargs_bellan2012, "ion": Particle("He +2"), "theta": 0 * u.deg},
            ),
            #
            # z_mean defaults to 1
            (
                {**_kwargs_bellan2012, "ion": Particle("He"), "theta": 0 * u.deg},
                {**_kwargs_bellan2012, "ion": Particle("He+"), "theta": 0 * u.deg},
            ),
        ],
    )
    def test_z_mean_override(self, kwargs, expected):
        """Test overriding behavior of kw 'z_mean'."""
        ws = two_fluid(**kwargs)
        ws_expected = two_fluid(**expected)

        for mode in ws:
            assert np.isclose(ws[mode], ws_expected[mode], atol=0, rtol=1.7e-4)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({**_kwargs_bellan2012, "theta": 0 * u.deg}, {"shape": ()}),
            (
                {
                    **_kwargs_bellan2012,
                    "theta": 0 * u.deg,
                    "k": [1, 2, 3] * u.rad / u.m,
                },
                {"shape": (3,)},
            ),
            (
                {
                    **_kwargs_bellan2012,
                    "theta": [10, 20, 30, 40, 50] * u.deg,
                    "k": [1, 2, 3] * u.rad / u.m,
                },
                {"shape": (3, 5)},
            ),
            (
                {**_kwargs_bellan2012, "theta": [10, 20, 30, 40, 50] * u.deg},
                {"shape": (5,)},
            ),
        ],
    )
    def test_return_structure(self, kwargs, expected):
        """Test the structure of the returned values."""
        ws = two_fluid(**kwargs)

        assert isinstance(ws, dict)
        assert len({"acoustic_mode", "alfven_mode", "fast_mode"} - set(ws.keys())) == 0

        for mode, val in ws.items():
            assert isinstance(val, u.Quantity)
            assert val.unit == u.rad / u.s
            assert val.shape == expected["shape"]
