"""Tests for the hirose dispersion solution."""

import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.numerical.hirose_ import hirose
from plasmapy.formulary import speeds
from plasmapy.particles import Particle
from plasmapy.utils.exceptions import PhysicsWarning


class TestHirose:
    _kwargs_single_valued = {
        "k": 0.01 * u.rad / u.m,
        "theta": 88 * u.deg,
        "n_i": 5 * u.cm ** -3,
        "B": 2.2e-8 * u.T,
        "T_e": 1.6e6 * u.K,
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
            ({**_kwargs_single_valued, "n_i": [5e6, 6e6] * u.m ** -3}, ValueError),
            ({**_kwargs_single_valued, "n_i": -5e6 * u.m ** -3}, ValueError),
            ({**_kwargs_single_valued, "n_i": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "T_e": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "T_e": [1.4e6, 1.7e6] * u.K}, ValueError),
            ({**_kwargs_single_valued, "T_e": -10 * u.eV}, ValueError),
            ({**_kwargs_single_valued, "T_e": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "T_i": "anything"}, TypeError),
            ({**_kwargs_single_valued, "theta": np.ones((3, 2)) * u.deg}, ValueError),
            ({**_kwargs_single_valued, "theta": 5 * u.eV}, u.UnitTypeError),
            ({**_kwargs_single_valued, "gamma_e": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "gamma_i": "wrong type"}, TypeError),
        ],
    )
    def test_raises(self, kwargs, _error):
        """Test scenarios that raise an `Exception`."""
        with pytest.raises(_error):
            hirose(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, _warning",
        [
            # w/w_ci << 1 PhysicsWarning
            (
                {
                    "k": 0.01 * u.rad / u.m,
                    "theta": 88 * u.deg,
                    "n_i": 0.05 * u.cm ** -3,
                    "B": 2.2e-8 * u.T,
                    "T_e": 1.6e6 * u.K,
                    "ion": Particle("p+"),
                },
                PhysicsWarning,
            ),
        ],
    )
    def test_warning(self, kwargs, _warning):
        """Test scenarios that raise a `Warning`."""
        with pytest.warns(_warning):
            hirose(**kwargs)

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
        ws = hirose(**kwargs)
        ws_expected = hirose(**expected)

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
        ws = hirose(**kwargs)

        assert isinstance(ws, dict)
        assert len({"acoustic_mode", "alfven_mode", "fast_mode"} - set(ws.keys())) == 0

        for mode, val in ws.items():
            assert isinstance(val, u.Quantity)
            assert val.unit == u.rad / u.s
            assert val.shape == expected["shape"]
