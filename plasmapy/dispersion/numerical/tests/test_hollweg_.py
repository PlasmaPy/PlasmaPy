"""Tests for the hollweg dispersion solution."""

import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.numerical.hollweg_ import hollweg
from plasmapy.formulary import parameters as pfp
from plasmapy.particles import Particle
from plasmapy.utils.exceptions import PhysicsWarning


class TestHollweg:
    _kwargs_single_valued = {
        # Values may need to be changed
        "k": 0.01 * u.rad / u.m,
        "theta": 88 * u.deg,
        "n_i": 5 * u.cm ** -3,
        "B": 2.2e-8 * u.T,
        "T_e": 1.6e6 * u.K,
        "T_i": 4.0e5 * u.K,
        "ion": Particle("p+"),
    }
    
    _kwargs_hollweg1999 = {
        "k": 10e-2 * u.rad / u.m,
        "theta": 88 * u.deg,
        "n_i": 5 * u.cm ** -3,
        "B": 6.98e-8 * u.T,
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
            hollweg(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, _warning",
        [
            # check the low-frequency limit (w<<w_ci)
            (
                {
                    "k": 0.01 * u.rad / u.m,
                    "theta": 88 * u.deg,
                    "n_i": 5 * u.cm ** -3,
                    "B": 2.2e-8 * u.T,
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
    
   """ 
   @pytest.mark.parametrize(
        "kwargs, expected",
        [
            
        ]
    )
    def test_on_hollweg1999_vals(self, kwargs, expected):
        """
        Test calculated values based on Figure 2 of Hollweg 1999
        (DOI: https://doi.org/10.1029/1998JA900132).
        """
        beta = []
        B_vals = [6.97178e-7, 6.971e-8, 2.205e-8, 1.097e-8]
        for x in range(0,4):
            # Need different cs and va values for different betas. Just with B changed.
            cs = pfp.cs_(kwargs["T_e"], kwargs["T_i"], kwargs["ion"])
            va = pfp.va_(kwargs["B"], kwargs["n_i"], ion=kwargs["ion"])
            beta[x] = (cs / va).value ** 2
    """
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
