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
