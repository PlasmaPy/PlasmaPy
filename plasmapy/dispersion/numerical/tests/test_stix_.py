"""Test functionality of Stix in `plasmapy.dispersion.numerical.stix_`."""
import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.numerical.stix_ import stix
from plasmapy.particles import Particle


class TestStix:
    _kwargs_single_valued = {
        "B": 8.3e-9 * u.T,
        "w": 0.001 * u.rad / u.s,
        "ions": [Particle("He+"), Particle("H+")],
        "n_i": [4.0e5, 2.0e5] * u.m ** -3,
        "theta": 30 * u.deg,
    }

    @pytest.mark.parametrize(
        "kwargs, _error",
        [
            ({**_kwargs_single_valued, "B": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "B": [8e-9, 8.5e-9] * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": -1 * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": 5 * u.m}, u.UnitTypeError),
            ({**_kwargs_single_valued, "w": -1.0 * u.rad / u.s}, ValueError),
            ({**_kwargs_single_valued, "w": 5 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "ions": {"not": "a particle"}}, TypeError),
            ({**_kwargs_single_valued, "n_i": "wrong type"}, TypeError),
            (
                {**_kwargs_single_valued, "n_i": 6 * u.m / u.s},
                u.UnitTypeError,
            ),
            ({**_kwargs_single_valued, "theta": 5 * u.eV}, u.UnitTypeError),
        ],
    )
    def test_raises(self, kwargs, _error):
        with pytest.raises(_error):
            stix(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({**_kwargs_single_valued, "w": 0 * u.rad / u.s}, {"shape": 1}),
            (
                {**_kwargs_single_valued, "w": [10] * u.rad / u.s},
                {"shape": 1},
            ),
            (
                {**_kwargs_single_valued, "w": [10, 20, 30] * u.rad / u.s},
                {"shape": 3},
            ),
            ({**_kwargs_single_valued, "ions": ["He+", "e-"]}, {"shape": 2}),
            (
                {
                    **_kwargs_single_valued,
                    "ions": ["He+"],
                    "n_i": [1] * u.m ** -3,
                },
                {"shape": 1},
            ),
            ({**_kwargs_single_valued, "ions": ["He+", "e-"]}, {"shape": 2}),
            (
                {
                    **_kwargs_single_valued,
                    "ions": ["He+", "H+"],
                    "n_i": [1, 2] * u.m ** -3,
                },
                {"shape": 2},
            ),
        ],
    )
    def test_return_structure(self, kwargs, expected):

        w = stix(**kwargs)

        assert isinstance(w, dict)

        for key in w.keys():
            for val in w[key]:
                assert isinstance(val, u.Quantity)
                assert val.unit == u.rad / u.m

        assert np.shape(w) == expected["shape"]

    @pytest.mark.parametrize(
        "kwargs, _warning",
        [
            ({**_kwargs_single_valued, "w": 0 * u.s / u.s}, u.UnitTypeError),
        ],
    )
    def test_warns(self, kwargs, _warning):
        with pytest.warns(_warning):
            stix(**kwargs)
