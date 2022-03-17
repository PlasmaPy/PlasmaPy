"""Test functinality of Stix in `plasmapy.dispersion.numerical.stix_`."""
import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.numerical.stix_ import stix


class TestStix:
    _kwargs_single_valued = {
        "B": 8.3e-9 * u.T,
        "k": 0.001 * u.rad / u.m,
        "ions": ["e-", "H+"],
        "omega_ions": [4.0e5, 2.0e5] * u.rad / u.s,
        "theta": 30 * u.deg,
    }

    @pytest.mark.parametrize(
        "kwargs, _error",
        [
            ({**_kwargs_single_valued, "B": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "B": [8e-9, 8.5e-9] * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": -1 * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": 5 * u.m}, u.UnitTypeError),
            ({**_kwargs_single_valued, "k": -1.0 * u.rad / u.m}, ValueError),
            ({**_kwargs_single_valued, "k": 5 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "ions": {"not": "a particle"}}, TypeError),
            ({**_kwargs_single_valued, "ions": "e-"}, ValueError),
            ({**_kwargs_single_valued, "omega_ions": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "omega_ions": 6 * u.rad / u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "theta": np.ones((3, 2)) * u.deg}, ValueError),
            ({**_kwargs_single_valued, "theta": 5 * u.eV}, u.UnitTypeError),
        ],
    )
    def test_raises(self, kwargs, _error):
        with pytest.raises(_error):
            stix(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({**_kwargs_single_valued, "k": 0 * u.rad / u.m}, {"shape": ()}),
            (
                {**_kwargs_single_valued, "k": [10] * u.rad / u.m},
                {"shape": (1,)},
            ),
            (
                {**_kwargs_single_valued, "k": [10, 20, 30] * u.rad / u.m},
                {"shape": (3,)},
            ),
            ({**_kwargs_single_valued, "ions": "He+"}, {"shape": ()}),
            (
                {
                    **_kwargs_single_valued,
                    "ions": "He+",
                    "omega_ions": [1] * u.rad / u.s,
                },
                {"shape": (1,)},
            ),
            ({**_kwargs_single_valued, "ions": "He+"}, {"shape": ()}),
            (
                {
                    **_kwargs_single_valued,
                    "ions": ["He+", "H+"],
                    "omega_ions": [1, 2] * u.rad / u.s,
                },
                {"shape": (2,)},
            ),
        ],
    )
    def test_return_structure(self, kwargs, expected):

        w = stix(**kwargs)

        assert isinstance(w, dict)

        for val in w.values():
            assert isinstance(val, u.Quantity)
            assert val.unit == u.rad / u.s
            assert val.shape == expected["shape"]

    @pytest.mark.parametrize(
        "kwargs, _warning",
        [
            ({**_kwargs_single_valued, "k": 0 * u.rad / u.m}, u.UnitTypeError),
        ],
    )
    def test_warns(self, kwargs, _warning):
        with pytest.warns(_warning):
            stix(**kwargs)


print("Out")
