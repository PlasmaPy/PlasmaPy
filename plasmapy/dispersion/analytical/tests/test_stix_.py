"""Test functionality of Stix in `plasmapy.dispersion.numerical.stix_`."""
import numpy as np
import pytest

from astropy import units as u
from astropy.constants.si import c

from plasmapy.dispersion.analytical.stix_ import stix
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError

c_unitless = c.value


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
            ({**_kwargs_single_valued, "w": [-1, 2] * u.rad / u.s}, ValueError),
            ({**_kwargs_single_valued, "w": np.ones((2, 2)) * u.rad / u.s}, ValueError),
            ({**_kwargs_single_valued, "w": 5 * u.s}, u.UnitTypeError),
            (
                {**_kwargs_single_valued, "ions": {"not": "a particle"}},
                InvalidParticleError,
            ),
            ({**_kwargs_single_valued, "n_i": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "n_i": 6 * u.m / u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "theta": 5 * u.eV}, u.UnitTypeError),
            (
                {**_kwargs_single_valued, "theta": np.ones((2, 2)) * u.rad},
                TypeError,
            ),
            ({**_kwargs_single_valued, "ions": Particle("e-")}, ValueError),
            ({**_kwargs_single_valued, "n_i": [4, 2, 3] * u.m ** -3}, ValueError),
            (
                {**_kwargs_single_valued, "n_i": np.ones((2, 2)) * u.m ** -3},
                ValueError,
            ),
        ],
    )
    def test_raises(self, kwargs, _error):
        with pytest.raises(_error):
            stix(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({**_kwargs_single_valued, "w": 2 * u.rad / u.s}, {"shape": (4,)}),
            (
                {**_kwargs_single_valued, "w": [10] * u.rad / u.s},
                {"shape": (4,)},
            ),
            (
                {**_kwargs_single_valued, "w": [10, 20, 30] * u.rad / u.s},
                {"shape": (3, 4)},
            ),
            ({**_kwargs_single_valued, "ions": ["He+", "H+"]}, {"shape": (4,)}),
            (
                {
                    **_kwargs_single_valued,
                    "ions": ["He+"],
                    "n_i": [1] * u.m ** -3,
                },
                {"shape": (4,)},
            ),
            ({**_kwargs_single_valued, "ions": ["He+", "H+"]}, {"shape": (4,)}),
            (
                {
                    **_kwargs_single_valued,
                    "ions": ["He+", "H+"],
                    "n_i": [1, 2] * u.m ** -3,
                },
                {"shape": (4,)},
            ),
            ({**_kwargs_single_valued, "theta": [10, 20, 30]}, {"shape": (3, 4)}),
            (
                {
                    **_kwargs_single_valued,
                    "w": [10, 20],
                    "theta": 10 * u.rad,
                },
                {"shape": (2, 4)},
            ),
        ],
    )
    def test_return_structure(self, kwargs, expected):
        k = stix(**kwargs)

        assert isinstance(k, u.Quantity)
        assert np.shape(k) == expected["shape"]
        assert k.unit == u.rad / u.m

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            (
                {"theta": 0 * u.deg, "gamma": 1000, "beta": 1000, "mu": 1000.5},
                {
                    "n_1": 31.638591944805466,
                    "n_2": 1414.5674255763668,
                },
            ),
            (
                {"theta": np.pi / 2 * u.deg, "gamma": 1000, "beta": 1000, "mu": 1000.5},
                {
                    "n_1": 44.73253849828242,
                    "n_2": 0,
                },
            ),
        ],
    )
    def test_vals(self, kwargs, expected):

        gamma = kwargs["gamma"]
        beta = kwargs["beta"]
        mu = kwargs["mu"]

        R = (
            1
            - (gamma * beta ** 2) / (beta + 1)
            + (gamma * beta ** 2) / (beta - (1 / mu))
        )
        L = (
            1
            + (gamma * beta ** 2) / (beta - 1)
            + (gamma * beta ** 2) / (beta + (1 / mu))
        )
        P = 1 - gamma * beta ** 2 - gamma * mu * beta ** 2

        S = 0.5 * (R + L)

        if kwargs["theta"].value == 0:
            n_1 = np.sqrt(R)
            n_2 = np.sqrt(L)
        elif kwargs["theta"].value == np.pi / 2:
            n_1 = np.sqrt(R * L / S)
            n_2 = np.sqrt(P)

        assert np.isclose(n_1, expected["n_1"])
        assert np.isclose(n_2, expected["n_2"])
