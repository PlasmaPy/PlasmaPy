"""Test functionality of timescales in `plasmapy.formulary.collisions.timescales`."""
import numpy as np
import pytest

from astropy import units as u
from astropy.constants.si import c

from plasmapy.formulary.collisions.timescales import Hellinger
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError


class TestTimescales:
    version = 2009
    args_2009 = {
        "T": 8 * u.K,
        "n_i": 4.0e5 * u.m**-3,
        "ions": [Particle("He+"), Particle("H+")],
        "par_speeds": [30 * u.m / u.s, 50 * u.m / u.s],
    }

    _kwargs_2009 = {"version": version, "kwargs": args_2009}

    @pytest.mark.parametrize(
        "version, kwargs, _error",
        [
            (
                    "wrong type", _kwargs_2009["kwargs"], TypeError,
            ),
            (
                    2000, _kwargs_2009["kwargs"], ValueError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "T": "wrong type"},
                TypeError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "T": [8e-9, 8.5e-9] * u.K},
                ValueError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "T": -1 * u.K},
                ValueError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "n_i": "wrong type"},
                TypeError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "n_i": 6 * u.m / u.s},
                u.UnitTypeError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "n_i": [4, 2, 3] * u.m**-3},
                ValueError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "n_i": np.ones((2, 2)) * u.m**-3},
                ValueError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "ions": {"not": "a particle"}},
                InvalidParticleError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "par_speeds": [10, 40, 60] * u.m / u.s},
                ValueError,
            ),
            (
                _kwargs_2009["version"],
                {**_kwargs_2009["kwargs"], "par_speeds": [10, 40] * u.rad / u.s},
                u.UnitTypeError,
            ),

        ],
    )
    def test_raises(self, version, kwargs, _error):
        with pytest.raises(_error):
            Hellinger(kwargs, version)
