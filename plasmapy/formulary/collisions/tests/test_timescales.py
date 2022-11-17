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

    _kwargs_single_valued = {"version": version, "kwargs": args_2009}

    @pytest.mark.parametrize(
        "version, kwargs, _error",
        [
            (
                _kwargs_single_valued["version"],
                {**_kwargs_single_valued["kwargs"], "T": "wrong type"},
                TypeError,
            ),
            (
                _kwargs_single_valued["version"],
                {**_kwargs_single_valued["kwargs"], "T": [8e-9, 8.5e-9] * u.K},
                ValueError,
            ),
            (
                _kwargs_single_valued["version"],
                {**_kwargs_single_valued["kwargs"], "T": -1 * u.K},
                ValueError,
            ),
            (
                _kwargs_single_valued["version"],
                {**_kwargs_single_valued["kwargs"], "n_i": "wrong type"},
                TypeError,
            ),
            (
                _kwargs_single_valued["version"],
                {**_kwargs_single_valued["kwargs"], "n_i": 6 * u.m / u.s},
                u.UnitTypeError,
            ),
            (
                _kwargs_single_valued["version"],
                {**_kwargs_single_valued["kwargs"], "n_i": [4, 2, 3] * u.m**-3},
                ValueError,
            ),
            (
                _kwargs_single_valued["version"],
                {**_kwargs_single_valued["kwargs"], "n_i": np.ones((2, 2)) * u.m**-3},
                ValueError,
            ),
            (
                _kwargs_single_valued["version"],
                {**_kwargs_single_valued["kwargs"], "ions": {"not": "a particle"}},
                InvalidParticleError,
            ),
        ],
    )
    def test_raises(self, version, kwargs, _error):
        with pytest.raises(_error):
            Hellinger(kwargs, version)
