"""Test functionality of `plasmapy.formulary.collisions.helio.collisional_analysis`."""
import pytest

from astropy import units as u
from astropy.constants.si import c

from plasmapy.formulary.collisions.helio.collisional_analysis import (
    thermalization_ratio,
)
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError

c_si_unitless = c.value


class Testcollisional_thermalizastion:
    _kwargs_single_valued = {
        "r_0": [0.1, 0.1, 0.1] * u.au,
        "r_n": [1.0, 1.0, 1.0] * u.au,
        "n_1": [1, 2, 3] * u.cm**-3,
        "n_2": [0.5, 0.6, 0.4] * u.cm**-3,
        "v_1": [400, 500, 600] * u.m / u.s,
        "T_1": [50, 55, 60] * 10**3 * u.K,
        "T_2": [30, 40, 60] * 10**3 * u.K,
        "ions": [Particle("p+"), Particle("He-4++")],
        "n_step": 1000,
        "density_scale": -1.8,
        "velocity_scale": -0.2,
        "temperature_scale": -0.77,
    }

    @pytest.mark.parametrize(
        "kwargs, _error",
        [
            ({**_kwargs_single_valued, "n_1": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "n_1": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "n_2": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "n_2": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "v_1": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "v_1": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "T_1": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "T_1": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "T_2": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "T_2": 2 * u.s}, u.UnitTypeError),
            ({**_kwargs_single_valued, "ions": [Particle("p+"), "He"]}, ValueError),
            (
                {**_kwargs_single_valued, "ions": [Particle("p+"), "not a particle"]},
                InvalidParticleError,
            ),
            ({**_kwargs_single_valued, "ions": [Particle("p+")]}, ValueError),
            ({**_kwargs_single_valued, "n_step": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "n_step": 4.5}, TypeError),
            ({**_kwargs_single_valued, "density_scale": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "velocity_scale": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "temperature_scale": "wrong type"}, TypeError),
        ],
    )
    def test_raises(self, kwargs, _error):
        with pytest.raises(_error):
            thermalization_ratio(**kwargs)
