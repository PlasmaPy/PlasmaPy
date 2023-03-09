"""Test functionality of Stix in `plasmapy.formulary.collisions.collisional_analysis`."""
import numpy as np
import pytest

from astropy import units as u
from astropy.constants.si import c

from plasmapy.formulary.collisions.tools.collisional_analysis import collisional_thermalization
from plasmapy.particles import Particle, ParticleList
from plasmapy.particles.exceptions import InvalidParticleError

c_si_unitless = c.value


class Testcollisional_thermalizastion:
    _kwargs_single_valued = {
        "n_a": [1, 2, 3] * u.cm ** -3,
        "n_b": [0.5, 0.6, 0.4] * u.cm ** -3,
        "v_a": [400, 500, 600] * u.m / u.s,
        "T_a": [50, 55, 60] * 10**3 * u.K,
        "T_b": [30, 40, 60] * 10**3 * u.K,
        "ions": [Particle("p+"), Particle("He-4++")],
        "r_0": [0.1, 0.1, 0.1] * u.au,
        "r_n": [1.0, 1.0, 1.0] * u.au,
    }

    @pytest.mark.parametrize(
        "kwargs, _error",
        [
            ({**_kwargs_single_valued, "B": "wrong type"}, TypeError),
        ],
    )
    def test_raises(self, kwargs, _error):
        with pytest.raises(_error):
            collisional_thermalization(**kwargs)


