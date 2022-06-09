"""Test functionality of Stix in `plasmapy.dispersion.numerical.kinetic_alfven_`."""
import numpy as np
import pytest

from astropy import units as u
from astropy.constants.si import c

from plasmapy.dispersion.numerical.kinetic_alfven_ import kinetic_alfven
from plasmapy.formulary import gyrofrequency, plasma_frequency
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError

c_si_unitless = c.value


class TestKinetic_Alfven:
    _kwargs_single_valued = {
        "B": 8.3e-9 * u.T,
        "ion": Particle("p+"),
        "k": np.logspace(-7, -2, 2) * u.rad / u.m,
        "n_i": 5 * u.m ** -3,
        "T_e": 1.6e6 * u.K,
        "T_i": 4.0e5 * u.K,
        "theta": 30 * u.deg,
        "gamma_e": 3,
        "gamma_i": 3,
        "z_mean": 1,
    }

    @pytest.mark.parametrize(
        "kwargs, _error",
        [
            ({**_kwargs_single_valued, "B": "wrong type"}, TypeError),
            ({**_kwargs_single_valued, "B": [1e-9, 2e-9, 3e-9] * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": -1 * u.T}, ValueError),
            ({**_kwargs_single_valued, "B": 5 * u.m}, u.UnitTypeError),
            ({**_kwargs_single_valued, "ion": "not a particle"}, InvalidParticleError),
            ({**_kwargs_single_valued, "ion": "e-"}, ValueError),
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
            kinetic_alfven(**kwargs)