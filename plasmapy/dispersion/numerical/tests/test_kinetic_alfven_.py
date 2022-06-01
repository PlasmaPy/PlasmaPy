"""Test functionality of Stix in `plasmapy.dispersion.numerical.stix_`."""
import numpy as np
import pytest

from astropy import units as u
from astropy.constants.si import c

from plasmapy.dispersion.numerical.kinetic_alfven_ import kinetic_alfven
from plasmapy.formulary import gyrofrequency, plasma_frequency
from plasmapy.particles import Particle
from plasmapy.particles.exceptions import InvalidParticleError

c_si_unitless = c.value


class TestStix:
    _kwargs_single_valued = {
        "k": np.logspace(-7, -2, 2) * u.rad / u.m,
        "theta": 30 * u.deg,
        "B": 8.3e-9 * u.T,
        "n_i": 5 * u.m ** -3,
        "T_e": 1.6e6 * u.K,
        "T_i": 4.0e5 * u.K,
        "ion": Particle("p+"),
    }
