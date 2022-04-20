from astropy import units as u

from plasmapy.dispersion.numerical.stix_ import stix
from plasmapy.particles import Particle

inputs = {
    "B": 8.3e-9 * u.T,
    "w": 0.001 * u.rad / u.s,
    "ions": [Particle("H+"), Particle("He+")],
    "n_i": [4.0e5, 2.0e5] * u.m ** -3,
    "theta": 30 * u.deg,
}
w = stix(**inputs)
print(w)
