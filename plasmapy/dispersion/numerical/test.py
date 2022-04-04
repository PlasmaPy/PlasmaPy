from astropy import units as u

from plasmapy.dispersion.numerical.stix_ import stix
from plasmapy.particles import Particle

inputs = {
    "B": 5e-9 * u.T,
    "k": 0.001 * u.rad / u.m,
    "ions": Particle("H+"),
    "n_i": 1.0 * u.m**-3,
    "theta": 30 * u.deg,
}
w = stix(**inputs)
print(w[0.001])
