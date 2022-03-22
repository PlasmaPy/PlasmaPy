from astropy import units as u

from plasmapy.dispersion.numerical.stix_ import stix
from plasmapy.particles import Particle

inputs = {
    "B": 5e-9 * u.T,
    "k": 0.001 * u.rad / u.m,
    "species": Particle("H+"),
    "omega_species": 1.0 * u.rad / u.s,
    "theta": 30 * u.deg,
}
w = stix(**inputs)
print(w[0.001])
