import numpy as np
from astropy import units as u
from numerical import kinetic_alfven_ as kinetic_alfven

inputs = {
    "k": np.logspace(-7, -2, 2) * u.rad / u.m,
    "theta": 30 * u.deg,
    "B": 8.3e-9 * u.T,
    "n_i": 5 * u.m ** -3,
    "T_e": 1.6e6 * u.K,
    "T_i": 4.0e5 * u.K,
    "ion": Particle("p+"),
}
omegas = kinetic_alfven(**inputs)
print(omegas)
