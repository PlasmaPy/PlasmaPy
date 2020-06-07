import numpy as np
import astropy.units as u
from plasmapy.formulary import chemical_potential

ne = u.Quantity(np.logspace(15,23, 100), u.m**-3)
T = 11000 * u.K

mu = [chemical_potential(N, T) for N in ne]

import matplotlib.pyplot as plt
plt.semilogx(ne, mu, "o-")
