"""
1D Maxwellian distribution function
===================================

We import the usual modules, and the hero of this notebook,
the Maxwellian 1D distribution:
"""


import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from plasmapy.constants import (m_e, k_B)

from plasmapy.theory.physics.distribution import Maxwellian_1D


############################################################
# Given we'll be plotting, import astropy's quantity support:


from astropy.visualization import quantity_support
quantity_support()

############################################################
# As a first example, let's get the probability density of
# finding an electron with a speed of 1 m/s if we have a
# plasma at a temperature of 30 000 K:


p_dens = Maxwellian_1D(v=1 * u.m / u.s,
                       T=30000 * u.K,
                       particle='e',
                       v_drift=0 * u.m / u.s)
print(p_dens)

############################################################
# Note the units! Integrated over speed, this will give us a
# probability. Let's test that for a bunch of particles:

T = 3e4 * u.K
dv = 10 * u.m / u.s
v = np.arange(-5e6, 5e6, 10) * u.m / u.s

############################################################
# Check that the integral over all speeds is 1
# (the particle has to be somewhere):

for particle in ['p', 'e']:
    pdf = Maxwellian_1D(v, T=T, particle=particle)
    integral = (pdf).sum() * dv
    print(f"Integral value for {particle}: {integral}")
    plt.plot(v, pdf, label=particle)
plt.legend()


############################################################
# The standard deviation of this distribution should give us back the
# temperature:


std = np.sqrt((Maxwellian_1D(v, T=T, particle='e') * v ** 2 * dv).sum())
T_theo = (std ** 2 / k_B * m_e).to(u.K)

print('T from standard deviation:', T_theo)
print('Initial T:', T)
