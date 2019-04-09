"""
Analysing ITER parameters
=========================

Let's try to look at ITER plasma conditions using the `parameters` subpackage.
"""

from astropy import units as u
from plasmapy.theory import parameters
import matplotlib.pyplot as plt
import numpy as np

######################################################
# The radius of electric field shielding clouds, also known as the Debye length,
# would be

electron_temperature = 8.8 * u.keV
electron_concentration = 10.1e19 / u.m**3
print(parameters.Debye_length(electron_temperature, electron_concentration))

############################################################
# Note that we can also neglect the unit for the concentration, as
# 1/m^3 is the a standard unit for this kind of Quantity:

print(parameters.Debye_length(electron_temperature, 10.1e19))

############################################################
# Assuming the magnetic field as 5.3 Teslas (which is the value at the major
# radius):

B = 5.3 * u.T

print(parameters.gyrofrequency(B, particle='e'))

print(parameters.gyroradius(B, T_i=electron_temperature, particle='e'))

######################################################################
# The electron inertial length would be
print(parameters.inertial_length(electron_concentration, particle='e'))

######################################################################
# In these conditions, they should reach thermal velocities of about
print(parameters.thermal_speed(T=electron_temperature, particle='e'))

######################################################################
# And the Langmuir wave plasma frequency should be on the order of
print(parameters.plasma_frequency(electron_concentration))

############################################################
# Let's try to recreate some plots and get a feel for some of these quantities.

n_e = np.logspace(4, 30, 100) / u.m**3
plt.plot(n_e, parameters.plasma_frequency(n_e))
plt.scatter(
    electron_concentration,
    parameters.plasma_frequency(electron_concentration))
plt.xlabel("Electron Concentration (m^-3)")
plt.ylabel("Langmuir Wave Plasma Frequency (rad/s)")
plt.show()
