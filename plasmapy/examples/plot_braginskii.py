# coding: utf-8
"""
Braginskii coefficients
=========================

A short example of how to calculate classical transport coefficients
from Bragiński's theory.
"""

from astropy import units as u
from plasmapy import physics
import matplotlib.pyplot as plt
import numpy as np

#####################################################
# We'll use some sample ITER data, without much regard for whether
# the regime is even fit for classical transport theory:

electron_temperature = 8.8 * u.keV
electron_concentration = 10.1e19 / u.m**3

ion_temperature = 8.0 * u.keV
ion_concentration = electron_concentration
ion_particle = 'D+'  # a crude approximation

######################################################
# We now make the default classical_transport object:

braginskii = physics.transport.classical_transport(electron_temperature,
                                                   electron_concentration,
                                                   ion_temperature,
                                                   ion_concentration,
                                                   ion_particle)

######################################################
# These variables are calculated during initialization and can be
# referred to straight away:

print(braginskii.coulomb_log_ei)
print(braginskii.coulomb_log_ii)
print(braginskii.hall_e)
print(braginskii.hall_i)

######################################################
# These quantities are not calculated during initialization and can be
# referred to via methods. To signify the need to calculate them, we
# call them via (). They could be made to act like variables via @Property,
# but I'm not sure that's wise.

print(braginskii.resistivity())
print(braginskii.thermoelectric_conductivity())
