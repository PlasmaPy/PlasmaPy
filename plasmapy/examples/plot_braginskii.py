"""
Braginskii coefficients
=========================

A short example of how to calculate classical transport coefficients
from Bragiński's theory.
"""

from astropy import units as u
from plasmapy.theory.transport import ClassicalTransport

#####################################################
# We'll use some sample ITER data, without much regard for whether
# the regime is even fit for classical transport theory:

thermal_energy_per_electron = 8.8 * u.keV
electron_concentration = 10.1e19 / u.m**3

thermal_energy_per_ion = 8.0 * u.keV
ion_concentration = electron_concentration
ion_particle = 'D+'  # a crude approximation

######################################################
# We now make the default ClassicalTransport object:

braginskii = ClassicalTransport(thermal_energy_per_electron,
                                electron_concentration,
                                thermal_energy_per_ion,
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
# call them via ().

print(braginskii.resistivity())
print(braginskii.thermoelectric_conductivity())
print(braginskii.electron_thermal_conductivity())
print(braginskii.ion_thermal_conductivity())

######################################################
# They also change with magnetization:

mag_braginskii = ClassicalTransport(thermal_energy_per_electron,
                                    electron_concentration,
                                    thermal_energy_per_ion,
                                    ion_concentration,
                                    ion_particle,
                                    B = 0.1 * u.T)

print(mag_braginskii.resistivity())
print(mag_braginskii.thermoelectric_conductivity())
print(mag_braginskii.electron_thermal_conductivity())
print(mag_braginskii.ion_thermal_conductivity())

######################################################
# They also change with direction with respect to the magnetic field. Here,
# we choose to print out, as arrays, the (parallel, perpendicular,
# and cross) directions. Take a look at the docs to `ClassicalTransport`
# for more information on these.

all_direction_braginskii = ClassicalTransport(thermal_energy_per_electron,
                                    electron_concentration,
                                    thermal_energy_per_ion,
                                    ion_concentration,
                                    ion_particle,
                                    B = 0.1 * u.T,
                                    field_orientation = 'all')

print(all_direction_braginskii.resistivity())
print(all_direction_braginskii.thermoelectric_conductivity())
print(all_direction_braginskii.electron_thermal_conductivity())
print(all_direction_braginskii.ion_thermal_conductivity())

######################################################
# The viscosities return arrays:

print(braginskii.electron_viscosity())
print(mag_braginskii.electron_viscosity())
print(braginskii.ion_viscosity())
print(mag_braginskii.ion_viscosity())
