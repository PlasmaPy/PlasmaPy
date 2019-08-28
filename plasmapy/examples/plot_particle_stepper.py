"""
Particle stepper
================

An example of PlasmaPy's particle stepper class, currently in need of a rewrite
for speed.
"""


import numpy as np
from astropy import units as u
from plasmapy.classes import Plasma
from plasmapy.simulation import ParticleTracker
from plasmapy.physics.parameters import gyrofrequency

############################################################
# Initialize a plasma. This will be a source of electric and magnetic
# fields for our particles to move in.

plasma = Plasma(domain_x=np.linspace(-1, 1, 10) * u.m,
                domain_y=np.linspace(-1, 1, 10) * u.m,
                domain_z=np.linspace(-1, 1, 10) * u.m)

############################################################
# Initialize the fields. We'll take B in the x direction
# and E in the y direction, which gets us an E cross B drift
# in the z direction.

B0 = 4 * u.T
plasma.magnetic_field[0, :, :, :] = np.ones((10, 10, 10)) * B0

E0 = 2 * u.V / u.m
plasma.electric_field[1, :, :, :] = np.ones((10, 10, 10)) * E0

############################################################
# Calculate the timestep. We'll take one proton `p`, take its gyrofrequency, invert that
# to get to the gyroperiod, and resolve that into 10 steps for higher accuracy.

freq = gyrofrequency(B0, 'p').to(u.Hz, equivalencies = u.dimensionless_angles())
gyroperiod = (1/freq).to(u.s)
steps_to_gyroperiod = 10
timestep = gyroperiod / steps_to_gyroperiod

############################################################
# Initialize the trajectory calculation.

number_steps = steps_to_gyroperiod * int(2 * np.pi)
particle_trajectory = ParticleTracker(plasma, 'p', 1, 1, timestep, number_steps)

############################################################
# We still have to initialize the particle's velocity. We'll limit ourselves to
# one in the x direction, parallel to the magnetic field B -
# that way, it won't turn in the z direction.

particle_trajectory.v[0][0] = 1 * (u.m / u.s)

############################################################
# Run the pusher and plot the trajectory versus time.

particle_trajectory.run()
particle_trajectory.plot_time_trajectories()

############################################################
# Plot the shape of the trajectory in 3D.

particle_trajectory.plot_trajectories()

############################################################
# As a test, we calculate the mean velocity in the z direction from the
# velocity and position

vmean = particle_trajectory.velocity_history[:, :, 2].mean()
print(f"The calculated drift velocity is {vmean:.4f} to compare with the"
      f"theoretical E0/B0 = {E0/B0:.4f}")

############################################################
# and from position:
Vdrift = particle_trajectory.position_history[-1, 0, 2] / (particle_trajectory.NT * particle_trajectory.dt)
print(f"The calculated drift velocity from position is {Vdrift:.4f}")
