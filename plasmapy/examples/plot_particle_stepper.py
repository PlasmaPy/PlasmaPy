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

############################################################
# Initialize a plasma. This will be a source of electric and magnetic
# fields for our particles to move in.

plasma = Plasma(domain_x=np.linspace(-1, 1, 10) * u.m,
                domain_y=np.linspace(-1, 1, 10) * u.m,
                domain_z=np.linspace(-1, 1, 10) * u.m)

############################################################
# Initialize the fields. We'll take $\vec{B}$ in the $\hat{x}$ direction
# and $E$ in the $\hat{y}$ direction, which gets us an $E \times B$ drift
# in $\hat{z}$.

B0 = 4 * u.T
plasma.magnetic_field[0, :, :, :] = np.ones((10, 10, 10)) * B0

E0 = 2 * u.V / u.m
plasma.electric_field[1, :, :, :] = np.ones((10, 10, 10)) * E0

############################################################
# Initialize the particle. We'll take one proton `p` with a timestep of
# $10^{-4}s$ and run it for 40 iterations.

particle = ParticleTracker(plasma, 'p', 1, 1, 1e-4 * u.s, 40)

############################################################
# We still have to initialize the particle's velocity. We'll limit ourselves to
# one in the $\hat{x}$ direction, parallel to the magnetic field $\vec{B}$ -
# that way, it won't turn in the $\hat{z}$ direction.

particle.v[0][0] = 1 * (u.m / u.s)

############################################################
# Run the pusher and plot the trajectory versus time.

particle.run()
particle.plot_time_trajectories()

############################################################
# Plot the shape of the trajectory in 3D.

particle.plot_trajectories()

############################################################
# As a test, we calculate the mean velocity in the z direction from the
# velocity and position

vmean = particle.velocity_history[:, :, 2].mean()
print(f"The calculated drift velocity is {vmean:.4f} to compare with the"
      f"theoretical E0/B0 = {E0/B0:.4f}")

############################################################
# and from position:
Vdrift = particle.position_history[-1, 0, 2] / (particle.NT * particle.dt)
normdrift = Vdrift
print(f"The calculated drift velocity from position is {normdrift:.4f}")
