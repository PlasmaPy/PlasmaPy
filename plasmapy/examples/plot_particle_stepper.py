"""
Particle stepper
================

An example of using PlasmaPy's particle stepper class.
"""


import numpy as np
from astropy import units as u
from plasmapy.classes.sources.analyticalplasma import AnalyticalPlasma
from plasmapy.simulation import ParticleTracker
from plasmapy.formulary import gyrofrequency

############################################################
# Initialize a plasma. This will be a source of electric and magnetic
# fields for our particles to move in.
# We'll take the magnetic field in the x direction
# and the electric field in the y direction, which gets us an E cross B drift
# in the z direction.


def magnetic_field(r):
    # TODO proper shape of array needs to be asserted
    return u.Quantity([4, 0, 0], u.T)

# precomputed for efficiency
E_unit = u.V / u.m
def electric_field(r):
    return u.Quantity([0, 2, 0], E_unit)

plasma = AnalyticalPlasma(magnetic_field, electric_field)

############################################################
# We'll now calculate the timestep. We'll take one proton `p`,
# take its gyrofrequency, invert that
# to get to the gyroperiod, and resolve that into 10 steps for higher accuracy.

freq = gyrofrequency(4 * u.T, 'p').to(u.Hz, equivalencies=u.dimensionless_angles())
gyroperiod = (1/freq).to(u.s)
steps_to_gyroperiod = 10
timestep = gyroperiod / steps_to_gyroperiod

############################################################
# Initialize the trajectory calculation.
# We still have to initialize the particle's velocity. We'll limit ourselves to
# one in the x direction, parallel to the magnetic field B -
# that way, it won't turn in the z direction.

number_steps = 5 * steps_to_gyroperiod * int(2 * np.pi)
trajectory = ParticleTracker(plasma, v = u.Quantity([[-1, 0, 0]] * u.m/u.s), particle_type = 'p')

############################################################
# Let's run the pusher and plot the trajectory versus time.
# We'll just show the y-z trajectories for clarity.

solution = trajectory.run(timestep, number_steps)
solution.plot_time_trajectories('yz')

############################################################
# Plot the shape of the trajectory in 3D.

solution.plot_trajectories()

############################################################
# If you have Mayavi, you can 

solution.visualize()


############################################################
# As a test, we calculate the mean velocity in the z direction from the
# velocity and position

vmean = solution.velocity_history[:, :, 2].mean()
print(f"The calculated drift velocity is {vmean:.4f} to compare with the "
      f"theoretical E0/B0 = {0.5 * u.m / u.s}")

############################################################
# and from position:
Vdrift = solution.position_history[-1, 0, 2] / solution.t.max()
print(f"The calculated drift velocity from position is {Vdrift:.4f}")

############################################################
# Supposing we wanted to examine the effect of the initial velocity in the x-y plane on the trajectory:
N = 20
np.random.seed(0)
v = np.zeros((N, 3))
v[:, :2] = np.random.normal(size=(N, 2))
trajectory = ParticleTracker(plasma, v = v * u.m / u.s, particle_type = 'p')
# we choose this as our example's thumbnail:
# sphinx_gallery_thumbnail_number = 3
solution = trajectory.run(timestep/100, number_steps*200)
solution.plot_trajectories(alpha=0.8)

############################################################
# Note how while each trajectory fans out in a different way,
# each one traverses the z direction in about the same time:

solution.plot_time_trajectories('z')
