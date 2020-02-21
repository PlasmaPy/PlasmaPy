#!/usr/bin/env python
# coding: utf-8

# # Particle movement in a simple model of a tokamak

# In[1]:


from plasmapy import simulation
from plasmapy.formulary import magnetostatics
from plasmapy.classes.sources import Coils
import astropy.units as u
import numpy as np

radius = 0.5 * u.m
main_current = 4 * u.kA
plasma_wire = magnetostatics.CircularWire(
    [0, 0, 1], u.Quantity((0, 0, 0), u.m), radius, main_current
)

# That's supposed to model just the plasma; let's add a few coils:

# In[2]:


n_coils = 8
currents = n_coils * [4 * u.kA]

coil_angles = np.linspace(0, 2 * np.pi, n_coils, endpoint=False)
coil_angles
minor_radius = 0.3 * u.m

coils = []
for i in range(n_coils):
    coil_angle = coil_angles[i]
    x = radius * np.cos(coil_angle)
    y = radius * np.sin(coil_angle)
    normal_angle = np.pi / 2 + coil_angle
    normal = u.Quantity([np.cos(normal_angle), np.sin(normal_angle), 0])
    center = u.Quantity([x, y, 0 * u.m])
    coil = magnetostatics.CircularWire(normal, center, minor_radius, currents[i])
    coils.append(coil)
all_currents = coils + [plasma_wire]
c = Coils(*all_currents)
# Note that a shortcut for this model can be quickly accessed as `Coils.toykamak`.

# In[3]:


x = u.Quantity([[0.6, 0, 0]], u.m)
v = u.Quantity([[0, 300, 0]], u.m / u.s)

solution = simulation.ParticleTracker(c, x, v, "p").run(1e-1 * u.s, dt=1e-5 * u.s)
solution


# In[4]:


solution.particle.sel(particle=list((0,)))


# In[5]:


solution.particletracker.animate("toykamak.mp4", nframes=1000, plasma=c)


# In[6]:


from IPython.display import Video

Video("toykamak.mp4")
