#!/usr/bin/env python
# coding: utf-8

# In[1]:


from plasmapy import simulation
from plasmapy.formulary import magnetostatics
from plasmapy.classes.sources import Coils

import astropy.units as u
import numpy as np

radius = 1 * u.m
main_current = 1 * u.A
plasma_wire = magnetostatics.CircularWire(
    [0, 0, 1], u.Quantity((0, 0, 0), u.m), radius, main_current
)


n_coils = 4
coil_angles = np.linspace(0, 2 * np.pi, n_coils, endpoint=False)


minor_radius = 0.3 * u.m
currents = u.Quantity(n_coils * [0.1], u.A)

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
x = u.Quantity([[0, 0, minor_radius.si.value / 2]], u.m)
v = u.Quantity([[0, 1000, 100]], u.m / u.s)
sim = simulation.ParticleTracker(c, x, v, "e")


def profile():
    sim.run(dt=1e-8 * u.s, nt=1e5)


profile()
