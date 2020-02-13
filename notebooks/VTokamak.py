#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

import astropy.units as u
import numpy as np
import pytest

from plasmapy import simulation
from plasmapy.formulary import magnetostatics
from plasmapy.classes.sources import Coils
import matplotlib.pyplot as plt
import pyvista

MINOR_RADIUS = 0.3 * u.m
RADIUS = 0.5 * u.m
MAIN_CURRENT = 4 * u.kA  # TOROIDAL

COIL_CURRENTS = 8 * [4 * u.kA]


coils = Coils.toykamak(MINOR_RADIUS, RADIUS, MAIN_CURRENT, COIL_CURRENTS)


x = u.Quantity([[0.6, 0, 0]], u.m)
v = u.Quantity([[0, 300, 0]], u.m / u.s)

sim_single = simulation.ParticleTracker(coils, x, v, "p")

for push in ["explicit_boris"]:
    print(push)
    solution = sim_single.run(1e-3 * u.s, 1e-5 * u.s, pusher=push)
    fig = pyvista.Plotter()
    solution.visualize(fig)
    fig.show()
