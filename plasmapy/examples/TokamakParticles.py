#!/usr/bin/env python
# coding: utf-8

# In[1]:


from plasmapy import simulation

from plasmapy.formulary import magnetostatics
from plasmapy.classes.sources import Coils

import astropy.units as u
import numpy as np
radius = 1 * u.m
main_current = 15 * u.MA
plasma_wire = magnetostatics.CircularWire([0, 0, 1], u.Quantity((0, 0, 0), u.m), radius, main_current)
plasma_wire


# In[2]:


n_coils = 8
coil_angles = np.linspace(0, 2*np.pi, n_coils, endpoint=False)
coil_angles


# In[3]:


minor_radius = 0.3 * u.m
currents = u.Quantity(n_coils * [10e6], u.A)
currents


# In[4]:


coils = []
for i in range(n_coils):
    coil_angle = coil_angles[i]
    x = radius * np.cos(coil_angle)
    y = radius * np.sin(coil_angle)
    normal_angle = np.pi/2 + coil_angle
    normal = u.Quantity([np.cos(normal_angle), np.sin(normal_angle), 0])
    center = u.Quantity([x, y, 0 * u.m])
    coil = magnetostatics.CircularWire(normal, center, minor_radius, currents[i])
    coils.append(coil)


# In[5]:


all_currents = coils + [plasma_wire]


# # TODO units in reprs for MagnetoStatics
# 
# # TODO MagnetoStatics.magnetic_field accepts no units, just numpy arrays!
# 
# # TODO mendeleev vs atomic - comparison

# In[6]:


plasma_wire.magnetic_field(center.value)


# In[ ]:


c = Coils(*all_currents)

sim = simulation.ParticleTracker(c, 'e', dt=1e-8 * u.s, nt=int(1e6))
sim._x[0][0] = 1 + minor_radius.si.value / 2 # * (u.m / u.s)
sim._v[0][1] = 10000 # * (u.m / u.s)
sim._v[0][2] = 100 # * (u.m / u.s)


sim.run() # this should return a Solution object or sth


# In[34]:


import mayavi
from mayavi import mlab
mlab.init_notebook()


# In[35]:


fig = mlab.figure()
c.visualize(fig)
sim.visualize(fig)
mlab.orientation_axes(figure=fig)
display(fig)


# In[36]:


sim2 = simulation.ParticleTracker(c, 'e', dt=1e-8 * u.s, nt=int(1e6))
sim2._x[0][0] = 1 + minor_radius.si.value / 2 # * (u.m / u.s)
sim2._v[0][0] = 0 # * (u.m / u.s)
sim2._v[0][1] = 10000 # * (u.m / u.s)
sim2._v[0][2] = 0 # * (u.m / u.s)
sim2.run() # this should return a Solution object or sth


# In[38]:


fig2 = mlab.figure()
c.visualize(fig2)
sim2.visualize(fig2)
mlab.orientation_axes(figure=fig2)
display(fig2)


# In[68]:


factor = 1.7
x = np.linspace(-factor * (radius + minor_radius), factor * (radius + minor_radius), 100)
z = np.linspace(-2*factor * minor_radius, 2*factor * minor_radius, 100)
X, Y, Z = np.meshgrid(x, x, z, indexing='ij')

Bval = np.zeros((3, *Z.shape))

from tqdm import auto as tqdm

for i, xi in tqdm.tqdm(enumerate(x.si.value), total=len(x)):
    for j, yi in enumerate(x.si.value):
        for k, zi in enumerate(z.si.value):
            pos = np.array([[xi, yi, zi]])
            field = c._interpolate_B(pos)
            Bval[:, i,j,k] = field

Bmag2 = np.sum(Bval**2, axis=0)


# In[70]:


fig = mlab.figure(size=(800, 600))
c.visualize(fig)

contours = mlab.contour3d(X.value, Y.value, Z.value, np.log10(Bmag2), figure=fig,
#                           contours = np.linspace(Bmag2.min(), Bmag2.max(), 50).tolist(),
#                           contours = np.logspace(np.log10(Bmag2.min()),
#                                                  np.log10(Bmag2.max()),
#                                                  num=10,
#                                                 ).tolist(),
                          contours = 15,
                          opacity = 0.3,
                         )
mlab.colorbar(contours, title="log10 |B|^2")
mlab.orientation_axes(figure=fig)
sim.visualize(fig2)
display(fig)

