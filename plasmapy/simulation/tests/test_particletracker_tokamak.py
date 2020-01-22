#!/usr/bin/env python
# coding: utf-8

import astropy.units as u
import numpy as np
import pytest

from plasmapy import simulation
from plasmapy.formulary import magnetostatics
from plasmapy.classes.sources import Coils

MINOR_RADIUS = 0.3 * u.m
RADIUS = 1 * u.m
MAIN_CURRENT = 15 * u.MA

COIL_CURRENT = 10 * u.MA

@pytest.fixture
def coils():
    n_coils = 1 # 8
    currents = u.Quantity(n_coils * [COIL_CURRENT])

    coil_angles = np.linspace(0, 2*np.pi, n_coils, endpoint=False)

    coils = []
    for i in range(n_coils):
        coil_angle = coil_angles[i]
        x = RADIUS * np.cos(coil_angle)
        y = RADIUS * np.sin(coil_angle)
        normal_angle = np.pi/2 + coil_angle
        normal = u.Quantity([np.cos(normal_angle), np.sin(normal_angle), 0])
        center = u.Quantity([x, y, 0 * u.m])
        coil = magnetostatics.CircularWire(normal, center, MINOR_RADIUS, currents[i])
        coils.append(coil)

    # plasma_wire = magnetostatics.CircularWire([0, 0, 1],
    #                                           u.Quantity((0, 0, 0), u.m),
    #                                           RADIUS, MAIN_CURRENT)
    # coils.append(plasma_wire)

    c = Coils(*coils)
    return c


@pytest.fixture
def sim1(coils):
    x = u.Quantity([[1 + MINOR_RADIUS.si.value / 2, 0, 0]],  u.m)
    v = u.Quantity([[0, 100, 10]], u.m / u.s)

    sim = simulation.ParticleTracker(coils, x, v, 'e-')
    return sim

def test_1(sim1):
    solution = sim1.run(1e-3 * u.s, 1e2)
    # c = sim1.plasma

    # from mayavi import mlab
    # fig = mlab.figure()
    # c.visualize(fig)
    # solution.visualize(fig)
    # mlab.orientation_axes(figure=fig)
    # mlab.show()

    assert abs(np.mean(solution.position_history[:,0,1])) < 4 * u.m  # should be about 5m for no B field
    assert 0.001 * u.m < abs(np.mean(solution.position_history[:,0,1])) < 0.01 * u.m
