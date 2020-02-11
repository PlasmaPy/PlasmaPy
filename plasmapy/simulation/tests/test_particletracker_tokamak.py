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
MAIN_CURRENT = 0 * u.MA

COIL_CURRENTS = 1 * [10 * u.MA]


@pytest.fixture
def coils():
    return Coils.toykamak(MINOR_RADIUS, RADIUS, MAIN_CURRENT, COIL_CURRENTS)


@pytest.fixture
def sim_single(coils):
    x = u.Quantity([[1 + MINOR_RADIUS.si.value / 2, 0, 0]], u.m)
    v = u.Quantity([[0, 100, 10]], u.m / u.s)

    sim = simulation.ParticleTracker(coils, x, v, "e-")
    return sim


@pytest.fixture
def sim_many(coils):
    N = 100
    x = u.Quantity(N * [[1 + MINOR_RADIUS.si.value / 2, 0, 0]], u.m)
    s = np.random.RandomState(0)
    v = u.Quantity(s.normal(size=(N, 3)), u.m / u.s)

    sim = simulation.ParticleTracker(coils, x, v, "e-")
    return sim


def test_1(sim_single, integrator_name):
    solution = sim_single.run(1 * u.s, 1e-3 * u.s, pusher=integrator_name)
    assert (
        abs(solution.data.position.sel(dimension="y").mean()).item() < 4
    )  # should be about 5m for no B field
    assert 0.001 < abs(solution.data.position.sel(dimension="y").mean()).item() < 0.01


def test_2(sim_many, integrator_name):
    solution = sim_many.run(1 * u.s, 1e-3 * u.s, pusher=integrator_name)
    assert (
        abs(solution.data.position.sel(dimension="y").mean()).item() < 4
    )  # should be about 5m for no B field
    assert 0.001 < abs(solution.data.position.sel(dimension="y").mean()).item() < 0.1
