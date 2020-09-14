"""
Tests for proton radiography functions
"""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.diagnostics import proton_radiography as prad


def test_coordinate_systems():
    """
    Check that specifying the same point in different coordinate systems
    ends up with identical source and detector vectors.

    """

    grid, E, B = prad.test_fields(mode="no fields")

    # Cartesian
    source = (-7.07 * u.mm, -7.07 * u.mm, 0 * u.mm)
    detector = (70.07 * u.mm, 70.07 * u.mm, 0 * u.mm)
    sim1 = prad.SimPrad(
        grid, E, B, source, detector, geometry="cartesian", verbose=False
    )

    # Cylindrical
    source = (-1 * u.cm, 45 * u.deg, 0 * u.mm)
    detector = (10 * u.cm, 45 * u.deg, 0 * u.mm)
    sim2 = prad.SimPrad(
        grid, E, B, source, detector, geometry="cylindrical", verbose=False
    )

    # In spherical
    source = (-0.01 * u.m, 90 * u.deg, 45 * u.deg)
    detector = (0.1 * u.m, 90 * u.deg, 45 * u.deg)
    sim3 = prad.SimPrad(
        grid, E, B, source, detector, geometry="spherical", verbose=False
    )

    assert np.allclose(sim1.source, sim2.source, atol=1e-2)
    assert np.allclose(sim2.source, sim3.source, atol=1e-2)
    assert np.allclose(sim1.detector, sim2.detector, atol=1e-2)
    assert np.allclose(sim2.detector, sim3.detector, atol=1e-2)


def test_regular_grid():
    """
    Run a simulation with a regular grid
    """
    grid, E, B = prad.test_fields(
        mode="electrostatic gaussian sphere",
        regular_grid=True,
        length=np.array([1, 1, 1]) * u.mm,
        num=(100, 100, 100),
    )

    source = (-10 * u.mm, 90 * u.deg, 45 * u.deg)
    detector = (100 * u.mm, 90 * u.deg, 45 * u.deg)
    sim = prad.SimPrad(grid, E, B, source, detector, geometry="spherical", verbose=True)

    sim.run(1e3, max_theta=np.pi / 6 * u.rad)
    hax, vax, values = sim.synthetic_radiograph()

    size = np.array([[-1, 1], [-1, 1]]) * 5e-2 * u.m
    bins = [200, 200]
    hax, vax, values = sim.synthetic_radiograph(bins=bins, size=size)

    sim.calc_ke()


def test_irregular_grid():
    """
    Run a simulation with an irregular grid
    """
    grid, E, B = prad.test_fields(mode="axial magnetic field", regular_grid=False)

    source = (-10 * u.mm, 90 * u.deg, 45 * u.deg)
    detector = (100 * u.mm, 90 * u.deg, 45 * u.deg)
    sim = prad.SimPrad(
        grid, E, B, source, detector, geometry="spherical", verbose=False
    )

    sim.run(1e3, max_theta=np.pi / 6 * u.rad)
    hax, vax, values = sim.synthetic_radiograph()


def test_other_test_fields():
    """
    Creates all test fields that aren't run in other tests'
    """

    grid, E, B = prad.test_fields(
        mode="electrostatic planar shock", num=(100, 100, 200)
    )
