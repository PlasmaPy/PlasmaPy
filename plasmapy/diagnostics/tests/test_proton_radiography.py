"""
Tests for proton radiography functions
"""

import astropy.units as u
import numpy as np
import pytest
import warnings

from plasmapy.diagnostics import proton_radiography as prad
from plasmapy.plasma import fields as fields


def test_coordinate_systems():
    """
    Check that specifying the same point in different coordinate systems
    ends up with identical source and detector vectors.

    """

    grid, E, B = fields.example_fields(model="no fields")

    # Cartesian
    source = (-7.07 * u.mm, -7.07 * u.mm, 0 * u.mm)
    detector = (70.07 * u.mm, 70.07 * u.mm, 0 * u.mm)
    sim1 = prad.SyntheticProtonRadiograph(
        grid, E, B, source, detector, geometry="cartesian", verbose=False
    )

    # Cylindrical
    source = (-1 * u.cm, 45 * u.deg, 0 * u.mm)
    detector = (10 * u.cm, 45 * u.deg, 0 * u.mm)
    sim2 = prad.SyntheticProtonRadiograph(
        grid, E, B, source, detector, geometry="cylindrical", verbose=False
    )

    # In spherical
    source = (-0.01 * u.m, 90 * u.deg, 45 * u.deg)
    detector = (0.1 * u.m, 90 * u.deg, 45 * u.deg)
    sim3 = prad.SyntheticProtonRadiograph(
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
    grid, E, B = fields.example_fields(
        model="electrostatic gaussian sphere",
        regular_grid=True,
        length=np.array([1, 1, 1]) * u.mm,
        num=(100, 100, 100),
    )

    source = (-10 * u.mm, 90 * u.deg, 45 * u.deg)
    detector = (100 * u.mm, 90 * u.deg, 45 * u.deg)
    sim = prad.SyntheticProtonRadiograph(
        grid, E, B, source, detector, geometry="spherical", verbose=True
    )

    sim.run(1e3, max_theta=np.pi / 12 * u.rad, field_weighting="volume averaged")
    hax, vax, values = sim.synthetic_radiograph()

    # Make a radiograph
    size = np.array([[-1, 1], [-1, 1]]) * 5e-2 * u.m
    bins = [200, 200]
    hax, vax, values = sim.synthetic_radiograph(bins=bins, size=size)

    # Make an OD radiograph
    hax, vax, values = sim.synthetic_radiograph(
        bins=bins, size=size, optical_density=True
    )

    sim.calc_ke()


def test_irregular_grid():
    """
    Run a simulation with an irregular grid
    """
    grid, E, B = fields.example_fields(
        model="axial magnetic field", num=50, regular_grid=False
    )

    source = (-10 * u.mm, 90 * u.deg, 45 * u.deg)
    detector = (100 * u.mm, 90 * u.deg, 45 * u.deg)
    sim = prad.SyntheticProtonRadiograph(
        grid, E, B, source, detector, geometry="spherical", verbose=False
    )

    sim.run(1e2, max_theta=np.pi / 12 * u.rad)
    hax, vax, values = sim.synthetic_radiograph()

    # Check that trying to run this simulation with volume-averaged fields
    # raises an exception
    with pytest.raises(ValueError):
        sim.run(1e1, max_theta=np.pi / 12 * u.rad, field_weighting="volume averaged")


def test_SyntheticProtonRadiograph_error_handling():
    """
    Intentionally raise a number of errors.
    """

    # INIT ERRORS

    grid, E, B = fields.example_fields(
        model="electrostatic gaussian sphere", num=(100, 100, 100)
    )
    source = (-10 * u.mm, 90 * u.deg, 45 * u.deg)
    detector = (100 * u.mm, 90 * u.deg, 45 * u.deg)

    # Check that an error is raised when an input grid has a nan or infty value
    E[0, 0, 0, :] = np.nan
    with pytest.raises(ValueError):
        sim = prad.SyntheticProtonRadiograph(
            grid, E, B, source, detector, geometry="spherical", verbose=False
        )
    E[0, 0, 0] = 0 * u.V / u.m  # Reset element for the rest of the tests

    E[0, 0, 0, :] = np.infty
    with pytest.raises(ValueError):
        sim = prad.SyntheticProtonRadiograph(
            grid, E, B, source, detector, geometry="spherical", verbose=False
        )
    E[0, 0, 0] = 0 * u.V / u.m  # Reset element for the rest of the tests

    # Set an edge field value to a large value
    E[0, 0, 0, :] = 0.2 * np.max(E)
    with pytest.warns(RuntimeWarning):
        sim = prad.SyntheticProtonRadiograph(
            grid, E, B, source, detector, geometry="spherical", verbose=False
        )
    E[0, 0, 0] = 0 * u.V / u.m  # Reset element for the rest of the tests

    # Raise error when source-to-detector vector doesn't pass through the
    # field grid
    source_bad = (10 * u.mm, -10 * u.mm, 0 * u.mm)
    detector_bad = (10 * u.mm, 100 * u.mm, 0 * u.mm)
    with pytest.raises(ValueError):
        sim = prad.SyntheticProtonRadiograph(
            grid, E, B, source_bad, detector_bad, geometry="cartesian", verbose=False
        )

    # RUNTIME ERRORS
    sim = prad.SyntheticProtonRadiograph(
        grid, E, B, source, detector, geometry="spherical", verbose=False
    )

    # Chose too large of a max_theta so that many particles miss the grid
    with pytest.warns(RuntimeWarning):
        sim.run(1e4, max_theta=0.99 * np.pi / 2 * u.rad)

    # SYNTHETIC RADIOGRAPH ERRORS
    sim.run(1e4, max_theta=np.pi / 10 * u.rad)

    # Choose a very small synthetic radiograph size that misses most of the
    # particles
    with pytest.warns(RuntimeWarning):
        size = np.array([[-1, 1], [-1, 1]]) * 1 * u.mm

        hax, vax, values = sim.synthetic_radiograph(size=size)
