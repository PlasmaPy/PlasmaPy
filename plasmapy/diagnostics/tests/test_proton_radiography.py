"""
Tests for proton radiography functions
"""

import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest
import warnings

from scipy.special import erf

from plasmapy.diagnostics import proton_radiography as prad
from plasmapy.plasma.grids import CartesianGrid

import matplotlib.pyplot as plt


def _test_grid(name, L=1 * u.mm, num=100, B0=10 * u.T):
    r"""
    Generates grids representing some common physical scenarios for testing
    and illustration. Valid example names are:
    * axially_magnetized_cylinder : A cylinder of radius L/4 magnetized in the
        Z-direction (like a solenoid, but without the fringe fields).
    * electrostatic_discontinuity : A discontinuity in the electric field at z=0
        with a radial gaussian profile in the xy plane.
    * electrostatic_gaussian_sphere : An electric field created by a sphere
        of potential of radius L/2 with a radial Gaussian distribution.
    Parameters
    ----------
    name : str
        Name of example to load (from list above)
    L : `~u.Quantity` (or array of three of the same)
        Length scale (or scales). -L and L are passed to the grid constructor
        as start and stop respectively. The default is 1 cm.
    num : int or list of three ints
        The number of points in each direction (or list of one for each dimension).
        Passed to the grid cosntructor as the num argument. The default is 100.
    Returns
    -------
    grid : CartesianGrid
        A CartesianGrid object containing quantity arrays representing
        the chosen example.
    """

    grid = CartesianGrid(-L, L, num=num)

    # If an array was provided to the constructor, reduce to a single
    # length scale now.
    if L.size > 1:
        L = np.max(L)

    if name == "empty":
        pass

    elif name == "constant_bz":
        Bz = np.ones(grid.shape) * B0
        grid.add_quantities(B_z=Bz)

    elif name == "constant_ex":
        Ex = np.ones(grid.shape) * 5e8 * u.V / u.m
        grid.add_quantities(E_x=Ex)

    elif name == "axially_magnetized_cylinder":
        a = L / 4
        radius = np.linalg.norm(grid.grid[..., 0:2] * grid.unit, axis=3)
        Bz = np.where(radius < a, 100 * u.T, 0 * u.T)
        grid.add_quantities(B_z=Bz)

    elif name == "electrostatic_discontinuity":
        a = L / 2
        delta = a / 120

        radius = np.linalg.norm(grid.grid[..., 0:2] * grid.unit, axis=3)
        z = grid.grids[2]

        potential = (1 - erf(z / delta)) * np.exp(-((radius / a) ** 2)) * u.V

        Ex, Ey, Ez = np.gradient(potential, grid.dax0, grid.dax1, grid.dax2)

        grid.add_quantities(E_x=Ex, E_y=Ey, E_z=Ez, phi=potential)

    elif name == "electrostatic_gaussian_sphere":
        a = L / 3
        b = L / 2
        radius = np.linalg.norm(grid.grid, axis=3)
        arg = (radius / a).to(u.dimensionless_unscaled)
        potential = np.exp(-(arg ** 2)) * u.V

        Ex, Ey, Ez = np.gradient(potential, grid.dax0, grid.dax1, grid.dax2)

        Ex = -np.where(radius < b, Ex, 0)
        Ey = -np.where(radius < b, Ey, 0)
        Ez = -np.where(radius < b, Ez, 0)

        grid.add_quantities(E_x=Ex, E_y=Ey, E_z=Ez, phi=potential)

    else:
        raise ValueError(
            "No example corresponding to the provided name " f"({name}) exists."
        )

    # If any of the following quantities are missing, add them as empty arrays
    req_quantities = ["E_x", "E_y", "E_z", "B_x", "B_y", "B_z"]

    for q in req_quantities:
        if q not in list(grid.ds.data_vars):
            unit = grid.recognized_quantities[q].unit
            arg = {q: np.zeros(grid.shape) * unit}
            grid.add_quantities(**arg)

    return grid


def run_1D_example(name):
    """
    Run a simulation through an example with parameters optimized to
    sum up to a lineout along x. The goal is to run a realtively fast
    sim with a quasi-1D field grid that can then be summed to get good
    enough statistics to use as a test.
    """
    grid = _test_grid(name, L=1 * u.mm, num=50)

    # Cartesian
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)

    # TODO: Catch warnings (test fields aren't well behaved at edges)

    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    sim.create_particles(1e4, 3 * u.MeV, max_theta=0.1 * u.deg)
    sim.run()

    size = np.array([[-1, 1], [-1, 1]]) * 10 * u.cm
    bins = [200, 60]
    hax, vax, values = sim.synthetic_radiograph(size=size, bins=bins)

    values = np.mean(values[:, 20:40], axis=1)

    return hax, values


def test_1D_deflections():
    # Check B-deflection
    hax, lineout = run_1D_example("constant_bz")
    loc = hax[np.argmax(lineout)]
    assert np.isclose(loc.si.value, 0.0165, 0.005)

    # Check E-deflection
    hax, lineout = run_1D_example("constant_ex")
    loc = hax[np.argmax(lineout)]
    assert np.isclose(loc.si.value, 0.0335, 0.005)


def test_coordinate_systems():
    """
    Check that specifying the same point in different coordinate systems
    ends up with identical source and detector vectors.
    """

    grid = _test_grid("empty")

    # Cartesian
    source = (-7.07 * u.mm, -7.07 * u.mm, 0 * u.mm)
    detector = (70.07 * u.mm, 70.07 * u.mm, 0 * u.mm)
    sim1 = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=True)

    # Cylindrical
    source = (-1 * u.cm, 45 * u.deg, 0 * u.mm)
    detector = (10 * u.cm, 45 * u.deg, 0 * u.mm)
    sim2 = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)

    # In spherical
    source = (-0.01 * u.m, 90 * u.deg, 45 * u.deg)
    detector = (0.1 * u.m, 90 * u.deg, 45 * u.deg)
    sim3 = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)

    assert np.allclose(sim1.source, sim2.source, atol=1e-2)
    assert np.allclose(sim2.source, sim3.source, atol=1e-2)
    assert np.allclose(sim1.detector, sim2.detector, atol=1e-2)
    assert np.allclose(sim2.detector, sim3.detector, atol=1e-2)


def test_input_validation():
    """
    Intentionally raise a number of errors.
    """

    # ************************************************************************
    # During initialization
    # ************************************************************************

    grid = _test_grid("electrostatic_gaussian_sphere")
    source = (-10 * u.mm, 90 * u.deg, 45 * u.deg)
    detector = (100 * u.mm, 90 * u.deg, 45 * u.deg)

    # Check that an error is raised when an input grid has a nan or infty value
    # First check NaN
    grid["E_x"][0, 0, 0] = np.nan * u.V / u.m
    with pytest.raises(ValueError):
        sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)

    grid["E_x"][0, 0, 0] = np.inf * u.V / u.m  # Reset element for the rest of the tests
    with pytest.raises(ValueError):
        sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    grid["E_x"][0, 0, 0] = 0 * u.V / u.m

    # Check what happens if a value is large realtive to the rest of the array
    grid["E_x"][0, 0, 0] = 0.5 * np.max(grid["E_x"])
    # with pytest.raises(ValueError):
    with pytest.warns(RuntimeWarning):
        sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    grid["E_x"][0, 0, 0] = 0 * u.V / u.m

    # Raise error when source-to-detector vector doesn't pass through the
    # field grid
    source_bad = (10 * u.mm, -10 * u.mm, 0 * u.mm)
    detector_bad = (10 * u.mm, 100 * u.mm, 0 * u.mm)
    with pytest.raises(ValueError):
        sim = prad.SyntheticProtonRadiograph(
            grid, source_bad, detector_bad, verbose=False
        )

    # Test raises warning when one (or more) of the required fields is missing
    grid_bad = CartesianGrid(-1 * u.mm, 1 * u.mm, num=50)
    with pytest.warns(RuntimeWarning):
        sim = prad.SyntheticProtonRadiograph(grid_bad, source, detector, verbose=True)

    # ************************************************************************
    # During create_particles
    # ************************************************************************
    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    sim.create_particles(1e3, 15 * u.MeV, max_theta=0.99 * np.pi / 2 * u.rad)

    # ************************************************************************
    # During runtime
    # ************************************************************************

    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    sim.create_particles(1e3, 15 * u.MeV)

    # Test an invalid field weighting keyword
    with pytest.raises(ValueError):
        sim.run(field_weighting="not a valid field weighting")

    # ************************************************************************
    # During runtime
    # ************************************************************************
    # SYNTHETIC RADIOGRAPH ERRORS
    sim.run()

    # Choose a very small synthetic radiograph size that misses most of the
    # particles
    with pytest.warns(RuntimeWarning):
        size = np.array([[-1, 1], [-1, 1]]) * 1 * u.mm
        hax, vax, values = sim.synthetic_radiograph(size=size)


def test_init():
    grid = _test_grid("electrostatic_gaussian_sphere", num=50)

    # Cartesian
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)

    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)

    # Test manually setting hdir and vdir
    hdir = np.array([1, 0, 1])
    sim = prad.SyntheticProtonRadiograph(
        grid, source, detector, verbose=False, detector_hdir=hdir
    )

    # Test special case hdir == [0,0,1]
    source = (0 * u.mm, 0 * u.mm, -10 * u.mm)
    detector = (0 * u.mm, 0 * u.mm, 200 * u.mm)
    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    assert all(sim.det_hdir == np.array([1, 0, 0]))

    # Test that hdir is calculated correctly if src-det axis is anti-parallel to z
    source = (0 * u.mm, 0 * u.mm, 10 * u.mm)
    detector = (0 * u.mm, 0 * u.mm, -200 * u.mm)
    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    assert all(sim.det_hdir == np.array([1, 0, 0]))


def test_create_particles():
    grid = _test_grid("electrostatic_gaussian_sphere", num=50)

    # Cartesian
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)

    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)

    sim.create_particles(
        1e3, 15 * u.MeV, max_theta=0.1 * u.rad, distribution="monte-carlo"
    )

    sim.create_particles(1e3, 15 * u.MeV, max_theta=0.1 * u.rad, distribution="uniform")

    # Test specifying particle
    charge = 3 * const.e.si
    mass = const.m_e.si
    sim.create_particles(1e3, 15 * u.MeV, particle="e")


def test_load_particles():

    grid = _test_grid("electrostatic_gaussian_sphere", num=50)

    # Cartesian
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)

    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    sim.create_particles(1e3, 15 * u.MeV, max_theta=0.1 * u.rad, distribution="uniform")

    # Test adding unequal numbers of particles
    x = np.zeros([100, 3]) * u.m
    v = np.ones([150, 3]) * u.m / u.s
    with pytest.raises(ValueError):
        sim.load_particles(x, v)

    # Test creating particles with explict keywords
    x = sim.x * u.m
    v = sim.v * u.m / u.s

    # Try setting particles going the wrong direction
    with pytest.warns(RuntimeWarning):
        sim.load_particles(x, -v)

    # Try specifying a larger ion (not a proton or electron)
    sim.load_particles(x, v, particle="C-12 +3")

    # Run the tracker to make sure everything works
    sim.run(field_weighting="nearest neighbor")


def test_run_options():
    grid = _test_grid("electrostatic_gaussian_sphere", num=50)

    # Cartesian
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)
    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)

    # Test that trying to call run() without creating particles
    # raises an exception
    with pytest.raises(ValueError):
        sim.run()

    sim.create_particles(1e4, 3 * u.MeV, max_theta=10 * u.deg)

    # Try running with nearest neighbor interpolator
    # Test manually setting a timestep
    sim.run(field_weighting="nearest neighbor", dt=1e-12 * u.s)

    # Test max_deflections
    sim.max_deflection

    # Test way too big of a max_theta
    sim.create_particles(1e4, 3 * u.MeV, max_theta=89 * u.deg)
    with pytest.warns(RuntimeWarning):
        sim.run(field_weighting="nearest neighbor", dt=1e-12 * u.s)

    # Test extreme deflections -> warns user
    # This requires instatiating a whole new example field with a really
    # big B-field
    grid = _test_grid("constant_bz", num=50, B0=500 * u.T)
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)
    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    sim.create_particles(1e4, 3 * u.MeV, max_theta=1 * u.deg)
    with pytest.warns(RuntimeWarning):
        sim.run(field_weighting="nearest neighbor", dt=1e-12 * u.s)


def test_synthetic_radiograph():

    # CREATE A RADIOGRAPH OBJECT
    grid = _test_grid("electrostatic_gaussian_sphere", num=50)
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)

    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    sim.create_particles(1e4, 3 * u.MeV, max_theta=10 * u.deg)
    sim.run(field_weighting="nearest neighbor")

    size = np.array([[-1, 1], [-1, 1]]) * 30 * u.cm
    bins = [200, 60]

    # Test size is None, default bins
    h, v, i = sim.synthetic_radiograph()

    # Test optical density
    h, v, i = sim.synthetic_radiograph(size=size, bins=bins, optical_density=True)



def test_insert_mesh():
    grid = _test_grid("empty", num=50)
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)

    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)



if __name__ == "__main__":
    test_coordinate_systems()
    test_input_validation()
    test_1D_deflections()
    test_init()
    test_create_particles()
    test_load_particles()
    test_run_options()
    test_synthetic_radiograph()
    pass
