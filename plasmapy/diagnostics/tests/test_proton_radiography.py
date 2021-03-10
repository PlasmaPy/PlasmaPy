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


def _test_grid(
    name,
    L=1 * u.mm,
    num=100,
    B0=10 * u.T,
    E0=5e8 * u.V / u.m,
    phi0=1.4e5 * u.V,
    a=None,
    b=None,
):
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

    E0, B0, phi0 : u.Quantities
        Scaling quantities used in the various examples

    a, b : u.Quantities
        Two length scales used in the various examples

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
        Ex = np.ones(grid.shape) * E0
        grid.add_quantities(E_x=Ex)

    elif name == "axially_magnetized_cylinder":
        if a is None:
            a = L / 4
        radius = np.linalg.norm(grid.grid[..., 0:2] * grid.unit, axis=3)
        Bz = np.where(radius < a, B0, 0 * u.T)
        grid.add_quantities(B_z=Bz)

    elif name == "electrostatic_discontinuity":
        if a is None:
            a = L / 2
        delta = a / 120

        radius = np.linalg.norm(grid.grid[..., 0:2] * grid.unit, axis=3)
        z = grid.grids[2]

        potential = (1 - erf(z / delta)) * np.exp(-((radius / a) ** 2)) * u.V

        Ex, Ey, Ez = np.gradient(potential, grid.dax0, grid.dax1, grid.dax2)

        grid.add_quantities(E_x=Ex, E_y=Ey, E_z=Ez, phi=potential)

    elif name == "electrostatic_gaussian_sphere":
        if a is None:
            a = L / 3
        if b is None:
            b = L / 2
        radius = np.linalg.norm(grid.grid, axis=3)
        arg = (radius / a).to(u.dimensionless_unscaled)
        potential = phi0 * np.exp(-(arg ** 2))

        Ex, Ey, Ez = np.gradient(potential, grid.dax0, grid.dax1, grid.dax2)

        Ex = np.where(radius < b, Ex, 0)
        Ey = np.where(radius < b, Ey, 0)
        Ez = np.where(radius < b, Ez, 0)

        grid.add_quantities(E_x=-Ex, E_y=-Ey, E_z=-Ez, phi=potential)

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

    # Expect warnings because these fields aren't well-behaved at the edges
    with pytest.warns(
        RuntimeWarning, match="Fields should go to zero at edges of grid to avoid "
    ):
        sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    sim.create_particles(1e4, 3 * u.MeV, max_theta=0.1 * u.deg)
    sim.run()

    size = np.array([[-1, 1], [-1, 1]]) * 10 * u.cm
    bins = [200, 60]
    hax, vax, values = sim.synthetic_radiograph(size=size, bins=bins)

    values = np.mean(values[:, 20:40], axis=1)

    return hax, values


def run_mesh_example(
    location=np.array([0, -2, 0]) * u.mm,
    extent=(2 * u.mm, 1.5 * u.mm),
    nwires=9,
    wire_diameter=20 * u.um,
    mesh_hdir=None,
    mesh_vdir=None,
    nparticles=1e4,
    problem="electrostatic_gaussian_sphere",
):
    """
    Takes all of the add_wire_mesh parameters and runs a standard example problem
    simulation using them.

    Returns the sim object for use in additional tests
    """

    grid = _test_grid(problem, num=100)
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)

    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)

    sim.add_wire_mesh(
        location,
        extent,
        nwires,
        wire_diameter,
        mesh_hdir=mesh_hdir,
        mesh_vdir=mesh_vdir,
    )

    sim.create_particles(nparticles, 3 * u.MeV, max_theta=10 * u.deg)
    sim.run(field_weighting="nearest neighbor")

    return sim


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
    Ex = grid["E_x"]
    Ex[0, 0, 0] = np.nan * u.V / u.m
    grid.add_quantities(E_x=Ex)
    with pytest.raises(ValueError):
        sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    Ex[0, 0, 0] = 0 * u.V / u.m

    Ex[0, 0, 0] = np.inf * u.V / u.m  # Reset element for the rest of the tests
    grid.add_quantities(E_x=Ex)
    with pytest.raises(ValueError):
        sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    Ex[0, 0, 0] = 0 * u.V / u.m

    # Check what happens if a value is large realtive to the rest of the array
    Ex[0, 0, 0] = 0.5 * np.max(Ex)
    grid.add_quantities(E_x=Ex)
    # with pytest.raises(ValueError):
    with pytest.warns(RuntimeWarning):
        sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    Ex[0, 0, 0] = 0 * u.V / u.m

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
    with pytest.warns(RuntimeWarning, match="is not specified for the provided grid."):
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
    with pytest.warns(
        RuntimeWarning, match="of the particles are shown on this synthetic radiograph."
    ):
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
    sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=True)

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
    with pytest.warns(RuntimeWarning, match="of " "particles entered the field grid"):
        sim.run(field_weighting="nearest neighbor", dt=1e-12 * u.s)

    # Test extreme deflections -> warns user
    # This requires instatiating a whole new example field with a really
    # big B-field
    grid = _test_grid("constant_bz", num=50, B0=250 * u.T)
    source = (0 * u.mm, -10 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 200 * u.mm, 0 * u.mm)

    # Expectwarnings because these fields aren't well-behaved at the edges
    with pytest.warns(
        RuntimeWarning, match="Fields should go to zero at edges of grid to avoid "
    ):
        sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)
    sim.create_particles(1e4, 3 * u.MeV, max_theta=0.1 * u.deg)
    with pytest.warns(
        RuntimeWarning,
        match="particles have been " "deflected away from the detector plane",
    ):
        sim.run(field_weighting="nearest neighbor", dt=1e-12 * u.s)
    # Calc max deflection: should be between 0 and pi/2
    # Note: that's only true because max_theta is very small
    # More generally, max_deflection can be a bit bigger than pi/2 for
    # particles that begin at an angle then deflect all the way around.
    assert 0 < sim.max_deflection.to(u.rad).value < np.pi / 2


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


def test_gaussian_sphere_analytical_comparison():
    """
    This test runs a known example problem and compares to a theoretical
    model for small deflections.

    Still under construction (comparing the actual form of the radiograph
    is possible but tricky to implement).
    """

    # The Gaussian sphere problem for small deflection potentials
    # is solved in Kugland2012relation, and the equations referenced
    # below are from that paper.
    # https://doi.org/10.1063/1.4750234

    a = (1 * u.mm / 3).to(u.mm).value
    phi0 = 1.4e5
    W = 15e6

    l = 10
    L = 200

    # Define and run the problem
    # Setting b to be much larger than the problem so that the field is not
    # cut off at the edges. This is required to be directly
    # comparable to the theoretical result.
    grid = _test_grid(
        "electrostatic_gaussian_sphere",
        num=100,
        phi0=phi0 * u.V,
        a=a * u.mm,
        b=20 * u.mm,
    )
    source = (0 * u.mm, -l * u.mm, 0 * u.mm)
    detector = (0 * u.mm, L * u.mm, 0 * u.mm)

    with pytest.warns(
        RuntimeWarning, match="Fields should go to zero at edges of grid to avoid "
    ):
        sim = prad.SyntheticProtonRadiograph(grid, source, detector, verbose=False)

    sim.create_particles(1e3, W * u.eV, max_theta=12 * u.deg)
    sim.run()

    size = np.array([[-1, 1], [-1, 1]]) * 4 * u.cm
    bins = [100, 100]
    h, v, i = sim.synthetic_radiograph(size=size, bins=bins)
    h = h.to(u.mm).value / sim.mag
    v = v.to(u.mm).value / sim.mag
    r0 = h

    # Calculate a lineout across the center of the plane (y=0)
    v0 = np.argmin(np.abs(v))

    line = np.mean(i[:, v0 - 6 : v0 + 6], axis=1)
    # Zero the edge of the radiograph
    line += -np.mean(line)
    line *= 1 / np.max(np.abs(line))

    # Calculate the theoretical deflection angles (Eq. 28)
    theory = phi0 / W * np.sqrt(np.pi) * (r0 / a) * np.exp(-((r0 / a) ** 2))

    max_deflection = np.max(np.abs(theory))
    mu = np.sqrt(np.pi) * (phi0 / W) * (l / a)

    # sim_mu = sim.max_deflection.to(u.rad).value*(l/a)

    # Calculate the theoretical inversion (Eq. 31 )
    theory_deflect = -2 * mu * (1 - (r0 / a) ** 2) * np.exp(-((r0 / a) ** 2))
    theory_deflect *= 1 / np.max(np.abs(theory_deflect))

    # Uncomment for debug
    """
    print(f"Theory max deflection: {max_deflection:.6f}")
    print(f"Theory mu: {mu:.3f}")
    print(f"Sim max deflection: {sim.max_deflection.to(u.rad).value:.6f}")
    print(f"Sim mu: {sim_mu:.3f}")

    import matplotlib.pyplot as plt
    print(f"Theory max deflection: {max_deflection:.6f}")
    print(f"Theory mu: {mu:.3f}")
    print(f"Sim max deflection: {sim.max_deflection.to(u.rad).value:.6f}")
    print(f"Sim mu: {sim_mu:.3f}")

    fig, ax = plt.subplots()
    ax.pcolormesh(h, v, i.T, shading='auto', cmap='Blues_r')
    ax.set_aspect('equal')

    fig, ax = plt.subplots()
    ax.plot(h, line )
    ax.plot(h, theory_deflect)
    """

    assert np.isclose(max_deflection, sim.max_deflection.to(u.rad).value, atol=1e-3)


def test_add_wire_mesh():

    # ************************************************************
    # Test various input configurations
    # ************************************************************

    # Test a circular mesh
    run_mesh_example(extent=1 * u.mm)

    # Test providng hdir
    run_mesh_example(mesh_hdir=np.array([0.5, 0, 0.5]))

    # Test providing hdir and vdir
    run_mesh_example(mesh_hdir=np.array([0.5, 0, 0.5]), mesh_vdir=np.array([0, 0.1, 1]))

    # ************************************************************
    # Test invalid inputs
    # ************************************************************

    # Test invalid extent (too many elements)
    with pytest.raises(ValueError):
        run_mesh_example(extent=(1 * u.mm, 2 * u.mm, 3 * u.mm))

    # Test wire mesh completely blocks all particles (in thise case because
    # the wire diameter is absurdely large)
    with pytest.raises(ValueError):
        run_mesh_example(wire_diameter=5 * u.mm)

    # Test if wire mesh is not between the source and object
    with pytest.raises(ValueError):
        run_mesh_example(location=np.array([0, 3, 0]) * u.mm)

    # ************************************************************
    # Test that mesh is the right size in the detector plane, and that
    # the wire spacing images correctly.
    # This is actually a good overall test of the whole proton radiography
    # particle tracing algorithm.
    # ************************************************************
    loc = np.array([0, -2, 0]) * u.mm
    extent = (1 * u.mm, 1 * u.mm)
    wire_diameter = 30 * u.um
    nwires = 9
    sim = run_mesh_example(
        problem="empty",
        nparticles=1e5,
        location=loc,
        extent=extent,
        wire_diameter=wire_diameter,
        nwires=nwires,
    )

    # Calculate the width that the grid SHOULD have on the image plane
    src_to_mesh = np.linalg.norm(loc.si.value - sim.source)
    mesh_to_det = np.linalg.norm(sim.detector - loc.si.value)
    mag = 1 + mesh_to_det / src_to_mesh
    true_width = mag * extent[0].to(u.mm).value
    true_spacing = true_width / (nwires - 1)

    # Create a synthetic radiograph
    size = np.array([[-1, 1], [-1, 1]]) * 2 * u.cm
    bins = [100, 50]
    # Expect a warning because many particles are off the radiograph
    # (Chose max_theta so corners are covered)
    with pytest.warns(RuntimeWarning):
        h, v, i = sim.synthetic_radiograph(size=size, bins=bins)

    # Sum up the vertical direction
    line = np.sum(i, axis=1)

    # Determine the points that are on gridlines: where 1/line is above the
    # median by a lot
    ind = np.argwhere(1 / line > 2 * np.median(1 / line))
    hwhere = h.to(u.mm).value[ind]
    measured_width = np.max(hwhere) - np.min(hwhere)

    # Calculate the max spatial frequency (should be close to the grid spacing)
    dx = np.abs(size[0][1] - size[0][0]).to(u.mm).value / bins[0]
    fnyquist = int(bins[0] / 2)
    freqs = np.fft.fftfreq(h.size, d=dx)
    freqs = freqs[0:fnyquist]
    # Calculate the positive frequency power spectrum
    pspect = np.abs(np.fft.fft(1 / line)) ** 2
    pspect = pspect[0:fnyquist]
    pspect = np.where(np.abs(freqs) < 0.1, 0, pspect)  # Mask the low frequencies

    # Measured spacing is the inverse of the maximum spatial frequency
    measured_spacing = 1 / freqs[np.argmax(pspect)]

    # This test is somewhat tricky, so here's a matplotlib plot
    # that can be uncommented for debugging
    """
    fig, ax = plt.subplots(nrows=3, figsize=(4,15))
    ax[0].pcolormesh(h.to(u.mm).value, v.to(u.mm).value, i.T, cmap='Blues_r')
    ax[0].set_aspect('equal')
    ax[0].axvline(x=np.max(hwhere), color='red')
    ax[0].axvline(x=np.min(hwhere), color='red')

    ax[1].plot(h.to(u.mm).value, 1/line)
    ax[1].axhline(y=np.median(1/line))
    ax[1].axvline(x=np.max(hwhere), color='red')
    ax[1].axvline(x=np.min(hwhere), color='red')

    ax[2].plot(freqs, pspect)
    """

    # Verify that the edges of the mesh are imaged correctly
    assert np.isclose(measured_width, true_width, 1)

    # Verify that the spacing is correct by checking the FFT
    assert np.isclose(measured_spacing, true_spacing, 0.5)


if __name__ == "__main__":
    """
    test_coordinate_systems()
    test_input_validation()
    test_1D_deflections()
    test_init()
    test_create_particles()
    test_load_particles()
    test_run_options()
    test_synthetic_radiograph()
    test_add_wire_mesh()
    test_gaussian_sphere_analytical_comparison()
    """
    pass
