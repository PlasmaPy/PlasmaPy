"""
Tests for grids.py
"""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.plasma import grids as grids


def test_AbstractGrid():

    # Create grid with single u.Quantity args
    grid = grids.AbstractGrid(-1 * u.cm, 1 * u.cm, num=10)
    assert grid.is_uniform

    # Create grid with lists of  u.Quantity args
    grid = grids.AbstractGrid(
        [-1 * u.cm, -1 * u.cm, -1 * u.cm],
        [1 * u.cm, 1 * u.cm, 1 * u.cm],
        num=[10, 10, 10],
    )
    assert grid.is_uniform

    # Create grid with arrays of u.quantities
    grid = grids.AbstractGrid(
        np.array([-1] * 3) * u.cm, np.array([1] * 3) * u.cm, num=[10, 10, 10]
    )
    assert grid.is_uniform

    print(grid)

    array = grid.grid
    units = grid.units

    pts0, pts1, pts2 = grid.pts0, grid.pts1, grid.pts2

    # Test wrong number of positional arguments: 1 or more than 3
    with pytest.raises(TypeError):
        grid = grids.AbstractGrid(1 * u.cm, num=10)
    with pytest.raises(TypeError):
        grid = grids.AbstractGrid(-1 * u.cm, 1 * u.cm, 1 * u.cm, 1 * u.cm)

    # Test unequal lengths of arguments raises error
    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(-1 * u.m, [2 * u.m, 3 * u.m], num=10)

    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(
            np.random.randn(2, 5, 3) * u.m,
            np.random.randn(2, 5, 3) * u.m,
            np.random.randn(2, 5, 4) * u.m,
        )

    # Test incompatible units
    with pytest.raises(ValueError):
        grid = grids.AbstractGrid(1 * u.cm, 1 * u.eV, num=10)

    # Test adding a quantity
    q = np.random.randn(10, 10, 10) * u.T
    grid.add_quantities(B_x=q)

    # Test setting a subset of a quantity array
    grid["B_x"][0, 0, 0] = 21 * u.T
    assert grid["B_x"][0, 0, 0] == 21 * u.T

    # Test accessing a quantity using __getitem__ or directly
    Bx = grid.ds["B_x"]
    Bx = grid["B_x"]
    # Assert that the array returned is a u.Quantity
    assert isinstance(Bx, u.Quantity)
    # Assert that the array returned has the right shape
    assert Bx.shape == grid.shape

    # Test require_quantities
    # Test with a key that is there
    req_q = ["B_x"]
    grid.require_quantities(req_q, replace_with_zeros=False)
    req_q = ["B_x", "B_y"]
    # Test with a key that is not there, but can be replaced
    # Do not replace
    with pytest.raises(KeyError):
        grid.require_quantities(req_q, replace_with_zeros=False)
    # Do replace
    with pytest.warns(RuntimeWarning, match="This quantity will be assumed to be zero"):
        grid.require_quantities(req_q, replace_with_zeros=True)
    req_q = ["B_x", "B_y"]
    # Test with a key that is not there, but cannot be replaced because
    # it's not a recognized key
    req_q = ["B_x", "not_a_recognized_key"]
    with pytest.raises(KeyError):
        with pytest.warns(
            RuntimeWarning, match="This quantity will be assumed to be zero"
        ):
            grid.require_quantities(req_q, replace_with_zeros=True)

    # Test adding a quantity with wrong units
    q = np.random.randn(10, 10, 10) * u.kg
    with pytest.raises(ValueError):
        grid.add_quantities(B_x=q)

    # Testing adding a quantity with an unrecognized key name
    with pytest.warns(UserWarning):
        grid.add_quantities(not_a_recognized_key=q)

    # Test adding a quantity of incompatible size
    q = np.random.randn(5, 20, 5) * u.T
    with pytest.raises(ValueError):
        grid.add_quantities(B_x=q)

    # Test adding multiple quantities at once
    q = np.random.randn(10, 10, 10) * u.T
    grid.add_quantities(B_x=q, B_y=q, B_z=q)

    print(grid)


def test_CartesianGrid():

    grid = grids.CartesianGrid(
        np.array([-1, -1, -1]) * u.cm, np.array([1, 1, 1]) * u.cm, num=(10, 10, 10)
    )

    x_arr, y_arr, z_arr = grid.grids
    x_axis, y_axis, z_axis = grid.ax0, grid.ax1, grid.ax2
    d_x, d_y, d_z = grid.dax0, grid.dax1, grid.dax2
    is_uniform = grid.is_uniform
    shape = grid.shape
    unit = grid.units

    # Grid should be uniform
    assert grid.is_uniform

    # Test initializing with a provided grid and a quantity
    q = np.zeros(grid.shape)
    grid2 = grids.CartesianGrid(
        grid.grids[0], grid.grids[1], grid.grids[2], test_quantity=q
    )

    # Test that input with the wrong units will raise an exception
    L0 = [-1 * u.mm, 0 * u.rad, -1 * u.mm]
    L1 = [1 * u.mm, 2 * np.pi * u.rad, 1 * u.mm]
    with pytest.raises(ValueError):
        grid = grids.CartesianGrid(L0, L1, num=10)


def test_grid_methods():
    # ************ UNIFORM CARTESIAN ****************************
    grid = grids.CartesianGrid(
        np.array([-1, -1, -1]) * u.cm, np.array([1, 1, 1]) * u.cm, num=(10, 10, 10)
    )

    # Test on-grid
    pos = np.array([[0.1, -0.3, 0], [3, -0.3, 0]]) * u.cm
    out = grid.on_grid(pos)
    assert np.all(out == np.array([True, False]))

    # Test vector_intersects
    # This vector passes through the grid
    p1, p2 = np.array([0, -5, 0]) * u.cm, np.array([0, 5, 0]) * u.cm
    assert grid.vector_intersects(p1, p2)
    # Test going backwards yields the same result
    assert grid.vector_intersects(p2, p1)
    # This one doesn't
    p1, p2 = np.array([0, -5, 0]) * u.cm, np.array([0, -5, 10]) * u.cm
    assert not grid.vector_intersects(p1, p2)
    assert not grid.vector_intersects(p2, p1)

    # ************ NON-UNIFORM CARTESIAN ****************************

    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=10)

    pos = np.array([[0.1, -0.3, 0], [3, -0.3, 0]]) * u.cm
    out = grid.on_grid(pos)
    assert np.all(out == np.array([True, False]))


# **********************************************************************
# Uniform Cartesian grid tests
# **********************************************************************


@pytest.fixture
def uniform_cartesian_grid():
    """
    A `pytest` fixture that generates a CartesianGrid that spans
    -1 cm to 1 cm in all dimensions.  Three quantities are added to the
    grid:

    1. "x" which is the x position at each point in the grid [in cm]
    2. "y" which is the y position at each point in the grid [in cm]
    3. "z" which is the z position at each point in the grid [in cm]
    4. "rho" which is a mass density at each point in the grid [kg/m^-3]
    """
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=21)
    # Add some data to the grid
    grid.add_quantities(x=grid.grids[0])
    grid.add_quantities(y=grid.grids[1])
    grid.add_quantities(z=grid.grids[2])

    radius = np.sqrt(grid.pts0 ** 2 + grid.pts1 ** 2 + grid.pts2 ** 2)
    rho = radius.to(u.mm).value ** 4 * u.kg * u.m ** -3
    grid.add_quantities(rho=rho)

    return grid


@pytest.mark.parametrize(
    "pos,quantities,expected",
    [  # Test one point
        (np.array([0.1, -0.3, 0]) * u.cm, ["x"], np.array([0.1]) * u.cm),
        # Test two points and two quantities
        (
            np.array([[0.1, -0.3, 0], [0.2, 0.5, 0.2]]) * u.cm,
            ["x", "y"],
            np.array([[0.1, 0.2], [-0.3, 0.5]]) * u.cm,
        ),
        # Test an out of bounds point
        (np.array([2, -0.3, 0]) * u.cm, ["x"], np.array([np.nan]) * u.cm),
    ],
)
def test_uniform_cartesian_NN_interp(pos, quantities, expected, uniform_cartesian_grid):
    """
    Test that the uniform Cartesian NN interpolator returns the correct values
    at various points to within the grid tolerance.

    """
    pout = uniform_cartesian_grid.nearest_neighbor_interpolator(pos, *quantities)

    # Should be correct to within dx/2, so 0.6 leaves some room
    assert np.allclose(
        pout, expected, atol=0.6 * uniform_cartesian_grid.dax0, equal_nan=True
    )


@pytest.mark.parametrize(
    "pos,quantities,error",
    [  # Quantity not in
        (np.array([0.1, -0.3, 0]) * u.cm, ["not_a_quantity"], KeyError),
    ],
)
def test_uniform_cartesian_NN_interp_errors(
    pos, quantities, error, uniform_cartesian_grid
):
    """
    Test that the uniform cartesian NN interpolator returns the expected
    errors.

    """
    with pytest.raises(error):
        uniform_cartesian_grid.nearest_neighbor_interpolator(pos, *quantities)


def test_uniform_cartesian_NN_interp_persistence(uniform_cartesian_grid):
    """
    Checks that the uniform Cartesian NN interpolator persistence feature
    performs correctly. Especially, this test ensures that changing the
    list of quantities while persistent==True doesn't crash the code.

    """
    pos = np.array([[0.1, -0.3, 0.1], [-0.5, 0, -0.6]]) * u.cm

    # Test with one set of quantities
    pout = uniform_cartesian_grid.nearest_neighbor_interpolator(
        pos, "x", "y", persistent=True
    )

    # Transpose of pos is required because pout is ordered first by quantity
    # (axis in this case) while pos is [N,3]
    assert np.allclose(pout, pos[:, [0, 1]].T, atol=0.6 * uniform_cartesian_grid.dax0)

    # Change quantities with persistent still True
    # Code should detect this and automatically run as
    # persistent==False for the first iteration to create the new
    # persistent arrays.
    pout = uniform_cartesian_grid.nearest_neighbor_interpolator(
        pos, "x", "z", persistent=True
    )

    assert np.allclose(pout, pos[:, [0, 2]].T, atol=0.6 * uniform_cartesian_grid.dax0)


# **********************************************************************
# Non-uniform Cartesian grid tests
# **********************************************************************


@pytest.fixture
def nonuniform_cartesian_grid():
    """
    A `pytest` fixture that generates a NonUniformCartesianGrid that spans
    -1 cm to 1 cm in all dimensions.  Three quantities are added to the
    grid:

    1. "x" which is the x position at each point in the grid [in cm]
    2. "y" which is the y position at each point in the grid [in cm]
    3. "z" which is the z position at each point in the grid [in cm]
    4. "rho" which is a mass density at each point in the grid [kg/m^-3]

    For testing purposes, we generate a special grid that is non-uniform
    but also has the following additional properties:
        - High resolution in x (for interpolation), low resolution in
         y and z to keep the array size down.
        - Points are created at the min and max of the range to ensure
          that out of bounds tests work correctly.

    """

    ax0 = np.sort(np.random.uniform(low=-1, high=1, size=100)) * u.cm
    ax0[0], ax0[-1] = -1 * u.cm, 1 * u.cm
    ax1 = np.linspace(-1, 1, num=5) * u.cm
    ax2 = np.linspace(-1, 1, num=5) * u.cm
    x, y, z = np.meshgrid(ax0, ax1, ax2, indexing="ij")

    # Create grid
    grid = grids.NonUniformCartesianGrid(x, y, z)
    # Add some data to the grid
    grid.add_quantities(x=x)
    grid.add_quantities(y=y)
    grid.add_quantities(z=z)

    radius = np.sqrt(grid.pts0 ** 2 + grid.pts1 ** 2 + grid.pts2 ** 2)
    rho = radius.to(u.mm).value ** 4 * u.kg * u.m ** -3
    grid.add_quantities(rho=rho)

    return grid


@pytest.mark.parametrize(
    "pos,quantities,expected",
    [  # Test one point
        (np.array([0.1, -0.3, 0]) * u.cm, ["x"], np.array([0.1]) * u.cm),
        # Test two points and two quantities
        # Same axis, since only x is dense
        (
            np.array([[0.1, -0.3, 0], [0.2, 0.5, 0.2]]) * u.cm,
            ["x", "x"],
            np.array([[0.1, 0.2], [0.1, 0.2]]) * u.cm,
        ),
        # Test an out of bounds point
        (np.array([2, -0.3, 0]) * u.cm, ["x"], np.array([np.nan]) * u.cm),
    ],
)
def test_nonuniform_cartesian_NN_interp(
    pos, quantities, expected, nonuniform_cartesian_grid
):
    """
    Test that the uniform Cartesian NN interpolator returns the correct values
    at various points to within the grid tolerance.
    """
    pout = nonuniform_cartesian_grid.nearest_neighbor_interpolator(pos, *quantities)

    # Determine the maximum grid spacing in x in order to set the tolerance
    # for this test
    dx_max = np.max(np.gradient(nonuniform_cartesian_grid.grid[:, 0]))
    
    assert np.allclose(pout, expected, atol=dx_max, equal_nan=True)


@pytest.mark.slow
def test_nonuniform_cartesian_nearest_neighbor_interpolator():
    """
    Note that this test is running on a very small grid, because otherwise it is
    very slow.

    """
    # Create a non-uniform grid
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=20)

    grid.add_quantities(x=grid.grids[0], y=grid.grids[1])

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x", persistent=True)
    assert np.allclose(pos[0], pout, atol=0.5)

    # Test out of bounds
    pos = np.array([-2, -0.3, 0]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x", persistent=True)
    assert np.isnan(pout)

    # Test persistence
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    pout1 = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=False)
    pout2 = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=True)
    assert np.allclose(pout1, pout2)


@pytest.mark.parametrize(
    "pos, what, expected",
    [
        (np.array([0.1, -0.3, 0.2]) * u.cm, ("x",), np.array([0.1]) * u.cm),
        (np.array([0.1, 0.25, 0.2]) * u.cm, ("x",), np.array([0.1]) * u.cm),
        (
            np.array(
                [
                    [0.1, -0.3, 0.2],
                    [-0.253, -0.1, 0.88],
                    [0.125, 0.11, -0.6],
                    [0.45, -0.16, 0.2],
                ]
            )
            * u.cm,
            ("x",),
            np.array([0.1, -0.253, 0.125, 0.45]) * u.cm,
        ),
        (
            np.array(
                [
                    [0.1, -0.3, 0.2],
                    [-0.253, -0.1, 0.88],
                    [0.125, 0.11, -0.6],
                    [0.45, -0.16, 0.2],
                ]
            )
            * u.cm,
            ("x", "y"),
            (
                np.array([0.1, -0.253, 0.125, 0.45]) * u.cm,
                np.array([-0.3, -0.1, 0.11, -0.16]) * u.cm,
            ),
        ),
    ],
)
def test_volume_averaged_interpolator_at_several_positions(
    pos, what, expected, uniform_cartesian_grid
):
    pout = uniform_cartesian_grid.volume_averaged_interpolator(pos, *what)
    if len(what) == 1:
        assert np.allclose(pout, expected)
    else:
        assert isinstance(pout, tuple)
        assert len(pout) == len(what)
        for ii in range(len(what)):
            assert np.allclose(pout[ii], expected[ii])


def test_volume_averaged_interpolator_missing_key(uniform_cartesian_grid):
    # Test quantity key not present in dataset
    pos = np.array([0.1, -0.3, 0.2]) * u.cm
    with pytest.raises(KeyError):
        uniform_cartesian_grid.volume_averaged_interpolator(pos, "B_x")


@pytest.mark.parametrize(
    "pos, nan_mask",
    [
        (np.array([-5.0, 0.0, 0.0]) * u.cm, None),
        (np.array([5.0, 0.0, 0.0]) * u.cm, None),
        (np.array([0.0, -1.2, 0.0]) * u.cm, None),
        (np.array([0.0, 1.5, 0.0]) * u.cm, None),
        (np.array([0.0, 0.0, -100.0]) * u.cm, None),
        (np.array([0.0, 0.0, 21.0]) * u.cm, None),
        (
            np.array(
                [[0.0, 0.0, 0.0], [-1.2, 0.0, 0.0], [0.9, 0.5, -0.2], [0.3, -2.0, 5.0]]
            )
            * u.cm,
            np.array([False, True, False, True]),
        ),
    ],
)
def test_volume_averaged_interpolator_handle_out_of_bounds(
    pos, nan_mask, uniform_cartesian_grid
):
    # Contains out-of-bounds values (must handle NaNs correctly)
    pout = uniform_cartesian_grid.volume_averaged_interpolator(pos, "x")
    if nan_mask is None:
        assert np.all(np.isnan(pout.value))
    else:
        assert np.all(np.isnan(pout.value[nan_mask]))
        assert np.all(~np.isnan(pout.value[~nan_mask]))


def test_volume_averaged_interpolator_persistence(uniform_cartesian_grid):
    # Try running with persistence
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    p1, p2 = uniform_cartesian_grid.volume_averaged_interpolator(
        pos, "x", "y", persistent=True
    )
    p1, p2 = uniform_cartesian_grid.volume_averaged_interpolator(
        pos, "x", "y", persistent=True
    )
    # Try changing the arg list, make sure it catches this and auto-reverts
    # to non-persistent interpolation in that case
    p1, p2 = uniform_cartesian_grid.volume_averaged_interpolator(
        pos, "x", persistent=True
    )
    assert p1.size == 1


def test_volume_averaged_interpolator_known_solutions():
    # Create a 4x4x4 test grid with positions -3, -1, 1, and 3 cm
    # Add a quantity that equals 0 when x=-3, 1 when x=-1,
    # 2 when x=1, and 3 when x= 3
    grid = grids.CartesianGrid(-3 * u.cm, 3 * u.cm, num=4)

    Bz = np.zeros([4, 4, 4]) * u.T
    Bz[1, :, :] = 1 * u.T
    Bz[2, :, :] = 2 * u.T
    Bz[3, :, :] = 3 * u.T
    grid.add_quantities(B_z=Bz)

    # Interpolate the following points:
    xpts = np.linspace(-4, 4, num=9) * u.cm
    pos = np.zeros([xpts.size, 3])
    pos[:, 0] = xpts

    pts = grid.volume_averaged_interpolator(pos, "B_z", persistent=True)

    assert np.allclose(
        pts.value, np.array([np.nan, 0, 0.5, 1, 1.5, 2, 2.5, 3, np.nan]), equal_nan=True
    )


def test_volume_averaged_interpolator_compare_NN_1D(uniform_cartesian_grid):
    # Create a low resolution test grid and check that the volume-avg
    # interpolator returns a higher resolution version
    npts = 150
    interp_pts = (
        np.array([np.linspace(-0.99, 1, num=npts), np.zeros(npts), np.zeros(npts)])
        * u.cm
    )
    interp_pts = np.moveaxis(interp_pts, 0, -1)

    interp_hax = interp_pts[:, 0].to(u.mm).value

    va_rho = uniform_cartesian_grid.volume_averaged_interpolator(interp_pts, "rho")
    nn_rho = uniform_cartesian_grid.nearest_neighbor_interpolator(interp_pts, "rho")

    a, b = np.argmin(np.abs(interp_hax + 9)), np.argmin(np.abs(interp_hax - 9))
    analytic = interp_hax ** 4
    va_error = np.sum(np.abs(analytic[a:b] - va_rho.value[a:b]))
    nn_error = np.sum(np.abs(analytic[a:b] - nn_rho.value[a:b]))

    # Assert that the volume averaged interpolator is more accurate than the
    # nearest neighbor interpolator
    assert va_error < nn_error


def test_volume_averaged_interpolator_compare_NN_3D(uniform_cartesian_grid):
    # Do the same computation as the NN_1D test but in 3D

    npts = 150
    interp_pts = (
        np.array(
            [
                np.linspace(-0.9, 0.9, num=npts),
                np.linspace(-0.9, 0.9, num=npts),
                np.linspace(-0.9, 0.9, num=npts),
            ]
        )
        * u.cm
    )
    interp_pts = np.moveaxis(interp_pts, 0, -1)

    xax = interp_pts[:, 0].to(u.mm).value
    yax = interp_pts[:, 1].to(u.mm).value
    zax = interp_pts[:, 2].to(u.mm).value
    analytic = np.sqrt(xax ** 2 + yax ** 2 + zax ** 2)

    va_rho = uniform_cartesian_grid.volume_averaged_interpolator(interp_pts, "rho")
    nn_rho = uniform_cartesian_grid.nearest_neighbor_interpolator(interp_pts, "rho")

    va_error = np.sum(np.abs(analytic - va_rho.value))
    nn_error = np.sum(np.abs(analytic - nn_rho.value))

    # Assert that the volume averaged interpolator is more accurate than the
    # nearest neighbor interpolator
    assert va_error < nn_error


def test_NonUniformCartesianGrid():
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=10)

    pts0, pts1, pts2 = grid.grids

    shape = grid.shape
    units = grid.units

    grid.add_quantities(x=pts0)
    print(grid)

    # Grid should be non-uniform
    assert not grid.is_uniform

    # Test assigning a quantity
    q1 = np.random.randn(10, 10, 10) * u.kg / u.cm ** 3
    grid.add_quantities(rho=q1)

    # Test grid resolution for non-uniform grids
    assert 0 < grid.grid_resolution < 2

    # Test volume interpolator not implemented yet
    pos = np.array([5, -0.3, 0]) * u.cm
    with pytest.raises(NotImplementedError):
        pout = grid.volume_averaged_interpolator(pos, "x")

    # Test that many properties are unavailable
    with pytest.raises(ValueError):
        grid.ax0
    with pytest.raises(ValueError):
        grid.ax1
    with pytest.raises(ValueError):
        grid.ax2
    with pytest.raises(ValueError):
        grid.dax0
    with pytest.raises(ValueError):
        grid.dax1
    with pytest.raises(ValueError):
        grid.dax2

    # Test that input with the wrong units will raise an exception
    L0 = [-1 * u.mm, 0 * u.rad, -1 * u.mm]
    L1 = [1 * u.mm, 2 * np.pi * u.rad, 1 * u.mm]
    with pytest.raises(ValueError):
        grid = grids.NonUniformCartesianGrid(L0, L1, num=10)


def debug_volume_avg_interpolator():
    """
    Plot the comparison of the nearest neighbor interpolator and volume
    averaged interpolator for `~plasmapy.plasma.grids.CartesianGrid`.
    """
    import matplotlib.pyplot as plt

    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=24)

    # add x and y positions to grid
    grid.add_quantities(x=grid.grids[0])
    grid.add_quantities(y=grid.grids[1])

    # add a mass density to grid
    radius = np.sqrt(grid.pts0 ** 2 + grid.pts1 ** 2 + grid.pts2 ** 2)
    rho = radius.to(u.mm).value ** 4 * u.kg * u.m ** -3
    grid.add_quantities(rho=rho)

    # Create a low resolution test grid
    npts = 150
    interp_pts = (
        np.array([np.linspace(-0.99, 1, num=npts), np.zeros(npts), np.zeros(npts)])
        * u.cm
    )
    interp_pts = np.moveaxis(interp_pts, 0, -1)

    interp_hax = interp_pts[:, 0].to(u.mm).value

    va_rho = grid.volume_averaged_interpolator(interp_pts, "rho")
    nn_rho = grid.nearest_neighbor_interpolator(interp_pts, "rho")

    analytic = interp_hax ** 4

    raw_hax = grid.ax0.to(u.mm).value
    half = int(25 / 2)
    raw_rho = grid["rho"][:, half, half]
    plt.plot(raw_hax, raw_rho, marker="*", label="Interpolated points")
    plt.plot(interp_hax, nn_rho, label="Nearest neighbor")
    plt.plot(interp_hax, va_rho, marker="o", label="Volume averaged")
    plt.plot(interp_hax, analytic, label="Analytic")
    plt.legend()
    plt.xlim(-11, 11)


def test_fast_nearest_neighbor_interpolate():
    """
    Confirms that the fast linear interpolator function is equivalent to
    the np.argmin(np.abs(x-y)) search method for ordered arrays
    """
    ax = 100 * np.linspace(0, 1, num=100)
    pos = np.random.random([300])
    # Make sure values outside the axis on either end are included
    pos[0] = -2
    pos[1] = 102
    # Make sure at least one value is closer to the top of an interval than
    # the bottom
    pos[2] = 3.9

    expected = np.abs(pos[:, None] - ax).argmin(axis=1)

    result = grids._fast_nearest_neighbor_interpolate(pos, ax)

    assert np.allclose(expected, result)


if __name__ == "__main__":
    # test_volume_averaged_interpolator_known_solutions()
    # test_fast_nearest_neighbor_interpolate()
    # test_uniform_cartesian_nearest_neighbor_interpolator()
    test_nonuniform_cartesian_nearest_neighbor_interpolator()
