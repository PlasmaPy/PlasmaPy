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


@pytest.mark.slow
def test_interpolate_indices():
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=25)

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    i = grid.interpolate_indices(pos)[0]
    # Assert that nearest grid cell was found
    pout = grid.grid[int(i[0]), int(i[1]), int(i[2])]
    assert np.allclose(pos, pout, atol=0.1)

    # One position, no units
    pos = np.array([0.1, -0.3, 0])
    i = grid.interpolate_indices(pos)[0]

    # Two positions
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    i = grid.interpolate_indices(pos)[0]

    # Contains out-of-bounds values (index array should contain NaNs)
    pos = np.array([5, -0.3, 0]) * u.cm
    i = grid.interpolate_indices(pos)[0]
    assert np.sum(np.isnan(i)) > 0

    # ***********************************************************************

    # Create a non-uniform grid
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=100)

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    i = grid.interpolate_indices(pos)[0]

    # Assert that nearest grid cell was found
    pout = grid.grid[int(i)]
    assert np.allclose(pos, pout, atol=0.5)


@pytest.mark.slow
def test_nearest_neighbor_interpolator():
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=25)
    # Add some data to the grid
    grid.add_quantities(x=grid.grids[0])
    grid.add_quantities(y=grid.grids[1])

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x")
    assert np.allclose(pos[0], pout, atol=0.1)

    # Test quantity key not present in dataset
    with pytest.raises(KeyError):
        pout = grid.nearest_neighbor_interpolator(pos, "B_x")

    # Two positions, two quantities
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y")

    # Contains out-of-bounds values (must handle NaNs correctly)
    pos = np.array([5, -0.3, 0]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x")
    assert np.allclose(pout, 0 * u.cm, atol=0.1)

    # Test persistence
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=True)
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=True)

    # ***********************************************************************

    # Create a non-uniform grid
    grid = grids.NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=100)

    grid.add_quantities(x=grid.grids[0], y=grid.grids[1])

    # One position
    pos = np.array([0.1, -0.3, 0]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x")
    assert np.allclose(pos[0], pout, atol=0.5)

    # Test persistence
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=True)
    pout = grid.nearest_neighbor_interpolator(pos, "x", "y", persistent=True)


@pytest.fixture
def example_grid():
    """
    A `pytest` fixture that generates a CartesianGrid that spans
    -1 cm to 1 cm in all dimensions.  Three quantities are added to the
    grid:

    1. "x" which is the x position at each point in the grid [in cm]
    2. "y" which is the y position at each point in the grid [in cm]
    3. "rho" which is a mass density at each point in the grid [kg/m^-3]
    """
    # Create grid
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=24)

    # Add some data to the grid
    grid.add_quantities(x=grid.grids[0])
    grid.add_quantities(y=grid.grids[1])

    radius = np.sqrt(grid.pts0 ** 2 + grid.pts1 ** 2 + grid.pts2 ** 2)
    rho = radius.to(u.mm).value ** 4 * u.kg * u.m ** -3
    grid.add_quantities(rho=rho)

    return grid


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
    pos, what, expected, example_grid
):
    pout = example_grid.volume_averaged_interpolator(pos, *what)
    if len(what) == 1:
        assert np.allclose(pout, expected)
    else:
        assert isinstance(pout, tuple)
        assert len(pout) == len(what)
        for ii in range(len(what)):
            assert np.allclose(pout[ii], expected[ii])


def test_volume_averaged_interpolator_missing_key(example_grid):
    # Test quantity key not present in dataset
    pos = np.array([0.1, -0.3, 0.2]) * u.cm
    with pytest.raises(KeyError):
        example_grid.volume_averaged_interpolator(pos, "B_x")


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
def test_volume_averaged_interpolator_handle_out_of_bounds(pos, nan_mask, example_grid):
    # Contains out-of-bounds values (must handle NaNs correctly)
    pout = example_grid.volume_averaged_interpolator(pos, "x")
    if nan_mask is None:
        assert np.all(np.isnan(pout.value))
    else:
        assert np.all(np.isnan(pout.value[nan_mask]))
        assert np.all(~np.isnan(pout.value[~nan_mask]))


def test_volume_averaged_interpolator_persistence(example_grid):
    # Try running with persistence
    pos = np.array([[0.1, -0.3, 0], [0.1, -0.3, 0]]) * u.cm
    p1, p2 = example_grid.volume_averaged_interpolator(pos, "x", "y", persistent=True)
    p1, p2 = example_grid.volume_averaged_interpolator(pos, "x", "y", persistent=True)
    # Try changing the arg list, make sure it catches this and auto-reverts
    # to non-persistent interpolation in that case
    p1, p2 = example_grid.volume_averaged_interpolator(pos, "x", persistent=True)
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


def test_volume_averaged_interpolator_compare_NN_1D(example_grid):
    # Create a low resolution test grid and check that the volume-avg
    # interpolator returns a higher resolution version
    npts = 150
    interp_pts = (
        np.array([np.linspace(-0.99, 1, num=npts), np.zeros(npts), np.zeros(npts)])
        * u.cm
    )
    interp_pts = np.moveaxis(interp_pts, 0, -1)

    interp_hax = interp_pts[:, 0].to(u.mm).value

    va_rho = example_grid.volume_averaged_interpolator(interp_pts, "rho")
    nn_rho = example_grid.nearest_neighbor_interpolator(interp_pts, "rho")

    a, b = np.argmin(np.abs(interp_hax + 9)), np.argmin(np.abs(interp_hax - 9))
    analytic = interp_hax ** 4
    va_error = np.sum(np.abs(analytic[a:b] - va_rho.value[a:b]))
    nn_error = np.sum(np.abs(analytic[a:b] - nn_rho.value[a:b]))

    # Assert that the volume averaged interpolator is more accurate than the
    # nearest neighbor interpolator
    assert va_error < nn_error


def test_volume_averaged_interpolator_compare_NN_3D(example_grid):
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

    va_rho = example_grid.volume_averaged_interpolator(interp_pts, "rho")
    nn_rho = example_grid.nearest_neighbor_interpolator(interp_pts, "rho")

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


if __name__ == "__main__":
    test_volume_averaged_interpolator_known_solutions()
