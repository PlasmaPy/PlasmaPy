"""
Tests for grids.py
"""

import astropy.units as u
import numpy as np
import pytest

from plasmapy.plasma import grids as grids

rs = np.random.RandomState(120921)


@pytest.fixture
def abstract_grid_uniform():
    """
    A `pytest` fixture that generates an abstract grid that spans
    -1 cm to 1 cm in all dimensions. The grid mesh would be a valid
    CartesianGrid. Three quantities are added to the grid:

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

    radius = np.sqrt(grid.pts0**2 + grid.pts1**2 + grid.pts2**2)
    rho = radius.to(u.mm).value ** 4 * u.kg * u.m**-3
    grid.add_quantities(rho=rho)

    return grid


@pytest.fixture
def abstract_grid_nonuniform():
    """
    A `pytest` fixture that generates an abstract grid that spans
    -1 cm to 1 cm in all dimensions. The grid mesh would be a
    NonUniformCartesian grid. Three quantities are added to the grid:

    1. "x" which is the x position at each point in the grid [in cm]
    2. "y" which is the y position at each point in the grid [in cm]
    3. "z" which is the z position at each point in the grid [in cm]
    4. "rho" which is a mass density at each point in the grid [kg/m^-3]
    """
    ax0 = np.array([-1, 2, 0, 0.5, 1])
    ax1 = np.array([-1, -0.5, 0, 0.5, 1])
    ax2 = np.array([-1, -0.5, 0, 0.5, 1])
    x, y, z = np.meshgrid(ax0, ax1, ax2)

    grid = grids.NonUniformCartesianGrid(x * u.cm, y * u.cm, z * u.cm)

    # Add some data to the grid
    grid.add_quantities(x=grid.grids[0])
    grid.add_quantities(y=grid.grids[1])
    grid.add_quantities(z=grid.grids[2])

    radius = np.sqrt(grid.pts0**2 + grid.pts1**2 + grid.pts2**2)
    rho = radius.to(u.mm).value ** 4 * u.kg * u.m**-3
    grid.add_quantities(rho=rho)

    return grid


create_args = [
    # Same start, stop and num for each axis
    ([-1 * u.cm, 1 * u.cm], {"num": 10}, (10, 10, 10), None),
    # Different start, stop, and num for each axis
    (
        [
            [-1 * u.cm, -2 * u.cm, -3 * u.cm],
            [1 * u.cm, 3 * u.cm, 2 * u.cm],
        ],
        {"num": [10, 5, 3]},
        (10, 5, 3),
        None,
    ),
    # Explicit arrays of points
    (
        np.meshgrid(*[np.linspace(-1, 1, num=5) * u.cm] * 3, indexing="ij"),
        {},
        (5, 5, 5),
        None,
    ),
    # Array of quantities instead of a list
    (
        [
            np.array([-1, -2, -3]) * u.cm,
            np.array([1, 3, 2]) * u.cm,
        ],
        {"num": [10, 5, 3]},
        (10, 5, 3),
        None,
    ),
    # Start is a list of a single element, stop is a single u.Quantity
    (
        [
            [
                -1 * u.cm,
            ],
            1 * u.cm,
        ],
        {"num": 10},
        (10, 10, 10),
        None,
    ),
    # num is a list of a single element
    ([-1 * u.cm, 1 * u.cm], {"num": [10]}, (10, 10, 10), None),
    # num is a list of three elements
    ([-1 * u.cm, 1 * u.cm], {"num": [10, 5, 2]}, (10, 5, 2), None),
    # start, stop, and num are tuples
    ([(-1 * u.cm), (1 * u.cm)], {"num": (10, 5, 2)}, (10, 5, 2), None),
    # Test wrong number of positional arguments: too few
    ([1 * u.cm], {"num": 10}, None, TypeError),
    # Test wrong number of positional arguments: too nmany
    ([1 * u.cm] * 4, {"num": 10}, None, TypeError),
    # Test unequal lengths of arguments raises error
    ([1 * u.cm, [2 * u.m, 3 * u.m]], {"num": 10}, None, TypeError),
    # Test arrays of points that are different shapes
    (
        [
            rs.randn(2, 5, 3) * u.m,
            rs.randn(2, 5, 3) * u.m,
            rs.randn(2, 5, 4) * u.m,
        ],
        {},
        None,
        ValueError,
    ),
    # Start and stop must be quantities
    (
        [
            [-1, -2, -3],
            [1, 3, 2],
        ],
        {"num": [10, 5, 3]},
        (10, 5, 3),
        TypeError,
    ),
    ([-1, 1], {"num": 10}, (10, 10, 10), TypeError),
    # Test incompatible grid units
    ([1 * u.cm, 1 * u.eV], {"num": 10}, None, ValueError),
    # Non-integer num
    ([-1 * u.cm, 1 * u.cm], {"num": 10.1}, (10, 10, 10), TypeError),
]


@pytest.mark.parametrize("args,kwargs,shape,error", create_args)
def test_AbstractGrid_creation(args, kwargs, shape, error):
    """
    Test the creation of AbstractGrids

    Use CartesianGrid as the test example

    """
    # If no exception is expected, create the grid and check its shape
    if error is None:
        grid = grids.CartesianGrid(*args, **kwargs)
        assert grid.shape == shape
    # If an exception is expected, verify that it is raised
    else:
        with pytest.raises(error):
            print(f"{args = }")
            print(f"{kwargs = }")
            grids.CartesianGrid(*args, **kwargs)


def test_print_summary(abstract_grid_uniform, abstract_grid_nonuniform):
    """
    Verify that both __str__ methods can be called without errors
    """
    print(abstract_grid_uniform)
    print(abstract_grid_nonuniform)

    # Test printing a grid with no quantities
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=3)
    print(grid)

    # Test printing a grid with unrecognized quantities
    grid = grids.CartesianGrid(-1 * u.cm, 1 * u.cm, num=3)
    grid.add_quantities(unrecognized_quantity=np.ones([3, 3, 3]) * u.T)
    print(grid)


abstract_attrs = [
    ("is_uniform", bool, None, True),
    ("shape", tuple, int, (21, 21, 21)),
    ("grid", u.Quantity, None, None),
    ("grids", tuple, u.Quantity, None),
    (
        "units",
        list,
        u.core.Unit,
        [
            u.cm,
        ]
        * 3,
    ),
    ("unit0", u.core.Unit, None, u.cm),
    ("unit1", u.core.Unit, None, u.cm),
    ("unit2", u.core.Unit, None, u.cm),
    ("unit", u.core.Unit, None, u.cm),
    ("pts0", u.Quantity, None, None),
    ("pts1", u.Quantity, None, None),
    ("pts2", u.Quantity, None, None),
    ("ax0", u.Quantity, None, np.linspace(-1, 1, num=21) * u.cm),
    ("ax1", u.Quantity, None, np.linspace(-1, 1, num=21) * u.cm),
    ("ax2", u.Quantity, None, np.linspace(-1, 1, num=21) * u.cm),
    ("dax0", u.Quantity, None, 0.1 * u.cm),
    ("dax1", u.Quantity, None, 0.1 * u.cm),
    ("dax2", u.Quantity, None, 0.1 * u.cm),
    ("si_scale_factors", list, float, [0.01, 0.01, 0.01]),
    ("_ax0_si", np.ndarray, None, 0.01 * np.linspace(-1, 1, num=21)),
    ("_ax1_si", np.ndarray, None, 0.01 * np.linspace(-1, 1, num=21)),
    ("_ax2_si", np.ndarray, None, 0.01 * np.linspace(-1, 1, num=21)),
    ("_dax0_si", float, None, 0.001),
    ("_dax1_si", float, None, 0.001),
    ("_dax2_si", float, None, 0.001),
]


@pytest.mark.parametrize("attr,type,type_in_iter,value", abstract_attrs)
def test_AbstractGrid_uniform_attributes(
    attr,
    type,
    type_in_iter,
    value,
    abstract_grid_uniform,
):
    """
    Tests that the attributes of AbstractGrid have the correct type and
    values for the fixture abstract_grid_uniform.
    """
    attr = getattr(abstract_grid_uniform, attr)
    assert isinstance(attr, type)

    # If the attribute is an iterable, check the type inside too
    if type_in_iter is not None:
        for elem in attr:
            isinstance(elem, type_in_iter)

    # If an expected value is given, verify the attribute matches
    if value is not None:
        if isinstance(value, np.ndarray):
            assert np.allclose(attr, value, rtol=0.1)
        elif isinstance(value, (float, int)):
            assert np.isclose(attr, value, rtol=0.1)
        else:
            assert attr == value


abstract_attrs = [
    ("is_uniform", bool, None, False),
    ("shape", tuple, int, (125,)),
]


@pytest.mark.parametrize("attr,type,type_in_iter,value", abstract_attrs)
def test_AbstractGrid_nonuniform_attributes(
    attr,
    type,
    type_in_iter,
    value,
    abstract_grid_nonuniform,
):
    """
    Tests that the attributes of AbstractGrid have the correct type and
    values for the fixture abstract_grid_uniform.
    """

    attr = getattr(abstract_grid_nonuniform, attr)
    assert isinstance(attr, type)

    # If the attribute is an iterable, check the type inside too
    if type_in_iter is not None:
        for elem in attr:
            isinstance(elem, type_in_iter)

    # If an expected value is given, verify the attribute matches
    if value is not None:
        if isinstance(value, np.ndarray):
            assert np.allclose(attr, value, rtol=0.1)
        elif isinstance(value, (float, int)):
            assert np.isclose(attr, value, rtol=0.1)
        else:
            assert attr == value


quantities = [
    # Test adding one quantity
    ("B_x", np.ones([21, 21, 21]) * u.T, None, None, None),
    # Quantity shape does not match grid shape
    ("B_x", np.ones([10, 10, 10]) * u.T, ValueError, None, None),
    # Adding quantity with units not matching recognized quantities
    ("B_x", np.ones([21, 21, 21]) * u.kg, ValueError, None, None),
    ("not_recognized_quantity", np.ones([21, 21, 21]) * u.kg, None, UserWarning, None),
]


@pytest.mark.skip("Not testable until cylindrical or spherical grids are implemented")
def test_unit_attribute_error_case():
    """
    Verify that the unit attribute raises an exception if the units on all
    axes are not the same.
    """
    grid = grids.AbstractGrid(
        [-1 * u.cm, 0 * u.rad, -2 * u.cm],
        [1 * u.cm, 2 * np.pi * u.rad, 2 * u.cm],
        num=5,
    )

    with pytest.raises(ValueError):
        grid.unit


@pytest.mark.parametrize("key,value,error,warning,match", quantities)
def test_AbstractGrid_add_quantities(
    abstract_grid_uniform, key, value, error, warning, match
):
    """
    Tests the add_quantities method of AbstractGrid
    """
    # If an error is expected, make sure it is raised
    if error is not None:
        with pytest.raises(error):
            abstract_grid_uniform.add_quantities(**{key: value})
    # If a warning is expected, make sure it occurs
    elif warning is not None:
        with pytest.warns(warning, match=match):
            abstract_grid_uniform.add_quantities(**{key: value})
    # Ensure that the quantity was correctly added and can be accessed
    else:
        abstract_grid_uniform.add_quantities(**{key: value})

        # Quantity is accessible and matches values in dataset
        assert np.all(abstract_grid_uniform[key] == value)
        assert np.all(abstract_grid_uniform.ds[key].data == value.value)
        # Quantity is correct type and unit
        assert isinstance(abstract_grid_uniform[key], u.Quantity)
        assert abstract_grid_uniform[key].unit == value.unit
        # Quantity has same shape as grid
        assert abstract_grid_uniform[key].shape == abstract_grid_uniform.shape


req_q = [
    # Requiring an existing keyword
    (["x"], False, None, None, None),
    # Requiring a keyword that does not exist raises an exception
    (["key_does_not_exist"], False, KeyError, None, None),
    # Requiring a recognized keyword that isn't defined
    # only raises a warning if replace_with_zeros is True, but
    # an error if replace_with_zeros is False
    (["E_x"], True, None, RuntimeWarning, "This quantity will be assumed to be zero"),
    (["E_x"], False, KeyError, None, None),
    # Cannot replace an unrecognized key with zeros
    # (because we don't know the units to use)
    (["key_does_not_exist"], True, KeyError, None, None),
]


@pytest.mark.parametrize("required,replace_with_zeros,error,warning,match", req_q)
def test_AbstractGrid_require_quantities(
    abstract_grid_uniform, required, replace_with_zeros, error, warning, match
):
    """
    Tests the AbstractGrid require_quantities method
    """
    # If an error is expected, make sure it is raised
    if error is not None:
        with pytest.raises(error):
            abstract_grid_uniform.require_quantities(
                required, replace_with_zeros=replace_with_zeros
            )
    # If a warning is expected, make sure it occurs
    elif warning is not None:
        with pytest.warns(warning, match=match):
            abstract_grid_uniform.require_quantities(
                required, replace_with_zeros=replace_with_zeros
            )
    # Ensure the quantities do exist
    else:
        abstract_grid_uniform.require_quantities(
            required, replace_with_zeros=replace_with_zeros
        )

        assert all(k in abstract_grid_uniform.quantities for k in required)


def test_AbstractGrid_indexing(abstract_grid_uniform):
    """
    Tests using indexing to directly get and set quantity array elements
    """
    # Test setting a subset of a quantity array
    abstract_grid_uniform["x"][0, 0, 0] = 2 * u.cm
    assert abstract_grid_uniform["x"][0, 0, 0] == 2 * u.cm


on_grid = [
    # Test with two points: one on and one off
    (
        "uniform",
        np.array([[0.1, -0.3, 0], [3, -0.3, 0]]) * u.cm,
        np.array([True, False]),
    ),
    (
        "nonuniform",
        np.array([[0.1, -0.3, 0], [3, -0.3, 0]]) * u.cm,
        np.array([True, False]),
    ),
]


@pytest.mark.parametrize("fixture,pos,result", on_grid)
def test_AbstractGrid_on_grid(
    abstract_grid_uniform, abstract_grid_nonuniform, fixture, pos, result
):
    # Select one of the grid fixtures
    if fixture == "uniform":
        grid = abstract_grid_uniform
    else:
        grid = abstract_grid_nonuniform

    out = grid.on_grid(pos)
    assert np.all(out == result)


vector_intersect = [
    # This vector goes through the grid
    ("uniform", np.array([0, -5, 0]) * u.cm, np.array([0, 5, 0]) * u.cm, True),
    # This one doesn't
    ("uniform", np.array([0, -5, 0]) * u.cm, np.array([0, -5, 10]) * u.cm, False),
    # Nonuniform grid: This vector goes through the grid
    ("nonuniform", np.array([0, -5, 0]) * u.cm, np.array([0, 5, 0]) * u.cm, True),
    # Nonuniform grid: This one doesn't
    ("nonuniform", np.array([0, -5, 0]) * u.cm, np.array([0, -5, 10]) * u.cm, False),
]


@pytest.mark.parametrize("fixture,p1,p2,result", vector_intersect)
def test_AbstractGrid_vector_intersects(
    abstract_grid_uniform, abstract_grid_nonuniform, fixture, p1, p2, result
):
    # Select one of the grid fixtures
    if fixture == "uniform":
        grid = abstract_grid_uniform
    else:
        grid = abstract_grid_nonuniform

    assert grid.vector_intersects(p1, p2) == result
    # Test going backwards yields the same result
    assert grid.vector_intersects(p2, p1) == result


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

    radius = np.sqrt(grid.pts0**2 + grid.pts1**2 + grid.pts2**2)
    rho = radius.to(u.mm).value ** 4 * u.kg * u.m**-3
    grid.add_quantities(rho=rho)

    return grid


create_args_uniform_cartesian = [
    # Same start, stop and num for each axis
    ([-1 * u.cm, 1 * u.cm], {"num": 10}, (10, 10, 10), None),
    # Incompatible units raises an exception
    (
        [[-1 * u.cm, 0 * u.rad, -2 * u.cm], [1 * u.cm, 2 * np.pi * u.rad, 2 * u.cm]],
        {"num": 3},
        (3, 3, 3),
        ValueError,
    ),
]


@pytest.mark.parametrize("args,kwargs,shape,error", create_args_uniform_cartesian)
def test_CartesianGrid_creation(args, kwargs, shape, error):
    # If no exception is expected, create the grid and check its shape
    if error is None:
        grid = grids.CartesianGrid(*args, **kwargs)
        assert grid.shape == shape
    # If an exception is expected, verify that it is raised
    else:
        with pytest.raises(error):
            grid = grids.CartesianGrid(*args, **kwargs)


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

    ax0 = np.sort(rs.uniform(low=-1, high=1, size=100)) * u.cm
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

    radius = np.sqrt(grid.pts0**2 + grid.pts1**2 + grid.pts2**2)
    rho = radius.to(u.mm).value ** 4 * u.kg * u.m**-3
    grid.add_quantities(rho=rho)

    return grid


create_args_nonuniform_cartesian = [
    # Same start, stop and num for each axis
    ([-1 * u.cm, 1 * u.cm], {"num": 3}, (27,), None),
    # Incompatible units raises an exception
    (
        [[-1 * u.cm, 0 * u.rad, -2 * u.cm], [1 * u.cm, 2 * np.pi * u.rad, 2 * u.cm]],
        {"num": 3},
        (27,),
        ValueError,
    ),
]


@pytest.mark.parametrize("args,kwargs,shape,error", create_args_nonuniform_cartesian)
def test_NonUniformCartesianGrid_creation(args, kwargs, shape, error):
    # If no exception is expected, create the grid and check its shape
    if error is None:
        grid = grids.NonUniformCartesianGrid(*args, **kwargs)
        assert grid.shape == shape
    # If an exception is expected, verify that it is raised
    else:
        with pytest.raises(error):
            grid = grids.NonUniformCartesianGrid(*args, **kwargs)


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
    analytic = interp_hax**4
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
    analytic = np.sqrt(xax**2 + yax**2 + zax**2)

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
    q1 = rs.randn(10, 10, 10) * u.kg / u.cm**3
    grid.add_quantities(rho=q1)

    # Test grid resolution for non-uniform grids
    assert 0 * u.cm < grid.grid_resolution < 2 * u.cm

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
    radius = np.sqrt(grid.pts0**2 + grid.pts1**2 + grid.pts2**2)
    rho = radius.to(u.mm).value ** 4 * u.kg * u.m**-3
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

    analytic = interp_hax**4

    raw_hax = grid.ax0.to(u.mm).value
    half = 12
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
    # Seed random number generator for repeatability
    pos = rs.random([300])
    # Make sure values outside the axis on either end are included
    pos[0] = -2
    pos[1] = 102
    # Make sure at least one value is closer to the top of an interval than
    # the bottom
    pos[2] = 3.9

    expected = np.abs(pos[:, None] - ax).argmin(axis=1)

    result = grids._fast_nearest_neighbor_interpolate(pos, ax)

    assert np.allclose(expected, result)
