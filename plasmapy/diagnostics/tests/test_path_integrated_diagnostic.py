import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pytest

from plasmapy.diagnostics.path_integrated_diagnostic import (
    Interferometer,
    LineIntegratedDiagnostic,
    LineIntegrateScalarQuantities,
    PathIntegratedDiagnostic,
)
from plasmapy.plasma.grids import CartesianGrid, NonUniformCartesianGrid


@pytest.fixture
def grid():
    ax = np.linspace(-2, 2, num=200) * u.mm
    xarr, yarr, zarr = np.meshgrid(ax, ax, ax, indexing="ij")
    radius = np.sqrt(xarr ** 2 + yarr ** 2)
    field = np.where(radius < 1 * u.mm, 1, 0) * u.kg / u.m ** 3
    grid = CartesianGrid(xarr, yarr, zarr)
    grid.add_quantities(rho=field)
    return grid


@pytest.fixture
def nonuniform_cartesian_grid():
    # Make this grid intentionally narrow in Z so the line-integral has
    # lots of particles per xy bin
    grid = NonUniformCartesianGrid(-1 * u.cm, 1 * u.cm, num=(100, 100, 10))
    xarr, yarr, zarr = grid.grids
    radius = np.sqrt(xarr ** 2 + yarr ** 2)
    field = np.where(radius < 1 * u.mm, 1, 0) * u.kg / u.m ** 3
    grid.add_quantities(rho=field)
    return grid


def test_abstract_line_integrated_diagnostic(grid):
    source = (0 * u.mm, 0 * u.mm, -5 * u.mm)
    detector = (0 * u.mm, 0 * u.mm, 5 * u.mm)

    # Test that you can't instantiate the base class
    with pytest.raises(TypeError):
        obj = LineIntegratedDiagnostic(grid, source, detector)


def test_integrate_scalar_quantities(grid):
    source = (0 * u.mm, 0 * u.mm, -5 * u.mm)
    detector = (0 * u.mm, 0 * u.mm, 5 * u.mm)

    # Test that you can't calculate a line integral on the abstract base class
    with pytest.raises(ValueError):
        obj = LineIntegrateScalarQuantities(grid, source, detector, "B_x")


@pytest.mark.parametrize("gridname", ["grid", "nonuniform_cartesian_grid"])
def test_constant_cylinder(gridname, request):
    # Load the grid fixture
    grid = request.getfixturevalue(gridname)

    source = (0 * u.mm, -5 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 5 * u.mm, 0 * u.mm)

    obj = LineIntegrateScalarQuantities(grid, source, detector, "rho", verbose=True)

    # Test that line-integrating with collimated = False yields
    size = np.array([[-2, 2], [-2, 2]]) * u.mm
    hax, vax, integral = obj._line_integral(
        size=size, bins=[100, 100], collimated=True, num=100
    )
    hax = hax.to(u.mm).value
    vax = vax.to(u.mm).value
    integral = integral.to(u.kg / u.m ** 2).value

    line = np.mean(integral, axis=1)

    """
    theory = np.where(np.abs(hax) < 0.999, 2 * np.sqrt(np.abs(1 - hax ** 2)) * 1e-3, 0)
    plt.plot(hax, line)
    plt.plot(hax, theory)
    plt.show()
    """

    height = np.max(line)
    wi = np.argmax(np.where(line > 0.00010, 1, 0))
    halfwidth = np.abs(hax[wi])

    assert np.isclose(height, 0.002, atol=0.01)
    assert np.isclose(halfwidth, 1, atol=0.25)


def test_non_collimated(grid):
    """
    Tests line-integration with a point source by checking that the
    width of the feature scales with the magnification.
    """
    source = (0 * u.mm, -5 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 10 * u.mm, 0 * u.mm)
    obj = LineIntegrateScalarQuantities(grid, source, detector, "rho", verbose=False)

    # Run the same test but with a point source (not collimated)

    size = np.array([[-1, 1], [-1, 1]]) * 5 * u.mm
    hax, vax, integral = obj._line_integral(
        size=size, bins=[100, 100], collimated=False, num=100
    )

    # Scaling hax/mag gives hax in the object plane
    # That will now agree with the theory
    hax = hax.to(u.mm).value / obj.mag
    vax = vax.to(u.mm).value

    integral = integral.to(u.kg / u.m ** 2).value
    line = np.mean(integral, axis=1)

    theory = np.where(np.abs(hax) < 0.999, 2 * np.sqrt(np.abs(1 - hax ** 2)) * 1e-3, 0)

    """
    import matplotlib.pyplot as plt
    plt.plot(hax, line)
    plt.plot(hax, theory)
    plt.show()
    """

    assert np.allclose(line, theory, atol=5e-4)


def test_constant_box(grid):
    x, y, z = grid.grid.T
    t1 = np.where(np.abs(x) < 1 * u.mm, 1, 0)
    t2 = np.where(np.abs(z) < 1 * u.mm, 1, 0)
    field = np.where(t1 * t2 != 0, 1, 0) * u.kg / u.m ** 3
    grid.add_quantities(rho=field)

    source = (0 * u.mm, -5 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 5 * u.mm, 0 * u.mm)
    obj = LineIntegrateScalarQuantities(grid, source, detector, "rho", verbose=False)

    # Test that line-integrating with collimated = False yields

    size = np.array([[-2, 2], [-2, 2]]) * u.mm
    hax, vax, integral = obj._line_integral(
        size=size, bins=[100, 100], collimated=True, num=100
    )
    hax = hax.to(u.mm).value
    vax = vax.to(u.mm).value

    line = np.mean(integral, axis=1)

    # The value should be the density in the rectangular region times
    # it's width, which is 2 mm
    val = (1 * u.kg / u.m ** 3 * 2 * u.mm).to(integral.unit)
    theory = np.where(np.abs(hax) < 1, val, 0 * integral.unit)

    """
    plt.plot(hax, line)
    plt.plot(hax, theory)
    plt.show()
    """

    assert np.allclose(line, theory, atol=2e-4)


def test_interferogram_sphere(grid):
    x, y, z = grid.grid.T
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    n_e = np.where(r < 1 * u.mm, 1, 0) * 3e19 / u.cm ** 3
    grid.add_quantities(n_e=n_e)

    source = (0 * u.mm, -5 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 5 * u.mm, 0 * u.mm)
    obj = Interferometer(grid, source, detector, verbose=False)

    size = np.array([[-2, 2], [-2, 2]]) * u.mm
    bins = [350, 350]

    hax, vax, phase = obj.evaluate(
        1.14e15 * u.Hz,
        size=size,
        bins=bins,
        num=100,
        unwrapped=False,
    )

    hax, vax, phase = obj.evaluate(
        1.14e15 * u.Hz,
        size=size,
        bins=bins,
        num=100,
        unwrapped=True,
    )
