import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from plasmapy.diagnostics.line_integrated_diagnostic import (
    Interferometer,
    LineIntegrateScalarQuantities,
)
from plasmapy.plasma.grids import CartesianGrid

# TODO: Write a test that shows that this works when the integration axis is tilted?


def test_constant_cylinder():
    ax = np.linspace(-2, 2, num=200) * u.mm
    xarr, yarr, zarr = np.meshgrid(ax, ax, ax, indexing="ij")

    radius = np.sqrt(xarr ** 2 + yarr ** 2)

    field = np.where(radius < 1 * u.mm, 1, 0) * u.kg / u.m ** 3

    grid = CartesianGrid(xarr, yarr, zarr)
    grid.add_quantities(rho=field)

    source = (0 * u.mm, -5 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 5 * u.mm, 0 * u.mm)
    obj = LineIntegrateScalarQuantities(grid, source, detector, "rho", verbose=False)

    # Test that line-integrating with collimated = False yields

    size = np.array([[-2, 2], [-2, 2]]) * u.mm
    hax, vax, integral = obj.line_integral(
        size=size, bins=[100, 100], collimated=True, num=100
    )
    hax = hax.to(u.mm).value
    vax = vax.to(u.mm).value
    integral = integral.to(u.kg / u.m ** 2).value

    line = np.mean(integral, axis=1)

    theory = np.where(np.abs(hax) < 0.999, 2 * np.sqrt(np.abs(1 - hax ** 2)) * 1e-3, 0)

    """
    plt.plot(hax, line)
    plt.plot(hax, theory)
    plt.show()
    """

    assert np.allclose(line, theory, atol=2e-4)


def test_constant_box():
    ax = np.linspace(-2, 2, num=200) * u.mm
    xarr, yarr, zarr = np.meshgrid(ax, ax, ax, indexing="ij")

    t1 = np.where(np.abs(xarr) < 1 * u.mm, 1, 0)
    t2 = np.where(np.abs(zarr) < 1 * u.mm, 1, 0)

    field = np.where(t1 * t2 != 0, 1, 0) * u.kg / u.m ** 3

    grid = CartesianGrid(xarr, yarr, zarr)
    grid.add_quantities(rho=field)

    source = (0 * u.mm, -5 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 5 * u.mm, 0 * u.mm)
    obj = LineIntegrateScalarQuantities(grid, source, detector, "rho", verbose=False)

    # Test that line-integrating with collimated = False yields

    size = np.array([[-2, 2], [-2, 2]]) * u.mm
    hax, vax, integral = obj.line_integral(
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


def test_interferogram_sphere():
    ax = np.linspace(-2, 2, num=200) * u.mm
    xarr, yarr, zarr = np.meshgrid(ax, ax, ax, indexing="ij")
    r = np.sqrt(xarr ** 2 + yarr ** 2 + zarr ** 2)

    n_e = np.where(r < 1 * u.mm, 1, 0) * 3e19 / u.cm ** 3

    grid = CartesianGrid(xarr, yarr, zarr)
    grid.add_quantities(n_e=n_e)

    source = (0 * u.mm, -5 * u.mm, 0 * u.mm)
    detector = (0 * u.mm, 5 * u.mm, 0 * u.mm)
    obj = Interferometer(grid, source, detector, verbose=False)

    size = np.array([[-2, 2], [-2, 2]]) * u.mm
    bins = [350, 350]

    hax, vax, phase = obj.Interferogram(
        1.14e15 * u.Hz, size=size, bins=bins, num=100, interference=True
    )

    hax = hax.to(u.mm).value
    vax = vax.to(u.mm).value

    fig, ax = plt.subplots()
    ax.pcolormesh(hax, vax, phase.T, cmap="binary", shading="auto")
    ax.set_aspect("equal")
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    plt.show()


if __name__ == "__main__":
    # test_constant_cylinder()
    # test_constant_box()
    test_interferogram_sphere()
