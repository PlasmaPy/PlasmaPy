import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from plasmapy.diagnostics.line_integrated_diagnostic import (
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
    obj = LineIntegrateScalarQuantities(grid, source, detector, "rho")

    # Test that line-integrating with collimated = False yields

    size = np.array([[-2, 2], [-2, 2]]) * u.mm
    hax, vax, integral = obj.line_integral(
        size=size, bins=[100, 100], collimated=True, num=100
    )
    hax = hax.to(u.mm).value
    vax = vax.to(u.mm).value

    line = np.mean(integral, axis=1)
    line *= 1 / np.max(line)

    # TODO: Actually show that this integral should be sqrt(cos(x)) ??
    theory = np.where(np.abs(hax) < 1, np.sqrt(np.abs(np.cos(hax * np.pi / 2))), 0)

    """
    plt.plot(hax, line)
    plt.plot(hax, theory)
    plt.show()
    """

    assert np.allclose(line, theory, atol=0.1)


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
    obj = LineIntegrateScalarQuantities(grid, source, detector, "rho")

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

    assert np.allclose(line, theory, rtol=0.05)


if __name__ == "__main__":
    # test_constant_cylinder()
    test_constant_box()
