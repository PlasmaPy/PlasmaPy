
import astropy.units as u
import numpy as np

from scipy.special import erf

from plasmapy.diagnostics.proton_radiography import SyntheticProtonRadiograph
from plasmapy.plasma.grids import CartesianGrid

import matplotlib.pyplot as plt


def example_grid(name, L=1 * u.mm, num=100):
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


    if name == 'constant_bz':
        Bz = np.ones(grid.shape)*10*u.T
        grid.add_quantities(B_z=Bz)

    elif name == 'constant_ex':
        Ex = np.ones(grid.shape)*5e8*u.V/u.m
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
        radius = np.linalg.norm(grid.grid * grid.unit, axis=3)
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
    req_quantities = ['E_x', 'E_y', 'E_z', 'B_x', 'B_y', 'B_z']

    for q in req_quantities:
        if q not in list(grid.ds.data_vars):
            unit = grid.recognized_quantities[q].unit
            arg = {q:np.zeros(grid.shape) * unit}
            grid.add_quantities(**arg)

    return grid


def define_problem(grid='constant_ex'):

    grid = example_grid(grid)

    source = source = (0*u.mm,  -100*u.mm, 0*u.mm)
    detector = ( 0*u.mm, 100*u.mm,  0*u.mm)

    size = np.array([[-1,1],[-1 ,1]])*20*u.cm

    bins=[100,100]

    return grid, source, detector, size, bins


def test_linear_prad():

    grid, source, detector, size, bins = define_problem()
    radio = SyntheticProtonRadiograph(grid, source, detector, verbose=True)


    hax, vax, intensity = radio.linear_radiograph(size=size, bins=bins,
                                                    collimated=False, num=300,
                                                    proton_energy = 15*u.MeV)

    fig, ax = plt.subplots()
    ax.pcolormesh(hax.to(u.mm).value/radio.mag,
                  vax.to(u.mm).value/radio.mag,
                  intensity.T, cmap='Blues')
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")


def test_particle_tracker_prad():
    grid, source, detector, size, bins = define_problem()
    radio = SyntheticProtonRadiograph(grid, source, detector, verbose=True)

    radio.run(1e4, max_theta=0.9 * np.pi / 2 * u.rad, proton_energy = 15*u.MeV)



    hax, vax, intensity = radio.linear_radiograph(size=size, bins=bins,
                                                    )

    fig, ax = plt.subplots()
    ax.pcolormesh(hax.to(u.mm).value/radio_e.mag,
                  vax.to(u.mm).value/radio_e.mag,
                  intensity.T, cmap='Blues')
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")



if __name__ == '__main__':
    test_linear_prad()
    #test_particle_tracker_prad()


