import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
import numpy as np

from plasmapy import simulation
from plasmapy.classes.sources import Coils
import xarray


def test_saving_loading(tmp_path):
    MINOR_RADIUS = 0.3 * u.m
    RADIUS = 0.5 * u.m
    MAIN_CURRENT = 4 * u.kA  # TOROIDAL

    COIL_CURRENTS = []

    coils = Coils.toykamak(MINOR_RADIUS, RADIUS, MAIN_CURRENT, COIL_CURRENTS)
    x = u.Quantity([[0.6, 0, 0]], u.m)
    v = u.Quantity([[0, 300, 0]], u.m / u.s)

    sim_single = simulation.ParticleTracker(coils, x, v, "p")

    solution = sim_single.run(1e-4 * u.s, 1e-5 * u.s, pusher="explicit_boris")

    filename = tmp_path / "dataset.nc"
    solution.to_netcdf(filename)
    loaded_solution = xarray.open_dataset(filename)
    xarray.testing.assert_identical(solution, loaded_solution)
    for key in solution:
        assert_quantity_allclose(solution[key], loaded_solution[key])
