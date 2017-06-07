import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from pathlib import Path
from ..plasma import Plasma


def gaussian(x, mean=0.0, std=1.0, amp=1.0):
    """Simple function to return a Gaussian distribution"""
    if isinstance(x, list):
        x = np.array(x)
    power = -((x - mean) ** 2.0) / (2.0 * (std ** 2.0))
    f = amp * np.exp(power)
    if amp == 1:
        f = f / max(f)
    return f


def test_mhd_waves():
    """
    """

    # Define simulation grid and coordinates
    print('- Initiating Simulation... ', end='')
    waves = Plasma(domain_x=np.linspace(-0.5, 0.5, 128)*u.m,
                   domain_y=np.linspace(-0.5, 0.5, 128)*u.m,
                   domain_z=np.linspace(0, 1, 1)*u.m)
    grid = waves.domain_shape
    x, y, z = waves.grid
    r = np.sqrt(x**2 + y**2)
    print('Success')

    # Define initial parameter values - only in perturbed component
    print('- Setting initial conditions... ', end='')
    bfield = np.zeros((3, *grid)) * u.T
    density = (gaussian(r.value, std=0.06, amp=0.9) + 0.1) * u.kg / u.m**3
    energy = (gaussian(r.value, std=0.06, amp=0.9) + 0.1) * u.J / u.m**3

    waves.density = density
    waves.energy = energy
    waves.magnetic_field = bfield
    print('Success')

    for max_i in range(0, 5001, 1000):
        print('- Running Simulation... ', end='')
        waves.simulate(max_its=max_i)
        print("Simulation complete")

        fig, ax = plt.subplots()
        plt.imshow(waves.density.value, cmap='plasma')
        plt.colorbar()
        savedir = Path("/home/drew/PlasmaPy/plasmapy/tests/test_output")
        if not savedir.exists():
            savedir.mkdir()
        plt.savefig(str(savedir/"mhd_waves_{}".format(max_i)))
        plt.close()
