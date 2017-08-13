import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from pathlib import Path
import datetime

from plasmapy import Plasma

savedir = Path("~/PlasmaPy/plasmapy/tests/test_output").expanduser()
if not savedir.exists():
    savedir.mkdir(parents=True)


def gaussian(x, mean=0.0, std=1.0, amp=1.0):
    """Simple function to return a Gaussian distribution"""
    if isinstance(x, list):
        x = np.array(x)
    power = -((x - mean) ** 2.0) / (2.0 * (std ** 2.0))
    f = amp * np.exp(power)
    if amp == 1:
        f = f / max(f)
    return f


def riemann_shock():
    """
    Run simulation for the Riemann shock tube problem.
    """

    print('riemann_shock()')

    # Define simulation grid and coordinates
    print('- Initiating Simulation... ', end='')
    riemann = Plasma(domain_x=np.linspace(0, 1, 128)*u.m,
                     domain_y=[0]*u.m, domain_z=[0]*u.m,
                     gamma=1.4)
    grid = riemann.domain_shape
    x = riemann.grid[0]
    print('Success')

    # Define initial parameter values - only in perturbed component
    print('- Setting initial conditions... ', end='')
    density = np.zeros(grid) * u.kg / u.m**3
    energy = np.zeros(grid) * u.J / u.m**3

    discont = 0.5 * u.m
    density[x < discont] = 1.0 * density.unit
    density[x > discont] = 0.125 * density.unit
    energy[x < discont] = 2.5 * energy.unit
    energy[x > discont] = 0.25 * energy.unit

    # Set simulation parameters to values defined above
    riemann.density = density
    riemann.energy = energy
    print('Success')

    t0 = datetime.datetime.now()
    for maxt in np.arange(0, 0.21, 0.02):
        print('- Running Simulation... ', end='')
        riemann.simulate(max_time=maxt*u.s)
        print("Simulation complete")

        fig, ax = plt.subplots(2, 2)
        ax = ax.flatten()
        ax[0].plot(x, riemann.density)
        ax[0].set_ylim(0, 1.2)
        ax[1].plot(x, riemann.velocity[0, :, ...])
        ax[1].set_ylim(0, 1.1)
        ax[2].plot(x, riemann.energy)
        ax[3].plot(x, riemann.pressure)
        ax[3].set_ylim(0, 1.2)
        plt.savefig(
            str(savedir/"riemann_shock_{:.4f}".format(maxt)).replace('.', '_'))
        plt.close()

    riemann.simulate(max_time=0.2*u.s)
    t1 = datetime.datetime.now()
    dt = t1 - t0

    print(f'Simulation begun {t0}, ended {t1}, duration {dt}')
    params = [riemann.density, riemann.velocity[0, :, ...],
              riemann.pressure, riemann.energy]
    labels = [r'$\rho$', r'$v_x$', r'$p$', r'$e$']

    fig, axes = plt.subplots(2, 2, figsize=(12, 6))
    axes = axes.flatten()
    for i, ax in enumerate(axes):
        ax.plot(x, params[i])
        ax.set_ylabel(labels[i])
    fig.savefig(str(savedir/"riemann_shock_final"))
    plt.close(fig)


def test_mhd_waves():
    """
    """

    # Define simulation grid and coordinates
    print('- Initiating Simulation... ', end='')
    waves = Plasma(domain_x=np.linspace(-0.5, 0.5, 128)*u.m,
                   domain_y=np.linspace(-0.5, 0.5, 128)*u.m,
                   domain_z=np.linspace(0, 1, 1)*u.m)
    # TODO Make these vaiable names make more sense
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
        IM = ax.imshow(waves.density.value, cmap='plasma')
        plt.colorbar(IM)
        fig.savefig(str(savedir/f"mhd_waves_{max_i}"))
        plt.close(fig)


if __name__ == '__main__':
    riemann_shock()
