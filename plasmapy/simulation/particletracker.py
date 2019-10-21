"""
Class representing a group of particles moving in a plasma's electric and
magnetic fields.
"""

import numpy as np
from plasmapy import atomic
from astropy import constants
from astropy import units as u
import numba
from plasmapy.utils import PhysicsError, check_quantity
import tqdm

__all__ = [
    "ParticleTracker",
]


class ParticleTracker:
    """
    Object representing a species of particles: ions, electrons, or simply
    a group of particles with a particular initial velocity distribution.

    Parameters
    ----------
    plasma : `Plasma`
        plasma from which fields can be pulled
    type : str
        particle type. See `plasmapy.atomic.atomic` for suitable arguments.
        The default is a proton.
    n_particles : int
        number of macroparticles. The default is a single particle.
    dt : `astropy.units.Quantity`
        length of timestep
    nt : int
        number of timesteps

    Attributes
    ----------
    _x : `np.ndarray`
    x : `astropy.units.Quantity`
    _v : `np.ndarray`
    v : `astropy.units.Quantity`
        Current position and velocity, without and with units. Shape (n, 3).
    _position_history : `np.ndarray`
    position_history : `astropy.units.Quantity`
    _velocity_history : `np.ndarray`
    velocity_history : `astropy.units.Quantity`
        History of position and velocity. Shape (nt, n, 3).
    q : `astropy.units.Quantity`
    m : `astropy.units.Quantity`
        Charge and mass of particle.
    kinetic_energy
        calculated from `v`, as in, current velocity.
    kinetic_energy_history
        calculated from `velocity_history`.

    Examples
    ----------
    See `plasmapy/examples/plot_particle_stepper.ipynb`.

    """

    @atomic.particle_input
    @check_quantity()
    def __init__(self,
                 plasma,
                 particle_type: atomic.Particle = 'p',
                 n_particles: int = 1,
                 dt: u.s = np.inf * u.s,
                 nt: int = np.inf,
                 ):

        if np.isinf(dt) and np.isinf(nt):  # coverage: ignore
            raise ValueError("Both dt and nt are infinite.")

        self.q = particle_type.charge
        self.m = particle_type.mass
        self.N = int(n_particles)
        self.NT = int(nt)
        self.name = particle_type.element_name

        self.plasma = plasma

        self._dt = dt.si.value
        self._t = np.arange(nt) * self._dt

        self._x = np.zeros((n_particles, 3), dtype=float)
        self._v = np.zeros((n_particles, 3), dtype=float)

        self._position_history = np.zeros((self.NT, *self.x.shape),
                                          dtype=float)
        self._velocity_history = np.zeros((self.NT, *self.v.shape),
                                          dtype=float)
        self._hqmdt = (self.q / self.m / 2 * dt).si.value
        self._check_field_size()

    def _check_field_size(self):
        b = self.plasma.interpolate_B(self.x)
        e = self.plasma.interpolate_E(self.x)
        if b.shape != self._x.shape:
            raise ValueError(
                f"""Invalid shape {b.shape} for the magnetic field array!
                `plasma.interpolate_B` must return an array of shape (N, 3),
                where N is the number of particles in the simulation, currently {N}."""
            )
        if e.shape != self._x.shape:
            raise ValueError(
                f"""Invalid shape {e.shape} for the electric field array!
                `plasma.interpolate_E` must return an array of shape (N, 3),
                where N is the number of particles in the simulation, currently {N}."""
            )


    # TODO: find way to clean up the lines below!
    @property
    def x(self):
        return self._x * u.m

    @x.setter
    def x(self, value):
        self._x = value.si.value

    @property
    def v(self):
        return self._v * u.m / u.s

    @v.setter
    def v(self, value):
        self._v = value.si.value

    @property
    def dt(self):
        return self._dt * u.s

    @property
    def t(self):
        return self._t * u.s

    @property
    def position_history(self):
        return self._position_history * u.m

    @property
    def velocity_history(self):
        return self._velocity_history * u.m / u.s

    @property
    def kinetic_energy_history(self):
        r"""
        Calculates the kinetic energy history for each particle.

        Returns
        --------
        ~astropy.units.Quantity
            Array of kinetic energies, shape (nt, n).
        """
        return (self.velocity_history ** 2).sum(axis=-1) * self.m / 2

    def boris_push(self, init=False):
        r"""
        Implements the Boris algorithm for moving particles and updating their
        velocities.

        Arguments
        ----------
        init : bool (optional)
            If `True`, does not change the particle positions and sets dt
            to -dt/2.

        Notes
        ----------
        The Boris algorithm is the standard energy conserving algorithm for
        particle movement in plasma physics. See [1]_ for more details.

        Conceptually, the algorithm has three phases:

        1. Add half the impulse from electric field.
        2. Rotate the particle velocity about the direction of the magnetic
           field.
        3. Add the second half of the impulse from the electric field.

        This ends up causing the magnetic field action to be properly
        "centered" in time, and the algorithm conserves energy.

        References
        ----------
        .. [1] C. K. Birdsall, A. B. Langdon, "Plasma Physics via Computer
               Simulation", 2004, p. 58-63
        """
        b = self.plasma.interpolate_B(self.x).si.value
        e = self.plasma.interpolate_E(self.x).si.value
        if init:
            self._boris_push(self._x,
                             self._v,
                             b, e, -0.5 * self._hqmdt, -0.5*self._dt)
            self._x = self._x - self._v * 0.5 * self._dt
        else:
            self._boris_push(self._x,
                             self._v,
                             b, e, self._hqmdt, self._dt)

    @staticmethod
    @numba.njit()
    def _boris_push(x, v, b, e, hqmdt, dt):
        for i in numba.prange(len(x)):
            # add first half of electric impulse
            vminus = v[i] + hqmdt * e[i]

            # rotate to add magnetic field
            t = -b[i] * hqmdt
            s = 2 * t / (1 + (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]))
            cross_result = np.cross(vminus, t)
            vprime = vminus + cross_result
            cross_result_2 = np.cross(vprime, s)
            vplus = vminus + cross_result_2

            # add second half of electric impulse
            v[i] = vplus + e[i] * hqmdt
            x[i] += v[i] * dt

    def run(self):
        r"""
        Runs a simulation instance.
        """
        with np.errstate(all='raise'):
            self.boris_push(init=True)
            self._position_history[0] = self._x
            self._velocity_history[0] = self._v
            for i in tqdm.trange(1, self.NT):
                self.boris_push()
                self._position_history[i] = self._x
                self._velocity_history[i] = self._v

    def __repr__(self, *args, **kwargs):
        return f"Species(q={self.q:.4e},m={self.m:.4e},N={self.N}," \
               f"name=\"{self.name}\",NT={self.NT})"

    def __str__(self):  # coverage: ignore
        return f"{self.N} {self.name} with " \
               f"q = {self.q:.2e}, m = {self.m:.2e}, " \
               f"{self.saved_iterations} saved history " \
               f"steps over {self.NT} iterations"

    def plot_trajectories(self, *args, **kwargs):  # coverage: ignore
        r"""Draws trajectory history."""
        from astropy.visualization import quantity_support
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        quantity_support()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for p_index in range(self.N):
            r = self.position_history[:, p_index]
            x, y, z = r.T
            ax.plot(x, y, z, *args, **kwargs)
        ax.set_title(self.name)
        ax.set_xlabel("$x$ position")
        ax.set_ylabel("$y$ position")
        ax.set_zlabel("$z$ position")
        plt.show()

    def plot_time_trajectories(self, plot="xyz"):  # coverage: ignore
        r"""
        Draws position history versus time.

        Parameters
        ----------
        plot : str (optional)
            Enable plotting of position component x, y, z for each of these
            letters included in `plot`.
        """
        from astropy.visualization import quantity_support
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        quantity_support()
        fig, ax = plt.subplots()
        for p_index in range(self.N):
            r = self.position_history[:, p_index]
            x, y, z = r.T
            if "x" in plot:
                ax.plot(self.t, x, label=f"x_{p_index}")
            if "y" in plot:
                ax.plot(self.t, y, label=f"y_{p_index}")
            if "z" in plot:
                ax.plot(self.t, z, label=f"z_{p_index}")
        ax.set_title(self.name)
        ax.legend(loc='best')
        ax.grid()
        plt.show()

    def test_kinetic_energy(self):
        r"""Test conservation of kinetic energy."""
        conservation = np.allclose(self.kinetic_energy_history,
                                   self.kinetic_energy_history.mean(),
                                   atol=3 * self.kinetic_energy_history.std())
        if not conservation:
            try:
                from astropy.visualization import quantity_support
                import matplotlib.pyplot as plt
                from mpl_toolkits.mplot3d import Axes3D

                quantity_support()
                fig, ax = plt.subplots()
                difference = self.kinetic_energy_history - self.kinetic_energy_history[0]
                ax.plot(difference)
                plt.show()
            except ImportError:
                pass
            raise PhysicsError("Kinetic energy is not conserved!")
