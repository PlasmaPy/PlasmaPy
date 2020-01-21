"""
Class representing a group of particles moving in a plasma's electric and
magnetic fields.
"""

import numpy as np
from astropy import units as u
import numba
import tqdm.auto

from plasmapy import atomic
from plasmapy.utils.decorators import check_units
from plasmapy.utils import PhysicsError

__all__ = [
    "ParticleTracker",
    "ParticleTrackerSolution",
]

@numba.njit()
def _boris_push(x, v, b, e, hqmdt, dt):
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

class ParticleTrackerSolution:
    """
    _position_history : `np.ndarray`
    position_history : `astropy.units.Quantity`
    _velocity_history : `np.ndarray`
    velocity_history : `astropy.units.Quantity`
        History of position and velocity. Shape (nt, n, 3).
    kinetic_energy_history
        calculated from `velocity_history`.
    """
    def __init__(self, x, v, NT, dt):
        self.init_x = x.copy()
        self.init_v = v.copy()
        assert x.shape == v.shape
        self.N, dims = x.shape
        assert dims == 3
        self._position_history = np.zeros((NT, *x.shape),
                                          dtype=float)
        self._velocity_history = np.zeros((NT, *v.shape),
                                          dtype=float)
        self._dt = dt.si.value
        self._t = np.arange(NT) * self._dt

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
        # ax.set_title(self.name)
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
        # ax.set_title(self.name)
        ax.legend(loc='best')
        ax.grid()
        ax.set_xlabel(f"Time $t$ [{u.s}]")
        ax.set_ylabel(f"Position [{u.m}]")
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

    def visualize(self, figure = None):
        import mayavi
        from mayavi import mlab
        if figure is None:
            fig = mlab.figure()
        else:
            fig = figure
        x, y, z = self.position_history[:,0,:].T   # FIXME
        trajectory = mlab.plot3d(x,y,z, self.t, figure=fig, line_width=1e-13, representation='surface')
        mlab.colorbar(trajectory, title="Trajectory - Time", orientation="vertical")
        if figure is None:
            mlab.show()

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


class ParticleTracker:
    """
    A group of particles moving in a plasma's electromagnetic field.

    Parameters
    ----------
    plasma : `Plasma`
        plasma from which fields can be pulled
    type : str
        particle type. See `plasmapy.atomic.atomic` for suitable arguments.
        The default is a proton.
    x: `astropy.units.Quantity`
    v: `astropy.units.Quantity`
        Initial conditions for position and velocity, with a shape of (N, 3),
        N being the number of particles in your simulation.

    Attributes
    ----------
    _x : `np.ndarray`
    _v : `np.ndarray`
        Current position and velocity, without units. Shape (N, 3).
    q : `astropy.units.Quantity`
    m : `astropy.units.Quantity`
        Charge and mass of particle.

    Examples
    ----------
    See `plasmapy/examples/plot_particle_stepper.ipynb`.

    """

    @atomic.particle_input
    @check_units()
    def __init__(self,
                 plasma,
                 x: u.m,
                 v: u.m/u.s,
                 particle_type: atomic.Particle = 'p',
                 ):

        self.q = particle_type.charge
        self.m = particle_type.mass
        self.name = particle_type.particle

        self.plasma = plasma

        self._x = x.si.value
        self._v = v.si.value
        assert v.shape == x.shape
        self.N, dims = x.shape
        assert dims == 3

        self._check_field_size()

    def _check_field_size(self):
        b = self.plasma.interpolate_B(self.x)
        e = self.plasma.interpolate_E(self.x)
        if b.shape != self._x.shape:
            raise ValueError(
                f"""Invalid shape {b.shape} for the magnetic field array!
                `plasma.interpolate_B` must return an array of shape (N, 3),
                where N is the number of particles in the simulation, currently {self.N}."""
            )
        if e.shape != self._x.shape:
            raise ValueError(
                f"""Invalid shape {e.shape} for the electric field array!
                `plasma.interpolate_E` must return an array of shape (N, 3),
                where N is the number of particles in the simulation, currently {self.N}."""
            )


    # TODO: find way to clean up the lines below!
    @property
    def x(self):
        return u.Quantity(self._x, u.m, copy = False)

    # @check_units() # TODO
    @x.setter
    def x(self, value: u.m):
        self._x = value.si.value

    @property
    def v(self):
        return u.Quantity(self._v, u.m / u.s, copy = False)

    # @check_units()
    @v.setter
    def v(self, value: u.m / u.s):
        self._v = value.si.value

    @check_units()
    def run(self, dt: u.s, nt: int):
        r"""
        Runs a simulation instance.
         dt: u.s = np.inf * u.s,
         nt: int = np.inf,
        """
        if np.isinf(dt) and np.isinf(nt):  # coverage: ignore
            raise ValueError("Both dt and nt are infinite.")

        _hqmdt = (self.q / self.m / 2 * dt).si.value
        _dt = dt.si.value

        nt = int(nt)

        _x = self._x.copy()
        _v = self._v.copy()

        solution = ParticleTrackerSolution(self.x, self.v, nt, dt)

        with np.errstate(all='raise'):
            b = self.plasma._interpolate_B(_x)
            e = self.plasma.interpolate_E(_x).si.value
            _boris_push(_x, _v, b, e, -0.5 * _hqmdt, -0.5*_dt)

            _x = _x - _v * 0.5 * _dt

            solution._position_history[0] = _x
            solution._velocity_history[0] = _v
            for i in tqdm.auto.trange(1, nt):
                b = self.plasma._interpolate_B(_x)
                e = self.plasma.interpolate_E(_x).si.value
                _boris_push(_x, _v, b, e, _hqmdt, _dt)
                solution._position_history[i] = _x
                solution._velocity_history[i] = _v
        return solution

    def __repr__(self, *args, **kwargs):
        return f"Species(q={self.q:.4e},m={self.m:.4e},N={self.N}," \
               f"name=\"{self.name}\""

    def __str__(self):  # coverage: ignore
        return f"{self.N} {self.name} with " \
               f"q = {self.q:.2e}, m = {self.m:.2e}"
