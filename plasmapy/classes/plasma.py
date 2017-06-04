"""
plasmapy.classes.plasma
=======================

Defines the core Plasma class used by PlasmaPy to represent plasma properties.
"""

import numpy as np
import astropy.units as u
from astropy.utils.console import ProgressBar
from .simulation import MHDSimulation, dot
from ..constants import mu0


class Plasma:
    """Core class for describing and calculating plasma parameters.

    Attributes
    ----------
    x : `astropy.units.Quantity`
        x-coordinates within the plasma domain. Equal to `domain_x`.
    y : `astropy.units.Quantity`
        y-coordinates within the plasma domain. Equal to `domain_y`.
    z : `astropy.units.Quantity`
        z-coordinates within the plasma domain. Equal to `domain_z`.
    grid : `astropy.units.Quantity`
        (3, x, y, z) array containing the values of each coordinate
        at every point in the domain.
    domain_shape : tuple
        Shape of the plasma domain.
    density : `astropy.units.Quantity`
        (x, y, z) array of mass density at every point in the domain.
    momentum : `astropy.units.Quantity`
        (3, x, y, z) array of the momentum vector at every point in the domain.
    pressure : `astropy.units.Quantity`
        (x, y, z) array of pressure at every point in the domain.
    magnetic_field : `astropy.units.Quantity`
        (3, x, y, z) array of the magnetic field vector at every point in the
        domain.

    Parameters
    ----------
    domain_x : `astropy.units.Quantity`
        1D array of x-coordinates for the plasma domain. Must have units
        convertable to length.
    domain_y : `astropy.units.Quantity`
        1D array of y-coordinates for the plasma domain. Must have units
        convertable to length.
    domain_z : `astropy.units.Quantity`
        1D array of z-coordinates for the plasma domain. Must have units
        convertable to length.
    """
    @u.quantity_input(domain_x=u.m, domain_y=u.m, domain_z=u.m)
    def __init__(self, domain_x, domain_y, domain_z, gamma=5/3):
        # Define domain sizes
        self.x = domain_x
        self.y = domain_y
        self.z = domain_z

        x, y, z = self.x.si.value, self.y.si.value, self.z.si.value
        self.grid = np.meshgrid(x, y, z, indexing='ij') * u.m
        self.grid = np.squeeze(self.grid)
        self.domain_shape = []
        for length in (len(x), len(y), len(z)):
            if length > 1:
                self.domain_shape.append(length)
        self.domain_shape = tuple(self.domain_shape)
        print(self.domain_shape)
        self.gamma = gamma

        # Initiate core plasma variables
        self._density = np.zeros(self.domain_shape) * u.kg / u.m**3
        self._momentum = np.zeros((3, *self.domain_shape)) * u.kg \
            / (u.m**2 * u.s)
        self._energy = np.zeros(self.domain_shape) * u.J / u.m**3
        self._magnetic_field = np.zeros((3, *self.domain_shape)) * u.T

        # Collect core variables into a list for usefulness
        self.core_variables

        # Connect a simulation object for simulating
        self.simulation_physics = MHDSimulation(self)

    """
    Define getters and setters for variables.
    """
    # ==== Core variables ====
    # Density
    @property
    def density(self):
        return self._density

    @density.setter
    @u.quantity_input
    def density(self, density: u.kg/u.m**3):
        """Sets the simulation's density profile to the specified array.
        Other arrays which depend on the density values, such as the kinetic
        pressure, are then redefined automatically.

        Parameters
        ----------

        density : ndarray
            Array of density values. Shape and size must be equal to those of
            the simulation grid.
            Must have units of density.

        """

        assert density.shape == self.domain_shape, """
            Specified density array shape {} does not match simulation grid {}.
            """.format(density.shape, self.domain_shape)
        self._density = density

    # Momentum
    @property
    def momentum(self):
        return self._momentum

    @momentum.setter
    @u.quantity_input
    def momentum(self, momentum: u.kg/(u.m**2 * u.s)):
        """Sets the simulation's momentum profile to the specified array.
        Other arrays which depend on the velocity values, such as the kinetic
        pressure,
        are then redefined automatically.

        Parameters
        ----------

        momentum : ndarray
            Array of momentum vectors. Shape must be (3, x, [y, z]), where x,
            y and z are the dimensions of the simulation grid.
            Note that a full 3D vector is necessary even if the simulation has
            fewer than 3 dimensions.

        """
        assert momentum.shape == (3, *self.domain_shape), """
            Specified density array shape {} does not match simulation grid {}.
            """.format(momentum.shape, (3, *self.domain_shape))
        self._momentum = momentum

    # Internal energy
    @property
    def energy(self):
        return self._energy

    @energy.setter
    @u.quantity_input
    def energy(self, energy: u.J/u.m**3):
        """Sets the simulation's total energy density profile to the specified array.
        Other arrays which depend on the energy values, such as the kinetic
        pressure, are then redefined automatically.

        Parameters
        ----------

        energy : ndarray
            Array of energy values. Shape must be (x, y, z), where x, y, and z
            are the grid sizes of the simulation in the x, y, and z dimensions.
            Must have units of energy.

        """

        assert energy.shape == self.domain_shape, """
            Specified density array shape {} does not match simulation grid {}.
            """.format(energy.shape, self.domain_shape)
        self._energy = energy

    # Magnetic field
    @property
    def magnetic_field(self):
        return self._magnetic_field

    @magnetic_field.setter
    @u.quantity_input
    def magnetic_field(self, magnetic_field: u.Tesla):
        """
        Sets the simulation's magnetic field profile to the specified array.
        Other arrays which depend on the magnetic field, such as the magnetic
        pressure, are then redefined automatically.

        Parameters
        ----------

        magnetic_field : ndarray
            Array of magnetic field values. Shape must be (3, x, [y, z]),
            where x, y, and z are the grid sizes of the simulation in the x, y,
            and z dimensions.
            Note that a full 3D vector is necessary even if the simulation has
            fewer than 3 dimensions.

        """

        assert magnetic_field.shape == (3, *self.domain_shape), """
            Specified density array shape {} does not match simulation grid {}.
            """.format(magnetic_field.shape, (3, *self.domain_shape))
        self._magnetic_field = magnetic_field

    @property
    def core_variables(self):
        """Returns an up-to-date list of the core variables used in the calculations.
        """
        return [self.density, self.momentum, self.energy, self.magnetic_field]

    # ==== Derived variables ====
    # Velocity
    @property
    def velocity(self):
        """Returns the velocity profile of the simulation, as calculated from the
        momentum and total density.
        """

        return self.momentum / self.density

    @velocity.setter
    @u.quantity_input
    def velocity(self, velocity: u.m / u.s):
        """Defines the velocity throughout the simulation, and automatically
        updates the momentum based on the current density values.

        Parameters
        ----------

        velocity : ndarray
            Array of velocity vectors with shape (3, x, [y, z]) where x, y and
            z are the spatial grid sizes of the simulation.
            Note that a full 3D vector is required even if the simulation is
            run for fewer than 3 dimensions.
            Must have units of velocity.

        """
        assert velocity.shape == (3, *self.domain_shape), """Specified velocity
            array shape does not match simulation grid."""
        self.momentum = velocity * self.density

    @property
    def pressure(self):
        """Sets the simulation's kinetic pressure profile to the specified array.

        The kinetic pressure is defined as:

        .. math::

            p = (\\gamma - 1) (e_0 - \\frac{\\rho\\textbf{v}^2}{2})

        """
        v = self.velocity

        return (self.gamma - 1) \
            * (self.energy - ((self.density * dot(v, v)) / 2))

    @property
    def sound_speed(self):
        """Calculate the sound speed everywhere in the domain based on the pressure,
        density and adiabatic index:

        .. math::

            c_s = \\sqrt{\\frac{\\gamma p}{\\rho}}

        """

        return np.sqrt((self.gamma * self.pressure) / self.density)

    @property
    def magnetic_field_strength(self):
        """
        """

        return np.sqrt(dot(self.magnetic_field, self.magnetic_field))

    @property
    def alfven_speed(self):
        """
        """

        return self.magnetic_field_strength / np.sqrt(mu0 * self.density)

    @u.quantity_input(max_time=u.s)
    def simulate(self, max_its=np.inf, max_time=np.inf*u.s):
        """Simulates the plasma as set up, either for the given number of
        iterations or until the simulation reaches the given time.

        Parameters
        ----------
        max_its : int
            Tells the simulation to run for a set number of iterations.

        max_time : astropy.units.Quantity
            Maximum total (in-simulation) time to allow the simulation to run.
            Must have units of time.

        Examples
        --------
        >>> # Run a simulation for exactly one thousand iterations.
        >>> myplasma.simulate(max_time=1000)

        >>> # Run a simulation for up to half an hour of simulation time.
        >>> myplasma.simulate(max_time=30*u.minute)
        """
        if np.isinf(max_its) and np.isinf(max_time.value):
            raise ValueError("Either max_time or max_its must be set.")

        physics = self.simulation_physics
        dt = physics.dt

        if np.isinf(max_time):
            pb = ProgressBar(max_its)
        else:
            pb = ProgressBar(int(max_time / dt))

        with pb as bar:
            while (physics.current_iteration < max_its
                   and physics.current_time < max_time):
                physics.time_stepper()
                bar.update()
