"""
Module containing the definition for the general particle tracker.
"""

__all__ = [
    "AbstractDiskSaveRoutine",
    "AbstractIntervalSaveRoutine",
    "AbstractMemorySaveRoutine",
    "AbstractSaveRoutine",
    "AbstractStopCondition",
    "DiskIntervalSaveRoutine",
    "MemoryIntervalSaveRoutine",
    "ParticleTracker",
    "TimeElapsedStopCondition",
]

import astropy.units as u
import collections
import h5py
import numpy as np
import sys
import warnings

from abc import ABC, abstractmethod
from collections.abc import Iterable
from pathlib import Path
from tqdm import tqdm
from typing import Optional, Union

from plasmapy.particles import Particle, particle_input
from plasmapy.plasma.grids import AbstractGrid
from plasmapy.plasma.plasma_base import BasePlasma
from plasmapy.simulation.particle_integrators import boris_push


class AbstractStopCondition(ABC):
    """Abstract base class containing the necessary methods for a ParticleTracker stopping condition."""

    @property
    def particle_tracker(self):
        """Return the `ParticleTracker` object for this stop condition."""
        return self._particle_tracker

    @particle_tracker.setter
    def particle_tracker(self, particle_tracker):
        self._particle_tracker = particle_tracker

    @abstractmethod
    def require_uniform_dt(self):
        """Return whether or not this stop condition requires a uniform dt to be specified."""
        ...

    @abstractmethod
    def description(self):
        """Return a small string describing the relevant quantity shown on the meter."""
        ...

    @abstractmethod
    def units(self):
        """Return the units of `total`."""
        ...

    @abstractmethod
    def is_finished(self):
        """Return `True` if the simulation has finished."""
        ...

    @abstractmethod
    def progress(self):
        """Return a number representing the progress of the simulation (compared to total).
        This number represents the numerator of the completion percentage.
        """
        ...

    @abstractmethod
    def total(self):
        """Return a number representing the total number of steps in a simulation.
        This number represents the denominator of the completion percentage.
        """
        ...


class TimeElapsedStopCondition(AbstractStopCondition):
    """Stop condition corresponding to the elapsed time of a ParticleTracker."""

    def __init__(self, stop_time: u.Quantity):
        self.stop_time = stop_time.to(u.s).value

    @property
    def require_uniform_dt(self):
        """The elapsed time stop condition requires a uniform time step
        to stop all particles at the same time.
        """
        return True

    @property
    def description(self):
        """The time elapsed stop condition depends on elapsed time,
        therefore the relevant quantity is time remaining.
        """
        return "Time remaining"

    @property
    def units(self):
        """The units for the time elapsed condition have the units of seconds."""

        return "seconds"

    @property
    def is_finished(self):
        """Conclude the simulation if all particles have been tracked over the specified stop time."""
        return self._particle_tracker.time >= self.stop_time

    @property
    def progress(self):
        """Return the current time step of the simulation."""
        return self._particle_tracker.time

    @property
    def total(self):  # noqa: ARG002
        """Return the total amount of time over which the particles are tracked."""
        return self.stop_time


class AbstractSaveRoutine(ABC):
    # TODO: Change this docstring to use the |ParticleTracker| alias
    """Abstract base class containing the necessary methods for a
    `~plasmapy.simulation.particle_tracker.ParticleTracker` save routine.

    The save routine class is responsible for defining the conditions and hooks
    for saving.

    Parameters
    ----------
    output_directory : `~pathlib.Path`
        Output for objects that are saved to disk


    Notes
    -----
    After every push, the `post_push_hook` method is called with the
    respective `~plasmapy.simulation.particle_tracker.ParticleTracker` object passed as a parameter.
    Then, the hook calls `save_now` to determine whether or not the simulation state should be saved.
    If it is determined that the simulation state should be saved, the save routine may call either
    `save_to_disk` or `save_to_memory` depending on the routine implemented.
    """

    @abstractmethod
    def require_uniform_dt(self):
        """Return whether or not this save routine requires a uniform dt to be specified."""
        ...

    @abstractmethod
    def save_now(self, particle_tracker):
        """Determine whether or not to save on the current push step."""
        ...

    @abstractmethod
    def save(self, particle_tracker):
        """The abstract method for saving the current state of the |ParticleTracker|."""
        ...

    def post_push_hook(self, particle_tracker):
        """Function called after a push step.

        This function is responsible for handling two steps of save routine, namely:
            - Deciding whether or not to save on the current time step
            - How the simulation data is saved (i.e. to disk or memory)

        """
        if self.save_now(particle_tracker):
            self.save(particle_tracker)


class AbstractIntervalSaveRoutine(AbstractSaveRoutine, ABC):
    """Abstract class describing a save routine that saves every given interval."""

    def __init__(self, interval: u.Quantity):
        self.save_interval = interval.to(u.s).value
        self.time_of_last_save = 0

    @property
    def require_uniform_dt(self):
        """Save output only makes sense for uniform time steps,
        therefore this routine requires a uniform time step.
        """
        return True

    def save_now(self, particle_tracker):
        """Save at every interval given in instantiation."""
        saves_per_step = np.floor(self.save_interval / particle_tracker.dt)

        # If the user is requesting multiple saves per step then raise a ValueError
        if saves_per_step == 0:
            raise ValueError(
                "You must specify a time step smaller than the save interval!"
            )

        if particle_tracker.time - self.time_of_last_save >= self.save_interval:
            self.time_of_last_save = particle_tracker.time

            return True
        else:
            return False


class AbstractDiskSaveRoutine(AbstractSaveRoutine, ABC):
    """Abstract save routine corresponding to writing a hdf5 file to disk."""

    def __init__(self, output_directory: Path):
        self.output_directory = output_directory

    def save(self, particle_tracker):
        """Save a hdf5 file containing simulation positions and velocities."""

        path = self.output_directory / f"{particle_tracker.iteration_number}.hdf5"

        with h5py.File(path, "w") as output_file:
            output_file["x"] = particle_tracker.x
            output_file["v"] = particle_tracker.v


class AbstractMemorySaveRoutine(AbstractSaveRoutine, ABC):
    """Abstract save routine corresponding to saving the state of the tracker to memory."""

    def __init__(self):
        self.x_all = []
        self.v_all = []

    def save(self, particle_tracker):
        """Append simulation positions and velocities to save routine object."""
        self.x_all.append(np.copy(particle_tracker.x))
        self.v_all.append(np.copy(particle_tracker.v))


class DiskIntervalSaveRoutine(AbstractDiskSaveRoutine, AbstractIntervalSaveRoutine):
    """Save routine corresponding to saving a hdf5 file every given interval."""

    def __init__(self, output_directory, interval: u.Quantity):
        AbstractDiskSaveRoutine.__init__(self, output_directory)
        AbstractIntervalSaveRoutine.__init__(self, interval)


class MemoryIntervalSaveRoutine(AbstractMemorySaveRoutine, AbstractIntervalSaveRoutine):
    """Save the state of the tracker every given interval."""

    def __init__(self, interval: u.Quantity):
        AbstractMemorySaveRoutine.__init__(self)
        AbstractIntervalSaveRoutine.__init__(self, interval)


class ParticleTracker:
    """
    General particle tracer.

    """

    def __init__(
        self,
        grids: Union[AbstractGrid, Iterable[AbstractGrid]],
        req_quantities=None,
        verbose=True,
    ):
        # self.grid is the grid object
        if isinstance(grids, AbstractGrid):
            self.grids = [
                grids,
            ]
        elif isinstance(grids, collections.abc.Iterable):
            self.grids = grids
        elif isinstance(grids, BasePlasma):
            raise TypeError(
                "It appears you may be trying to access an older version of the ParticleTracker class."
                "This class has been deprecated."
                "Please revert to PlasmaPy version 2023.5.1 to use this version of ParticleTracker."
            )
        else:
            raise TypeError("Type of argument `grids` not recognized.")

        # self.grid_arr is the grid positions in si units. This is created here
        # so that it isn't continuously called later
        self.grids_arr = [grid.grid.to(u.m).value for grid in self.grids]

        self.verbose = verbose

        # This flag records whether the simulation has been run
        self._has_run = False

        # *********************************************************************
        # Validate required fields
        # *********************************************************************

        for grid in self.grids:
            grid.require_quantities(req_quantities, replace_with_zeros=True)

            for rq in req_quantities:
                # Check that there are no infinite values
                if not np.isfinite(grid[rq].value).all():
                    raise ValueError(
                        f"Input arrays must be finite: {rq} contains "
                        "either NaN or infinite values."
                    )

                # Check that the max values on the edges of the arrays are
                # small relative to the maximum values on that grid
                #
                # Array must be dimensionless to re-assemble it into an array
                # of max values like this
                arr = np.abs(grid[rq]).value
                edge_max = np.max(
                    np.array(
                        [
                            np.max(a)
                            for a in (
                                arr[0, :, :],
                                arr[-1, :, :],
                                arr[:, 0, :],
                                arr[:, -1, :],
                                arr[:, :, 0],
                                arr[:, :, -1],
                            )
                        ]
                    )
                )

                if edge_max > 1e-3 * np.max(arr):
                    unit = grid.recognized_quantities[rq].unit
                    warnings.warn(
                        "Fields should go to zero at edges of grid to avoid "
                        f"non-physical effects, but a value of {edge_max:.2E} {unit} was "
                        f"found on the edge of the {rq} array. Consider applying a "
                        "envelope function to force the fields at the edge to go to "
                        "zero.",
                        RuntimeWarning,
                    )

    @property
    def num_grids(self):  # noqa: D102
        return len(self.grids)

    def _log(self, msg):
        if self.verbose:
            print(msg)

    @particle_input
    def load_particles(
        self,
        x,
        v,
        particle: Particle = Particle("p+"),  # noqa: B008
    ):
        r"""
        Load arrays of particle positions and velocities.

        Parameters
        ----------
        x : `~astropy.units.Quantity`, shape (N,3)
            Positions for N particles

        v: `~astropy.units.Quantity`, shape (N,3)
            Velocities for N particles

        particle : |particle-like|, optional
            Representation of the particle species as either a |Particle| object
            or a string representation. The default particle is protons.

        distribution: str
            A keyword which determines how particles will be distributed
            in velocity space. Options are:

                - 'monte-carlo': velocities will be chosen randomly,
                    such that the flux per solid angle is uniform.

                - 'uniform': velocities will be distributed such that,
                   left unperturbed,they will form a uniform pattern
                   on the detection plane.

            Simulations run in the ``'uniform'`` mode will imprint a grid pattern
            on the image, but will well-sample the field grid with a
            smaller number of particles. The default is ``'monte-carlo'``.


        """
        # Raise an error if the run method has already been called.
        self._enforce_order()

        self.q = particle.charge.to(u.C).value
        self.m = particle.mass.to(u.kg).value

        if x.shape[0] != v.shape[0]:
            raise ValueError(
                "Provided x and v arrays have inconsistent numbers "
                " of particles "
                f"({x.shape[0]} and {v.shape[0]} respectively)."
            )
        else:
            self.nparticles = x.shape[0]

        self.x = x.to(u.m).value
        self.v = v.to(u.m / u.s).value

    def _stop_particles(self, particles_to_stop_mask):
        if len(particles_to_stop_mask) != self.x.shape[0]:
            raise ValueError(
                f"Expected mask of size {self.x.shape[0]}, got {len(particles_to_stop_mask)}"
            )

        self.v[particles_to_stop_mask] = np.NaN

    def _remove_particles(self, particles_to_remove_mask):
        if len(particles_to_remove_mask) != self.x.shape[0]:
            raise ValueError(
                f"Expected mask of size {self.x.shape[0]}, got {len(particles_to_remove_mask)}"
            )

        self.x[particles_to_remove_mask] = np.NaN
        self.v[particles_to_remove_mask] = np.NaN

    # *************************************************************************
    # Run/push loop methods
    # *************************************************************************

    def _adaptive_dt(self, Ex, Ey, Ez, Bx, By, Bz):  # noqa: ARG002
        r"""
        Calculate the appropriate dt for each grid based on a number of
        considerations
        including the local grid resolution (ds) and the gyroperiod of the
        particles in the current fields.
        """
        # If dt was explicitly set, skip the rest of this function
        if self.dt.size == 1:
            return self.dt

        # candidate timesteps includes one per grid (based on the grid resolution)
        # plus additional candidates based on the field at each particle
        candidates = np.ones([self.nparticles, self.num_grids + 1]) * np.inf

        # Compute the timestep indicated by the grid resolution
        ds = np.array([grid.grid_resolution.to(u.m).value for grid in self.grids])
        gridstep = 0.5 * (ds / self.vmax)

        # Wherever a particle is on a grid, include that grid's gridstep
        # in the list of candidate timesteps
        for i, _grid in enumerate(self.grids):  # noqa: B007
            candidates[:, i] = np.where(self.on_grid[:, i] > 0, gridstep[i], np.inf)

        # If not, compute a number of possible timesteps
        # Compute the cyclotron gyroperiod
        Bmag = np.max(np.sqrt(Bx**2 + By**2 + Bz**2)).to(u.T).value
        # Compute the gyroperiod
        if Bmag == 0:
            gyroperiod = np.inf
        else:
            gyroperiod = 2 * np.pi * self.m / (self.q * np.max(Bmag))

        candidates[:, self.num_grids] = gyroperiod / 12

        # TODO: introduce a minimum timestep based on electric fields too!

        # Enforce limits on dt
        candidates = np.clip(candidates, self.dt[0], self.dt[1])

        # dt is the min of all the candidates for each particle
        # a separate dt is returned for each particle
        dt = np.min(candidates, axis=-1)

        # dt should never actually be infinite, so replace any infinities
        # with the largest gridstep
        return np.where(dt == np.inf, np.max(gridstep), dt)

    def _push(self):
        r"""
        Advance particles using an implementation of the time-centered
        Boris algorithm.
        """
        # Get a list of positions (input for interpolator)
        tracked_mask = self._tracked_particle_mask()

        self.iteration_number += 1

        # nparticles_tracked may fluctuate with stopping
        # all particles are therefore used regardless of whether or not they reach the grid
        pos_all = self.x
        pos_tracked = pos_all[tracked_mask]

        vel_all = self.v
        vel_tracked = vel_all[tracked_mask]

        # Update the list of particles on and off the grid
        # shape [nparticles, ngrids]
        self.on_grid = np.array([grid.on_grid(pos_all * u.m) for grid in self.grids]).T

        # entered_grid is zero at the end if a particle has never
        # entered any grid
        self.entered_grid += np.sum(self.on_grid, axis=-1)

        Ex = np.zeros(self.nparticles_tracked) * u.V / u.m
        Ey = np.zeros(self.nparticles_tracked) * u.V / u.m
        Ez = np.zeros(self.nparticles_tracked) * u.V / u.m
        Bx = np.zeros(self.nparticles_tracked) * u.T
        By = np.zeros(self.nparticles_tracked) * u.T
        Bz = np.zeros(self.nparticles_tracked) * u.T
        for grid in self.grids:
            # Estimate the E and B fields for each particle
            # Note that this interpolation step is BY FAR the slowest part of the push
            # loop. Any speed improvements will have to come from here.
            if self.field_weighting == "volume averaged":
                _Ex, _Ey, _Ez, _Bx, _By, _Bz = grid.volume_averaged_interpolator(
                    pos_tracked * u.m,
                    "E_x",
                    "E_y",
                    "E_z",
                    "B_x",
                    "B_y",
                    "B_z",
                    persistent=True,
                )
            elif self.field_weighting == "nearest neighbor":
                _Ex, _Ey, _Ez, _Bx, _By, _Bz = grid.nearest_neighbor_interpolator(
                    pos_tracked * u.m,
                    "E_x",
                    "E_y",
                    "E_z",
                    "B_x",
                    "B_y",
                    "B_z",
                    persistent=True,
                )

            # Interpret any NaN values (points off the grid) as zero
            # Do this before adding to the totals, because 0 + nan = nan
            _Ex = np.nan_to_num(_Ex, nan=0.0 * u.V / u.m)
            _Ey = np.nan_to_num(_Ey, nan=0.0 * u.V / u.m)
            _Ez = np.nan_to_num(_Ez, nan=0.0 * u.V / u.m)
            _Bx = np.nan_to_num(_Bx, nan=0.0 * u.T)
            _By = np.nan_to_num(_By, nan=0.0 * u.T)
            _Bz = np.nan_to_num(_Bz, nan=0.0 * u.T)

            # Add the values interpolated for this grid to the totals
            Ex += _Ex
            Ey += _Ey
            Ez += _Ez
            Bx += _Bx
            By += _By
            Bz += _Bz

        # Create arrays of E and B as required by push algorithm
        E = np.array(
            [Ex.to(u.V / u.m).value, Ey.to(u.V / u.m).value, Ez.to(u.V / u.m).value]
        )
        E = np.moveaxis(E, 0, -1)
        B = np.array([Bx.to(u.T).value, By.to(u.T).value, Bz.to(u.T).value])
        B = np.moveaxis(B, 0, -1)

        # Calculate the adaptive timestep from the fields currently experienced
        # by the particles
        # If user sets dt explicitly, that's handled in _adaptive_dt
        dt = self._adaptive_dt(Ex, Ey, Ez, Bx, By, Bz)

        # TODO: Test v/c and implement relativistic Boris push when required
        # vc = np.max(v)/_c

        if self.is_uniform_time:
            # If a uniform timestep is specified, use that
            self.time += dt
        else:
            # If dt is not a scalar (i.e. uniform time), make sure it can be multiplied by a
            # [nparticles, 3] shape field array
            dt = dt[tracked_mask, np.newaxis]

            # Increment the tracked particles' time by dt
            self.time[tracked_mask] += dt

        self.x[tracked_mask], self.v[tracked_mask] = boris_push(
            pos_tracked, vel_tracked, B, E, self.q, self.m, dt, inplace=False
        )

    @property
    def on_any_grid(self):
        """
        Binary array for each particle indicating whether it is currently
        on ANY grid.
        """
        return np.sum(self.on_grid, axis=-1) > 0

    @property
    def vmax(self):
        """
        Calculate the maximum velocity.
        Used for determining the grid crossing maximum timestep.
        """
        tracked_mask = self._tracked_particle_mask()

        return np.max(np.linalg.norm(self.v[tracked_mask], axis=-1))

    def _validate_run_inputs(
        self, stop_condition, save_routine, dt, field_weighting: str
    ):
        if not isinstance(stop_condition, AbstractStopCondition):
            raise TypeError("Please specify a valid stop condition.")

        if not isinstance(save_routine, AbstractSaveRoutine):
            raise TypeError("Please specify a valid save routine")

        # Will the simulation follow a uniform time step?
        # This must be the case for certain stopping conditions and saving routines
        self.is_uniform_time = isinstance(dt, u.Quantity) and isinstance(
            dt.value, float
        )

        # Raise a ValueError if dt is required by stop condition or save routine but not specified
        if (
            stop_condition.require_uniform_dt or save_routine.require_uniform_dt
        ) and not self.is_uniform_time:
            raise ValueError(
                "Please specify a uniform time step to use the simulation with this configuration!"
            )

        # Load and validate inputs
        field_weightings = ["volume averaged", "nearest neighbor"]
        if field_weighting in field_weightings:
            self.field_weighting = field_weighting
        else:
            raise ValueError(
                f"{field_weighting} is not a valid option for ",
                "field_weighting. Valid choices are",
                f"{field_weightings}",
            )

    def _tracked_particle_mask(self):
        """
        Calculates a boolean mask corresponding to particles that have not been stopped or removed.
        """
        return ~np.logical_or(np.isnan(self.x[:, 0]), np.isnan(self.v[:, 0]))

    @property
    def nparticles_tracked(self):
        """Return the number of particles that don't have NaN position or velocity."""
        return self._tracked_particle_mask().sum()

    def run(
        self,
        stop_condition: AbstractStopCondition,
        save_routine: Optional[AbstractSaveRoutine] = None,
        dt=None,
        field_weighting="volume averaged",
    ):
        r"""
        Runs a particle-tracing simulation.
        Timesteps are adaptively calculated based on the
        local grid resolution of the particles and the electric and magnetic
        fields they are experiencing.

        Parameters
        ----------
        dt : `~astropy.units.Quantity`, optional
            An explicitly set timestep in units convertible to seconds.
            Setting this optional keyword overrules the adaptive time step
            capability and forces the use of this timestep throughout. If a tuple
            of timesteps is provided, the adaptive timestep will be clamped
            between the first and second values.

        field_weighting : str
            String that selects the field weighting algorithm used to determine
            what fields are felt by the particles. Options are:

            * 'nearest neighbor': Particles are assigned the fields on
                the grid vertex closest to them.
            * 'volume averaged' : The fields experienced by a particle are a
                volume-average of the eight grid points surrounding them.

            The default is 'volume averaged'.

        stop_condition : callable
            A function responsible for determining when the simulation has reached
            a suitable termination condition. The function is passed an instance
            of `~plasmapy.simulation.particle_tracker.ParticleTracker`. The function
            must return a tuple of a boolean, corresponding to the completion status of
            the simulation, and a number, representing the progress of the simulation.

        Returns
        -------
        None

        """

        # Validate inputs to the run function
        # Sets is_uniform_time property
        self._validate_run_inputs(stop_condition, save_routine, dt, field_weighting)
        self._enforce_particle_creation()

        # By default, set dt as an infinite range (auto dt with no restrictions)
        self.dt = np.array([0.0, np.inf]) * u.s if dt is None else dt
        self.dt = (self.dt).to(u.s).value

        # Keep track of how many push steps have occurred for trajectory tracing
        self.iteration_number = 0

        self.time = 0 if self.is_uniform_time else np.zeros((self.nparticles, 1))
        # Create flags for tracking when particles during the simulation
        # on_grid -> zero if the particle is off grid, 1
        # shape [nparticles, ngrids]
        self.on_grid = np.zeros([self.nparticles, self.num_grids])

        # Entered grid -> non-zero if particle EVER entered a grid
        self.entered_grid = np.zeros([self.nparticles])

        # Update the `particle_tracker` attribute so that the stop condition can be used
        stop_condition.particle_tracker = self

        # Initialize a "progress bar" (really more of a meter)
        # Setting sys.stdout lets this play nicely with regular print()
        pbar = tqdm(
            initial=0,
            total=stop_condition.total,
            disable=not self.verbose,
            desc=stop_condition.description,
            unit=stop_condition.units,
            bar_format="{l_bar}{bar}{n:.1e}/{total:.1e} {unit}",  # noqa: FS003
            file=sys.stdout,
        )

        # Push the particles until the stop condition is satisfied
        # (no more particles on the simulation grid)
        is_finished = False
        while not is_finished:
            is_finished = stop_condition.is_finished
            progress = min(stop_condition.progress, stop_condition.total)

            pbar.n = progress
            pbar.last_print_n = progress
            pbar.update(0)

            self._push()

            if save_routine is not None:
                save_routine.post_push_hook(self)

        pbar.close()

        # Log a summary of the run

        self._log("Run completed")

        # Simulation has not run, because creating new particles changes the simulation
        self._has_run = True

    def _enforce_particle_creation(self):
        """Ensure the array position array `x` has been populated."""

        # Check to make sure particles have already been generated
        if not hasattr(self, "x"):
            raise ValueError(
                "Either the create_particles or load_particles method must be "
                "called before running the particle tracing algorithm."
            )

    def _enforce_order(self):
        r"""
        The `Tracker` methods could give strange results if setup methods
        are used again after the simulation has run. This method
        raises an error if the simulation has already been run.

        """

        if self._has_run:
            raise RuntimeError(
                "Modifying the `Tracker` object after running the "
                "simulation is not supported. Create a new `Tracker` "
                "object for a new simulation."
            )
