"""
Module containing the definition for the general particle tracker.
"""

__all__ = [
    "AbstractIntervalSaveRoutine",
    "AbstractSaveRoutine",
    "AbstractStopCondition",
    "DiskIntervalSaveRoutine",
    "MemoryIntervalSaveRoutine",
    "NoParticlesOnGridsStoppingCondition",
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
    def tracker(self):
        """Return the `ParticleTracker` object for this stop condition."""
        return self._particle_tracker

    @tracker.setter
    def tracker(self, particle_tracker):
        self._particle_tracker = particle_tracker

    @property
    @abstractmethod
    def require_synchronized_dt(self):
        """Return whether or not this stop condition requires a synchronized time step."""
        ...

    @property
    @abstractmethod
    def progress_description(self):
        """Return a small string describing the relevant quantity shown on the meter."""
        ...

    @property
    @abstractmethod
    def units(self):
        """Return the units of `total`."""
        ...

    @property
    @abstractmethod
    def is_finished(self):
        """Return `True` if the simulation has finished."""
        ...

    @property
    @abstractmethod
    def progress(self):
        """Return a number representing the progress of the simulation (compared to total).
        This number represents the numerator of the completion percentage.
        """
        ...

    @property
    @abstractmethod
    def total(self):
        """Return a number representing the total number of steps in a simulation.
        This number represents the denominator of the completion percentage.
        """
        ...


class TimeElapsedStopCondition(AbstractStopCondition):
    """Stop condition corresponding to the elapsed time of a ParticleTracker."""

    def __init__(self, stop_time: u.Quantity):
        self._particle_tracker = None
        self.stop_time = stop_time.to(u.s).value

    @property
    def require_synchronized_dt(self):
        """The elapsed time stop condition requires a synchronized step
        to stop all particles at the same time.
        """
        return True

    @property
    def progress_description(self):
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
        return int(self._particle_tracker.time)

    @property
    def total(self):
        """Return the total amount of time over which the particles are tracked."""
        return self.stop_time


class NoParticlesOnGridsStoppingCondition(AbstractStopCondition):
    """Stopping condition corresponding to stopping the simulation when all particles have exited the grid."""

    def __init__(self):
        self._particle_tracker = None

    @property
    def require_synchronized_dt(self):
        """The no field stopping condition does not require a synchronized time step."""
        return False

    @property
    def progress_description(self):
        """The progress meter is described in terms of the fraction of particles still on the grid."""
        return "Number of particles still on grid"

    @property
    def units(self):
        """The progress meter is described in terms of the fraction of particles still on the grid."""
        return "Particles"

    @property
    def is_finished(self):
        """The simulation is finished when no more particles are on any grids."""
        is_not_on_grid = self._particle_tracker.on_grid == 0

        return is_not_on_grid.all() and self._particle_tracker.iteration_number > 0

    @property
    def progress(self):
        """The progress of the simulation is measured by how many particles are no longer on a grid."""
        is_not_on_grid = self._particle_tracker.on_grid == 0

        return is_not_on_grid.sum()

    @property
    def total(self):
        """The progress of the simulation is measured against the total number of particles in the simulation."""
        return self._particle_tracker.nparticles


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

    def __init__(self, output_directory: Optional[Path] = None):
        self.output_directory = output_directory

        self.x_all = []
        self.v_all = []

        self._particle_tracker = None

    @property
    def tracker(self):
        """Return the `ParticleTracker` object for this stop condition."""
        return self._particle_tracker

    @tracker.setter
    def tracker(self, particle_tracker):
        self._particle_tracker = particle_tracker

    @property
    @abstractmethod
    def require_synchronized_dt(self):
        """Return whether or not this save routine requires a synchronized time step."""
        ...

    @property
    @abstractmethod
    def save_now(self):
        """Determine whether or not to save on the current push step."""
        ...

    def save(self):
        """Save the current state of the simulation to disk or memory based on whether the output directory was set."""

        if self.output_directory is not None:
            self.save_to_disk()
        else:
            self.save_to_memory()

    def save_to_disk(self):
        """Save a hdf5 file containing simulation positions and velocities."""

        if self.output_directory is None:
            raise ValueError(
                "Please pass an `output_directory` parameter during instantiation to use disk save "
                "routines!"
            )

        path = self.output_directory / f"{self._particle_tracker.iteration_number}.hdf5"

        with h5py.File(path, "w") as output_file:
            output_file["x"] = self._particle_tracker.x
            output_file["v"] = self._particle_tracker.v

    def save_to_memory(self):
        """Append simulation positions and velocities to save routine object."""
        self.x_all.append(np.copy(self._particle_tracker.x))
        self.v_all.append(np.copy(self._particle_tracker.v))

    def post_push_hook(self, force_save=False):
        """Function called after a push step.

        This function is responsible for handling two steps of save routine, namely:
            - Deciding whether or not to save on the current time step
            - How the simulation data is saved (i.e. to disk or memory)

        """
        if self.save_now or force_save:
            self.save()


class AbstractIntervalSaveRoutine(AbstractSaveRoutine, ABC):
    """Abstract class describing a save routine that saves every given interval."""

    def __init__(self, interval: u.Quantity, **kwargs):
        super().__init__(**kwargs)
        self.save_interval = interval.to(u.s).value
        self.time_of_last_save = 0

    @property
    def require_synchronized_dt(self):
        """Save output only makes sense for synchronized time steps."""
        return True

    @property
    def save_now(self):
        """Save at every interval given in instantiation."""

        if self._particle_tracker.time - self.time_of_last_save >= self.save_interval:
            self.time_of_last_save = self._particle_tracker.time

            return True
        else:
            return False


class DiskIntervalSaveRoutine(AbstractIntervalSaveRoutine):
    """Save routine corresponding to saving a hdf5 file every given interval."""

    def __init__(self, interval: u.Quantity, output_directory: Path):
        super().__init__(interval, output_directory=output_directory)


class MemoryIntervalSaveRoutine(AbstractIntervalSaveRoutine):
    """Save the state of the tracker every given interval."""

    def __init__(self, interval: u.Quantity):
        super().__init__(interval)


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

        self._is_adaptive_time_step = False
        self._is_synchronized_time_step = False

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
        if not self._is_adaptive_time_step:
            return self.dt

        # candidate time steps includes one per grid (based on the grid resolution)
        # plus additional candidates based on the field at each particle
        candidates = np.ones([self.nparticles, self.num_grids + 1]) * np.inf

        # Compute the time step indicated by the grid resolution
        ds = np.array([grid.grid_resolution.to(u.m).value for grid in self.grids])
        gridstep = 0.5 * (ds / self.vmax)

        # Wherever a particle is on a grid, include that grid's grid step
        # in the list of candidate time steps
        for i, _grid in enumerate(self.grids):  # noqa: B007
            candidates[:, i] = np.where(self.on_grid[:, i] > 0, gridstep[i], np.inf)

        # If not, compute a number of possible time steps
        # Compute the cyclotron gyroperiod
        Bmag = np.max(np.sqrt(Bx**2 + By**2 + Bz**2)).to(u.T).value
        # Compute the gyroperiod
        if Bmag == 0:
            gyroperiod = np.inf
        else:
            gyroperiod = 2 * np.pi * self.m / (self.q * np.max(Bmag))

        candidates[:, self.num_grids] = gyroperiod / 12

        # TODO: introduce a minimum time step based on electric fields too!

        # Enforce limits on dt
        candidates = np.clip(candidates, self.dt_range[0], self.dt_range[1])

        if not self._is_synchronized_time_step:
            # dt is the min of all the candidates for each particle
            # a separate dt is returned for each particle
            dt = np.min(candidates, axis=-1)

            # dt should never actually be infinite, so replace any infinities
            # with the largest gridstep
            dt[dt == np.inf] = np.max(gridstep)
        else:
            # a single value for dt is returned
            # this is the time step used for all particles
            dt = np.min(candidates)

        return dt

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

        # Calculate the adaptive time step from the fields currently experienced
        # by the particles
        # If user sets dt explicitly, that's handled in _adaptive_dt
        dt = self._adaptive_dt(Ex, Ey, Ez, Bx, By, Bz)

        # TODO: Test v/c and implement relativistic Boris push when required
        # vc = np.max(v)/_c

        # Make sure the time step can be multiplied by a [nparticles, 3] shape field array
        if dt.size > 1:
            dt = dt[tracked_mask, np.newaxis]

            # Increment the tracked particles' time by dt
            self.time[tracked_mask] += dt
        else:
            self.time += dt

        self.x[tracked_mask], self.v[tracked_mask] = boris_push(
            pos_tracked, vel_tracked, B, E, self.q, self.m, dt, inplace=False
        )

        self.dt = dt

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
        Used for determining the grid crossing maximum time step.
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

        require_synchronized_time = (
            stop_condition.require_synchronized_dt or save_routine.require_synchronized_dt
        )

        # Will the simulation follow a synchronized time step?
        # This must be the case for certain stopping conditions and saving routines
        if isinstance(dt, u.Quantity):
            if isinstance(dt.value, np.ndarray):
                # If an array is specified for the time step, a synchronized time step is implied if all
                # the entries are equal
                self._is_synchronized_time_step = np.all(dt.value[0] == dt.value[:])
            else:
                self._is_synchronized_time_step = True
        elif dt is None:
            self._is_synchronized_time_step = True
            self._is_adaptive_time_step = True

        # Raise a ValueError if a synchronized dt is required by stop condition or save routine but one is not given
        # This is only the case if an array with differing entries is specified for dt
        if require_synchronized_time and not self._is_synchronized_time_step:
            raise ValueError(
                "Please specify a synchronized time step to use the simulation with this configuration!"
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

    @property
    def is_adaptive_time_step(self):
        """Return whether the simulation is calculating an adaptive time step or using the user-provided time step."""
        return self._is_adaptive_time_step

    @property
    def is_synchronized_time_step(self):
        """Return whether or not the simulation is applying the same time step across all particles."""
        return self._is_synchronized_time_step

    def run(
        self,
        stop_condition: AbstractStopCondition,
        save_routine: Optional[AbstractSaveRoutine] = None,
        dt=None,
        dt_range=None,
        field_weighting="volume averaged",
    ):
        r"""
        Runs a particle-tracing simulation.
        Time steps are adaptively calculated based on the
        local grid resolution of the particles and the electric and magnetic
        fields they are experiencing.

        Parameters
        ----------
        stop_condition : `~plasmapy.simulation.particle_tracker.AbstractStopCondition`
            A class responsible for determining when the simulation has reached
            a suitable termination condition. The class is passed an instance
            of `~plasmapy.simulation.particle_tracker.ParticleTracker`. See `~plasmapy.simulation.particle_tracker.AbstractStopCondition`
            for more details.

        save_routine: `~plasmapy.simulation.particle_tracker.AbstractSaveRoutine`
            A class responsible for determining when the state of the simulation
            is saved. See `~plasmapy.simulation.particle_tracker.AbstractSaveRoutine` for more details.

        dt : `~astropy.units.Quantity`, optional
            An explicitly set time step in units convertible to seconds.
            Setting this optional keyword overrules the adaptive time step
            capability and forces the use of this time step throughout. If a tuple
            of time steps is provided, the adaptive time step will be clamped
            between the first and second values.

        field_weighting : str
            String that selects the field weighting algorithm used to determine
            what fields are felt by the particles. Options are:

            * 'nearest neighbor': Particles are assigned the fields on
                the grid vertex closest to them.
            * 'volume averaged' : The fields experienced by a particle are a
                volume-average of the eight grid points surrounding them.

            The default is 'volume averaged'.

        Returns
        -------
        None

        """

        dt_range = [0, np.inf] * u.s if dt_range is None else dt_range
        self.dt_range = dt_range.to(u.s).value

        # Validate inputs to the run function
        # Sets is_synchronized_time property
        self._validate_run_inputs(stop_condition, save_routine, dt, field_weighting)
        self._enforce_particle_creation()

        self.dt = dt.to(u.s).value if dt is not None else None

        # Keep track of how many push steps have occurred for trajectory tracing
        self.iteration_number = 0

        self.time = np.zeros((self.nparticles, 1)) if not self.is_synchronized_time_step else 0
        # Create flags for tracking when particles during the simulation
        # on_grid -> zero if the particle is off grid, 1
        # shape [nparticles, ngrids]
        self.on_grid = np.zeros([self.nparticles, self.num_grids])

        # Entered grid -> non-zero if particle EVER entered a grid
        self.entered_grid = np.zeros([self.nparticles])

        # Update the `tracker` attribute so that the stop condition & save routine can be used
        stop_condition.tracker = self
        save_routine.tracker = self

        # Initialize a "progress bar" (really more of a meter)
        # Setting sys.stdout lets this play nicely with regular print()
        pbar = tqdm(
            initial=0,
            total=stop_condition.total,
            disable=not self.verbose,
            desc=stop_condition.progress_description,
            unit=stop_condition.units,
            bar_format="{l_bar}{bar}{n:.1e}/{total:.1e} {unit}",  # noqa: FS003
            file=sys.stdout,
        )

        # Push the particles until the stop condition is satisfied
        is_finished = False
        while not (is_finished or self.nparticles_tracked == 0):
            is_finished = stop_condition.is_finished
            progress = min(stop_condition.progress, stop_condition.total)

            pbar.n = progress
            pbar.last_print_n = progress
            pbar.update(0)

            self._push()

            if save_routine is not None:
                save_routine.post_push_hook()

        save_routine.post_push_hook(force_save=True)
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
