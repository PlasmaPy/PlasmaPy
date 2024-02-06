"""
Module containing the definition for the general particle tracker.
"""

__all__ = [
    "AbstractSaveRoutine",
    "AbstractTerminationCondition",
    "DoNotSaveSaveRoutine",
    "IntervalSaveRoutine",
    "NoParticlesOnGridsTerminationCondition",
    "ParticleTracker",
    "TimeElapsedTerminationCondition",
]

import collections
import sys
import warnings
from abc import ABC, abstractmethod
from collections.abc import Iterable
from pathlib import Path
from typing import Optional, Union

import astropy.units as u
import h5py
import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from plasmapy.particles import Particle, particle_input
from plasmapy.plasma.grids import AbstractGrid
from plasmapy.plasma.plasma_base import BasePlasma
from plasmapy.simulation.particle_integrators import boris_push


class AbstractTerminationCondition(ABC):
    """Abstract base class containing the necessary methods for a ParticleTracker termination condition."""

    @property
    def tracker(self) -> Optional["ParticleTracker"]:
        """Return the `ParticleTracker` object for this termination condition."""
        return self._particle_tracker

    @tracker.setter
    def tracker(self, particle_tracker: "ParticleTracker") -> None:
        self._particle_tracker: ParticleTracker = particle_tracker

    @property
    @abstractmethod
    def require_synchronized_dt(self) -> bool:
        """Return if this termination condition requires a synchronized time step."""
        ...

    @property
    @abstractmethod
    def progress_description(self) -> str:
        """Return a small string describing the relevant quantity shown on the meter."""
        ...

    @property
    @abstractmethod
    def units_string(self) -> str:
        """Return a string representation of the units of `total`."""
        ...

    @property
    @abstractmethod
    def is_finished(self) -> bool:
        """Return `True` if the simulation has finished."""
        ...

    @property
    @abstractmethod
    def progress(self) -> float:
        """Return a number representing the progress of the simulation (compared to total).
        This number represents the numerator of the completion percentage.
        """
        ...

    @property
    @abstractmethod
    def total(self) -> float:
        """Return a number representing the total number of steps in a simulation.
        This number represents the denominator of the completion percentage.
        """
        ...


class TimeElapsedTerminationCondition(AbstractTerminationCondition):
    """Termination condition corresponding to the elapsed time of a ParticleTracker."""

    def __init__(self, termination_time: u.Quantity) -> None:
        self.termination_time: float = termination_time.to(u.s).value

    @property
    def require_synchronized_dt(self) -> bool:
        """The elapsed time termination condition requires a synchronized step
        to stop all particles at the same time.
        """
        return True

    @property
    def progress_description(self) -> str:
        """The time elapsed termination condition depends on elapsed time,
        therefore the relevant quantity is time remaining.
        """
        return "Time remaining"

    @property
    def units_string(self) -> str:
        """The units for the time elapsed condition have the units of seconds."""

        return "seconds"

    @property
    def is_finished(self) -> bool:
        """Conclude the simulation if all particles have been tracked over the specified termination time."""

        return bool(float(self.tracker.time) >= self.termination_time)

    @property
    def progress(self) -> float:
        """Return the current time step of the simulation."""
        return float(self.tracker.time)

    @property
    def total(self) -> float:
        """Return the total amount of time over which the particles are tracked."""
        return self.termination_time


class NoParticlesOnGridsTerminationCondition(AbstractTerminationCondition):
    """Termination condition corresponding to stopping the simulation when all particles have exited the grid."""

    def __init__(self) -> None:
        pass

    @property
    def require_synchronized_dt(self) -> bool:
        """The no field termination condition does not require a synchronized time step."""
        return False

    @property
    def progress_description(self) -> str:
        """The progress meter is described in terms of the fraction of particles still on the grid."""
        return "Number of particles still on grid"

    @property
    def units_string(self) -> str:
        """The progress meter is described in terms of the fraction of particles still on the grid."""
        return "Particles"

    @property
    def is_finished(self) -> bool:
        """The simulation is finished when no more particles are on any grids."""
        is_not_on_grid = self.tracker.on_any_grid == False  # noqa: E712

        return is_not_on_grid.all() and self.tracker.iteration_number > 0

    @property
    def progress(self) -> float:
        """The progress of the simulation is measured by how many particles are no longer on a grid."""
        is_not_on_grid: NDArray[np.bool_] = self.tracker.on_grid == 0

        return float(is_not_on_grid.sum())

    @property
    def total(self) -> float:
        """The progress of the simulation is measured against the total number of particles in the simulation."""
        return float(self.tracker.nparticles)


class AbstractSaveRoutine(ABC):
    """Abstract base class containing the necessary methods for a
    |ParticleTracker| save routine.

    The save routine class is responsible for defining the conditions and hooks
    for saving.

    Parameters
    ----------
    output_directory : `~pathlib.Path`, optional
        Output for objects that are saved to disk. If a directory is not specified
        then a memory save routine is used.


    Notes
    -----
    After every push, the `post_push_hook` method is called with the
    respective `~plasmapy.simulation.particle_tracker.ParticleTracker` object passed as a parameter.
    Then, the hook calls `save_now` to determine whether or not the simulation state should be saved.
    """

    def __init__(self, output_directory: Optional[Path] = None) -> None:
        self.output_directory = output_directory

        self.x_all = []
        self.v_all = []

        self._particle_tracker: Optional[ParticleTracker] = None

    @property
    def tracker(self) -> Optional["ParticleTracker"]:
        """Return the `ParticleTracker` object for this stop condition."""
        return self._particle_tracker

    @tracker.setter
    def tracker(self, particle_tracker: "ParticleTracker") -> None:
        self._particle_tracker = particle_tracker

    @property
    @abstractmethod
    def require_synchronized_dt(self) -> bool:
        """Return if this save routine requires a synchronized time step."""
        ...

    @property
    @abstractmethod
    def save_now(self) -> bool:
        """Determine if to save on the current push step."""
        ...

    def save(self) -> None:
        """Save the current state of the simulation to disk or memory based on whether the output directory was set."""

        if self.output_directory is not None:
            self._save_to_disk()
        else:
            self._save_to_memory()

    def _save_to_disk(self) -> None:
        """Save a hdf5 file containing simulation positions and velocities."""

        path = self.output_directory / f"{self.tracker.iteration_number}.hdf5"

        with h5py.File(path, "w") as output_file:
            output_file["x"] = self.tracker.x
            output_file["v"] = self.tracker.v

    def _save_to_memory(self) -> None:
        """Append simulation positions and velocities to save routine object."""
        self.x_all.append(np.copy(self._particle_tracker.x))
        self.v_all.append(np.copy(self._particle_tracker.v))

    def post_push_hook(self, force_save=False) -> None:
        """Function called after a push step.

        This function is responsible for handling two steps of save routine, namely:
            - Deciding to save on the current time step
            - How the simulation data is saved (i.e. to disk or memory)

        """
        if self.save_now or force_save:
            self.save()


class DoNotSaveSaveRoutine(AbstractSaveRoutine):
    """The default save routine for the |ParticleTracker| class.

    This save routine is a placeholder and will not save the state of the particle tracker.
    """

    def __init__(self) -> None:
        super().__init__()

    @property
    def require_synchronized_dt(self) -> bool:
        """The do not save save routine does not require a synchronized time step."""
        return False

    @property
    def save_now(self) -> bool:
        """The do not save save routine will never save by definition."""
        return False


class IntervalSaveRoutine(AbstractSaveRoutine):
    """Abstract class describing a save routine that saves every given interval."""

    def __init__(self, interval: u.Quantity, **kwargs) -> None:
        super().__init__(**kwargs)
        self.t_all: list[float] = []

        self.save_interval: float = interval.to(u.s).value
        self.time_of_last_save: float = 0

    @property
    def require_synchronized_dt(self) -> bool:
        """Save output only makes sense for synchronized time steps."""
        return True

    @property
    def save_now(self) -> bool:
        """Save at every interval given in instantiation."""

        return bool(self.tracker.time - self.time_of_last_save >= self.save_interval)

    def save(self) -> None:
        """Save the current state of the simulation.
        Sets the time of last save attribute and log the timestamp.
        """
        super().save()

        self.time_of_last_save = self.tracker.time
        self.t_all.append(self.tracker.time)

    def results(self) -> tuple[u.Quantity, u.Quantity, u.Quantity]:
        """Return the results of the simulation.
        The quantities returned are the times, positions, and velocities, respectively.
        """

        return (
            np.asarray(self.t_all) * u.s,
            np.asarray(self.x_all) * u.m,
            np.asarray(self.v_all) * u.m / u.s,
        )


class ParticleTracker:
    r"""A particle tracker for particles in electric and magnetic fields without inter-particle interactions.

    Particles are instantiated and pushed through a grid of provided E and
    B fields using the Boris push algorithm. These fields are specified as
    part of a grid which are then interpolated to determine the local field
    acting on each particle.

    The time step used in the push routine can be specified, or an adaptive
    time step will be determined based off the gyroradius of the particle.
    Some save routines involve time stamping the location and velocities of
    particles at fixed intervals. In order for this data to be coherent, it
    is required that all particles follow the same time step. This is
    referred to as a synchronized time step. If no time step is specified
    and the provided save routine does not require a synchronized time step,
    then an adaptive time step is calculated independently for each particle.

    The simulation will push particles through the provided fields until a
    condition is met. This termination condition is provided as an instance
    of the `~plasmapy.simulation.particle_tracker.AbstractTerminationCondition`
    class as arguments to the simulation constructor. The results of a simulation
    can be exported by specifying an instance of the `~plasmapy.simulation.particle_tracker.AbstractSaveRoutine`
    class to the ``run`` method.

    Parameters
    ----------
    grids : An instance of `~plasmapy.plasma.grids.AbstractGrid`
        A Grid object or list of grid objects containing the required quantities.
        The list of required quantities varies depending on other keywords.

    termination_condition : `~plasmapy.simulation.particle_tracker.AbstractTerminationCondition`
        An instance of `~plasmapy.simulation.particle_tracker.AbstractTerminationCondition` which determines when the simulation has finished.
        See `~plasmapy.simulation.particle_tracker.AbstractTerminationCondition` for more details.

    save_routine : `~plasmapy.simulation.particle_tracker.AbstractSaveRoutine`, optional
        An instance of `~plasmapy.simulation.particle_tracker.AbstractSaveRoutine` which determines which
        time steps of the simulation to save. The default is `~plasmapy.simulation.particle_tracker.DoNotSaveSaveRoutine`.
        See `~plasmapy.simulation.particle_tracker.AbstractSaveRoutine` for more details.

    dt : `~astropy.units.Quantity`, optional
        An explicitly set time step in units convertible to seconds.
        Setting this optional keyword overrules the adaptive time step
        capability and forces the use of this time step throughout.

    dt_range : tuple of shape (2,) of `~astropy.units.Quantity`, optional
        If specified, the calculated adaptive time step will be clamped
        between the first and second values.

    field_weighting : str
        String that selects the field weighting algorithm used to determine
        what fields are felt by the particles. Options are:

        * 'nearest neighbor': Particles are assigned the fields on
            the grid vertex closest to them.
        * 'volume averaged' : The fields experienced by a particle are a
            volume-average of the eight grid points surrounding them.

        The default is 'volume averaged'.

    req_quantities : `list` of str, optional
        A list of quantity keys required to be specified on the Grid object.
        The base particle pushing simulation requires the quantities [E_x, E_y, E_z, B_x, B_y, B_z].
        If any of these quantities are missing, a warning will be given and that
        quantity will be assumed to be zero everywhere. This keyword is for specifying quantities in
        addition to these six. The default is None.

    verbose : bool, optional
        If true, updates on the status of the program will be printed
        into the standard output while running. The default is True.

    Notes
    -----
    We adopt the convention of ``NaN`` values to represent various states of a particle.

    If the particle's position and velocity are not ``NaN``, the particle is being tracked and evolved.
    If the particle's position is not ``NaN``, but the velocity is ``NaN``, the particle has been stopped
    (i.e. it is still in the simulation but is no longer evolved.)
    If both the particle's position and velocity are set to ``NaN``, then the particle has been removed from the simulation.
    """

    def __init__(
        self,
        grids: Union[AbstractGrid, Iterable[AbstractGrid]],
        termination_condition: Union[AbstractTerminationCondition, None] = None,
        save_routine: Union[AbstractSaveRoutine, None] = None,
        dt=None,
        dt_range=None,
        field_weighting="volume averaged",
        req_quantities=None,
        verbose=True,
    ) -> None:
        # self.grid is the grid object
        if isinstance(grids, AbstractGrid):
            self.grids = [
                grids,
            ]
        elif isinstance(grids, collections.abc.Iterable):
            self.grids = grids
        else:
            self.grids = None

        # Errors for unsupported grid types are raised in the validate constructor inputs method

        # Instantiate the "do not save" save routine if no save routine was specified
        if save_routine is None:
            save_routine = DoNotSaveSaveRoutine()

        # Validate inputs to the run function
        self._validate_constructor_inputs(
            grids, termination_condition, save_routine, field_weighting
        )

        self._set_time_step_attributes(dt, termination_condition, save_routine)

        if dt_range is not None and not self._is_adaptive_time_step:
            raise ValueError(
                "Specifying a time step range is only possible for an adaptive time step."
            )

        self.verbose = verbose

        # This flag records whether the simulation has been run
        self._has_run = False

        # Raise a ValueError if a synchronized dt is required by termination condition or save routine but one is
        # not given. This is only the case if an array with differing entries is specified for dt
        if self._require_synchronized_time and not self._is_synchronized_time_step:
            raise ValueError(
                "Please specify a synchronized time step to use the simulation with this configuration!"
            )

        self._preprocess_grids(req_quantities)

        # self.grid_arr is the grid positions in si units. This is created here
        # so that it isn't continuously called later
        self.grids_arr = [grid.grid.to(u.m).value for grid in self.grids]

        self.dt = dt.to(u.s).value if dt is not None else None

        dt_range = [0, np.inf] * u.s if dt_range is None else dt_range
        self.dt_range = dt_range.to(u.s).value

        # Update the `tracker` attribute so that the stop condition & save routine can be used
        termination_condition.tracker = self

        save_routine.tracker = self

        self.termination_condition = termination_condition
        self.save_routine = save_routine

    def _set_time_step_attributes(
        self, dt, termination_condition, save_routine
    ) -> None:
        """Determines whether the simulation will follow a synchronized or adaptive time step."""

        self._require_synchronized_time = (
            termination_condition.require_synchronized_dt
            or (save_routine is not None and save_routine.require_synchronized_dt)
        )

        if isinstance(dt, u.Quantity):
            if isinstance(dt.value, np.ndarray):
                # If an array is specified for the time step, a synchronized time step is implied if all
                # the entries are equal
                self._is_synchronized_time_step = bool(
                    np.all(dt.value[0] == dt.value[:])
                )
            else:
                self._is_synchronized_time_step = True

            self._is_adaptive_time_step = False
        elif dt is None:
            self._is_synchronized_time_step = self._require_synchronized_time
            self._is_adaptive_time_step = True

        if self._is_adaptive_time_step:
            # Initialize default values for time steps per gyroperiod and Courant parameter
            self.setup_adaptive_time_step()

    def setup_adaptive_time_step(
        self,
        time_steps_per_gyroperiod: Optional[int] = 12,
        Courant_parameter: Optional[float] = 0.5,
    ) -> None:
        """Set parameters for the adaptive time step candidates.

        Parameters
        ----------
        time_steps_per_gyroperiod : int, optional
            The minimum ratio of the particle gyroperiod to the timestep. Higher numbers
            correspond to higher temporal resolution. The default is twelve.

        Courant_parameter : float, optional
            The Courant parameter is the minimum ratio of the timestep to the grid crossing time,
            grid cell length / particle velocity. Lower Courant numbers correspond to higher temporal resolution.


        Notes
        -----
        Two candidates are calculated for the adaptive time step: a time step based on the gyroradius
        of the particle and a time step based on the resolution of the grid. The candidate associated
        with the gyroradius of the particle takes a ``time_steps_per_gyroperiod`` parameter that specifies
        how many times the orbit of a gyrating particles will be subdivided. The other candidate,
        associated with the spatial resolution of the grid object, calculates a time step using the time
        it would take the fastest particle to cross some fraction of a grid cell length. This fraction is the Courant number.
        """

        if not self._is_adaptive_time_step:
            raise ValueError(
                "The setup adaptive time step method only applies to adaptive time steps!"
            )

        self._steps_per_gyroperiod = time_steps_per_gyroperiod
        self._Courant_parameter = Courant_parameter

    def _validate_constructor_inputs(
        self, grids, termination_condition, save_routine, field_weighting: str
    ) -> None:
        """
        Ensure the specified termination condition and save routine are actually
        a termination routine class and save routine, respectively. This function also
        sets the `_is_synchronized_time_step` and `_is_adaptive_time_step` attributes.
        """

        if isinstance(grids, BasePlasma):
            raise TypeError(
                "It appears you may be trying to access an older version of the ParticleTracker class."
                "This class has been deprecated."
                "Please revert to PlasmaPy version 2023.5.1 to use this version of ParticleTracker."
            )
        # The constructor did not recognize the provided grid object
        elif self.grids is None:
            raise TypeError("Type of argument `grids` not recognized.")

        if not isinstance(termination_condition, AbstractTerminationCondition):
            raise TypeError("Please specify a valid termination condition.")

        if not isinstance(save_routine, AbstractSaveRoutine):
            raise TypeError("Please specify a valid save routine.")

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

    def _preprocess_grids(self, additional_required_quantities) -> None:
        """Add required quantities to grid objects.

        Grids lacking the required quantities will be filled with zeros.
        """

        # Some quantities are necessary for the particle tracker to function regardless of other configurations
        required_quantities = {"E_x", "E_y", "E_z", "B_x", "B_y", "B_z"}

        if additional_required_quantities is not None:
            # Add additional required quantities based off simulation configuration
            required_quantities.update(additional_required_quantities)

        for grid in self.grids:
            grid.require_quantities(required_quantities, replace_with_zeros=True)

            for rq in required_quantities:
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
                        "Quantities should go to zero at edges of grid to avoid "
                        f"non-physical effects, but a value of {edge_max:.2E} {unit} was "
                        f"found on the edge of the {rq} array. Consider applying a "
                        "envelope function to force the quantities at the edge to go to "
                        "zero.",
                        RuntimeWarning,
                    )

    @property
    def num_grids(self) -> int:
        """The number of grids specified at instantiation."""
        return len(self.grids)

    def _log(self, msg) -> None:
        if self.verbose:
            print(msg)  # noqa: T201

    @particle_input
    def load_particles(
        self,
        x,
        v,
        particle: Particle = Particle("p+"),  # noqa: B008
    ) -> None:
        r"""
        Load arrays of particle positions and velocities.

        Parameters
        ----------
        x : `~astropy.units.Quantity`, shape (N,3)
            Positions for N particles

        v : `~astropy.units.Quantity`, shape (N,3)
            Velocities for N particles

        particle : |particle-like|, optional
            Representation of the particle species as either a |Particle| object
            or a string representation. The default particle is protons.
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
            self.nparticles: int = x.shape[0]

        self.x = x.to(u.m).value
        self.v = v.to(u.m / u.s).value

    def run(self) -> None:
        r"""
        Runs a particle-tracing simulation.
        Time steps are adaptively calculated based on the
        local grid resolution of the particles and the electric and magnetic
        fields they are experiencing.

        Returns
        -------
        None

        """

        self._enforce_particle_creation()

        # Keep track of how many push steps have occurred for trajectory tracing
        self.iteration_number = 0

        self.time: Union[NDArray[np.float64], float] = (
            np.zeros((self.nparticles, 1)) if not self.is_synchronized_time_step else 0
        )
        # Create flags for tracking when particles during the simulation
        # on_grid -> zero if the particle is off grid, 1
        # shape [nparticles, ngrids]
        self.on_grid: NDArray[np.bool_] = np.zeros(
            [self.nparticles, self.num_grids]
        ).astype(np.bool_)

        # Entered grid -> non-zero if particle EVER entered a grid
        self.entered_grid: NDArray[np.bool_] = np.zeros([self.nparticles]).astype(
            np.bool_
        )

        # Initialize a "progress bar" (really more of a meter)
        # Setting sys.stdout lets this play nicely with regular print()
        pbar = tqdm(
            initial=0,
            total=self.termination_condition.total,
            disable=not self.verbose,
            desc=self.termination_condition.progress_description,
            unit=self.termination_condition.units_string,
            bar_format="{l_bar}{bar}{n:.1e}/{total:.1e} {unit}",
            file=sys.stdout,
        )

        # Push the particles until the termination condition is satisfied
        is_finished = False
        while not (is_finished or self.nparticles_tracked == 0):
            is_finished = self.termination_condition.is_finished
            progress = min(
                self.termination_condition.progress, self.termination_condition.total
            )

            pbar.n = progress
            pbar.last_print_n = progress
            pbar.update(0)

            self._push()

            if self.save_routine is not None:
                self.save_routine.post_push_hook()

        if self.save_routine is not None:
            self.save_routine.post_push_hook(force_save=True)

        pbar.close()

        # Log a summary of the run

        self._log("Run completed")

        # Simulation has not run, because creating new particles changes the simulation
        self._has_run = True

    def _stop_particles(self, particles_to_stop_mask) -> None:
        """Stop tracking the particles specified by the stop mask.

        This is represented by setting the particle's velocity to NaN.
        """

        if len(particles_to_stop_mask) != self.x.shape[0]:
            raise ValueError(
                f"Expected mask of size {self.x.shape[0]}, got {len(particles_to_stop_mask)}"
            )

        self.v[particles_to_stop_mask] = np.nan

    def _remove_particles(self, particles_to_remove_mask) -> None:
        """Remove the specified particles from the simulation.

        For the sake of keeping consistent array lengths, the position and velocities of the
        removed particles are set to NaN.
        """

        if len(particles_to_remove_mask) != self.x.shape[0]:
            raise ValueError(
                f"Expected mask of size {self.x.shape[0]}, got {len(particles_to_remove_mask)}"
            )

        self.x[particles_to_remove_mask] = np.nan
        self.v[particles_to_remove_mask] = np.nan

    # *************************************************************************
    # Run/push loop methods
    # *************************************************************************

    def _adaptive_dt(self, Ex, Ey, Ez, Bx, By, Bz) -> Union[NDArray[np.float64], float]:  # noqa: ARG002
        r"""
        Calculate the appropriate dt for each grid based on a number of
        considerations
        including the local grid resolution (ds) and the gyroperiod of the
        particles in the current fields.
        """

        # candidate time steps includes one per grid (based on the grid resolution)
        # plus additional candidates based on the field at each particle
        candidates = np.ones([self.nparticles, self.num_grids + 1]) * np.inf

        # Compute the time step indicated by the grid resolution
        ds = np.array([grid.grid_resolution.to(u.m).value for grid in self.grids])
        gridstep = self._Courant_parameter * (ds / self.vmax)

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

        # Subdivide the gyroperiod into a provided number of steps
        # Use the result as the candidate associated with gyration in B field
        candidates[:, self.num_grids] = gyroperiod / self._steps_per_gyroperiod

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

    def _push(self) -> None:
        r"""
        Advance particles using an implementation of the time-centered
        Boris algorithm.
        """
        # Get a list of positions (input for interpolator)
        tracked_mask = self._tracked_particle_mask()

        self.iteration_number += 1

        # nparticles_tracked may fluctuate with stopping
        # all particles are therefore used regardless of if they reach the grid
        pos_all = self.x
        pos_tracked = pos_all[tracked_mask]

        vel_all = self.v
        vel_tracked = vel_all[tracked_mask]

        # Update the list of particles on and off the grid
        # shape [nparticles, ngrids]
        self.on_grid = (
            np.array([grid.on_grid(pos_all * u.m) for grid in self.grids])
            .astype(np.bool_)
            .T
        )

        # entered_grid is zero at the end if a particle has never
        # entered any grid
        self.entered_grid += np.sum(self.on_grid, axis=-1).astype(np.bool_)

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
        if self._is_adaptive_time_step:
            dt = self._adaptive_dt(Ex, Ey, Ez, Bx, By, Bz)
        else:
            dt = self.dt

        # TODO: Test v/c and implement relativistic Boris push when required
        # vc = np.max(v)/_c

        # Make sure the time step can be multiplied by a [nparticles, 3] shape field array
        if isinstance(dt, np.ndarray) and dt.size > 1:
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
    def on_any_grid(self) -> NDArray[np.bool_]:
        """
        Binary array for each particle indicating whether it is currently
        on ANY grid.
        """
        return np.sum(self.on_grid, axis=-1) > 0

    @property
    def vmax(self) -> float:
        """The maximum velocity of any particle in the simulation.

        This quantity is used for determining the grid crossing maximum time step.
        """
        tracked_mask = self._tracked_particle_mask()

        return float(np.max(np.linalg.norm(self.v[tracked_mask], axis=-1)))

    def _tracked_particle_mask(self) -> NDArray[np.bool_]:
        """
        Calculates a boolean mask corresponding to particles that have not been stopped or removed.
        """
        return ~np.logical_or(np.isnan(self.x[:, 0]), np.isnan(self.v[:, 0]))

    @property
    def nparticles_tracked(self) -> int:
        """Return the number of particles that don't have NaN position or velocity."""
        return int(self._tracked_particle_mask().sum())

    @property
    def is_adaptive_time_step(self) -> bool:
        """Return whether the simulation is calculating an adaptive time step or using the user-provided time step."""
        return self._is_adaptive_time_step

    @property
    def is_synchronized_time_step(self) -> bool:
        """Return if the simulation is applying the same time step across all particles."""
        return self._is_synchronized_time_step

    def _enforce_particle_creation(self) -> None:
        """Ensure the array position array `x` has been populated."""

        # Check to make sure particles have already been generated
        if not hasattr(self, "x"):
            raise ValueError(
                "Either the create_particles or load_particles method must be "
                "called before running the particle tracing algorithm."
            )

    def _enforce_order(self) -> None:
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
