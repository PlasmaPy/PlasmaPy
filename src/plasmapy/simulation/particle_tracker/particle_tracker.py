"""
Module containing the definition for the general particle tracker.
"""

__all__ = [
    "ParticleTracker",
]

import collections
import sys
import warnings
from collections.abc import Iterable
from functools import cached_property
from typing import Literal

import astropy.constants as const
import astropy.units as u
import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from plasmapy.formulary.collisions.misc import Bethe_stopping_lite
from plasmapy.particles import Particle, particle_input
from plasmapy.particles.atomic import stopping_power
from plasmapy.plasma.grids import AbstractGrid
from plasmapy.plasma.plasma_base import BasePlasma
from plasmapy.simulation.particle_integrators import (
    AbstractIntegrator,
    RelativisticBorisIntegrator,
)
from plasmapy.simulation.particle_tracker.save_routines import (
    AbstractSaveRoutine,
    DoNotSaveSaveRoutine,
)
from plasmapy.simulation.particle_tracker.termination_conditions import (
    AbstractTerminationCondition,
)
from plasmapy.utils.exceptions import PhysicsWarning, RelativityWarning

_c = const.c
_m_p = const.m_p


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
    of the `~plasmapy.simulation.particle_tracker.termination_conditions.AbstractTerminationCondition`
    class as arguments to the simulation constructor. The results of a simulation
    can be exported by specifying an instance of the `~plasmapy.simulation.particle_tracker.save_routines.AbstractSaveRoutine`
    class to the ``run`` method.

    Parameters
    ----------
    grids : An instance of `~plasmapy.plasma.grids.AbstractGrid`
        A Grid object or list of grid objects containing the required quantities.
        The list of required quantities varies depending on other keywords.

    termination_condition : `~plasmapy.simulation.particle_tracker.termination_conditions.AbstractTerminationCondition`
        An subclass of `~plasmapy.simulation.particle_tracker.termination_conditions.AbstractTerminationCondition` which determines when the simulation has finished.
        See `~plasmapy.simulation.particle_tracker.termination_conditions.AbstractTerminationCondition` for more details.

    save_routine : `~plasmapy.simulation.particle_tracker.save_routines.AbstractSaveRoutine`, optional
        An subclass of `~plasmapy.simulation.particle_tracker.save_routines.AbstractSaveRoutine` which determines which
        time steps of the simulation to save. The default is `~plasmapy.simulation.particle_tracker.save_routines.DoNotSaveSaveRoutine`.
        See `~plasmapy.simulation.particle_tracker.save_routines.AbstractSaveRoutine` for more details.

    particle_integrator : `~plasmapy.simulation.particle_integrators.AbstractIntegrator`, optional
        An subclass of `~plasmapy.simulation.particle_integrators.AbstractIntegrator` that is responsible for implementing the push behavior
        of the simulation when provided the electric and magnetic fields. The default value is set to `~plasmapy.simulation.particle_integrators.RelativisticBorisIntegrator`.
        See `~plasmapy.simulation.particle_integrators.AbstractIntegrator` for more information on how to implement custom push routines.

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

    verbose : bool, optional
        If true, updates on the status of the program will be printed
        into the standard output while running. The default is True.

    Warns
    -----
    `~plasmapy.utils.exceptions.RelativityWarning`
        The provided integrator does not account for relativistic corrections
        and particles have reached a non-negligible fraction of the speed of
        light.

    Notes
    -----
    We adopt the convention of ``NaN`` values to represent various states of a particle.

    If the particle's position and velocity are not ``NaN``, the particle is being tracked and evolved.
    If the particle's position is not ``NaN``, but the velocity is ``NaN``, the particle has been stopped
    (i.e. it is still in the simulation but is no longer evolved.)
    If both the particle's position and velocity are set to ``NaN``, then the particle has been removed from the simulation.

    Example
    -----
    >>> from plasmapy.particles import Particle
    >>> from plasmapy.plasma.grids import CartesianGrid
    >>> from plasmapy.simulation.particle_tracker.particle_tracker import (
    ...     ParticleTracker,
    ... )
    >>> from plasmapy.simulation.particle_tracker.save_routines import (
    ...     IntervalSaveRoutine,
    ... )
    >>> from plasmapy.simulation.particle_tracker.termination_conditions import (
    ...     TimeElapsedTerminationCondition,
    ... )
    >>> import astropy.units as u
    >>> import numpy as np
    >>> example_particle = Particle("p+")
    >>> grid = CartesianGrid(-1e6 * u.m, 1e6 * u.m, num=2)
    >>> grid_shape = (2, 2, 2)
    >>> Bz = np.full(grid_shape, 1) * u.T
    >>> grid.add_quantities(B_z=Bz)
    >>> x = [[0, 0, 0]] * u.m
    >>> v = [[1, 0, 0]] * u.m / u.s
    >>> termination_condition = TimeElapsedTerminationCondition(6.28 * u.s)
    >>> save_routine = IntervalSaveRoutine(6.28 * u.s / 10)
    >>> simulation = ParticleTracker(
    ...     grid,
    ...     termination_condition,
    ...     dt=1e-2 * u.s,
    ...     save_routine=save_routine,
    ...     field_weighting="nearest neighbor",
    ...     verbose=False,
    ... )
    >>> simulation.load_particles(x, v, example_particle)
    >>> simulation.run()
    >>> print(
    ...     simulation.save_routine.results["time"][-1],
    ...     simulation.save_routine.results["x"][-1],
    ... )
    6.2999999999... s [[-1.73302...e-08  1.31539...e-05  0.00000...e+00]] m
    """

    def __init__(
        self,
        grids: AbstractGrid | Iterable[AbstractGrid],
        termination_condition: AbstractTerminationCondition | None = None,
        save_routine: AbstractSaveRoutine | None = None,
        particle_integrator: type[AbstractIntegrator] | None = None,
        dt=None,
        dt_range=None,
        field_weighting: str = "volume averaged",
        verbose: bool = True,
    ) -> None:
        # Verbose flag controls whether or not output is printed to stdout
        self.verbose = verbose

        # This flag records whether the simulation has been run
        self._has_run = False

        # Should the tracker update particle energies after every time step to
        # reflect stopping?
        self._do_stopping = False

        # Instantiate the integrator object for use in the _push() method
        self._integrator: AbstractIntegrator = (
            RelativisticBorisIntegrator()
            if not particle_integrator
            else particle_integrator()
        )

        self._raised_relativity_warning = False

        # Quantities required for tracking - if they are not defined on a grid
        # they will be assumed to be zero
        self._required_quantities: list[str] = []

        # self.grid is the grid object
        self.grids = self._grid_factory(grids)

        self._validate_grids()

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

        # Raise a ValueError if a synchronized dt is required by termination condition or save routine but one is
        # not given. This is only the case if an array with differing entries is specified for dt
        if self._require_synchronized_time and not self._is_synchronized_time_step:
            raise ValueError(
                "Please specify a synchronized time step to use the simulation with this configuration!"
            )

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

    @staticmethod
    def _grid_factory(grids):
        """
        Take the user provided argument for grids and convert it into the proper type.
        """

        if isinstance(grids, AbstractGrid):
            return [
                grids,
            ]
        elif isinstance(grids, collections.abc.Iterable):
            return grids
        else:
            return None

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

        sync = (
            "a synchronized" if self._is_synchronized_time_step else "an asynchronous"
        )
        adaptive = "adaptive" if self._is_adaptive_time_step else "fixed-value"
        self._log(f"Tracker setup using {sync} {adaptive} timestep setup")

    def setup_adaptive_time_step(
        self,
        time_steps_per_gyroperiod: int | None = 12,
        Courant_parameter: float | None = 0.5,
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
        a termination routine class and save routine, respectively.
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
        particle: Particle,
    ) -> None:
        r"""
        Load arrays of particle positions and velocities.

        Parameters
        ----------
        x : `~astropy.units.Quantity`, shape (N,3)
            Positions for N particles

        v : `~astropy.units.Quantity`, shape (N,3)
            Velocities for N particles

        particle : |particle-like|
            Representation of the particle species as either a |Particle| object
            or a string representation.
        """
        # Raise an error if the run method has already been called.
        self._enforce_order()

        self.q = particle.charge.to(u.C).value
        self.m = particle.mass.to(u.kg).value
        self._particle = particle

        if self.q != 0:
            self._required_quantities += ["E_x", "E_y", "E_z", "B_x", "B_y", "B_z"]

        if x.shape[0] != v.shape[0]:
            raise ValueError(
                "Provided x and v arrays have inconsistent numbers "
                " of particles "
                f"({x.shape[0]} and {v.shape[0]} respectively)."
            )
        else:
            self.num_particles: int = x.shape[0]

        self.x = x.to(u.m).value
        self.v = v.to(u.m / u.s).value

    def _validate_stopping_inputs(
        self,
        method: Literal["NIST", "Bethe"],
        materials: list[str | None] | None = None,
        I: list[u.Quantity[u.J] | None] | None = None,  # noqa: E741
    ) -> bool:
        r"""
        Validate inputs to the `add_stopping` method. Raises errors if the
        proper keyword arguments are not provided for a given method.
        """

        match method:
            case "NIST":
                if materials is None or len(materials) != self.num_grids:
                    raise ValueError(
                        "Please provide an array of length ngrids for the materials."
                    )

                for i, grid in enumerate(self.grids):
                    if materials[i] is not None and "rho" not in grid.quantities:
                        raise ValueError(
                            f"Material {materials[i]} was provided for stopping on grid {i},"
                            " but quantity ``rho`` is not defined on that grid."
                        )

            case "Bethe":
                if I is None or len(I) != self.num_grids:
                    raise ValueError(
                        "Please provide an array of length ngrids for the mean excitation energy."
                    )

                for i, grid in enumerate(self.grids):
                    if I[i] is not None and "n_e" not in grid.quantities:
                        raise ValueError(
                            f"I={I[i]} was provided for stopping on grid {i},"
                            " but quantity ``n_e`` is not defined on that grid."
                        )

        return True

    def add_stopping(
        self,
        method: Literal["NIST", "Bethe"],
        materials: list[str | None] | None = None,
        I: list[u.Quantity[u.J] | None] | None = None,  # noqa: E741
    ):
        r"""
        Enable particle stopping using experimental stopping powers.

        Interpolation of stopping powers is conducted using data from the NIST
        PSTAR database. This information is combined with the mass density
        quantity provided in the grids to calculate the energy loss over the
        distance travelled during a timestep.

        Parameters
        ----------
        method : str
            Sets the method for calculating particle stopping powers. Valid inputs
            are 'Bethe' and 'NIST'.

        materials : list[str|None]
            A list of materials (one per grid) required when using the NIST
            stopping power table method. Grids with no stopping material
            should be set to None.

        I : list[ `~astropy.units.Quantity` | None]
            The mean excitation energy ``I`` in the Bethe stopping model,
            only required when using the 'Bethe' stopping method. In the
            Bethe model, ``I`` completely describes the material, and is
            approximately equal to

            .. math:: I = 10 eV * A

            where A is the atomic number of the atoms in the material.

        """
        # TODO: Add a reference somewhere to the valid material strings?

        # Check inputs for user error and raise respective exceptions/warnings if
        # necessary.

        self._validate_stopping_inputs(method, materials, I)

        match method:
            case "NIST":
                self._required_quantities += ["rho"]
                stopping_power_interpolators = [
                    stopping_power(self._particle, material, return_interpolator=True)
                    if material is not None
                    else None
                    for material in materials
                ]

            case "Bethe":
                self._required_quantities += ["n_e"]
                self._raised_energy_warning = False

                # These functions are used to represent that the mean excitation energy
                # does not change over space for a given grid.
                def wrapped_Bethe_stopping(I_grid):
                    def inner_Bethe_stopping(v, n_e):
                        return Bethe_stopping_lite(
                            I_grid, n_e, v, self._particle.charge_number
                        )

                    return inner_Bethe_stopping

                stopping_power_interpolators = [
                    wrapped_Bethe_stopping(I_grid.si.value)
                    if I_grid is not None
                    else None
                    for I_grid in I
                ]

            case _:
                raise ValueError(
                    f"Please provide one of 'NIST' or 'Bethe' for the method keyword. (Got: {method})"
                )

        self._do_stopping = True
        self._stopping_method = method
        self._stopping_power_interpolators = stopping_power_interpolators

        self._log(f"Stopping module activated using the {method} method")

    def _validate_grids(self) -> None:
        """Validate input grids."""
        for i, grid in enumerate(self.grids):
            for q in grid.quantities:
                # Check that there are no infinite values
                if not np.isfinite(grid[q].value).all():
                    raise ValueError(
                        f"Input arrays must be finite: {q} contains "
                        "either NaN or infinite values."
                    )

                # Check that the max values on the edges of the arrays are
                # small relative to the maximum values on that grid
                #
                # Array must be dimensionless to re-assemble it into an array
                # of max values like this
                arr = np.abs(grid[q]).value
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
                    unit = grid.recognized_quantities()[q].unit
                    warnings.warn(
                        "Quantities should go to zero at edges of grid to avoid "
                        f"non-physical effects, but a value of {edge_max:.2E} {unit} was "
                        f"found on the edge of the {q} array of grid {i}. Consider applying a "
                        "envelope function to force the quantities at the edge to go to "
                        "zero.",
                        RuntimeWarning,
                    )

    def _setup_for_interpolator(self) -> None:
        """Prepare grids for interpolation.

        This will be called at the beginning of run, since the required
        quantities depends on how the Tracker has been set up.
        """
        # Assemble an array of the resolutions of all of the grids
        self._grid_resolutions = np.array(
            [grid.grid_resolution.to(u.m).value for grid in self.grids]
        )

        # Determine which quantities should be interpolated from each grid
        # This set contains all quantities defined on any grid
        self._interpolated_quantities_any_grid: set[str] = set()
        # This list of sets contains all quantities defined on each grid
        self._interpolated_quantities_per_grid: list[str] = []

        for i, grid in enumerate(self.grids):
            quantities = set(self._required_quantities).intersection(grid.quantities)
            self._interpolated_quantities_per_grid.append(quantities)
            self._interpolated_quantities_any_grid = (
                self._interpolated_quantities_any_grid.union(quantities)
            )
            self._log(f"On grid {i}, interpolating: {quantities}")

        # Construct a dictionary to store the interpolation results
        #  Each quantity is initialized as a zeros array with its respective units
        # Arrays are ``num_particles`` sized so they don't need to be recreated
        # when the number of tracked particles changes
        self._total_grid_values = {
            field_name: np.zeros(self.num_particles)
            * AbstractGrid.recognized_quantities()[field_name].unit
            for field_name in self._interpolated_quantities_any_grid
        }

        # These variables indicate whether any E or B fields exist
        # on any of the grids
        self._E_on_grids = any(
            k in self._total_grid_values for k in ["E_x", "E_y", "E_z"]
        )
        self._B_on_grids = any(
            k in self._total_grid_values for k in ["B_x", "B_y", "B_z"]
        )

        # Make sure the complete set of E, B arrays are available
        for key in ["E_x", "E_y", "E_z"]:
            if key not in self._total_grid_values:
                self._total_grid_values[key] = 0.0 * u.V / u.m
        for key in ["B_x", "B_y", "B_z"]:
            if key not in self._total_grid_values:
                self._total_grid_values[key] = 0.0 * u.T

        # Initialize some variables for interpolator
        self._E = np.zeros((self.num_particles, 3))
        self._B = np.zeros((self.num_particles, 3))

    def run(self) -> None:
        r"""
        Runs a particle-tracing simulation.
        Time steps are adaptively calculated based on the local grid resolution
        of the particles and the electric and magnetic fields they are
        experiencing.

        Returns
        -------
        None

        """

        self._enforce_particle_creation()

        self._setup_for_interpolator()

        # Keep track of how many push steps have occurred for trajectory tracing
        # This number is independent of the current "time" of the simulation
        self.iteration_number = 0

        # The time state of a simulation with synchronized time step can be described
        # by a single number. Otherwise, a time value is required for each particle.
        self.time: NDArray[np.float64] | float = (
            np.zeros((self.num_particles, 1))
            if not self.is_synchronized_time_step
            else 0
        )

        self.on_grid: NDArray[np.bool_] = np.zeros(
            [self.num_particles, self.num_grids]
        ).astype(np.bool_)

        # non-zero if particle EVER entered ANY grid
        self.ever_entered_any_grid: NDArray[np.bool_] = np.zeros(
            [self.num_particles]
        ).astype(np.bool_)

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
        # or the number of particles being evolved is zero
        is_finished = False
        while not (is_finished or self.num_particles_tracked == 0):
            is_finished = self.termination_condition.is_finished
            progress = min(
                self.termination_condition.progress, self.termination_condition.total
            )

            pbar.n = progress
            pbar.last_print_n = progress
            pbar.update(0)

            self._push()

            if self.verbose:
                desc = (
                    f"Iter. {self.iteration_number},  "
                    f"Avg. Pos.=({np.nanmean(self.x[:, 0]):.1e},{np.nanmean(self.x[:, 1]):.1e},{np.nanmean(self.x[:, 2]):.1e}) m, "
                    f"Avg. Vel.=({np.nanmean(self.v[:, 0]):.1e},{np.nanmean(self.v[:, 1]):.1e},{np.nanmean(self.v[:, 2]):.1e}) m/s, "
                    f"dt range = ({np.min(self.dt):.1e},{np.max(self.dt):.1e}) s, "
                    f"{self.termination_condition.progress_description}"
                )
                pbar.set_description(desc)

            # The state of a step is saved after each time step by calling `post_push_hook`
            # The save routine may choose to do nothing with this information
            if self.save_routine is not None:
                self.save_routine.post_push_hook()

        # Simulation has finished running
        self._has_run = True

        # Force save of the final state of the simulation if a save routine is
        # provided
        if self.save_routine is not None:
            self.save_routine.save()

        pbar.close()

        self._log("Run completed")

    @property
    def num_entered(self):
        """Count the number of particles that have entered the grids.
        This number is calculated by summing the number of non-zero entries in the
        entered grid array.
        """

        return (self.ever_entered_any_grid > 0).sum()

    @property
    def fract_entered(self):
        """The fraction of particles that have entered the grid.
        The denominator of this fraction is based off the number of tracked
        particles, and therefore does not include stopped or removed particles.
        """
        return (
            self.num_entered - self.num_particles_removed - self.num_particles_stopped
        ) / self.num_particles_tracked

    def _stop_particles(self, particles_to_stop_mask) -> None:
        """Stop tracking the particles specified by the stop mask.

        This is represented by setting the particle's velocity to NaN.
        """

        if len(particles_to_stop_mask) != self.x.shape[0]:
            raise ValueError(
                f"Expected mask of size {self.x.shape[0]}, got {len(particles_to_stop_mask)}"
            )

        self.v[particles_to_stop_mask] = np.nan

        # Reset the cache to update the particle masks
        self._reset_cache()

    def _remove_particles(self, particles_to_remove_mask) -> None:
        """Remove the specified particles from the simulation.

        For the sake of keeping consistent array lengths, the position and
        velocities of the removed particles are set to NaN.
        """

        if len(particles_to_remove_mask) != self.x.shape[0]:
            raise ValueError(
                f"Expected mask of size {self.x.shape[0]}, got {len(particles_to_remove_mask)}"
            )

        self.x[particles_to_remove_mask] = np.nan
        self.v[particles_to_remove_mask] = np.nan

        # Reset the cache to update the particle masks
        self._reset_cache()

    # *************************************************************************
    # Run/push loop methods
    # *************************************************************************

    def _adaptive_dt(self) -> NDArray[np.float64] | float:
        r"""
        Calculate the appropriate dt for each grid based on a number of
        considerations including the local grid resolution (ds) and the
        gyroperiod of the particles in the current fields.
        """

        # candidate time steps includes one per grid (based on the grid resolution)
        # plus additional _dt_candidates based on the field at each particle
        candidates = np.full(
            (self.num_particles, self.num_grids + 1), fill_value=np.inf
        )

        # Compute the time step indicated by the grid resolution
        _gridstep = self._Courant_parameter * np.abs(self._grid_resolutions / self.vmax)
        _max_gridstep = np.max(_gridstep)
        _min_gridstep = np.min(_gridstep)

        # Wherever a particle is on a grid, include that grid's grid step
        # in the list of candidate time steps. If the particle is on no grid,
        # give it the grid step of the highest resolution grid
        for i, _grid in enumerate(self.grids):
            candidates[:, i] = np.where(
                self.particles_on_grid[:, i] > 0, _gridstep[i], _min_gridstep
            )

        # If not, compute a number of possible time steps

        # Compute the cyclotron gyroperiod for each particle
        if self._B_on_grids:
            Bmag = np.linalg.norm(self._B, axis=-1)
            mask = Bmag != 0
            gyroperiod = np.full(Bmag.shape, fill_value=np.inf)
            # TODO: Replace with formulary gyrofrequency lite function once available
            gyroperiod[mask] = 2 * np.pi * self.m / (np.abs(self.q) * Bmag[mask])

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
            dt[dt == np.inf] = _max_gridstep
        else:
            # a single value for dt is returned
            # this is the time step used for all particles
            dt = np.min(candidates)

        return dt

    def _interpolate_grid(self) -> None:
        # Get a list of positions (input for interpolator)
        pos_tracked = self.x[self._tracked_particle_mask]

        # Zero out the array of results
        self._total_grid_values = {
            field_name: required_quantity * 0
            for field_name, required_quantity in self._total_grid_values.items()
        }

        for i, grid in enumerate(self.grids):
            match self.field_weighting:
                case "volume averaged":
                    interpolation_method = grid.volume_averaged_interpolator
                case "nearest neighbor":
                    interpolation_method = grid.nearest_neighbor_interpolator

            # Use the keys of `total_grid_values` as input quantity strings to the interpolator
            grid_values = interpolation_method(
                pos_tracked * u.m,
                *self._interpolated_quantities_per_grid[i],
                persistent=True,
            )

            # If only a single quantity is defined on the grid, the interpolator does not
            # wrap it in a tuple, so we do so here
            if not isinstance(grid_values, tuple):
                grid_values = (grid_values,)

            # Iterate through the interpolated fields and add them to the running sum
            # NaN values are zeroed
            for grid_value, field_name in zip(
                grid_values,
                self._interpolated_quantities_per_grid[i],
                strict=True,
            ):
                self._total_grid_values[field_name][self._tracked_particle_mask] += (
                    np.nan_to_num(
                        grid_value,
                        0.0 * AbstractGrid.recognized_quantities()[field_name].unit,
                    )
                )
        self._E[:, 0] = self._total_grid_values["E_x"].to(u.V / u.m).value
        self._E[:, 1] = self._total_grid_values["E_y"].to(u.V / u.m).value
        self._E[:, 2] = self._total_grid_values["E_z"].to(u.V / u.m).value

        self._B[:, 0] = self._total_grid_values["B_x"].to(u.T).value
        self._B[:, 1] = self._total_grid_values["B_y"].to(u.T).value
        self._B[:, 2] = self._total_grid_values["B_z"].to(u.T).value

    def _update_time(self):
        r"""
        Calculate an appropriate time step for the simulation with respect to
        user configuration. Returns the calculated time step.
        """
        # Calculate the adaptive time step from the fields currently experienced
        # by the particles
        # If user sets dt explicitly, that's handled in _adaptive_dt

        dt = self._adaptive_dt() if self._is_adaptive_time_step else self.dt

        # Make sure the time step can be multiplied by a [num_particles, 3] shape field array
        if isinstance(dt, np.ndarray) and dt.size > 1:
            dt = dt[self._tracked_particle_mask, np.newaxis]
            self.time[self._tracked_particle_mask] += dt
        else:
            self.time += dt

        return dt

    def _update_position(self) -> None:
        r"""
        Update the positions and velocities of the simulated particles using the
        integrator provided at instantiation.
        """
        x_results, v_results = self._integrator.push(
            self.x[self._tracked_particle_mask],
            self.v[self._tracked_particle_mask],
            self._B[self._tracked_particle_mask],
            self._E[self._tracked_particle_mask],
            self.q,
            self.m,
            self.dt,
        )

        self.x[self._tracked_particle_mask], self.v[self._tracked_particle_mask] = (
            x_results,
            v_results,
        )

        # Reset cached properties since particles may have moved off-grid
        self._reset_cache()

    def _update_velocity_stopping(self) -> None:
        r"""
        Apply stopping to the simulated particles using the provided stopping
        routine. The stopping is applied to the simulation by calculating the
        new energy values of the particles, and then updating the particles'
        velocity to match these energies.
        """

        current_speeds = np.linalg.norm(
            self.v[self._tracked_particle_mask], axis=-1, keepdims=True
        )
        velocity_unit_vectors = np.multiply(
            1 / current_speeds, self.v[self._tracked_particle_mask]
        )
        dx = np.multiply(current_speeds, self.dt)

        stopping_power = np.zeros((self.num_particles_tracked, 1))
        relevant_kinetic_energy = (
            self._particle_kinetic_energy[self._tracked_particle_mask] * u.J
        )

        # Apply all previously created stopping power interpolators
        # These are created prior to run, for each grid where they apply
        # in add_stopping()
        match self._stopping_method:
            case "NIST":
                for cs in self._stopping_power_interpolators:
                    if cs is not None:
                        stopping_power += cs(relevant_kinetic_energy).si.value

                energy_loss_per_length = np.multiply(
                    stopping_power,
                    self._total_grid_values["rho"].si.value[
                        self._tracked_particle_mask, np.newaxis
                    ],
                )
            case "Bethe":
                for cs in self._stopping_power_interpolators:
                    if cs is not None:
                        interpolation_result = cs(
                            current_speeds,
                            self._total_grid_values["n_e"].si.value[
                                self._tracked_particle_mask, np.newaxis
                            ],
                        )

                        stopping_power += interpolation_result

                energy_loss_per_length = stopping_power

                if (
                    not self._raised_energy_warning
                    and np.min(energy_loss_per_length * u.J).to(u.MeV).value < 1
                ):
                    self._raised_energy_warning = True

                    warnings.warn(
                        "The Bethe model is only valid for high energy particles. Consider using"
                        "NIST stopping if you require accurate stopping powers at lower energies.",
                        category=PhysicsWarning,
                    )

        dE = -np.multiply(energy_loss_per_length, dx)

        # Update the velocities of the particles using the new energy values
        # TODO: again, figure out how to differentiate relativistic and classical cases
        E = self._particle_kinetic_energy[self._tracked_particle_mask] + dE

        particles_to_be_stopped_mask = np.full(
            shape=self._tracked_particle_mask.shape, fill_value=False
        )
        tracked_particles_to_be_stopped_mask = (
            E < 0
        ).flatten()  # A subset of the tracked particles!
        # Of the tracked particles, stop the ones indicated by the subset mask
        particles_to_be_stopped_mask[self._tracked_particle_mask] = (
            tracked_particles_to_be_stopped_mask
        )

        # Eliminate negative energies before calculating new speeds
        E = np.where(E < 0, 0, E)
        new_speeds = np.sqrt(2 * E / self.m)
        self.v[self._tracked_particle_mask] = np.multiply(
            new_speeds, velocity_unit_vectors
        )

        # Stop particles, which resets the cache
        self._stop_particles(particles_to_be_stopped_mask)

    def _push(self) -> None:
        r"""
        Advance particles using an implementation of the time-centered
        Boris algorithm.
        """
        self.iteration_number += 1

        # Interpolate fields at particle positions
        self._interpolate_grid()

        # Calculate an appropriate timestep (uniform, synchronized)
        self.dt = self._update_time()

        # Update the position and velocities of the particles using timestep
        # calculations as well as the magnitude of E and B fields
        self._update_position()

        if not self._integrator.is_relativistic and not self._raised_relativity_warning:
            beta_max = self.vmax / const.c.si.value

            if beta_max >= 0.001:
                warnings.warn(
                    f"Particles have reached {beta_max}% of the speed of light. Consider using a relativistic integrator for more accurate results.",
                    RelativityWarning,
                )

                self._raised_relativity_warning = True

        # Update velocities to reflect stopping
        if self._do_stopping:
            self._update_velocity_stopping()

        # Update this array, which will have zero elements at the end
        # only for particles that never entered any grid
        self.ever_entered_any_grid += np.sum(self.particles_on_grid, axis=-1).astype(
            np.bool_
        )

    def _reset_cache(self) -> None:
        """
        Reset the cached properties.

        Called after any method that updates the particle positions or
        velocities. This includes the ``_update_positions`` and
        the ``_stop_particles`` and ``_remove_particles`` methods,
        which set some positions or velocities to NaN.
        """
        properties = [
            "particles_on_grid",
            "on_any_grid",
            "vmax",
            "_particle_kinetic_energy",
            "_tracked_particle_mask",
            "num_particles_tracked",
            "_stopped_particle_mask",
            "num_particles_stopped",
            "_removed_particle_mask",
            "num_particles_removed",
        ]

        for prop in properties:
            if prop in self.__dict__:
                del self.__dict__[prop]

    @cached_property
    def particles_on_grid(self):
        r"""
        Returns a boolean mask of shape [ngrids, num_particles] corresponding to
        whether or not the particle is on the associated grid.
        """

        all_particles = np.array([grid.on_grid(self.x * u.m) for grid in self.grids]).T
        all_particles[~self._tracked_particle_mask] = False

        return all_particles

    @cached_property
    def on_any_grid(self) -> NDArray[np.bool_]:
        """
        Binary array for each particle indicating whether it is currently
        on ANY grid.
        """
        return np.sum(self.particles_on_grid, axis=-1) > 0

    @cached_property
    def vmax(self) -> float:
        """The maximum velocity of any particle in the simulation.

        This quantity is used for determining the grid crossing maximum time step.
        """
        return float(
            np.max(np.linalg.norm(self.v[self._tracked_particle_mask], axis=-1))
        )

    @cached_property
    def _particle_kinetic_energy(self):
        r"""
        Return the non-relativistic kinetic energy of the particles.
        """

        # TODO: how should the relativistic case be handled?
        return 0.5 * self.m * np.square(np.linalg.norm(self.v, axis=-1, keepdims=True))

    @cached_property
    def _tracked_particle_mask(self) -> NDArray[np.bool_]:
        """
        Calculates a boolean mask corresponding to particles that have not been stopped or removed.
        """
        # See Class docstring for definition of `stopped` and `removed`
        return ~np.logical_or(np.isnan(self.x[:, 0]), np.isnan(self.v[:, 0]))

    @cached_property
    def num_particles_tracked(self) -> int:
        """Return the number of particles currently being tracked.
        That is, they do not have NaN position or velocity.
        """
        return int(self._tracked_particle_mask.sum())

    @cached_property
    def _stopped_particle_mask(self) -> NDArray[np.bool_]:
        """
        Calculates a boolean mask corresponding to particles that have not been stopped or removed.
        """
        # See Class docstring for definition of `stopped` and `removed`
        return np.logical_and(~np.isnan(self.x[:, 0]), np.isnan(self.v[:, 0]))

    @cached_property
    def num_particles_stopped(self) -> int:
        """Return the number of particles currently being tracked.
        That is, they do not have NaN position or velocity.
        """
        return int(self._stopped_particle_mask.sum())

    @cached_property
    def _removed_particle_mask(self) -> NDArray[np.bool_]:
        """
        Calculates a boolean mask corresponding to particles that have not been stopped or removed.
        """
        # See Class docstring for definition of `stopped` and `removed`
        return np.logical_and(np.isnan(self.x[:, 0]), np.isnan(self.v[:, 0]))

    @cached_property
    def num_particles_removed(self) -> int:
        """Return the number of particles currently being tracked.
        That is, they do not have NaN position or velocity.
        """
        return int(self._removed_particle_mask.sum())

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
