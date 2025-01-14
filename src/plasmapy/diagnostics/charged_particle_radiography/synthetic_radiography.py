"""
Routines for the analysis of proton radiographs. These routines can be broadly
classified as either creating synthetic radiographs from prescribed fields or
methods of 'inverting' experimentally created radiographs to reconstruct the
original fields (under some set of assumptions).
"""

__all__ = ["Tracker", "synthetic_radiograph"]

import warnings
from collections.abc import Iterable
from pathlib import Path
from typing import Literal

import astropy.constants as const
import astropy.units as u
import h5py
import numpy as np

from plasmapy import particles
from plasmapy.formulary.mathematics import rot_a_to_b
from plasmapy.particles import Particle
from plasmapy.plasma.grids import AbstractGrid
from plasmapy.simulation.particle_tracker.particle_tracker import ParticleTracker
from plasmapy.simulation.particle_tracker.save_routines import SaveOnceOnCompletion
from plasmapy.simulation.particle_tracker.termination_conditions import (
    AllParticlesOffGridTerminationCondition,
)


def _coerce_to_cartesian_si(pos):
    """
    Takes a tuple of `astropy.units.Quantity` values representing a position
    in space in either Cartesian, cylindrical, or spherical coordinates, and
    returns a numpy array representing the same point in Cartesian
    coordinates and units of meters.
    """
    # Auto-detect geometry based on units
    geo_units = [x.unit for x in pos]
    if geo_units[2].is_equivalent(u.rad):
        geometry = "spherical"
    elif geo_units[1].is_equivalent(u.rad):
        geometry = "cylindrical"
    else:
        geometry = "cartesian"

    # Convert geometrical inputs between coordinates systems
    pos_out = np.zeros(3)
    if geometry == "cartesian":
        x, y, z = pos
        pos_out[0] = x.to(u.m).value
        pos_out[1] = y.to(u.m).value
        pos_out[2] = z.to(u.m).value

    elif geometry == "cylindrical":
        r, t, z = pos
        r = r.to(u.m)
        t = t.to(u.rad).value
        z = z.to(u.m)
        pos_out[0] = (r * np.cos(t)).to(u.m).value
        pos_out[1] = (r * np.sin(t)).to(u.m).value
        pos_out[2] = z.to(u.m).value

    elif geometry == "spherical":
        r, t, p = pos
        r = r.to(u.m)
        t = t.to(u.rad).value
        p = p.to(u.rad).value

        pos_out[0] = (r * np.sin(t) * np.cos(p)).to(u.m).value
        pos_out[1] = (r * np.sin(t) * np.sin(p)).to(u.m).value
        pos_out[2] = (r * np.cos(t)).to(u.m).value

    return pos_out


class _SyntheticRadiographySaveRoutine(SaveOnceOnCompletion):
    def __init__(
        self, output_directory: Path | None = None, output_basename: str = "output"
    ) -> None:
        super().__init__(
            output_directory=output_directory, output_basename=output_basename
        )

        self._quantities = {
            "source": (u.m, "attribute"),
            "detector": (u.m, "attribute"),
            "mag": (u.m, "attribute"),
            "max_deflection": (None, "attribute"),
            "num_particles": (None, "attribute"),
            "x": (u.m, "dataset"),
            "y": (u.m, "dataset"),
            "v": (u.m / u.s, "dataset"),
            "x0": (u.m, "dataset"),
            "y0": (u.m, "dataset"),
            "v0": (u.m, "dataset"),
        }

    def save(self) -> None:
        result_dictionary = self._particle_tracker.results_dict

        if self.output_directory is None:
            return

        output_file_path = self.output_directory / Path(f"{self.output_basename}.h5")

        with h5py.File(output_file_path, "w") as output_file:
            for key, (_units, data_type) in self._quantities.items():
                match data_type:
                    case "attribute":
                        output_file.attrs.create(key, result_dictionary[key])
                    case "dataset":
                        output_file.create_dataset(key, data=result_dictionary[key])


class Tracker(ParticleTracker):
    r"""
    Represents a charged particle radiography experiment with simulated or
    calculated E and B fields given at positions defined by a grid of spatial
    coordinates. The particle source and detector plane are defined by vectors
    from the origin of the grid.

    Parameters
    ----------
    grids : `~plasmapy.plasma.grids.AbstractGrid` or subclass thereof, or list
        of same.
        A Grid object or list of grid objects containing the required
        quantities [E_x, E_y, E_z, B_x, B_y, B_z].
        If any of these quantities are missing, a warning will be given and that
        quantity will be assumed to be zero everywhere. Grids with large values at
        the edges can cause numerical artifacts.
        The `~plasmapy.plasma.grids.CartesianGrid.soften_edges` method provides
        one way of smoothing this discontinuity out.

    source : `~astropy.units.Quantity`, shape (3)
        A vector pointing from the origin of the grid to the location
        of the particle source. This vector will be interpreted as
        being in either Cartesian, cylindrical, or spherical coordinates
        based on its units. Valid geometries are:

        * Cartesian (x,y,z) : (meters, meters, meters)
        * Cylindrical (r, theta, z) : (meters, radians, meters)
        * Spherical (r, theta, phi) : (meters, radians, radians)

        In spherical coordinates theta is the polar angle.

    detector : `~astropy.units.Quantity`, shape (3)
        A vector pointing from the origin of the grid to the center
        of the detector plane. The vector from the source point to this
        point defines the normal vector of the detector plane. This vector
        can also be specified in Cartesian, cylindrical, or spherical
        coordinates (see the ``source`` keyword).

    dt : `~astropy.units.Quantity`, optional
        An explicitly set time step in units convertible to seconds.
        Setting this optional keyword overrules the adaptive time step
        capability and forces the use of this time step throughout.

    dt_range : tuple of shape (2,) of `~astropy.units.Quantity`, optional
        If specified, the calculated adaptive time step will be clamped
        between the first and second values.

    field_weighting : str, default: ``"volume averaged"``
        String that selects the field weighting algorithm used to determine
        what fields are felt by the particles. Options are:

        * ``"nearest neighbor"``:
            Particles are assigned the fields on the grid vertex closest to
            them.
        * ``"volume averaged"``:
            The fields experienced by a particle are a volume-average of the
            eight grid points surrounding them.

    detector_hdir : `numpy.ndarray`, shape (3), optional
        A unit vector (in Cartesian coordinates) defining the horizontal
        direction on the detector plane. By default, the horizontal axis in the
        detector plane is defined to be perpendicular to both the
        source-to-detector vector and the z-axis (unless the source-to-detector
        axis is parallel to the z axis, in which case the horizontal axis is
        the x-axis).

        The detector vertical axis is then defined
        to be orthogonal to both the source-to-detector vector and the
        detector horizontal axis.

    output_directory : `~pathlib.Path`, optional
        Directory for objects that are saved to disk. If a directory is not
        specified then a memory save routine is used.

    output_basename : `str`, optional
        Optional base name for output files.

    fraction_exited_threshold : float, optional
        The fraction of particles that must leave the grids to terminate the
        simulation. This does not include particles that have never entered
        the grids. By default, this is set to ``0.999`` (corresponding to 99.9%).

    verbose : bool, optional
        If `True`, updates on the status of the program will be printed
        into the standard output while running.
    """

    def __init__(
        self,
        grids: AbstractGrid | Iterable[AbstractGrid],
        source: u.Quantity[u.m],
        detector: u.Quantity[u.m],
        dt=None,
        dt_range=None,
        field_weighting: Literal[
            "volume averaged", "nearest neighbor"
        ] = "volume averaged",
        detector_hdir=None,
        output_directory: Path | None = None,
        output_basename: str = "output",
        fraction_exited_threshold: float = 0.999,
        verbose: bool = True,
    ) -> None:
        # The synthetic radiography class handles logging, so we can disable logging for the particle tracker
        # The particle tracker class ensures that the provided grid argument has the proper type and
        # that the necessary grid quantities are created if they are not already specified
        save_routine = (
            _SyntheticRadiographySaveRoutine(output_directory, output_basename)
            if output_directory is not None
            else None
        )

        termination_condition = AllParticlesOffGridTerminationCondition(
            fraction_exited_threshold=fraction_exited_threshold
        )

        super().__init__(
            grids,
            termination_condition,
            save_routine,
            dt=dt,
            dt_range=dt_range,
            field_weighting=field_weighting,
            verbose=verbose,
        )

        # A list of wire meshes added to the grid with add_wire_mesh
        # Particles that would hit these meshes will be removed at runtime
        # by _apply_wire_mesh
        self.mesh_list = []

        # ************************************************************************
        # Setup the source and detector geometries
        # ************************************************************************

        self.source = _coerce_to_cartesian_si(source)
        self.detector = _coerce_to_cartesian_si(detector)
        self._log(f"Source: {self.source} m")
        self._log(f"Detector: {self.detector} m")

        # Calculate normal vectors (facing towards the grid origin) for both
        # the source and detector planes
        self.src_n = -self.source / np.linalg.norm(self.source)
        self.det_n = -self.detector / np.linalg.norm(self.detector)

        # Vector directly from source to detector
        self.src_det = self.detector - self.source

        # Magnification
        self.mag = 1 + (np.linalg.norm(self.detector) / np.linalg.norm(self.source))
        self._log(f"Magnification: {self.mag}")

        # Check that source-detector vector actually passes through the grid
        test = [
            grid.vector_intersects(self.source * u.m, self.detector * u.m)
            for grid in self.grids
        ]
        if not any(test):
            raise ValueError(
                "The vector between the source and the detector "
                "does not intersect the grid provided!"
            )

        # Determine the angle above which particles will not hit the grid
        # these particles can be ignored until the end of the simulation,
        # then immediately advanced to the detector grid with their original
        # velocities
        self.max_theta_hit_grid = self._max_theta_hit_grid()

        # *********************************************************************
        # Define the detector plane
        # *********************************************************************

        # Load or calculate the detector hdir
        if detector_hdir is not None:
            self.det_hdir = detector_hdir / np.linalg.norm(detector_hdir)
        else:
            self.det_hdir = self._default_detector_hdir()

        # Calculate the detector vdir
        ny = np.cross(self.det_hdir, self.det_n)
        self.det_vdir = -ny / np.linalg.norm(ny)

    def _default_detector_hdir(self):
        """
        Calculates the default horizontal unit vector for the detector plane
        (see __init__ description for details).
        """
        # Create unit vectors that define the detector plane
        # Define plane  horizontal axis
        if np.allclose(np.abs(self.det_n), np.array([0, 0, 1])):
            nx = np.array([1, 0, 0])
        else:
            nx = np.cross(np.array([0, 0, 1]), self.det_n)
        return nx / np.linalg.norm(nx)

    def _max_theta_hit_grid(self):
        r"""
        Using the grid and the source position, compute the maximum particle
        theta that will impact the grid. This value can be used to determine
        which particles are worth tracking.
        """
        theta = np.zeros([8, self.num_grids])

        for i, _grid in enumerate(self.grids):
            ind = 0
            for x in (0, -1):
                for y in (0, -1):
                    for z in (0, -1):
                        # Source to grid corner vector
                        vec = self.grids_arr[i][x, y, z, :] - self.source

                        # Calculate angle between vec and the source-to-detector
                        # axis, which is the central axis of the particle beam
                        theta[ind, i] = np.arccos(
                            np.dot(vec, self.src_det)
                            / np.linalg.norm(vec)
                            / np.linalg.norm(self.src_det)
                        )
                        ind += 1

        return np.max(theta)

    def _log(self, msg) -> None:
        if self.verbose:
            # We'll need to switch from print() to using logging library
            print(msg)  # noqa: T201

    # Define some constants so they don't get constantly re-evaluated
    _c = const.c.si.value

    # *************************************************************************
    # Create mesh
    # *************************************************************************

    def add_wire_mesh(
        self, location, extent, nwires, wire_diameter, mesh_hdir=None, mesh_vdir=None
    ):
        """
        Add a wire mesh grid between the particle source and the object grid
        that blocks particles whose paths intersect the wires.

        Parameters
        ----------
        location : `~astropy.units.Quantity`, shape (3)
            A vector pointing from the origin of the grid to the center of the
            mesh grid. This location must be between the source and the
            object grid.

            This vector will be interpreted as
            being in either cartesian, cylindrical, or spherical coordinates
            based on its units. Valid geometries are:

            * Cartesian (x,y,z) : (meters, meters, meters)
            * Cylindrical (r, theta, z) : (meters, radians, meters)
            * Spherical (r, theta, phi) : (meters, radians, radians)

            In spherical coordinates theta is the polar angle.

        extent : Tuple of 1 or 2 `~astropy.units.Quantity`
            The size of the mesh grid (in the mesh plane). If one value
            is provided, the mesh is circular and the value provided is
            interpreted as the diameter. If two values are provided, the
            mesh is rectangular and the values are interpreted as the
            width and height respectively.

        nwires : Tuple of 1 or 2 ints, or a single int
            The number of wires in the horizontal and vertical directions. If
            only one value is provided, the number in the two directions is
            assumed to be equal. Note that a wire will cross the center of the
            mesh only when nwires is odd.

        wire_diameter : `~astropy.units.Quantity`
            The diameter of the wires.

        mesh_hdir : `numpy.ndarray`, shape (3), optional
            A unit vector (in Cartesian coordinates) defining the horizontal
            direction on the mesh plane. Modifying this vector can rotate the
            mesh in the plane or tilt the mesh plane relative to the
            source-detector axis. By default, ``mesh_hdir`` is set equal to
            ``detector_hdir`` (see ``detector_hdir`` keyword in ``__init__``).

        mesh_vdir : `numpy.ndarray`, shape (3), optional
            A unit vector (in Cartesian coordinates) defining the vertical
            direction on the mesh plane. Modifying this vector can tilt the
            mesh relative to the source-detector axis. By default, ``mesh_vdir``
            is defined to be perpendicular to ``mesh_hdir`` and the detector
            plane normal (such that the mesh is parallel to the detector plane).

        Raises
        ------
        ValueError
            Raises a ValueError if the provided mesh location is not
            between the source and the object grid.

        """

        # Raise an error if the run method has already been called.
        self._enforce_order()

        location = _coerce_to_cartesian_si(location)
        wire_radius = wire_diameter.si.value / 2

        if not isinstance(extent, tuple):
            extent = (extent,)

        if len(extent) == 1:
            radius = 0.5 * extent[0].si.value
            width = extent[0].si.value
            height = extent[0].si.value
        elif len(extent) == 2:
            radius = None
            width = extent[0].si.value
            height = extent[1].si.value
        else:
            raise ValueError(
                "extent must be a tuple of 1 or 2 elements, but "
                f"{len(extent)} elements were provided."
            )

        if not isinstance(nwires, tuple):
            nwires = (nwires,)

        if len(nwires) != 2:
            nwires = (nwires[0], nwires[0])

        # If no hdir/vdir is specified, calculate a default value
        # If one is specified, make sure it is normalized
        if mesh_hdir is None:
            # Re-calculate the default here, in case the user
            # specified a different det_hdir
            mesh_hdir = self._default_detector_hdir()
        else:
            mesh_hdir = mesh_hdir / np.linalg.norm(mesh_hdir)

        if mesh_vdir is None:
            mesh_vdir = np.cross(mesh_hdir, self.det_n)
            mesh_vdir = -mesh_vdir / np.linalg.norm(mesh_vdir)
        else:
            mesh_vdir = mesh_vdir / np.linalg.norm(mesh_vdir)

        # Raise exception if mesh is AFTER the field grid
        if np.linalg.norm(location - self.source) > np.linalg.norm(self.source):
            raise ValueError(
                f"The specified mesh location, {location},"
                "is not between the source and the origin."
            )

        mesh_entry = {
            "location": location,
            "wire_radius": wire_radius,
            "radius": radius,
            "width": width,
            "height": height,
            "nwires": nwires,
            "mesh_hdir": mesh_hdir,
            "mesh_vdir": mesh_vdir,
        }

        self.mesh_list.append(mesh_entry)

    def _apply_wire_mesh(
        self,
        location=None,
        wire_radius=None,
        radius=None,
        width=None,
        height=None,
        nwires=None,
        mesh_hdir=None,
        mesh_vdir=None,
    ):
        """
        Apply wire meshes that were added to ``self.mesh_list``.
        """
        x = self._coast_to_plane(location, mesh_hdir, mesh_vdir)

        # Particle positions in 2D on the mesh plane
        xloc = np.dot(x - location, mesh_hdir)
        yloc = np.dot(x - location, mesh_vdir)

        # Create an array in which True indicates that a particle has hit
        # a wire and False indicates that it has not
        hit = np.zeros(self.num_particles, dtype=bool)

        # Mark particles that overlap vertical or horizontal position with
        # a wire
        h_centers = np.linspace(-width / 2, width / 2, num=nwires[0])
        for c in h_centers:
            hit |= np.isclose(xloc, c, atol=wire_radius)

        v_centers = np.linspace(-height / 2, height / 2, num=nwires[1])
        for c in v_centers:
            hit |= np.isclose(yloc, c, atol=wire_radius)

        # Put back any particles that are outside the mesh boundaries
        # First handle the case where the mesh is rectangular
        if radius is None:
            # Replace particles outside the x-boundary
            hit[
                np.logical_or(
                    xloc > np.max(h_centers) + wire_radius,
                    xloc < np.min(h_centers) - wire_radius,
                )
            ] = False
            # Replace particles outside the y-boundary
            hit[
                np.logical_or(
                    yloc > np.max(v_centers) + wire_radius,
                    yloc < np.min(v_centers) - wire_radius,
                )
            ] = False
        # Handle the case where the mesh is circular
        else:
            loc_rad = np.sqrt(xloc**2 + yloc**2)
            hit[loc_rad > radius] = False

            # In the case of a circular mesh, also create a round wire along
            # the outside edge
            hit[np.isclose(loc_rad, radius, atol=wire_radius)] = True

        # Identify the particles that have hit something, then remove them from
        # all of the arrays
        keep_these_particles = ~hit
        number_kept_particles = keep_these_particles.sum()
        nremoved = self.num_particles - number_kept_particles

        if self.num_particles - nremoved <= 0:
            raise ValueError(
                "The specified mesh is blocking all of the particles. "
                f"The wire diameter ({2*wire_radius}) may be too large."
            )

        self._stop_particles(~keep_these_particles)

    # *************************************************************************
    # Particle creation methods
    # *************************************************************************

    @staticmethod
    def _angles_monte_carlo(num_particles, max_theta, random_seed=None):
        """
        Generates angles for each particle randomly such that the flux
        per solid angle is uniform.
        """
        # Create a probability vector for the theta distribution
        # Theta must follow a sine distribution in order for the particle
        # flux per solid angle to be uniform.
        arg = np.linspace(0, max_theta, num=int(1e5))
        prob = np.sin(arg)
        prob *= 1 / np.sum(prob)

        # Create a numpy random number generator instance
        rng = np.random.default_rng(seed=random_seed)

        # Randomly choose theta's weighted with the sine probabilities
        theta = rng.choice(arg, size=num_particles, replace=True, p=prob)

        # Also generate a uniform phi distribution
        phi = rng.uniform(high=2 * np.pi, size=num_particles)

        return theta, phi

    @staticmethod
    def _angles_uniform(num_particles, max_theta):
        """
        Generates angles for each particle such that their velocities are
        uniformly distributed on a grid in theta and phi. This method
        requires that `num_particles` be a perfect square. If it is not,
        `num_particles` will be set as the largest perfect square smaller
        than the provided `num_particles`.
        """
        # Calculate the approximate square root
        n_per = np.floor(np.sqrt(num_particles)).astype(np.int32)

        # Create an imaginary grid positioned 1 unit from the source
        # and spanning max_theta at the corners
        extent = np.sin(max_theta) / np.sqrt(2)
        arr = np.linspace(-extent, extent, num=n_per)
        harr, varr = np.meshgrid(arr, arr, indexing="ij")

        # calculate the angles from the source for each point in
        # the grid.
        theta = np.arctan(np.sqrt(harr**2 + varr**2))
        phi = np.arctan2(varr, harr)

        return theta.flatten(), phi.flatten()

    @particles.particle_input
    def create_particles(
        self,
        num_particles,
        particle_energy,
        max_theta=None,
        particle: Particle = Particle("p+"),  # noqa: B008
        distribution: Literal["monte-carlo", "uniform"] = "monte-carlo",
        random_seed=None,
    ) -> None:
        r"""
        Generates the angular distributions about the Z-axis, then
        rotates those distributions to align with the source-to-detector axis.

        By default, particles are generated over almost the entire pi/2.
        However, if the detector is far from the source, many of these
        particles will never be observed. The max_theta keyword allows these
        extraneous particles to be neglected to focus computational resources
        on the particles who will actually hit the detector.

        Parameters
        ----------
        num_particles : integer
            The number of particles to include in the simulation. The default
            is 1e5.

        particle_energy : `~astropy.units.Quantity`
            The energy of the particle, in units convertible to eV.
            All particles are given the same energy.

        max_theta : `~astropy.units.Quantity`, optional
            The largest velocity vector angle (measured from the
            source-to-detector axis) for which particles should be generated.
            Decreasing this angle can eliminate particles that would never
            reach the detector region of interest. If no value is given, a
            guess will be made based on the size of the grid.
            Units must be convertible to radians.

        particle : ~plasmapy.particles.particle_class.Particle or string representation of same, optional
            Representation of the particle species as either a |Particle| object
            or a string representation. The default particle is protons.

        distribution: str
            A keyword which determines how particles will be distributed
            in velocity space. Options are:

                - 'monte-carlo': velocities will be chosen randomly,
                    such that the flux per solid angle is uniform.

                - 'uniform': velocities will be distributed such that,
                   left unperturbed,they will form a uniform pattern
                   on the detection plane. This method
                   requires that ``num_particles`` be a perfect square. If it is not,
                   ``num_particles`` will be set as the largest perfect square smaller
                   than the provided ``num_particles``.

            Simulations run in the ``'uniform'`` mode will imprint a grid pattern
            on the image, but will well-sample the field grid with a
            smaller number of particles. The default is ``'monte-carlo'``.

        random_seed : int, optional
            A random seed to be used when generating random particle
            distributions, e.g. with the ``monte-carlo`` distribution.
        """
        self._log("Creating Particles")

        # Raise an error if the run method has already been called.
        self._enforce_order()

        # Load inputs
        num_particles = int(num_particles)

        particle_energy = particle_energy.to(u.eV).value
        m = particle.mass.to(u.kg).value

        # If max_theta is not specified, make a guess based on the grid size
        if max_theta is None:
            max_theta = np.clip(1.5 * self.max_theta_hit_grid, 0.01, 0.99 * np.pi / 2)
        else:
            max_theta = max_theta.to(u.rad).value

        # Calculate the velocity corresponding to the particle energy
        ER = particle_energy * 1.6e-19 / (m * self._c**2)
        v0 = self._c * np.sqrt(1 - 1 / (ER + 1) ** 2)

        if distribution == "monte-carlo":
            theta, phi = self._angles_monte_carlo(
                num_particles, max_theta, random_seed=random_seed
            )
        elif distribution == "uniform":
            theta, phi = self._angles_uniform(num_particles, max_theta)

        # Adjust num_particles to reflex what the distribution function returned.
        # Some distributions will modify the number of particles to meet the
        # necessary criteria of the distribution.
        num_particles = theta.shape[0]  # TODO: make sure this works

        # Construct the velocity distribution around the z-axis
        v = np.zeros([num_particles, 3])
        v[:, 0] = v0 * np.sin(theta) * np.cos(phi)
        v[:, 1] = v0 * np.sin(theta) * np.sin(phi)
        v[:, 2] = v0 * np.cos(theta)

        # Calculate the rotation matrix that rotates the z-axis
        # onto the source-detector axis
        a = np.array([0, 0, 1])
        b = self.detector - self.source
        rot = rot_a_to_b(a, b)

        # Apply rotation matrix to calculated velocity distribution
        v = np.matmul(v, rot)

        # Place particles at the source
        x = np.tile(self.source, (num_particles, 1))

        # Call the underlying load method to ensure consistency with
        # other properties within the ParticleTracker
        self.load_particles(x * u.m, v * u.m / u.s, particle=particle)

    @particles.particle_input
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

        v: `~astropy.units.Quantity`, shape (N,3)
            Velocities for N particles

        particle : |particle-like|, optional
            Representation of the particle species as either a |Particle| object
            or a string representation. The default particle is protons.
        """
        # Load particles for particle tracker class
        super().load_particles(x, v, particle)

        # But also calculate geometry-dependent variables
        self.theta = np.arccos(
            np.inner(self.v, self.src_n) / np.linalg.norm(self.v, axis=-1)
        )

        n_wrong_way = np.sum(np.where(self.theta > np.pi / 2, 1, 0))
        if n_wrong_way > 1:
            warnings.warn(
                f"{100*n_wrong_way/self.num_particles:.2f}% of particles "
                "initialized are heading away from the grid. Check the "
                " orientation of the provided velocity vectors.",
                RuntimeWarning,
            )

    # *************************************************************************
    # Run/push loop methods
    # *************************************************************************

    def _coast_to_grid(self) -> None:
        r"""
        Coasts all particles to the timestep when the first particle should
        be entering the grid. Doing in this in one step (rather than pushing
        the particles through zero fields) saves computation time.
        """
        # Distance from the source to the nearest point on any grid
        dist = min(
            np.min(np.linalg.norm(arr - self.source, axis=3)) for arr in self.grids_arr
        )

        tracked_mask = self._tracked_particle_mask

        # Find speed of each particle towards the grid.
        vmax = np.dot(self.v[tracked_mask], self.src_n)

        # Time for each particle to reach the grid
        t = dist / vmax

        # Coast the particles to the advanced position
        self.x[tracked_mask] = (
            self.x[tracked_mask] + self.v[tracked_mask] * t[:, np.newaxis]
        )

    def _coast_to_plane(self, center, hdir, vdir, x=None, mask=None):
        """
        Calculates the positions where the current trajectories of each
        particle impact a plane, described by the plane's center and
        horizontal and vertical unit vectors.

        Returns an [num_particles, 3] array of the particle positions in the plane

        By default this function does not alter self.x. The optional keyword
        x can be used to pass in an output array that will used to hold
        the positions in the plane. This can be used to directly update self.x
        as follows:

        self._coast_to_plane(self.detector, self.det_hdir, self.det_vdir, x
                             = self.x)

        Parameters
        ----------
        center : float
            The center of the plane towards which the particles are advanced

        hdir : `numpy.ndarray`, shape (3)
            A unit vector (in Cartesian coordinates) defining the horizontal
            direction of the plane.

        vdir : `numpy.ndarray`, shape (3)
            A unit vector (in Cartesian coordinates) defining the vertical
            direction of the plane.

        x : `numpy.ndarray`, shape (num_particles), optional
            The array to which the resulting particle positions are stored.
            By default, the current position array will be used.

        mask : `numpy.ndarray`, shape (num_particles), optional
            A boolean mask representing the particles to perform the coasting
            operation. By default, only the tracked particles (i.e. those that
            are going to hit the grids) will be coasted.
        """

        normal = np.cross(hdir, vdir)

        if mask is None:
            mask = self._tracked_particle_mask

        # Calculate the time required to evolve each particle into the
        # plane
        t = np.inner(center[np.newaxis, :] - self.x[mask], normal) / np.inner(
            self.v[mask], normal
        )

        # Calculate particle positions in the plane
        if x is None:
            # If no output array is provided, preallocate
            x = np.copy(self.x)

        x[mask] = self.x[mask] + self.v[mask] * t[:, np.newaxis]

        # Check that all points are now in the plane
        # (Eq. of a plane is nhat*x + d = 0)
        plane_eq = np.dot(x[mask] - center, normal)
        if not np.allclose(plane_eq, 0, atol=1e-6):
            raise ValueError("Coasting particles to plane failed.")

        return x

    def _remove_deflected_particles(self) -> None:
        r"""
        Removes any particles that have been deflected away from the detector
        plane (eg. those that will never hit the grid).
        """
        dist_remaining = np.dot(self.x, self.det_n) + np.linalg.norm(self.detector)

        v_towards_det = np.dot(self.v, -self.det_n)

        # If particles have not yet reached the detector plane and are moving
        # away from it, they will never reach the detector.
        # So, we can remove them from the arrays

        # Find the indices of all particles that we should keep:
        # i.e. those still moving towards the detector.
        ind = np.logical_not((v_towards_det < 0) & (dist_remaining > 0)).nonzero()[0]

        # Drop the other particles
        self._remove_particles((v_towards_det < 0) & (dist_remaining > 0))

        # Store the number of particles deflected
        self.fract_deflected = (self.num_particles - ind.size) / self.num_particles

        # Warn the user if a large number of particles are being deflected
        if self.fract_deflected > 0.05:
            warnings.warn(
                f"{100*self.fract_deflected:.1f}% particles have been "
                "deflected away from the detector plane. The fields "
                "provided may be too high to successfully radiograph "
                "with this particle energy.",
                RuntimeWarning,
            )

    def run(self) -> None:
        r"""
        Runs a particle-tracing simulation.

        Timesteps are adaptively calculated based on the
        local grid resolution of the particles and the electric and magnetic
        fields they are experiencing. After all particles
        have left the grid, they are advanced to the
        detector plane where they can be used to construct a synthetic
        diagnostic image.

        Returns
        -------
        None
        """

        self._enforce_particle_creation()

        # If meshes have been added, apply them now
        for mesh in self.mesh_list:
            self._log("Applying meshes")
            self._apply_wire_mesh(**mesh)

        self._log("Coasting untracked particles to the detector plane")
        # These particles will not be pushed through the fields
        # but instead will be automatically advanced
        # to the detector plane
        theta_mask = self.theta >= self.max_theta_hit_grid

        # Take the intersection of the theta mask and the tracked particle mask
        # This must be done because some particles could have already been stopped
        # by _apply_wire_mesh
        self.x = self._coast_to_plane(
            self.detector,
            self.det_hdir,
            self.det_vdir,
            mask=theta_mask & self._tracked_particle_mask,
        )
        self._stop_particles(theta_mask)
        self.fract_tracked = self.num_particles_tracked / self.num_particles

        self._log("Generating null distribution in detector plane")
        # Store a copy of the initial velocity distribution in memory
        # This will be used later to calculate the maximum deflection
        self.v0 = np.copy(self.v)
        # Generate a null distribution of points (the result in the absence of
        # any fields) for statistical comparison
        self.x0 = self._coast_to_plane(self.detector, self.det_hdir, self.det_vdir)

        self._log("Advancing tracked particles to the start of the grid")
        # Advance the tracked particles to the near the start of the grid
        self._coast_to_grid()
        self.coasted_particles = np.copy(self.x)

        super().run()

        if self.num_entered < 0.1 * self.num_particles:
            warnings.warn(
                f"Only {100 * self.num_entered / self.num_particles:.2f}% of "
                "particles entered the field grid: consider "
                "decreasing the max_theta to increase this "
                "number.",
                RuntimeWarning,
            )

        # Remove particles that will never reach the detector
        self._remove_deflected_particles()

        # Advance the particles to the image plane
        self._coast_to_plane(self.detector, self.det_hdir, self.det_vdir, x=self.x)

        self.save_routine.save()

        # Log a summary of the run
        self._log("Run completed")

        self._log(f"Fraction of particles tracked: {self.fract_tracked:.1%}")

        self._log(
            "Fraction of tracked particles that entered the grid: "
            f"{self.fract_entered*100:.1f}%"
        )

        self._log(
            "Fraction of tracked particles deflected away from the "
            "detector plane: "
            f"{self.fract_deflected*100}%"
        )

    @property
    def max_deflection(self):
        """
        The maximum deflection experienced by one of the particles, determined
        by comparing their initial and final velocity vectors.

        This value can be used to determine the charged particle radiography regime
        using the dimensionless number defined by Kugland et al. 2012

        Returns
        -------
        max_deflection : float
            The maximum deflection in radians.
        """
        # Normalize the initial and final velocities
        v_norm = self.v / np.linalg.norm(self.v, axis=1, keepdims=True)
        v0_norm = self.v0 / np.linalg.norm(self.v0, axis=1, keepdims=True)
        # Compute the dot product
        proj = np.sum(v_norm * v0_norm, axis=1)
        # In case of numerical errors, make sure the output is within the domain of
        # arccos
        proj = np.where(proj > 1, 1, proj)
        max_deflection = np.nanmax(np.arccos(proj))

        return max_deflection * u.rad

    @property
    def results_dict(self):
        r"""
        A dictionary containing the results of the simulation.

        .. list-table:: Dictionary keys and descriptions.
           :width: 100%
           :widths: 1 1 4
           :header-rows: 1

           * - Key
             - Type
             - Description
           * - ``"source"``
             - `~numpy.ndarray`
             - The source location vector, in meters.
           * - ``"detector"``
             - `~numpy.ndarray`
             - The detector location vector, in meters.
           * - ``"mag"``
             - `float`
             - The system magnification.
           * - ``"num_particles"``
             - `int`
             - Number of particles in the simulation.
           * - ``"max_deflection"``
             - `~numpy.ndarray`
             - The maximum deflection experienced by a particle in the
               simulation, in radians.
           * - ``"x"``
             - `~numpy.ndarray`, [``num_particles``,]
             - The x-coordinate location where each particle hit the
               detector plane, in meters.
           * - ``"y"``
             - `~numpy.ndarray`, [``num_particles``,]
             - The y-coordinate location where each particle hit the
               detector plane, in meters.
           * - ``"v"``
             - `~numpy.ndarray`, [``num_particles``, 3]
             - The velocity of each particle when it hits the detector
               plane, in meters per second. The velocity is in a
               coordinate system relative to the detector plane. The
               components are [normal, horizontal, vertical] relative
               to the detector plane coordinates.
           * - ``"x0"``
             - `~numpy.ndarray`, [``num_particles``,]
             - The x-coordinate location where each particle would have
               hit the detector plane if the grid fields were zero, in
               meters. Useful for calculating the source profile.
           * - ``"y0"``
             - `~numpy.ndarray`, [``num_particles``,]
             - The y-coordinate location where each particle would have
               hit the detector plane if the grid fields were zero, in
               meters. Useful for calculating the source profile.
           * - ``"v0"``
             - `~numpy.ndarray`, [``num_particles``, 3]
             - The velocity of each particle when it hit the detector
               plan if the grid fields were zero, in meters per second.
               The velocity is in a coordinate system relative to the
               detector plane. The components are [normal, horizontal,
               vertical] relative to the detector plane coordinates.
        """

        if not self._has_run:
            raise RuntimeError(
                "The simulation must be run before a results "
                "dictionary can be created."
            )

        # Determine locations of points in the detector plane using unit
        # vectors
        xloc = np.dot(self.x - self.detector, self.det_hdir)
        yloc = np.dot(self.x - self.detector, self.det_vdir)

        x0loc = np.dot(self.x0 - self.detector, self.det_hdir)
        y0loc = np.dot(self.x0 - self.detector, self.det_vdir)

        # Determine the velocity components of each particle in the
        # coordinate frame of the detector, eg. the components of v are
        # now [normal, det_hdir, det_vdir]
        v = np.zeros(self.v.shape)
        v[:, 0] = np.dot(self.v, self.det_n)
        v[:, 1] = np.dot(self.v, self.det_hdir)
        v[:, 2] = np.dot(self.v, self.det_vdir)

        v0 = np.zeros(self.v.shape)
        v0[:, 0] = np.dot(self.v0, self.det_n)
        v0[:, 1] = np.dot(self.v0, self.det_hdir)
        v0[:, 2] = np.dot(self.v0, self.det_vdir)

        return {
            "source": self.source,
            "detector": self.detector,
            "mag": self.mag,
            "num_particles": self.num_particles,
            "max_deflection": self.max_deflection.to(u.rad).value,
            "x": xloc,
            "y": yloc,
            "v": v,
            "x0": x0loc,
            "y0": y0loc,
            "v0": v0,
        }


# *************************************************************************
# Synthetic diagnostic methods (creating output)
# *************************************************************************


def synthetic_radiograph(obj, size=None, bins=None, ignore_grid: bool = False):  # noqa: C901, PLR0912
    r"""
    Calculate a "synthetic radiograph" (particle count histogram in the
    image plane).

    .. |Tracker| replace:: `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.Tracker`
    .. |results_dict| replace:: `~plasmapy.diagnostics.charged_particle_radiography.synthetic_radiography.Tracker.results_dict`

    Parameters
    ----------
    obj: `dict` or `~pathlib.Path` or |Tracker|
        Either a |Tracker|
        object that has been run, a dictionary equivalent to
        |results_dict|, or path to a saved output file
        from a |Tracker| object (HDF5 file).

    size : `~astropy.units.Quantity`, shape ``(2, 2)``, optional
        The size of the detector array, specified as the minimum
        and maximum values included in both the horizontal and vertical
        directions in the detector plane coordinates. Shape is
        ``((hmin, hmax), (vmin, vmax))``. If not specified, the size will be
        set to include all particles on the detector. Units must be convertible
        to meters.

    bins : array of integers, shape ``(2)``
        The number of bins in each direction in the format
        ``(hbins, vbins)``.  The default is ``(200, 200)``.

    ignore_grid: `bool`
        If `True`, returns the intensity in the image plane in the absence
        of simulated fields.

    Returns
    -------
    hax : `~astropy.units.Quantity` array shape ``(hbins,)``
        The horizontal axis of the synthetic radiograph in meters.

    vax : `~astropy.units.Quantity` array shape ``(vbins, )``
        The vertical axis of the synthetic radiograph in meters.

    intensity : `~numpy.ndarray`, shape ``(hbins, vbins)``
        The number of particles counted in each bin of the histogram.


    Notes
    -----
    This function ignores any particles that are stopped or removed before
    reaching the detector plane.

    """

    # condition `obj` input
    if isinstance(obj, Tracker):
        # results_dict raises an error if the simulation has not been run.
        results_dict = obj.results_dict

    elif isinstance(obj, dict):
        results_dict = obj

    elif isinstance(obj, str | Path):
        results_dict = {}
        obj = Path(obj)
        # Create a dictionary of all of the datasets and attributes in the save file
        # Equivalent to |results_dict|
        with h5py.File(obj, "r") as f:
            for key in f:
                results_dict[key] = f[key][...]
            for key in f.attrs:
                results_dict[key] = f.attrs[key][...]
    else:
        raise TypeError(
            f"Expected type `Path`, `dict` or {Tracker} for argument `obj`, but "
            f"got type {type(obj)}."
        )

    if bins is None:
        bins = [200, 200]

    # Note that, at the end of the simulation, all particles were moved
    # into the image plane.

    # If ignore_grid is True, use the predicted positions in the absence of
    # simulated fields
    if ignore_grid:
        xloc = results_dict["x0"]
        yloc = results_dict["y0"]
        v = results_dict["v0"][:, 0]
    else:
        xloc = results_dict["x"]
        yloc = results_dict["y"]
        v = results_dict["v"][:, 0]

    if size is None:
        # If a detector size is not given, choose a size based on the
        # particle positions
        w = np.max([np.nanmax(np.abs(xloc)), np.nanmax(np.abs(yloc))])
        size = np.array([[-w, w], [-w, w]]) * u.m
    elif not isinstance(size, u.Quantity):
        raise TypeError(
            "Argument `size` must be an astropy.units.Quantity object with "
            "units convertible to meters."
        )
    elif not size.unit.is_equivalent(u.m):
        raise ValueError("Argument `size` must have units convertible to meters.")
    elif size.shape != (2, 2):
        raise ValueError(
            f"Argument `size` must have shape (2, 2), but got {size.shape}."
        )

    # Exclude NaN positions (deleted particles) and velocities
    # (stopped particles)
    nan_mask = ~np.isnan(xloc) * ~np.isnan(yloc) * ~np.isnan(v)
    sanitized_xloc = xloc[nan_mask]
    sanitized_yloc = yloc[nan_mask]

    # Generate the histogram
    intensity, h, v = np.histogram2d(
        sanitized_xloc, sanitized_yloc, range=size.to(u.m).value, bins=bins
    )

    # h, v are the bin edges: compute the centers to produce arrays
    # of the right length (then trim off the extra point)
    h = ((h + np.roll(h, -1)) / 2)[:-1]
    v = ((v + np.roll(v, -1)) / 2)[:-1]

    # Throw a warning if < 50% of the particles are included on the
    # histogram
    percentage = np.sum(intensity) / results_dict["num_particles"]
    if percentage < 0.5:
        warnings.warn(
            f"Only {percentage:.2%} of the particles are shown "
            "on this synthetic radiograph. Consider increasing "
            "the size to include more.",
            RuntimeWarning,
        )

    return h * u.m, v * u.m, intensity
