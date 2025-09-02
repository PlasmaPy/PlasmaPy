"""
Definition of the material class and helper functions. The module also
includes type definitions for stopping and ranging data.
"""

from collections.abc import Callable
from functools import cached_property
from typing import TypedDict

import astropy.units as u
import numpy as np
from scipy.interpolate import CubicSpline

from plasmapy.particles.particle_class import Particle

__all__ = [
    "MaterialSTARData",
    "Material",
]

from plasmapy.utils.decorators import validate_quantities


class MaterialSTARData(TypedDict):
    """Type definition for formatting data related to beam-target interactions."""

    incident_particle: Particle

    energies: u.Quantity[u.MeV]
    stopping_power: u.Quantity[
        u.MeV * u.cm**2 / u.g
    ]  # Units of energy per areal density
    CSDA_range: u.Quantity[u.cm]
    projected_range: u.Quantity[u.cm]


def _construct_log_interpolator(
    x0: u.Quantity, y0: u.Quantity
) -> Callable[[u.Quantity], u.Quantity]:
    """
    Apply a logarithmic transformation to input and output data points,
    allowing for approximation of an exponential relation using
    a cubic spline.
    """
    x_unit = x0.unit
    y_unit = y0.unit

    log_interpolator = CubicSpline(
        x=np.log(x0.value),
        y=np.log(y0.value),
    )

    def interpolator(x_unsanitized: u.Quantity[x_unit]) -> u.Quantity[y_unit]:
        x = x_unsanitized.to(x_unit).value
        log_x = np.log(x)

        return np.exp(log_interpolator(log_x)) * y_unit

    return interpolator


class Material:
    """A representation of a compound material.

    Parameters
    ----------
    A : `~astropy.units.Quantity`
        An array of the atomic weights of the elements in the compound in units
        convertible to atomic mass units.

    Z : `~astropy.units.Quantity`
        An array of the atomic numbers of the elements making up the compound.
        Does not have any units.

    fractional_weights : `~astropy.units.Quantity`
        An array containing the fractional composition of the compound. For
        example, if the value at index 0 has a value of 0.33, that element would
        make up 33% of the compound by weight.

    density : `~astropy.units.Quantity`
        The density of the material in units convertible to kilograms per cubic
        meter.

    name : str, optional
        The name of the material.

    STAR_data : `~plasmapy.plasma.material.MaterialSTARData`, optional
        The stopping and ranging (STAR) data of the material. An iterable of
        these structures can be provided if multiple incident particle species
        are of interest.
    """

    @validate_quantities(A={"equivalencies": u.molar_mass_amu()})
    def __init__(
        self,
        A: u.Quantity[u.kg / u.mol],
        Z: u.Quantity[u.dimensionless_unscaled],
        fractional_weights: u.Quantity[
            u.dimensionless_unscaled
        ],  # TODO: Make sure these sum to unity
        density: u.Quantity[u.kg / u.m**3],
        name: str | None = None,
        STAR_data: list[MaterialSTARData]
        | MaterialSTARData
        | None = None,  # User can provide STAR data for multiple incident species
    ):
        self.A = A
        self.Z = Z
        self.fractional_weights = fractional_weights
        self.density = density
        self.name = name

        self.interpolators = {}

        if STAR_data is not None:
            self._process_STAR_data(STAR_data)

    def _process_STAR_data(
        self, STAR_data: list[MaterialSTARData] | MaterialSTARData | None
    ):
        match type(STAR_data):
            case type(None):
                return
            case type({}):
                # TODO: Make sure each of the provided incident particle specie is unique

                for incident_STAR_data in STAR_data:
                    self._initialize_interpolators(incident_STAR_data)
            case type(MaterialSTARData()):
                self._initialize_interpolators(STAR_data)
            case _:
                raise TypeError(
                    f"Expected `MaterialSTARData` (or iterable equivalent) for `STAR_data`"
                    f"(Got unknown type `{type(STAR_data)}`)"
                )

    def _initialize_interpolators(self, STAR_data: MaterialSTARData):
        interpolators = {}
        interpolators.stopping_power_interpolator = _construct_log_interpolator(
            STAR_data["energies"],
            STAR_data["stopping_power"],
        )
        interpolators.range_interpolator = _construct_log_interpolator(
            STAR_data["energies"],
            STAR_data["projected_range"],
        )
        interpolators.transmitted_energy_interpolator = _construct_log_interpolator(
            STAR_data["projected_range"],
            STAR_data["energies"],
        )

        # TODO: Ensure it's clear that these interpolators correspond to the particle
        #  specie associated with the provided `STAR_data`
        self.interpolators[STAR_data["incident_particle"]] = interpolators

    def _interpolator_check(self, interpolator: str):
        """
        Check to see if the user provided material data at instantiation. If
        so, return it.
        """
        if not hasattr(self, f"_{interpolator}"):
            raise RuntimeError(
                f"Attempted to access `{interpolator}` without material stopping and "
                "range data. Please provide a value for `STAR_data` at instantiation to "
                "use this method."
            )

        return getattr(self, f"_{interpolator}")

    @cached_property
    def stopping_power_interpolator(self):
        """
        The stopping power of the particles in this material at a given kinetic
        energy.
        """
        return self._interpolator_check("stopping_power")

    @cached_property
    def range_interpolator(self):
        """
        The mean range of the particles in this material at a provided kinetic
        energy.
        """

        return self._interpolator_check("range_interpolator")

    @cached_property
    def transmitted_energy_interpolator(self):
        """
        The mean kinetic energy of the particles that make it through the
        provided target thickness.
        """
        return self._interpolator_check("transmitted_energy_interpolator")
