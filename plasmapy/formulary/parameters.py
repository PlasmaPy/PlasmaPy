"""Functions to calculate fundamental plasma parameters."""

__all__ = [
    "magnetic_energy_density",
    "magnetic_pressure",
]
__aliases__ = [
    "pmag_",
    "ub_",
]
__lite_funcs__ = []

import astropy.units as u

from astropy.constants.si import e, eps0, k_B, mu0

from plasmapy.utils.decorators import deprecated, validate_quantities
from plasmapy.utils.exceptions import PlasmaPyFutureWarning

from plasmapy.formulary import dimensionless, frequencies, lengths, misc, speeds  # noqa


__all__ += (
    dimensionless.__all__.copy()
    + frequencies.__all__.copy()
    + lengths.__all__.copy()
    + misc.__all__.copy()
    + speeds.__all__.copy()
    + __aliases__
    + __lite_funcs__
)

__aliases__ += (
    dimensionless.__aliases__.copy()
    + frequencies.__aliases__.copy()
    + lengths.__aliases__.copy()
    + misc.__aliases__.copy()
    + speeds.__aliases__.copy()
)

__lite_funcs__ += frequencies.__lite_funcs__.copy() + speeds.__lite_funcs__.copy()

# remove from __all__ and __aliases__ added functionality names that are not
# actually in this file
for name in (
    "beta",
    "betaH_",
    "Mag_Reynolds",
    "quantum_theta",
    "Re_",
    "Reynolds_number",
    "Rm_",
):  # coverage: ignore
    try:
        __aliases__.remove(name)
    except ValueError:
        # name was not in __aliases__
        pass

    try:
        __all__.remove(name)
    except ValueError:
        # name was not in __all__
        pass

e_si_unitless = e.value
eps0_si_unitless = eps0.value
k_B_si_unitless = k_B.value

funcs_to_deprecate_wrap = [  # (module_name, func_name)
    ("dimensionless", "Debye_number"),
    ("dimensionless", "nD_"),
    ("dimensionless", "Hall_parameter"),
    ("dimensionless", "betaH_"),
    ("frequencies", "gyrofrequency"),
    ("frequencies", "oc_"),
    ("frequencies", "wc_"),
    ("frequencies", "plasma_frequency"),
    ("frequencies", "plasma_frequency_lite"),
    ("frequencies", "wp_"),
    ("frequencies", "lower_hybrid_frequency"),
    ("frequencies", "wlh_"),
    ("frequencies", "upper_hybrid_frequency"),
    ("frequencies", "wuh_"),
    ("lengths", "Debye_length"),
    ("lengths", "lambdaD_"),
    ("lengths", "gyroradius"),
    ("lengths", "rc_"),
    ("lengths", "rhoc_"),
    ("lengths", "inertial_length"),
    ("lengths", "cwp_"),
    ("misc", "_grab_charge"),
    ("misc", "Bohm_diffusion"),
    ("misc", "DB_"),
    ("misc", "mass_density"),
    ("misc", "rho_"),
    ("misc", "thermal_pressure"),
    ("misc", "pth_"),
    ("speeds", "Alfven_speed"),
    ("speeds", "va_"),
    ("speeds", "ion_sound_speed"),
    ("speeds", "cs_"),
    ("speeds", "thermal_speed"),
    ("speeds", "thermal_speed_coefficients"),
    ("speeds", "thermal_speed_lite"),
    ("speeds", "vth_"),
    ("speeds", "kappa_thermal_speed"),
    ("speeds", "vth_kappa_"),
]
for modname, name in funcs_to_deprecate_wrap:
    globals()[name] = deprecated(
        since="0.7.0",
        warning_type=PlasmaPyFutureWarning,
        message=(
            f"The {name}() function has been moved to "
            f"plasmapy.formulary.{modname}.  Update your import to get "
            f"rid of this warning.  The 'plasmapy.formulary.parameters' module "
            f"will be officially removed in release v0.9.0."
        ),
    )(getattr(globals()[f"{modname}"], name))

del modname, name


@validate_quantities
def magnetic_pressure(B: u.T) -> u.Pa:
    r"""
    Calculate the magnetic pressure.

    **Aliases:** `pmag_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field in units convertible to tesla.

    Returns
    -------
    p_B : `~astropy.units.Quantity`
        The magnetic pressure in units in pascals (newtons per square meter).

    Raises
    ------
    `TypeError`
        If the input is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the input is not in units convertible to tesla.

    `ValueError`
        If the magnetic field strength is not a real number between
        :math:`±∞`\ .

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed.

    Notes
    -----
    The magnetic pressure is given by:

    .. math::
        p_B = \frac{B^2}{2 μ_0}

    The motivation behind having two separate functions for magnetic
    pressure and magnetic energy density is that it allows greater
    insight into the physics that are being considered by the user and
    thus more readable code.

    See Also
    --------
    magnetic_energy_density : returns an equivalent `~astropy.units.Quantity`,
        except in units of joules per cubic meter.

    Examples
    --------
    >>> from astropy import units as u
    >>> magnetic_pressure(0.1*u.T).to(u.Pa)
    <Quantity 3978.87... Pa>

    """
    return (B ** 2) / (2 * mu0)


pmag_ = magnetic_pressure
"""Alias to `~plasmapy.formulary.parameters.magnetic_pressure`."""


@validate_quantities
def magnetic_energy_density(B: u.T) -> u.J / u.m ** 3:
    r"""
    Calculate the magnetic energy density.

    **Aliases:** `ub_`

    Parameters
    ----------
    B : `~astropy.units.Quantity`
        The magnetic field in units convertible to tesla.

    Returns
    -------
    E_B : `~astropy.units.Quantity`
        The magnetic energy density in units of joules per cubic meter.

    Raises
    ------
    `TypeError`
        If the input is not a `~astropy.units.Quantity`.

    `~astropy.units.UnitConversionError`
        If the input is not in units convertible to tesla.

    `ValueError`
        If the magnetic field strength does not have an appropriate.
        value.

    Warns
    -----
    : `~astropy.units.UnitsWarning`
        If units are not provided, SI units are assumed

    Notes
    -----
    The magnetic energy density is given by:

    .. math::
        E_B = \frac{B^2}{2 μ_0}

    The motivation behind having two separate functions for magnetic
    pressure and magnetic energy density is that it allows greater
    insight into the physics that are being considered by the user and
    thus more readable code.

    See Also
    --------
    magnetic_pressure : Returns an equivalent `~astropy.units.Quantity`,
        except in units of pascals.

    Examples
    --------
    >>> from astropy import units as u
    >>> magnetic_energy_density(0.1*u.T)
    <Quantity 3978.87... J / m3>

    """
    return magnetic_pressure(B)


ub_ = magnetic_energy_density
"""Alias to `~plasmapy.formulary.parameters.magnetic_energy_density`."""
