"""
.. attention::

   The `plasmapy.formulary.parameters` module is deprecated and will be
   removed in the release of ``v0.9.0``.  Please update your code to
   import directly from `plasmapy.formulary` or the specific submodule
   the functionality is defined in.  For example,

   .. code-block:: python

      from plasmapy.formulary import Alfven_speed

      # or

      from plasmapy.formulary.speeds import Alfven_speed

Functions to calculate fundamental plasma parameters.
"""

__all__ = []

import contextlib

from plasmapy.utils.decorators import deprecated
from plasmapy.utils.exceptions import PlasmaPyFutureWarning

from plasmapy.formulary import dimensionless, frequencies, lengths, misc, speeds  # noqa


__all__ += (
    dimensionless.__all__.copy()
    + frequencies.__all__.copy()
    + lengths.__all__.copy()
    + misc.__all__.copy()
    + speeds.__all__.copy()
)

__aliases__ = (
    dimensionless.__aliases__.copy()
    + frequencies.__aliases__.copy()
    + lengths.__aliases__.copy()
    + misc.__aliases__.copy()
    + speeds.__aliases__.copy()
)

__lite_funcs__ = frequencies.__lite_funcs__.copy() + speeds.__lite_funcs__.copy()

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
    "Lundquist_number",
):  # coverage: ignore
    with contextlib.suppress(ValueError):  # ValueError if name was not in __aliases__
        __aliases__.remove(name)

    with contextlib.suppress(ValueError):  # ValueError if name was not in __all__
        __all__.remove(name)

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
    ("misc", "magnetic_energy_density"),
    ("misc", "ub_"),
    ("misc", "magnetic_pressure"),
    ("misc", "pmag_"),
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

del funcs_to_deprecate_wrap, modname, name
