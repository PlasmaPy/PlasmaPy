"""
Module containing functionality focused on the plasma dispersion function
:math:`Z(ζ)`.

.. deprecated::

    The location of this module has been moved, and will be removed in a
    subsequent release. Use
    `plasmapy.dispersion.dispersion_functions` instead.
"""
__all__ = ["plasma_dispersion_func", "plasma_dispersion_func_deriv"]

import warnings
from typing import Union

import astropy.units as u
import numpy as np

from plasmapy.dispersion import dispersion_functions
from plasmapy.utils.decorators import bind_lite_func, preserve_signature
from plasmapy.utils.exceptions import PlasmaPyFutureWarning


@preserve_signature
def plasma_dispersion_func_lite(zeta):
    r"""
    The :term:`lite-function` version of
    `~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func`.

    .. deprecated::

        The lite function for ``plasma_dispersion_func`` has been
        deprecated and will be removed in a future release. Use
        `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func`
        instead.
    """

    warnings.warn(
        (
            "The location of plasma_dispersion_func_lite has changed. Use "
            "plasmapy.dispersion.dispersion_functions.plasma_dispersion_func.lite"
            "instead."
        ),
        PlasmaPyFutureWarning,
    )

    return dispersion_functions.plasma_dispersion_func(zeta)  # coverage: ignore


@bind_lite_func(plasma_dispersion_func_lite)
def plasma_dispersion_func(
    zeta: Union[complex, np.ndarray, u.Quantity],
) -> Union[complex, np.ndarray, u.Quantity]:
    r"""
    Calculate the plasma dispersion function.

    .. deprecated::

       The location of this function has moved. Use
       `plasmapy.dispersion.dispersion_functions.plasma_dispersion_func`
       instead. This function will be removed from the old location in
       a subsequent release.
    """
    warnings.warn(
        (
            "The location of plasma_dispersion_func has changed. Use "
            "plasmapy.dispersion.dispersion_functions.plasma_dispersion_func"
            "instead."
        ),
        PlasmaPyFutureWarning,
    )

    return dispersion_functions.plasma_dispersion_func(zeta)  # coverage: ignore


@preserve_signature
def plasma_dispersion_func_deriv_lite(zeta):
    r"""
    The :term:`lite-function` version of
    `~plasmapy.dispersion.dispersionfunction.plasma_dispersion_func_deriv`.

    .. deprecated::

        The lite function for plasma_dispersion_func_deriv has been
        deprecated and will be removed in a future release. Use
        `~plasmapy.dispersion.dispersion_functions.plasma_dispersion_func_deriv`
        instead.
    """
    warnings.warn(
        (
            "The location of plasma_dispersion_func_deriv_lite has changed. Use"
            "`plasmapy.dispersion.dispersion_functions.plasma_dispersion_func_deriv_lite`"
            "instead."
        ),
        PlasmaPyFutureWarning,
    )

    return dispersion_functions.plasma_dispersion_func_deriv(zeta)  # coverage: ignore


@bind_lite_func(plasma_dispersion_func_deriv_lite)
def plasma_dispersion_func_deriv(
    zeta: Union[complex, np.ndarray, u.Quantity],
) -> Union[complex, np.ndarray, u.Quantity]:
    r"""
    Calculate the derivative of the plasma dispersion function.

    .. deprecated::

        The location of this function has moved. Use
        `plasmapy.dispersion.dispersion_functions.plasma_dispersion_func_deriv`
        instead. This function will be removed from the old location in
        a subsequent release.
    """

    warnings.warn(
        (
            "The location of plasma_dispersion_func_deriv has been moved. Use "
            "plasmapy.dispersion.dispersion_functions.plasma_dispersion_func_deriv"
            "instead."
        ),
        PlasmaPyFutureWarning,
    )

    return dispersion_functions.plasma_dispersion_func_deriv(zeta)  # coverage: ignore
