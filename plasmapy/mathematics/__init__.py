"""This module gathers highly theoretical mathematical formulas
relevant to plasma physics. Usually, those are used somewhere else in
the code but were deemed general enough for the mathematical apparatus
to be abstracted from the main function interface."""

from .mathematics import (plasma_dispersion_func,
                          plasma_dispersion_func_deriv,
                          Fermi_integral)
