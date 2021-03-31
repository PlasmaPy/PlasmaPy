"""Utilities to help with testing."""

__all__ = [
    "assert_can_handle_nparray",
]

import astropy.tests.helper as astrohelper
import astropy.units as u
import inspect
import numpy as np
import warnings

from plasmapy.tests.helpers.exceptions import (
    InvalidTestError,
    UnexpectedExceptionFail,
    UnexpectedResultFail,
)
from plasmapy.utils.code_repr import call_string
from plasmapy.utils.exceptions import PlasmaPyWarning


def assert_can_handle_nparray(
    function_to_test, insert_some_nans=None, insert_all_nans=None, kwargs=None,
):
    """
    Test for ability to handle numpy array quantities.

    Parameters
    ----------
    function_to_test
        The function to be tested for ability to handle numpy array quantities.
        Arguments are automatically given a vector input based on their
        variable name. Current args that are interpreted as vectors are:
        `["T", "T_i", "T_e", "temperature"]`
        `["n", "n_i", "n_e", "density"]`
        `["B"]`
        `["V", "Vperp"]`
        `["coulomb_log"]`
        `["characteristic_length"]`

    insert_some_nans: `list`
        List of argument names in which to insert some np.nan values.
        These must be arguments that will be tested as vectors as listed
        above.

    insert_all_nans: `list`
        List of argument names to fill entirely with np.nan values.

    kwargs: `dict`
        Arguments to pass directly to the function in under test, in the
        normal kwargs python dictionary format.

    Raises
    ------
    ValueError
        If this function cannot interpret a parameter of function_to_test.

    Examples
    --------
    >>> from plasmapy.formulary.parameters import Alfven_speed, gyrofrequency
    >>> assert_can_handle_nparray(Alfven_speed)
    >>> assert_can_handle_nparray(gyrofrequency, kwargs={"signed": True})
    >>> assert_can_handle_nparray(gyrofrequency, kwargs={"signed": False})
    """

    if insert_some_nans is None:
        insert_some_nans = []

    if insert_all_nans is None:
        insert_all_nans = []

    if kwargs is None:
        kwargs = {}

    def _prepare_input(
        param_name, param_default, insert_some_nans, insert_all_nans, kwargs
    ):
        """
        Parse parameter names and set up values to input for 0d, 1d, and 2d array tests.
        """
        # first things first: let any passed in kwarg right through (VIP access)
        if param_name in kwargs.keys():
            return (kwargs[param_name],) * 4

        # else, if it's a recognized variable name, give it a reasonable unit and magnitude
        elif param_name in ["particle", "ion_particle", "ion"]:
            if not (param_default is inspect._empty or param_default is None):
                return (param_default,) * 4
            else:
                return ("p",) * 4
        elif param_name == "particles" or param_name == "species":
            if not (param_default is inspect._empty):
                return (param_default,) * 4
            else:
                return (("e", "p"),) * 4
        elif param_name in ["T", "T_i", "T_e", "temperature"]:
            unit = u.eV
            magnitude = 1.0
        elif param_name in ["n", "n_i", "n_e", "density"]:
            unit = u.m ** -3
            magnitude = 1e20
        elif param_name == "B":
            unit = u.G
            magnitude = 1e3
        elif param_name in ["V", "Vperp"]:
            unit = u.m / u.s
            magnitude = 1e5
        elif param_name == "coulomb_log":
            unit = 1.0
            magnitude = 1e1
        elif param_name == "characteristic_length":
            unit = u.m
            magnitude = 1.0
        elif param_name == "k":
            unit = u.m ** -1
            magnitude = 1.0

        # else, last resort, if it has a default argument, go with that:
        elif not (param_default is inspect._empty):
            return (param_default,) * 4

        else:
            raise ValueError(f"Unrecognized function input: {param_name}")

        # now knowing unit and magnitude, set up the 0d, 1d, 2d, and 3d arrays:
        input_data_3d = np.reshape(np.arange(1.0, 9.0, 1.0), (2, 2, 2))
        input_data_2d = np.reshape(np.arange(1.0, 5.0, 1.0), (2, 2))
        input_data_1d = np.arange(1.0, 5.0, 1.0)
        if param_name in insert_some_nans:
            input_data_3d[0, 0, 1] = np.nan
            input_data_3d[0, 1, 0] = np.nan
            input_data_2d[0, 1] = np.nan
            input_data_2d[1, 0] = np.nan
            input_data_1d[1] = np.nan
        elif param_name in insert_all_nans:
            input_data_3d = np.ones((2, 2, 2)) * np.nan
            input_data_2d = np.ones((2, 2)) * np.nan
            input_data_1d = np.ones(4) * np.nan
        input_data_3d *= magnitude
        input_data_3d *= unit
        input_data_2d *= magnitude
        input_data_2d *= unit
        input_data_1d *= magnitude
        input_data_1d *= unit
        input_data_0d = input_data_1d[3]
        return input_data_0d, input_data_1d, input_data_2d, input_data_3d

    # call _prepare_input to prepare 0d, 1d, and 2d sets of arguments for the function:
    function_sig = inspect.signature(function_to_test)
    function_params = function_sig.parameters
    args_0d = dict()
    args_1d = dict()
    args_2d = dict()
    args_3d = dict()
    param_names = [elm for elm in function_params.keys()]
    for idx, key in enumerate(function_params):
        args_0d[key], args_1d[key], args_2d[key], args_3d[key] = _prepare_input(
            param_names[idx],
            function_params[key].default,
            insert_some_nans,
            insert_all_nans,
            kwargs,
        )

    # call the function with the prepared argument sets:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=PlasmaPyWarning)
        result_0d = function_to_test(**args_0d)
        result_1d = function_to_test(**args_1d)
        result_2d = function_to_test(**args_2d)
        result_3d = function_to_test(**args_3d)

    # assert that the 1d, 2d, 3d versions get the same result (elementwise) as the 0d version:
    # (if the function returns multiple values, loop through and test each)
    try:
        scalar_testable = result_0d.value
    except AttributeError:
        scalar_testable = result_0d
    if np.isscalar(scalar_testable):
        astrohelper.assert_quantity_allclose(result_0d, result_1d[3])
        astrohelper.assert_quantity_allclose(result_0d, result_2d[1, 1])
        astrohelper.assert_quantity_allclose(result_0d, result_3d[0, 1, 1])
    else:
        for idx, res_0d in enumerate(result_0d):
            astrohelper.assert_quantity_allclose(res_0d, result_1d[idx][3])
            astrohelper.assert_quantity_allclose(res_0d, result_2d[idx][1, 1])
            astrohelper.assert_quantity_allclose(res_0d, result_3d[idx][0, 1, 1])
