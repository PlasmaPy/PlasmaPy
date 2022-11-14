"""
Created on Tue Aug  2 22:18:16 2022

@author: dschaffner
"""
import numpy as np

def check_sweep(
    variable: np.ndarray,
    strip_units: bool = True):

    """
    Function for checking that the delta array is properly
    formatted for analysis by `plasmapy.analysis.intermittency`.
    Parameters
    ----------
    variable: `numpy.ndarray`
        1D `numpy.ndarray` representing the value of delta(B)  *No units are assumed or
        checked at this time.*
    strip_units: `bool`
        (Default: `True`) If `True`, then the units on ``variable`` will be stripped if either are passed in as an Astropy
        `~astropy.units.Quantity`.
    Returns
    -------
    variable : `numpy.ndarray`
        Input argument ``variable`` after it goes through all of its checks
        and conditioning.
    Raises
    ------
    `TypeError`
        If the ``variable`` arrays are not instances of a
        `numpy.ndarray`.
    `ValueError`:
        If the ``variable`` arrays are not 1D.
    `ValueError`
        If either the ``variable`` array does not have a
        `numpy.dtype` of either `numpy.integer` or `numpy.floating`.
    """
    # -- examine variable array --
    # check type
    if isinstance(variable, np.ndarray):
        pass
    elif isinstance(variable, (list, tuple)):
        variable = np.array(variable)
    else:
        raise TypeError(
            f"Expected 1D numpy array for variable, but got {type(variable)}.",
        )

    # check array structure
    if not (
        np.issubdtype(variable.dtype, np.floating)
        or np.issubdtype(variable.dtype, np.integer)
    ):
        raise ValueError(
            f"Expected 1D numpy array of floats or integers for variable, but"
            f" got an array with dtype '{variable.dtype}'."
        )
    elif variable.ndim != 1:
        raise ValueError(
            f"Expected 1D numpy array for variable, but got array with "
            f"{variable.ndim} dimensions.",
        )

    #if isinstance(variable, u.Quantity) and strip_units:
        #variable = variable.value

    return variable


def generate_deltas_list(timeseries,dt,timestep):
    check_sweep(timeseries)
    deltas=[]
    indexstep = np.round(timestep/dt)
    total = len(timeseries)
    initial_index = 0
    final_index = int(indexstep)
    initial = timeseries[initial_index]
    final = timeseries[final_index]
    numloops = 0
    while initial_index < (total-(indexstep+1)):
        numloops+=1
        delta = final-initial
        deltas.append(delta)
        initial_index+=1
        final_index+=1
        initial = timeseries[initial_index]
        final = timeseries[final_index]
    return deltas
