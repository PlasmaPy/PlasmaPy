r"""
Collection of helper functions used to quantify intermittency in
turbulent plasma.

In laboratory plasma experiments, fluctuating magnetic fields are
often measured with magnetic pickup (B-dot) probes, which record
:math:`\dot{B}(t)`. Using Faraday's law and induced voltage, a
time integration step will give :math:`B(t)`. Both are useful
quantities when looking at intermittency in plasma; see
``plasmapy.analysis.magnetic_pickup_probe``.

Power spectra of these signals give a relative energy transfer rate
but miss the rare energetic bursts commonly developed in turbulent
plasmas. These intermittent fluctuations are hypothesized to reflect
physical coherent structures in plasma like current sheets and are
candidates for reconnection sites, which therefore gives credence to
measuring and quantifying intermittency in plasma through B-dot
probe data [1]. This module provides methods to calculate
intermittency from a B-dot probe dataset in the process outlined
below. This module is based on the analysis described in [1] and 
refrence implementation in [4].

The first step of the intermittency analysis is introducing a
timescale into the time series by taking differences of the data
points at different timesteps. Therefore, we construct a series of
:math:`\Delta \dot{B}_\tau` increments for a variety of timesteps
:math:`\tau` to examine small scale structure and isolated changes
at a scale :math:`\tau`. The following is the increment equation
(equation (13) of [1]):

.. math::

    \Delta \dot{B}_{\tau}(t) = \dot{B}(t + \tau) - \dot{B}(t)

A histogram of these increments, normalized to their root mean
square value and total integral of the histogram, is then
constructed and analyzed. A Gaussian function is fit to these
distributions and a noticeable transition from Gaussian structure at
large :math:`\tau` to non-Gaussian at small :math:`\tau` should be
realized. These non-Gaussian tails at small :math:`\tau` indicate
intermittency of the signal, with fat tails signifying large
fluctuations in the values of the time series and the timestep
indicating relative temporal size of intermittent events [1]. A
histogram plotting function is defined below, although it is
commented out as it really only provides a basic qualitative
overview of the data. A more insightful analysis comes from reducing
the probability distribution function to a quantitative analysis by
Flatness. We more adeptly can plot the level of intermittency as a
function of :math:`\tau` through Flatness (equation (14) of [1]):

.. math::

    F(\tau) = \frac{S_4(\tau)}{S_2(\tau)^{2}}

where :math:`S_p(\tau)` is the :math:`p`\ th order moment for the
distribution of increments :math:`\Delta \dot{B}` calculated in the
following (equation (15) of [1]):

.. math::

    S_p(\tau) = \int_{-x_{\max}}^{x_{\max}}
    P(\Delta \dot{B}) (\Delta \dot{B})^{p} \, dx
    = \langle (\Delta \dot{B}_{\tau}(t))^{p} \rangle

Evidently, the flatness can be equivalently calculated as a moment
of the PDF or directly through an average of the time series of
increments [1], and the methods below use this direct average. For 
refrence, an exact Gaussian distribution has a flatness equal to 3 
[1] [3].

We have separate methods for the structure function [2] and
flatness calculations, but in overview an increase in flatness
indicates an increase in the intermittency of the probability
density function.

References
----------
.. [1] D. A. Schaffner, V. S. Lukin, A. Wan, and M. R. Brown,
   "Turbulence analysis of an experimental flux-rope plasma",
   Plasma Phys. Control. Fusion **56**, 064003 (2014).
   https://doi.org/10.1088/0741-3335/56/6/064003

.. [2] D. A. Schaffner and M. R. Brown, "Multifractal and monofractal
   scaling in a laboratory magnetohydrodynamic turbulence experiment",
   Astrophys. J. **811**, 61 (2015).
   https://doi.org/10.1088/0004-637X/811/1/61

.. [3] U. Frisch, *Turbulence: The Legacy of A. N. Kolmogorov*
   (Cambridge University Press, 1995).

.. [4] D. A. Schaffner, Bryn Mawr Plasma Laboratory (BMPL)
   intermittency analysis scripts.
   https://github.com/dschaffner/BMPL/tree/master/InDevelopment
"""

__all__ = ["flatness", "increment", "structure_function"]

import numpy as np


def increment(timeseries, dt: float, timestep: float) -> np.ndarray:
    r"""
    Calculate the increments of a time series at a given timestep.

    Parameters
    ----------
    timeseries : 1D array 
        Time series array of raw B-dot measurements

    dt : float
        Sampling time of ``timeseries``, in the same units as
        ``timestep``.

    timestep : float
        Time step :math:`\tau` of the increments, in the same units
        as ``dt``.

    Returns
    -------
    `numpy.ndarray`
        Series of increments :math:`\Delta \dot{B}_\tau` of
        ``timeseries``, constructed at the given ``timestep`` tau. 

    Examples
    --------
    >>> from plasmapy.analysis.time_series.intermittency import increment
    >>> increment([1.0, 2.0, 4.0, 7.0, 11.0], dt=1, timestep=1)
    array([1., 2., 3., 4.])
    >>> increment([1.0, 2.0, 4.0, 7.0, 11.0], dt=1, timestep=2)
    array([3., 5., 7.])
    >>> increment([0.0, 1.0, 3.0, 6.0], dt=2.5e-9, timestep=5e-9)
    array([3., 5.])
    """
    # Ensure sampling time is positive 
    if dt <= 0:
        msg = "dt must be positive"
        raise ValueError(msg)

    # Ensure that input is all of type float array 
    timeseries = np.asarray(timeseries, dtype=float)

    # Ensure input array is 1D 
    if timeseries.ndim != 1:
        msg = "timeseries must be one dimensional"
        raise ValueError(msg)

    # Convert the physical time lag tau into a whole number of
    # samples. Uses round to find closest possible float, and then
    # converts to int, as slices in future require int indicies 
    indexstep = int(np.round(timestep / dt))

    # Lower bound prevents issue of lag shorter than half a sample. 
    # --> rounds indexstep to 0 and causes empty slice error
    # Upper bound prevents lag of our sample recording length 
    # --> returns empty array  
    if not 1 <= indexstep < len(timeseries):
        msg = "timestep/dt must round to a value between 1 and len(timeseries) - 1"
        raise ValueError(msg)

    # b(t + tau) - b(t) for every t at once. Eqn 13 [1] the
    # first slice is the series shifted forward by the lag b(t + tau)
    # second is the series with the last lag samples cut off, b(t)
    # their elementwise difference is the full set of increments. 
    return timeseries[indexstep:] - timeseries[:-indexstep]


# For reference, the original implementation of generating increments [4]
# built the same increments one at a time in a while-loop, which is much
# less efficient in practice. It also drops the final increment, since its 
# loop condition stops one iteration early
#
# def generate_deltas_list(timeseries,dt,timestep):
#     deltas=[]
#     indexstep = np.round(timestep/dt)
#     total = len(timeseries)
#     initial_index = 0
#     final_index = int(indexstep)
#     initial=timeseries[initial_index]
#     final=timeseries[final_index]
#     numloops = 0
#     while initial_index < (total-(indexstep+1)):
#         numloops+=1
#         delta = final-initial
#         deltas.append(delta)
#         initial_index+=1
#         final_index+=1
#         initial=timeseries[initial_index]
#         final=timeseries[final_index]
#     return deltas


def structure_function(timeseries, dt: float, timestep: float, moments=(1, 2, 3, 4, 5, 6)) -> np.ndarray:
    r"""
    Calculate structure functions of a time series at each timestep.

    Parameters
    ----------
    timeseries : 1D array
        Time series array of raw B-dot measurements

    dt : float
        Sampling time of ``timeseries``, in the same units as
        ``timesteps``.

    timestep : float
        Time step :math:`\tau` of the increments, in the same units
        as ``dt`` (see `increment`).

    orders : 1D array, default: ``(1, 2, 3, 4, 5, 6)``
        Orders :math:`p` of the structure functions. 

    Returns
    -------
    `numpy.ndarray`
        Structure function values :math:`S_p(\tau)` of ``timeseries``
        at the given ``timestep``, one value per element of ``orders``. 
    """
    # Ensure orders is array of type float 
    moments = np.asarray(moments, dtype=float)

    # Ensure orders are positive 
    if np.any(moments <= 0):
        msg = "all elements of moments must be positive"
        raise ValueError(msg)

    # Increments of the series @ certain timestep tau.
    # Eqn 13 [1], timestep verified within increment function
    deltas = increment(timeseries, dt, timestep)

    # Center increments so S_p are central moments
    deltas = np.abs(deltas - np.mean(deltas))

    # S_p = <|delta|^p> --> direct time series average outlined 
    # in eqn 15 [1]. Returns array with each Structure function
    # as an element 
    return np.array([np.mean(deltas**p) for p in moments])

def flatness(timeseries, dt: float, timestep: float) -> float:
    r"""
    Calculate the flatness of the increments at a given timestep.

    Parameters
    ----------
    timeseries : 1D array
        Time series to array of raw B-dot measurements.

    dt : float
        Sampling time of ``timeseries``, in the same units as
        ``timesteps``.

    timestep : float
        Time step :math:`\tau` of the increments, in the same units
        as ``dt`` (see `increment`).

    Returns
    -------
    float
        Flatness :math:`F(\tau)` of the increments of ``timeseries``
        at the given ``timestep``.
    """
    # second and fourth order structure functions from set of increments
    s_2, s_4 = structure_function(timeseries, dt, timestep, orders=(2,4))

    # F = S_4 / S_2**2 (equation (14)); an exact Gaussian signal gives
    # F = 3, and values above 3 signal the fat tails of intermittency
    return float(s_4/s_2**2)


# The histogram, qualitative view of the increments is pasted below,
# kept from the original script this module is based off of. 
# It gives the qualitative fat-tail picture (figure 8 of [1]),
# while the functions above get the same moments exactly from the
# direct average of equation (15).
#
# def increment_histogram(timeseries, dt, timestep, nbins=300):
#     # generate list of deltas based on time interval
#     deltas = increment(timeseries, dt, timestep)
#     # convert deltas list to array
#     deltas = np.array(deltas)
#     # center deltas
#     deltas_mean = np.mean(deltas)
#     deltas = deltas-deltas_mean
#     # compute min and max of deltas (for histogramming)
#     deltas_min = np.floor(np.min(deltas))
#     deltas_max = np.ceil(np.max(deltas))
#     # calculate histogram
#     hist, bins = np.histogram(deltas,nbins,range=[deltas_min,deltas_max])
#     # shift bins to be centered
#     bins_centered = 0.5*(bins[1:]+bins[:-1])
#     # compute normalized S2 moment
#     S2 = np.sum(hist*bins_centered**2)/np.sum(hist)
#     S4 = np.sum(hist*bins_centered**4)/np.sum(hist)
#     # compute Flatness/Kurtosis
#     hist_flatness = S4/(S2)**2
#     return hist, bins_centered, hist_flatness
