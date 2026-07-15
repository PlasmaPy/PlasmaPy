"""Functionality for calculating magnetic field from magnetic pickup (B-dot) probe."""

__all__ = ["magnetic_field_calculation"]

import numpy as np
from scipy.integrate import cumulative_trapezoid


def magnetic_field_calculation(
    bdot_voltage: np.ndarray,
    time_array: np.ndarray,
    loop_area: float,
    num_loop: int = 1,
    gain: float = 1.0,
) -> np.ndarray:
    r"""
    Take the voltage output of a magnetic pickup probe (Bdot probe) and
    time base and return the associated array of the magnetic field as
    a function of time.

    Parameters
    ----------
    bdot_voltage: `numpy.ndarray`
        1-D array containing voltage fluctuations from a bdot probe
        (should be in units of volts). The array values are proportional
        to the time changing magnetic field via Faraday's law.

    time_array: `numpy.ndarray`
        1-D time series of the data collection (should be in seconds).

    loop_area: `float`
        The area through which the changing flux is measured (should be
        in square meters).

    num_loop: `int`
        The number of loops of the Bdot probe, each with corresponding
        ``loop_area``. (Default: 1)

    gain: `float`
        Any dimensionless gains that need to be applied. (Default: 1.0)

    Returns
    -------
    magnetic_field_array: `numpy.ndarray`
        The array containing the magnetic field (in tesla for SI
        inputs) as a function of time.

    Raises
    ------
    ValueError
        If ``bdot_voltage`` and ``time_array`` do not have equal sizes.

    Notes
    -----
    The integral form of Faraday's law states that the time change in
    flux through an area generates a voltage around a loop,

    .. math::

        V = -N\frac{d\Phi}{dt}

    where :math:`\Phi` is the magnetic flux

    .. math::

        \Phi = \int{BdA}

    and :math:`N` is the number of loops.

    A magnetic pickup probe works by providing a fixed area loop or
    loops using a length of conducting wire through which magnetic
    fields penetrate. As these fields change in time, a voltage is
    generated that can be measured by some kind of data acquisition
    device (such as an oscilloscope). We can also assume that the loop
    is small enough that the magnetic field through the loop is
    constant everywhere in the loop (and only changes in time). Given
    this, we can write

    .. math::

        V = -NA\frac{dB}{dt}

    where the voltage becomes proportional only to the change in
    magnetic field and the area of the loop. Magnetic field as a
    function of time can be recovered by time integrating
    :math:`-V(t)/(NAg)`, where :math:`g` is any dimensionless gain
    applied to the measured voltage.

    Examples
    --------
    >>> import numpy as np
    >>> bdot_voltage = np.array([1.0, 1.0, 1.0])
    >>> time_array = np.array([0.0, 0.5, 1.0])
    >>> magnetic_field_calculation(bdot_voltage, time_array, loop_area=1.0)
    array([ 0. , -0.5, -1. ])
    """
    if len(bdot_voltage) != len(time_array):
        msg = "sizes of time_array and bdot_voltage are not equal."
        raise ValueError(msg)

    scaled_voltage = -np.asarray(bdot_voltage) / (loop_area * num_loop * gain)

    return cumulative_trapezoid(scaled_voltage, time_array, initial=0)
