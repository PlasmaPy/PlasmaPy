__all__ = ["HarrisSheet"]

import astropy.constants as const
import astropy.units as u
import math
import numpy as np


class HarrisSheet:

    """
    Define a Harris Sheet Equilibrium.

    Parameters
    ----------
    B0 : `~astropy.units.Quantity` 
         magnetic field.

    delta : 'float'
        Delta is the thickness of the sheet.

    P0 : 'float'
        Plasma Pressure.

    Notes
    -----
    A current sheet is current limited to a surface.

    Examples
    --------
    >>> harris_sheet = HarrisSheet(delta = 3*u.m ,B0 = 2*u.T)
    >>> harris_sheet.magnetic_field(y = 5*u.m)
    <Quantity 1.8622... T>

    """

    def __init__(self, B0, delta, P0=0 * u.Pa):
        self.B0 = B0
        self.delta = delta
        self.P0 = P0

    def magnetic_field(self, y):
        r"""
        Compute the magnetic field along y = 0.

        This equation uses the asymptotic magnetic field strength along with y=0.

        Delta provides the thickness of the sheet.

        .. math::

            B_x(y) = B_0 \tanh \left( \frac{y}{\delta} \right)

        Parameters
        ----------

        y : `float`
          The axis of reference.

        """
        return self.B0 * math.tanh(y / self.delta)

    def current_density(self, y):
        r"""
        Compute the current density.

        .. math::

          B_0/(\delta\mu_0)(\mathrm{sech}(y/\delta)^2)

        Parameters
        ----------

        y : `float`
          The axis of reference.

        """
        return self.B0 / (self.delta * const.mu0) * (np.cosh((y / self.delta) ** -2))

    def plasma_pressure(self, y):
        r"""
        Compute plasma preassure.

        .. math::

            B_0/(\delta*mu_0)(\mathrm{sech}(y/\delta)^2)

        Parameters
        ----------

        y : `float`
            The axis of reference.

        """
        return (
            self.B0**2 / (2 * const.mu0) * (np.cosh(y / self.delta) ** -2) + self.P0
        )
