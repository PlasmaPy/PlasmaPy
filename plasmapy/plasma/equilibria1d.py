"""Functionality for representing one-dimensional equilibria."""

__all__ = ["HarrisSheet"]

import astropy.constants as const
import astropy.units as u
import numpy as np

from plasmapy.utils.decorators.validators import validate_quantities


class HarrisSheet:
    """
    Define a Harris Sheet Equilibrium.
    Magnetic field will be in the :math:`\pm x` direction and
    the current density will be in the :math:`\pm z` direction,
    :math:`\hat{x} \times \hat{y} = \hat{z}` coordinate system.
A Harris sheet is a 1D ideal MHD equilibrium. In resistive MHD if there is any resistivity, it won't be a true equilibrium since the resistivity will gradually smooth the profile out over time. A Harris sheet is often used as the initial condition for simulations of magnetic reconnection.

    Parameters
    ----------
    B0 : `~astropy.units.Quantity`
         Magnitude of magnetic field in the limit of :math:`y → ∞` in units
         convertible to teslas.

    delta : `~astropy.units.Quantity`
        The thickness of the current sheet in units convertible to meters.

    P0 : `~astropy.units.Quantity`
        The plasma pressure in the limit of :math:`y → ∞` in units
        convertible to pascals.

    Notes
    -----
    A current sheet is current limited to a surface.

    Examples
    --------
    >>> import astropy.units as u
    >>> harris_sheet = HarrisSheet(delta = 3 * u.m, B0 = 2 * u.T)
    >>> harris_sheet.magnetic_field(y = 5 * u.m)
    <Quantity 1.8622... T>
    """

    def __init__(self, B0, delta, P0=0 * u.Pa):
        self.B0 = B0
        self.delta = delta
        self.P0 = P0

    @validate_quantities
    def magnetic_field(self, y: u.m) -> u.T:
        r"""
        Compute the magnetic field.

        This equation uses the asymptotic magnetic field strength along with y=0.

        Delta provides the thickness of the sheet.

        .. math::

            B_x(y) = B_0 \tanh \left( \frac{y}{\delta} \right)

        Parameters
        ----------
        y : `~astropy.units.Quantity`
           Orthogonal distance from the current sheet center.

        """
        return self.B0 * np.tanh(u.rad * y / self.delta)

    @validate_quantities
    def current_density(self, y: u.m) -> u.A / u.m**2:
        r"""
        Compute the current density.

        .. math::

          B_0/(\delta\mu_0)(\mathrm{sech}^2(y/\delta))

        Parameters
        ----------
        y : `~astropy.units.Quantity`
          Orthogonal distance from the current sheet center.

        """
        return (
            -self.B0 / (self.delta * const.mu0) * np.cosh(u.rad * y / self.delta) ** -2
        )

    @validate_quantities
    def plasma_pressure(self, y: u.m) -> u.Pa:
        r"""
        Compute plasma pressure.

        .. math::

            B_0/(\delta*mu_0)(\mathrm{sech}^2(y/\delta))

        Parameters
        ----------
        y : `~astropy.units.Quantity`
          Orthogonal distance from the current sheet center.

        """
        return (
            self.B0**2 / (2 * const.mu0) * (np.cosh(u.rad * y / self.delta) ** -2)
            + self.P0
        )
