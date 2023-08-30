import numpy as np
import scipy.special


class ForceFreeFluxRope:
    """
    Define Analytical solution for lundquist solution.

    Parameters
    ----------
    B0 : `~astropy.units.Quantity`
        Magnetic field strength in units convertible to tesla.

    a : `~astropy.units.Quantity`
        Eigenvalue to make :math:`\mathbf{J} × \mathbf{B} = 0`, in units 
        convertible to inverse length.

    Notes
    -----
    Lundquist solutions describe the flow of a conducting fluid in an electric field


    """

    def __init__(self, B0, a):
        self.B0 = B0
        self.a = a

    def B_theta(self, r):
        r"""
        Compute the magnetic field in the theta direction.

        .. math::

            (self.B0*scipy.special.j1(self.a * r))

        Parameters
        ----------
        r : `~astropy.units.Quantity`
             A distance from the flux rope.
        """
        return self.B0 * scipy.special.j1(self.a * r)

    def B_z(self, r):
        r"""
        Compute the magnetic field in the z direction.

        .. math::

            self.B0*scipy.special.j0(self.a * r)

        Parameters
        ----------
        r : `float`
          distance from equilibria

        """
        return self.B0 * scipy.special.j0(self.a * r)

    def B_magnitude(self, r):
        r"""
        Compute the total magnetic field.

        .. math::

            np.sqrt(self.B_z(r)**2 + self.B_theta(r)**2)

        Parameters
        ----------
        r : `float`
            distance from equilibria

        """
        return np.sqrt(self.B_z(r) ** 2 + self.B_theta(r) ** 2)
