import numpy as np

# Define finite difference coefficients
# https://en.wikipedia.org/wiki/Finite_difference_coefficient
coeffs = {
    # Central difference
    'central':
        # 1st deriv, 2nd and 4th order accuracy
        [{2: [(-1, -0.5), (1, 0.5)],
          4: [(-2, 1/12), (-1, -2/3), (1, 2/3), (-1/12)]},
         # 2nd deriv, 2nd and 4th order accuracy
         {2: [(-1, 1), (0, -2), (1, 1)],
          4: [(-2, -1/12), (-1, 4/3), (0, -5/2), (1, 4/3), (2, -1/12)]},
         # 3rd deriv, 2nd and 4th order accuracy
         {2: [(-2, -1/2), (-1, 1), (1, -1), (2, 1/2)],
          4: [(-3, 1/8), (-2, -1), (-1, 13/8), (1, -13/8), (2, 1), (3, -1/8)]}],
    'forward':
        # 1st deriv, 1st - 3rd order accuracy
        [{1: [(0, -1), (1, 1)],
          2: [(0, -3/2), (1, 2), (2, -1/2)],
          3: [(0, -11/6), (1, 3), (2, -3/2), (3, 1/3)]},
         # 2nd deriv, 1st - 3rd order accuracy
         {1: [(0, 1), (1, -2), (2, 1)],
          2: [(0, 2), (1, -5), (2, 4), (3, -1)],
          3: [(0, 35/12), (1, -26/3), (2, 19/2), (3, -14/3), (11/12)]},
         # 3rd deriv, 1st - 3rd order accuracy
         {1: [(0, -1), (1, 3), (2, -3), (3, 1)],
          2: [(0, -5/2), (1, 9), (2, -12), (3, 7), (4, -3/2)],
          3: [(0, -17/4), (1, 71/4), (2, -59/2), (3, 49/2), (4, -41/4), (5, 7/4)]}],
    'backward':
        # 1st deriv, 1st and 2nd order accuracy
        [{1: [(0, 1), (-1, -1)],
          2: [(0, 3/2), (-1, -2), (-2, 1/2)]},
         # 2nd deriv, 1st - 3rd order accuracy
         {1: [(0, 1), (-1, -2), (-2, 1)],
          2: [(0, 2), (-1, -5), (-2, 4), (-3, -1)]},
         # 3rd deriv, 1st - 3rd order accuracy
         {1: [(0, 1), (-1, -3), (-2, 3), (-3, -1)],
          2: [(0, 5/2), (-1, -9), (-2, 12), (-3, -7), (-4, 3/2)]}]
}


def shift(f_shape, shift, axis, n_ghosts=4):
    """
    Return slice objects to correctly index padded array to get shifted array.
    """
    indices = [slice(None, None, None)] * len(f_shape)
    if shift-n_ghosts == 0:
        indices[axis] = slice(n_ghosts+shift, None, None)
    else:
        indices[axis] = slice(n_ghosts+shift, shift-n_ghosts, None)
    return indices


def set_boundaries(f):
    """
    Set values for boundaries of domain to constant value
    """
    newf = np.zeros(f.shape) * f.unit
    indices = [slice(1, -1, None)] * len(f.shape)
    newf[indices] = f[indices]

    return newf


class Solver():
    def __init__(self, dx, method='central', deriv=1, acc=2):
        """
        Create a callable instance to calculate derivatives of order `deriv`.
        Sets up the solver to use forward, backward or central differences for
        the specified derivative at accuracy corresponding to `acc`.

        Parameters
        ----------

        dx : tuple of floats
            Tuple of spatial step-sizes corresponding to each dimension of the
            simulation domain.
        method : str {'central' | 'forward' | 'backward'}
            Flavour of finite difference method to use.
        deriv : int {1 | 2 | 3}
            Order of the derivative to calculate with this Solver instance.
            Default is first derivative.
        acc : int
            Order of accuracy of derivative.
            Valid values are:
            - 2 or 4 for central difference;
            - 1, 2 or 3 for forward difference;
            - 1 or 2 for backward difference.
        """
        assert method in coeffs.keys(), """
            Invalid spatial solver method - must be one of {}
            """.format(coeffs.keys())

        # Define number of ghost cells on either side of the solution grid
        self.n_ghosts = 4

        # Define step-size and finite difference coefficients for spatial
        # differentiation
        self.dx = dx
        self.deriv = deriv
        self.method = method
        self.coeffs = coeffs[method][deriv-1][acc]

    def __call__(self, f, axis):
        """
        Calculates the derivative of function f with respect to the specified
        axis.

        Parameters
        ----------

        f : ndarray (x, y, z)
            Values of some function f(x, y, z) at every point in a 3D range.
        axis : int [0 | 1 | 2]
            Direction in which to calculate the derivative. 0, 1 and 2
            correspond to the x, y and z axes, respectively.
        """
        assert axis in [0, 1, 2], """
            Invalid axis identifier passed to spatial solver."""

        # Nest supplied array inside a larger array with n_ghosts ghost cells
        # either side in the direction of differentiation
        padding = [(0, 0)] * len(f.shape)
        if self.method == 'central':
            padding[axis] = (self.n_ghosts, self.n_ghosts)
        elif self.method == 'forward':
            padding[axis] = (0, self.n_ghosts*2)
        f = np.pad(f, padding, 'edge') * f.unit

        # Differentiate the array
        dfdx = sum([f[shift(f.shape, c[0], axis)] * c[1] for c in self.coeffs])\
            / (self.dx[axis] ** self.deriv)

        return dfdx
