"""
plasmapy.classes.simulation
============================

Classes and functionality for simulations.
"""

import numpy as np
import astropy.units as u
from .spatial_solvers import Solver
from ..constants import mu0


class MHDSimulation:
    r"""Physics class for magnetohydrodynamics.

    This class defines the MHD equations and implements time-stepping them:

    .. math::

       \frac{\partial \rho}{\partial t} + \nabla \cdot (\vec{v} \rho) = 0
       \frac{\partial (\rho \vec{v})}{\partial t} + \nabla \cdot (\vec{v} \rho \vec{v}) + \nabla p = 0
       \frac{\partial e}{\partial t} + \nabla \cdot (\vec{v} e + \vec{v} p) = 0

    for the fluid density, :math:`\rho`, momentum,
    :math:`\vec{m} = \vec{v} \rho`, energy :math:`e` and kinetic pressure
    :math:`p`. The pressure is a derived quantity defined as

    .. math::
       p = (\gamma - 1) (e - \frac{\rho \vec{v}^2}{2})

    Parameters
    ----------
    plasma : plasmapy.Plasma
        Plasma object describing variables that are being solved for, etc.
    gamma : float
        Value of the adiabatic index.

    Attributes
    ----------
    grid_size : tuple of ints
        Size of the simulation grid, as defined at initiation.
    gamma : float
        Adiabatic index for the simulation.
    """
    def __init__(self, plasma):
        """
        """
        self.dt = 1e-4 * u.s
        self.current_iteration = 0
        self.current_time = 0 * u.s
        self.plasma = plasma
        # Domain size
        # self.grid_size = grid_size

        # Physical parameters
        # self.gamma = gamma

        grids = (self.plasma.x.si, self.plasma.y.si, self.plasma.z.si)
        ranges = [grid for grid in grids if len(grid) > 1]
        stepsize = [range[1] - range[0] for range in ranges] * grids[0].unit
        self.solver = Solver(stepsize)

        # Collect equations into a nice easy-to-use list
        self.equations = [self._ddt_density, self._ddt_momentum,
                          self._ddt_energy, self._ddt_magfield]

    def time_stepper(self):
        """4th-order Runge-Kutta solver for stepping the simulation forward
        through time based on the equations defining the physics for the
        simulation.
        """
        half_dt = self.dt / 2
        kn = []
        derivs = []
        plasma = self.plasma
        orig_variables = [var.copy() for var in plasma.core_variables]

        for eq in self.equations:
            k1 = eq(self.current_time) * self.dt
            kn.append(k1)
            derivs.append(k1)
        for f, k1 in zip(plasma.core_variables, kn):
            np.add(f, k1/2, out=f)

        for i, (f, k1, eq) in enumerate(zip(plasma.core_variables, kn,
                                            self.equations)):
            k2 = eq(self.current_time + half_dt, f) * self.dt
            kn[i] = k2
            derivs[i] += (2 * k2)
        for f, f0, k2 in zip(plasma.core_variables, orig_variables, kn):
            np.add(f0, k2/2, out=f)

        for i, (f, k2, eq) in enumerate(zip(plasma.core_variables, kn,
                                            self.equations)):
            k3 = eq(self.current_time + half_dt, f) * self.dt
            kn[i] = k3
            derivs[i] += (2 * k3)
        for f, f0, k3 in zip(plasma.core_variables,
                             orig_variables, kn):
            np.add(f0, k3, out=f)

        for i, (f, k3, eq) in enumerate(zip(plasma.core_variables, kn,
                                            self.equations)):
            derivs[i] += eq(self.current_time+self.dt, f) * self.dt

        for f, f0, df in zip(plasma.core_variables,
                             orig_variables, derivs):
            np.add(f0, df/6, out=f)

        self.current_time += self.dt
        self.current_iteration += 1

    def _ddt_density(self, t, density=None):
        """
        """
        if not density:
            density = self.plasma.density
        return -div(self.plasma.velocity * density, self.solver)

    def _ddt_momentum(self, t, momentum=None):
        """
        """
        if not momentum:
            momentum = self.plasma.momentum
        v = self.plasma.velocity
        B = self.plasma.magnetic_field / np.sqrt(mu0)

        return (-grad(self.plasma.pressure, self.solver)
                - tensordiv(vdp(v, momentum) - vdp(B, B), self.solver))

    def _ddt_energy(self, t, energy=None):
        """
        """
        if not energy:
            energy = self.plasma.energy
        v = self.plasma.velocity
        B = self.plasma.magnetic_field / np.sqrt(mu0)

        return -div((v*energy) - (B * dot(B, v)) + (v*self.plasma.pressure),
                    self.solver)

    def _ddt_magfield(self, t, magfield=None):
        """
        """
        if not magfield:
            B = self.plasma.magnetic_field / np.sqrt(mu0)
        else:
            B = magfield / np.sqrt(mu0)
        v = self.plasma.velocity

        return -tensordiv(vdp(v, B) - vdp(B, v), self.solver) * np.sqrt(mu0)


def dot(vec1, vec2):
    r"""Calculates the dot product of two arrays of vector quantities.
    TODO: Replace this everywhere with the new NumPy way of doing this.

    Parameters
    ----------

    vec1, vec2 : array-like, shape=(3, x, [y, z])
        Arrays of vector values in a 1D, 2D or 3D domain.

    Returns
    -------

    scalar : ndarray, shape=(x, [y, z])
        3D grid of scalar values which are the dot products of specified
        vectors,

        .. math::

           a = \vec{v_1} \cdot \vec{v_2}
    """

    assert vec1.shape[0] == 3, "First argument provided is not a vector field"
    assert vec2.shape[0] == 3, "Second argument provided is not a vector field"
    assert vec1.shape == vec2.shape, """
        Shapes of vectors provided do not match: {}/{}
        """.format(vec1.shape, vec2.shape)

    product = np.sum(vec1 * vec2, axis=0)
    assert product.shape == vec1.shape[1:], """
        Result calculated has shape {}, should be {}
        """.format(product.shape, vec1.shape[1:])

    return product


def cross(vec1, vec2):
    r"""Calculates the cross product of two arrays of vector quantities.

    Parameters
    ----------

    vec1, vec2 : array-like, shape=(3, x, [y, z])
        Arrays of vector values in a 1D, 2D or 3D domain.

    Returns
    -------

    product : ndarray, shape=(3, x, [y, z])
        Vector field corresponding to the cross product of the specified
        vectors,

        .. math::

           \vec{a} = \vec{v_1} \times \vec{v_2}
    """

    assert vec1.shape[0] == 3, "First argument provided is not a vector field"
    assert vec2.shape[0] == 3, "Second argument provided is not a vector field"
    assert vec1.shape == vec2.shape, """
        Shapes of vectors provided do not match: {}/{}
        """.format(vec1.shape, vec2.shape)

    product = np.array((((vec1[1] * vec2[2]) - (vec1[2] * vec2[1])),
                        ((vec1[2] * vec2[0]) - (vec1[0] * vec2[2])),
                        ((vec1[0] * vec2[1]) - (vec1[1] * vec2[0]))))\
        * vec1.unit * vec2.unit
    assert product.shape == vec1.shape, """
        Result calculated has shape {}, should be {}
        """.format(product.shape, vec1.shape)

    return product


def grad(f, solver):
    r"""Calculates the gradient of a scalar field.

    Parameters
    ----------

    f : array-like, shape=(x, [y, z])
        1-, 2- or 3-dimensional scalar field.

    Returns
    -------

    gradient : astropy.units.Quantity, shape=(3, x, [y, z])
        Vector field corresponding to the gradient of the specified scalar
        field,

        .. math::

           \vec{a} = \nabla f
    """

    assert len(solver.dx) == len(f.shape), """
        Number of grid step sizes ({}) != to dimensionality of field ({})
        """.format(len(h), len(f.shape))

    gradient = np.zeros((3, *f.shape)) * f.unit / solver.dx.unit
    for dim in range(len(solver.dx)):
        gradient[dim] = solver(f, dim)

    return gradient


def div(vec, solver):
    r"""Calculates the divergence of a vector field.

    Parameters
    ----------

    vec : array-like, shape=(3, x, [y, z])
        3-dimensional vector field.

    Returns
    -------

    divergence : ndarray, shape=(x, [y, z])
        Scalar field of values corresponding to divergence of specified vector
        field,

        .. math::

           a = \nabla \cdot \vec{v}
    """

    assert vec.shape[0] == 3, "First argument provided is not a vector field"
    assert len(solver.dx) == len(vec.shape[1:]), """
        Number of grid step sizes ({}) != to dimensionality of field ({})
        """.format(len(solver.dx), len(vec.shape[1:]))

    dims = range(len(solver.dx))
    divergence = sum([solver(vec[i], i) for i in dims])
    assert divergence.shape == vec.shape[1:], """
        Output field has shape {}, should be {}
        """.format(divergence.shape, vec.shape[1:])

    return divergence


def curl(vec, solver):
    r"""Calculates the curl of a vector field.

    Parameters
    ----------

    vec : array-like, shape=(3, x, [y, z])
        3-dimensional vector field.

    Returns
    -------

    curl : ndarray, shape=(3, x, [y, z])
        Vector field corresponding to the curl of the input vector field.

        .. math::

           \vec{a} = \nabla \times \vec{v}
    """

    assert vec.shape[0] == 3, "First argument provided is not a vector field"
    assert len(solver.dx) == len(vec.shape[1:]), """
        Number of grid step sizes ({}) != to dimensionality of field ({})
        """.format(len(solver.dx), len(vec.shape[1:]))

    curl = np.zeros(vec.shape) * vec.unit / solver.dx.unit
    for k in range(3):
        for l in range(3):
            for m in range(3):
                try:
                    curl[k] += levi_civita3d(k, l, m) * solver(vec[m], l, solver.dx[l])
                except IndexError:
                    pass

    return curl


def vdp(vec1, vec2):
    r"""Calculate the Vector Direct Product of two vectors.

    Parameters
    ----------

    vec1, vec2 : array-like, shape=(3, x, [y, z])
        Arrays of vector values in a 1D, 2D or 3D domain.

    Returns
    -------

    tensor : ndarray, shape=(3, 3, x, [y, z])
        Tensor field resulting from direct product of the specified vectors.

        .. math::

           \textbf{A} = \vec{v_1} \vec{v_2}

    References
    ----------
    http://mathworld.wolfram.com/VectorDirectProduct.html

    """

    assert vec1.shape[0] == 3, "First argument provided is not a vector field"
    assert vec2.shape[0] == 3, "Second argument provided is not a vector field"
    assert vec1.shape == vec2.shape, """
        Shapes of vectors provided do not match: {}/{}
        """.format(vec1.shape, vec2.shape)

    tensor = vec1 * vec2.reshape(3, 1, *vec2.shape[1:])
    assert tensor.shape == (3, *vec1.shape), """
        Output field has shape {}, should be {}
        """.format(tensor.shape, (3, *vec1.shape))

    return tensor


def tensordiv(tensor, solver):
    r"""Calculates the divergence of a tensor field.

    Parameters
    ----------

    tensor : array-like, shape=(3, 3, x, [y, z])
        3-dimensional tensor field.

    Returns
    -------

    divergence : ndarray, shape=(3, x, [y, z])
        Vector field corresponding to the divergence of the input tensor field.

        .. math::

           \vec{a} = \nabla \cdot \textbf{T}
    """

    assert tensor.shape[:2] == (3, 3), """
        First argument provided is not a tensor field"""
    assert len(solver.dx) == len(tensor.shape[2:]), """
        Number of grid step sizes ({}) != to dimensionality of field ({})
        """.format(len(solver.dx), len(tensor.shape[2:]))

    dims = range(len(solver.dx))
    divergence = np.array([sum([solver(tensor[i, 0, ...], i) for i in dims]),
                           sum([solver(tensor[i, 1, ...], i) for i in dims]),
                           sum([solver(tensor[i, 2, ...], i) for i in dims])])\
        * tensor.unit / solver.dx.unit
    assert divergence.shape == tensor.shape[1:], """
        Output field has shape {}, should be {}
        """.format(divergence.shape, tensor.shape[1:])

    return divergence
