"""Functionality to find and analyze 3D magnetic null points."""

__all__ = [
    "Point",
    "NullPoint",
    "nullpoint_find",
    "trilinear_approx",
]

import numpy as np
import warnings

# Declare Constants & global variables
ATOL = 1e-10
MAX_DIVIDE = 10
global divide
divide = 0


class Point:
    """
    Abstract class for defining a point in 3D space.
    """

    def __init__(self, loc):
        self.loc = loc
        self.type = type

    def getLoc(self):
        r"""
        Returns the coordinates of the point object.
        """
        return self.loc


class NullPoint(Point):
    """
    A class for defining a null point in 3D space.
    """

    def __init__(self, null_loc, type):
        super().__init__(null_loc)
        self.type = type

    def getType(self):
        r"""
        Returns the type of the null point object.
        """
        return self.type

    def isEqual(self, point):
        r"""
        Returns True if two null point objects have the same coordinates.
        False otherwise.
        """
        return (
            np.isclose(self.getLoc()[0], point.getLoc()[0])
            and np.isclose(self.getLoc()[1], point.getLoc()[1])
            and np.isclose(self.getLoc()[2], point.getLoc()[2])
        )


def vector_space(
    x_arr=None,
    y_arr=None,
    z_arr=None,
    x_range=[0, 1],
    y_range=[0, 1],
    z_range=[0, 1],
    u_arr=None,
    v_arr=None,
    w_arr=None,
    func=(lambda x, y, z: [x, y, z]),
    precision=[0.05, 0.05, 0.05],
):
    r"""
    Returns a vector space in the form of a multi-dimensional array.

    Parameters
    ----------

    x_arr: array_like
        The array representing the coordinates in the x-dimension.
        If not given, then range values are used to construct a
        uniform array on that interval.

    y_arr: array_like
        The array representing the coordinates in the y-dimension.
        If not given, then range values are used to construct a
        uniform array on that interval.

    z_arr: array_like
        The array representing the coordinates in the z-dimension.
        If not given, then range values are used to construct a
        uniform array on that interval.

    x_range: array_like
        A 1 by 2 array containing the range of x-values for the vector spaces.
        If not given, the default interval [0,1] is assumed.

    y_range: array_like
        A 1 by 2 array containing the range of y-values for the vector spaces.
        If not given, the default interval [0,1] is assumed.

    z_range: array_like
        A 1 by 2 array containing the range of z-values for the vector spaces.
        If not given, the default interval [0,1] is assumed.

    u_arr: array_like
        A 3D array containing the x-component of the vector values for the vector
        space. If not given, the vector values are generated over the vector space
        using the function func.

    v_arr: array_like
        A 3D array containing the y-component of the vector values for the vector
        space. If not given, the vector values are generated over the vector space
        using the function func.

    w_arr: array_like
        A 3D array containing the z-component of the vector values for the vector
        space. If not given, the vector values are generated over the vector space
        using the function func.

    func: <class 'function'>
        A function that takes in 3 arguments, respectively representing a x, y, and z
        coordinate of a point and returns the vector value for that point in the form
        of a 1 by 3 array.

    precision: array_like
        A 1 by 3 array containing the approximate precision values for each dimension,
        in the case where uniform arrays are being used.
        The default value is [0.05, 0.05, 0.05].

    Returns
    -------
    ndarray
        A 1 by 3 array with
        the first element containing the coordinates.
        the second element containing the vector values.
        and the third element containing the delta values for each dimension.

    Raises
    ------
    This function does not raise any exceptions.

    Warns
    -----
    This function does not raise any warnings.


    Notes
    -----
    N/A
    """
    # Constructing the Meshgrid
    if (
        not isinstance(x_arr, type(None))
        and not isinstance(y_arr, type(None))
        and not isinstance(z_arr, type(None))
    ):
        x, y, z = np.meshgrid(
            x_arr,
            y_arr,
            z_arr,
            indexing="ij",
        )
        dx = np.diff(x_arr)
        dy = np.diff(y_arr)
        dz = np.diff(z_arr)
    else:
        x_den = int(np.around((x_range[1] - x_range[0]) / precision[0]) + 1)
        y_den = int(np.around((y_range[1] - y_range[0]) / precision[1]) + 1)
        z_den = int(np.around((z_range[1] - z_range[0]) / precision[2]) + 1)
        dx = np.diff(np.linspace(x_range[0], x_range[1], x_den))
        dy = np.diff(np.linspace(y_range[0], y_range[1], y_den))
        dz = np.diff(np.linspace(z_range[0], z_range[1], z_den))
        x, y, z = np.meshgrid(
            np.linspace(x_range[0], x_range[1], x_den),
            np.linspace(y_range[0], y_range[1], y_den),
            np.linspace(z_range[0], z_range[1], z_den),
            indexing="ij",
        )
    # Calculating the vector values
    if (
        not isinstance(u_arr, type(None))
        and not isinstance(v_arr, type(None))
        and not isinstance(w_arr, type(None))
    ):
        u = u_arr
        v = v_arr
        w = w_arr
    else:
        u, v, w = func(x, y, z)
    return np.array([x, y, z]), np.array([u, v, w]), np.array([dx, dy, dz])


def trilinear_coeff_cal(vspace, cell):
    r"""
    Return the coefficients for the trilinear approximation function.
    on a given grid cell in a given vector space.

    Parameters
    ----------

    vspace: array_like
        The vector space as constructed by the vector_space function which is
        A 1 by 3 array with the first element containing the coordinates,
        the second element containing the vector values,
        and the third element containing the delta values for each dimension.


    cell: array_like of integers
        A grid cell, represented by a 1 by 3 array
        of integers, which correspond to a grid cell
        in the vector space.


    Returns
    -------
    ndarray
        Returns a 1 by 3 array with
        the first element containing the coefficients for the trilinear approximation function
        for the x-component of the vector space,
        the second element containing the coefficients for the trilinear approximation function
        for the y-component of the vector space, and
        the third element containing the coefficients for the trilinear approximation function
        for the z-component of the vector space.


    Raises
    ------
    This function does not raise any exceptions.

    Warns
    -----
    This function does not raise any warnings.

    Notes
    -----
    N/A
    """
    u, v, w = vspace[1]
    deltax, deltay, deltaz = vspace[2]
    f000 = cell
    f001 = [cell[0], cell[1], cell[2] + 1]
    f010 = [cell[0], cell[1] + 1, cell[2]]
    f011 = [cell[0], cell[1] + 1, cell[2] + 1]
    f100 = [cell[0] + 1, cell[1], cell[2]]
    f101 = [cell[0] + 1, cell[1], cell[2] + 1]
    f110 = [cell[0] + 1, cell[1] + 1, cell[2]]
    f111 = [cell[0] + 1, cell[1] + 1, cell[2] + 1]
    x0 = float(vspace[0][0][f000[0]][f000[1]][f000[2]])
    x1 = float(x0 + deltax[cell[0]])
    y0 = float(vspace[0][1][f000[0]][f000[1]][f000[2]])
    y1 = float(y0 + deltay[cell[1]])
    z0 = float(vspace[0][2][f000[0]][f000[1]][f000[2]])
    z1 = float(z0 + deltaz[cell[2]])
    A = np.array(
        [
            [1, x0, y0, z0, x0 * y0, x0 * z0, y0 * z0, x0 * y0 * z0],
            [1, x1, y0, z0, x1 * y0, x1 * z0, y0 * z0, x1 * y0 * z0],
            [1, x0, y1, z0, x0 * y1, x0 * z0, y1 * z0, x0 * y1 * z0],
            [1, x1, y1, z0, x1 * y1, x1 * z0, y1 * z0, x1 * y1 * z0],
            [1, x0, y0, z1, x0 * y0, x0 * z1, y0 * z1, x0 * y0 * z1],
            [1, x1, y0, z1, x1 * y0, x1 * z1, y0 * z1, x1 * y0 * z1],
            [1, x0, y1, z1, x0 * y1, x0 * z1, y1 * z1, x0 * y1 * z1],
            [1, x1, y1, z1, x1 * y1, x1 * z1, y1 * z1, x1 * y1 * z1],
        ]
    )
    sx = np.array(
        [
            [u[f000[0]][f000[1]][f000[2]]],
            [u[f100[0]][f100[1]][f100[2]]],
            [u[f010[0]][f010[1]][f010[2]]],
            [u[f110[0]][f110[1]][f110[2]]],
            [u[f001[0]][f001[1]][f001[2]]],
            [u[f101[0]][f101[1]][f101[2]]],
            [u[f011[0]][f011[1]][f011[2]]],
            [u[f111[0]][f111[1]][f111[2]]],
        ]
    )
    sy = np.array(
        [
            [v[f000[0]][f000[1]][f000[2]]],
            [v[f100[0]][f100[1]][f100[2]]],
            [v[f010[0]][f010[1]][f010[2]]],
            [v[f110[0]][f110[1]][f110[2]]],
            [v[f001[0]][f001[1]][f001[2]]],
            [v[f101[0]][f101[1]][f101[2]]],
            [v[f011[0]][f011[1]][f011[2]]],
            [v[f111[0]][f111[1]][f111[2]]],
        ]
    )
    sz = np.array(
        [
            [w[f000[0]][f000[1]][f000[2]]],
            [w[f100[0]][f100[1]][f100[2]]],
            [w[f010[0]][f010[1]][f010[2]]],
            [w[f110[0]][f110[1]][f110[2]]],
            [w[f001[0]][f001[1]][f001[2]]],
            [w[f101[0]][f101[1]][f101[2]]],
            [w[f011[0]][f011[1]][f011[2]]],
            [w[f111[0]][f111[1]][f111[2]]],
        ]
    )

    ax, bx, cx, dx, ex, fx, gx, hx = np.linalg.solve(A, sx).reshape(1, 8)[0]
    ay, by, cy, dy, ey, fy, gy, hy = np.linalg.solve(A, sy).reshape(1, 8)[0]
    az, bz, cz, dz, ez, fz, gz, hz = np.linalg.solve(A, sz).reshape(1, 8)[0]

    return np.array(
        [
            [ax, bx, cx, dx, ex, fx, gx, hx],
            [ay, by, cy, dy, ey, fy, gy, hy],
            [az, bz, cz, dz, ez, fz, gz, hz],
        ]
    )


def trilinear_approx(vspace, cell):
    r"""
    Return a function whose input is a coordinate within a given grid cell
    and returns the trilinearly approximated vector value at that particular
    coordinate in that grid cell.

    Parameters
    ----------

    vspace: array_like
        The vector space as constructed by the vector_space function which is
        A 1 by 3 array with the first element containing the coordinates,
        the second element containing the vector values,
        and the third element containing the delta values for each dimension.

    cell: array_like of integers
        A grid cell, represented by a 1 by 3 array
        of integers, which correspond to a grid cell
        in the vector space.

    Returns
    -------
    <class 'function'>
        A function whose input is a coordinate within a given grid cell
        and returns the trilinearly approximated vector value at that particular
        coordinate in that grid cell.

    Raises
    ------
    This function does not raise any exceptions.

    Warns
    -----
    This function does not raise any warnings.


    Notes
    -----
    N/A
    """
    # Calculating coefficients
    ax, bx, cx, dx, ex, fx, gx, hx = trilinear_coeff_cal(vspace, cell)[0]
    ay, by, cy, dy, ey, fy, gy, hy = trilinear_coeff_cal(vspace, cell)[1]
    az, bz, cz, dz, ez, fz, gz, hz = trilinear_coeff_cal(vspace, cell)[2]

    def approx_func(xInput, yInput, zInput):
        Bx = (
            ax
            + bx * xInput
            + cx * yInput
            + dx * zInput
            + ex * xInput * yInput
            + fx * xInput * zInput
            + gx * yInput * zInput
            + hx * xInput * yInput * zInput
        )
        By = (
            ay
            + by * xInput
            + cy * yInput
            + dy * zInput
            + ey * xInput * yInput
            + fy * xInput * zInput
            + gy * yInput * zInput
            + hy * xInput * yInput * zInput
        )
        Bz = (
            az
            + bz * xInput
            + cz * yInput
            + dz * zInput
            + ez * xInput * yInput
            + fz * xInput * zInput
            + gz * yInput * zInput
            + hz * xInput * yInput * zInput
        )
        return np.array([Bx, By, Bz])

    return approx_func


def jacobian(vspace, cell):
    r"""
    Returns a function whose input is a coordinate within a given grid cell
    and returns the trilinearly approximated jacobian matrix for that particular
    coordinate in that grid cell.

    Parameters
    ----------

    vspace: array_like
        The vector space as constructed by the vector_space function which is
        A 1 by 3 array with the first element containing the coordinates,
        the second element containing the vector values,
        and the third element containing the delta values for each dimension.

    cell: array_like of integers
        A grid cell, represented by a 1 by 3 array
        of integers, which correspond to a grid cell
        in the vector space.


    Returns
    -------
    <class 'function'>
        A function whose input is a coordinate within a given grid cell
        and returns the trilinearly approximated jacobian matrix for that particular
        coordinate in that grid cell.


    Raises
    ------
    This function does not raise any exceptions.

    Warns
    -----
    This function does not raise any warnings.

    Notes
    -----
    N/A
    """
    # Calculating coefficients
    ax, bx, cx, dx, ex, fx, gx, hx = trilinear_coeff_cal(vspace, cell)[0]
    ay, by, cy, dy, ey, fy, gy, hy = trilinear_coeff_cal(vspace, cell)[1]
    az, bz, cz, dz, ez, fz, gz, hz = trilinear_coeff_cal(vspace, cell)[2]

    def jacobian_func(xInput, yInput, zInput):
        dBxdx = bx + ex * yInput + fx * zInput + hx * yInput * zInput
        dBxdy = cx + ex * xInput + gx * zInput + hx * xInput * yInput
        dBxdz = dx + fx * xInput + gx * yInput + hx * xInput * yInput
        dBydx = by + ey * yInput + fy * zInput + hy * yInput * zInput
        dBydy = cy + ey * xInput + gy * zInput + hy * xInput * yInput
        dBydz = dy + fy * xInput + gy * yInput + hy * xInput * yInput
        dBzdx = bz + ez * yInput + fz * zInput + hz * yInput * zInput
        dBzdy = cz + ez * xInput + gz * zInput + hz * xInput * yInput
        dBzdz = dz + fz * xInput + gz * yInput + hz * xInput * yInput
        jmatrix = np.array(
            [
                float(dBxdx),
                float(dBxdy),
                float(dBxdz),
                float(dBydx),
                float(dBydy),
                float(dBydz),
                float(dBzdx),
                float(dBzdy),
                float(dBzdz),
            ]
        ).reshape(3, 3)
        return jmatrix

    return jacobian_func


def reduction(vspace, cell):
    r"""
    Return a true or false based on weather
    a grid cell passes the reduction phase,
    meaning that they potentionally contain a null point.

    Parameters
    ----------

    vspace: array_like
        The vector space as constructed by the vector_space function which is
        A 1 by 3 array with the first element containing the coordinates,
        the second element containing the vector values,
        and the third element containing the delta values for each dimension.

    cell: array_like of integers
        A grid cell, represented by a 1 by 3 array
        of integers, which correspond to a grid cell
        in the vector space.

    Returns
    -------
    bool
        True if a grid cell passes the reduction phase.
        False, otherwise.

    Raises
    ------
    This function does not raise any exceptions.

    Warns
    -----
    This function does not raise any warnings.

    Notes
    -----
    Depending on the grid resolution, a cell containing more than
    one null point may not pass reduction, so it would not be detected.
    """
    u, v, w = vspace[1]
    f000 = cell
    f001 = [cell[0], cell[1], cell[2] + 1]
    f010 = [cell[0], cell[1] + 1, cell[2]]
    f011 = [cell[0], cell[1] + 1, cell[2] + 1]
    f100 = [cell[0] + 1, cell[1], cell[2]]
    f101 = [cell[0] + 1, cell[1], cell[2] + 1]
    f110 = [cell[0] + 1, cell[1] + 1, cell[2]]
    f111 = [cell[0] + 1, cell[1] + 1, cell[2] + 1]
    corners = [f000, f001, f010, f011, f100, f101, f110, f111]
    passX = False
    passY = False
    passZ = False
    # Check X-Component
    sign = np.sign(u[cell[0], cell[1], cell[2]])
    for point in corners:
        if (
            u[point[0]][point[1]][point[2]] == 0
            or np.sign(u[point[0]][point[1]][point[2]]) != sign
        ):
            passX = True

    # Check Y-Component
    sign = np.sign(v[cell[0], cell[1], cell[2]])
    for point in corners:
        if (
            v[point[0]][point[1]][point[2]] == 0
            or np.sign(v[point[0]][point[1]][point[2]]) != sign
        ):
            passY = True

    # Check Z-Component
    sign = np.sign(w[cell[0], cell[1], cell[2]])
    for point in corners:
        if (
            w[point[0]][point[1]][point[2]] == 0
            or np.sign(w[point[0]][point[1]][point[2]]) != sign
        ):
            passZ = True

    doesPassReduction = passX and passY and passZ

    return doesPassReduction


def bilinear_root(a1, b1, c1, d1, a2, b2, c2, d2):
    r"""
    Return the roots of a pair of bilinear equations of the following format.
    a1+b1x+c1y+d1xy=0
    a2+b2x+c2y+d2xy=0

    Parameters
    ----------
    a1: float
    b1: float
    c1: float
    d1: float
    a2: float
    b2: float
    c2: float
    d2: float

    Returns
    -------
    roots : array_like of floats
        A 1 by 2 array containing the two roots

    Raises
    ------
    This function does not raise any exceptions.

    Warns
    -----
    This function does not raise any warnings.

    Notes
    -----
    N/A
    """
    m1 = np.array([[a1, a2], [c1, c2]])
    m2 = np.array([[a1, a2], [d1, d2]])
    m3 = np.array([[b1, b2], [c1, c2]])
    m4 = np.array([[b1, b2], [d1, d2]])

    a = np.linalg.det(m4)
    b = np.linalg.det(m2) + np.linalg.det(m3)
    c = np.linalg.det(m1)

    if np.isclose(a, 0, atol=ATOL):
        if np.isclose(b, 0, atol=ATOL):
            return None, None
        else:
            x1 = (-1.0 * c) / b
            x2 = (-1.0 * c) / b
    else:
        if (b ** 2 - 4.0 * a * c) < 0:
            return None, None
        else:
            x1 = (-1.0 * b + (b ** 2 - 4.0 * a * c) ** 0.5) / (2.0 * a)
            x2 = (-1.0 * b - (b ** 2 - 4.0 * a * c) ** 0.5) / (2.0 * a)

    m1 = np.array([[a1, a2], [b1, b2]])
    m2 = np.array([[a1, a2], [d1, d2]])
    m3 = np.array([[b1, b2], [c1, c2]])
    m4 = np.array([[c1, c2], [d1, d2]])
    a = np.linalg.det(m4)
    b = np.linalg.det(m2) - np.linalg.det(m3)
    c = np.linalg.det(m1)

    if np.isclose(a, 0, atol=ATOL):
        if np.isclose(b, 0, atol=ATOL):
            return None, None
        else:
            y1 = (-1.0 * c) / b
            y2 = (-1.0 * c) / b
    else:
        if (b ** 2 - 4.0 * a * c) < 0:
            return None, None
        else:
            y1 = (-1.0 * b - (b ** 2 - 4.0 * a * c) ** 0.5) / (2.0 * a)
            y2 = (-1.0 * b + (b ** 2 - 4.0 * a * c) ** 0.5) / (2.0 * a)

    return [(x1, y1), (x2, y2)]


def trillinear_analysis(vspace, cell):
    r"""
    Return a true or false value based on weather
    a grid cell which has passed the reduction step,
    contains a null point, using trilinear analysis.

    Parameters
    ----------

    vspace: array_like
        The vector space as constructed by the vector_space function which is
        A 1 by 3 array with the first element containing the coordinates,
        the second element containing the vector values,
        and the third element containing the delta values for each dimension.

    cell: array_like of integers
        A grid cell, represented by a 1 by 3 array
        of integers, which correspond to a grid cell
        in the vector space.

    Returns
    -------
    bool
        True if a grid cell contains a nullpoint using trilinear analysis.
        False, otherwise.

    Raises
    ------
    This function does not raise any exceptions.

    Warns
    -----
    :'UserWarning'
        If there is a possible lack of grid resolution, so
        that a grid cell may contain more than one nullpoint.

    Notes
    -----
    N/A
    """

    # Helper Function
    def is_close(a, b):
        arr = np.isclose(a, b, atol=ATOL)
        if type(arr) == np.bool_:
            return arr
        res = True
        for b in arr:
            res = res and b
        return res

    # Critical Cell Corners
    f000 = cell
    f111 = [cell[0] + 1, cell[1] + 1, cell[2] + 1]

    # Calculating coefficients
    ax, bx, cx, dx, ex, fx, gx, hx = trilinear_coeff_cal(vspace, cell)[0]
    ay, by, cy, dy, ey, fy, gy, hy = trilinear_coeff_cal(vspace, cell)[1]
    az, bz, cz, dz, ez, fz, gz, hz = trilinear_coeff_cal(vspace, cell)[2]

    # Initial Position of the cell corner
    initial = np.array(
        [
            vspace[0][0][f000[0]][f000[1]][f000[2]],
            vspace[0][1][f000[0]][f000[1]][f000[2]],
            vspace[0][2][f000[0]][f000[1]][f000[2]],
        ]
    )

    # Arrays holding the endpoints
    BxByEndpoints = []
    BxBzEndpoints = []
    ByBzEndpoints = []

    # Front Surface
    yConst1 = vspace[0][1][f000[0]][f000[1]][
        f000[2]
    ]  # y-coordinate of the front surface
    # Bx=By=0 Curve Endpoint

    root1, root2 = bilinear_root(
        ax + cx * yConst1,
        bx + ex * yConst1,
        dx + gx * yConst1,
        fx + hx * yConst1,
        ay + cy * yConst1,
        by + ey * yConst1,
        dy + gy * yConst1,
        fy + hy * yConst1,
    )

    if root1 is not None:
        if is_close(root1, root2):
            BxByEndpoints.append((root1[0], yConst1, root1[1]))
        else:
            BxByEndpoints.append((root1[0], yConst1, root1[1]))
            BxByEndpoints.append((root2[0], yConst1, root2[1]))

    # Bx=BZ=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + cx * yConst1,
        bx + ex * yConst1,
        dx + gx * yConst1,
        fx + hx * yConst1,
        az + cz * yConst1,
        bz + ez * yConst1,
        dz + gz * yConst1,
        fz + hz * yConst1,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxBzEndpoints.append((root1[0], yConst1, root1[1]))
        else:
            BxBzEndpoints.append((root1[0], yConst1, root1[1]))
            BxBzEndpoints.append((root2[0], yConst1, root2[1]))
    # By=Bz=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ay + cy * yConst1,
        by + ey * yConst1,
        dy + gy * yConst1,
        fy + hy * yConst1,
        az + cz * yConst1,
        bz + ez * yConst1,
        dz + gz * yConst1,
        fz + hz * yConst1,
    )
    if root1 is not None:
        if is_close(root1, root2):
            ByBzEndpoints.append((root1[0], yConst1, root1[1]))
        else:
            ByBzEndpoints.append((root1[0], yConst1, root1[1]))
            ByBzEndpoints.append((root2[0], yConst1, root2[1]))

    # Back Surface
    yConst2 = vspace[0][1][f111[0]][f111[1]][
        f111[2]
    ]  # y-coordinate of the front surface
    # Bx=By=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + cx * yConst2,
        bx + ex * yConst2,
        dx + gx * yConst2,
        fx + hx * yConst2,
        ay + cy * yConst2,
        by + ey * yConst2,
        dy + gy * yConst2,
        fy + hy * yConst2,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxByEndpoints.append((root1[0], yConst2, root1[1]))
        else:
            BxByEndpoints.append((root1[0], yConst2, root1[1]))
            BxByEndpoints.append((root2[0], yConst2, root2[1]))
    # Bx=BZ=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + cx * yConst2,
        bx + ex * yConst2,
        dx + gx * yConst2,
        fx + hx * yConst2,
        az + cz * yConst2,
        bz + ez * yConst2,
        dz + gz * yConst2,
        fz + hz * yConst2,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxBzEndpoints.append((root1[0], yConst2, root1[1]))
        else:
            BxBzEndpoints.append((root1[0], yConst2, root1[1]))
            BxBzEndpoints.append((root2[0], yConst2, root2[1]))
    # By=Bz=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ay + cy * yConst2,
        by + ey * yConst2,
        dy + gy * yConst2,
        fy + hy * yConst2,
        az + cz * yConst2,
        bz + ez * yConst2,
        dz + gz * yConst2,
        fz + hz * yConst2,
    )
    if root1 is not None:
        if is_close(root1, root2):
            ByBzEndpoints.append((root1[0], yConst2, root1[1]))
        else:
            ByBzEndpoints.append((root1[0], yConst2, root1[1]))
            ByBzEndpoints.append((root2[0], yConst2, root2[1]))

    # Right Surface
    xConst1 = vspace[0][0][f111[0]][f111[1]][f111[2]]
    # Bx=By=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + bx * xConst1,
        cx + ex * xConst1,
        dx + fx * xConst1,
        gx + hx * xConst1,
        ay + by * xConst1,
        cy + ey * xConst1,
        dy + fy * xConst1,
        gy + hy * xConst1,
    )

    if root1 is not None:
        if is_close(root1, root2):
            BxByEndpoints.append((xConst1, root1[0], root1[1]))
        else:
            BxByEndpoints.append((xConst1, root1[0], root1[1]))
            BxByEndpoints.append((xConst1, root2[0], root2[1]))

    # Bx=BZ=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + bx * xConst1,
        cx + ex * xConst1,
        dx + fx * xConst1,
        gx + hx * xConst1,
        az + bz * xConst1,
        cz + ez * xConst1,
        dz + fz * xConst1,
        gz + hz * xConst1,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxBzEndpoints.append((xConst1, root1[0], root1[1]))
        else:
            BxBzEndpoints.append((xConst1, root1[0], root1[1]))
            BxBzEndpoints.append((xConst1, root2[0], root2[1]))

    # By=Bz=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ay + by * xConst1,
        cy + ey * xConst1,
        dy + fy * xConst1,
        gy + hy * xConst1,
        az + bz * xConst1,
        cz + ez * xConst1,
        dz + fz * xConst1,
        gz + hz * xConst1,
    )
    if root1 is not None:
        if is_close(root1, root2):
            ByBzEndpoints.append((xConst1, root1[0], root1[1]))
        else:
            ByBzEndpoints.append((xConst1, root1[0], root1[1]))
            ByBzEndpoints.append((xConst1, root2[0], root2[1]))

    # Left Surface
    xConst2 = vspace[0][0][f000[0]][f000[1]][f000[2]]
    # Bx=By=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + bx * xConst2,
        cx + ex * xConst2,
        dx + fx * xConst2,
        gx + hx * xConst2,
        ay + by * xConst2,
        cy + ey * xConst2,
        dy + fy * xConst2,
        gy + hy * xConst2,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxByEndpoints.append((xConst2, root1[0], root1[1]))
        else:
            BxByEndpoints.append((xConst2, root1[0], root1[1]))
            BxByEndpoints.append((xConst2, root2[0], root2[1]))
    # Bx=BZ=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + bx * xConst2,
        cx + ex * xConst2,
        dx + fx * xConst2,
        gx + hx * xConst2,
        az + bz * xConst2,
        cz + ez * xConst2,
        dz + fz * xConst2,
        gz + hz * xConst2,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxBzEndpoints.append((xConst2, root1[0], root1[1]))
        else:
            BxBzEndpoints.append((xConst2, root1[0], root1[1]))
            BxBzEndpoints.append((xConst2, root2[0], root2[1]))

    # By=Bz=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ay + by * xConst2,
        cy + ey * xConst2,
        dy + fy * xConst2,
        gy + hy * xConst2,
        az + bz * xConst2,
        cz + ez * xConst2,
        dz + fz * xConst2,
        gz + hz * xConst2,
    )
    if root1 is not None:
        if is_close(root1, root2):
            ByBzEndpoints.append((xConst2, root1[0], root1[1]))
        else:
            ByBzEndpoints.append((xConst2, root1[0], root1[1]))
            ByBzEndpoints.append((xConst2, root2[0], root2[1]))

    # Up Surface
    zConst1 = vspace[0][2][f111[0]][f111[1]][f111[2]]
    # Bx=By=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + dx * zConst1,
        bx + fx * zConst1,
        cx + gx * zConst1,
        ex + hx * zConst1,
        ay + dy * zConst1,
        by + fy * zConst1,
        cy + gy * zConst1,
        ey + hy * zConst1,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxByEndpoints.append((root1[0], root1[1], zConst1))
        else:
            BxByEndpoints.append((root1[0], root1[1], zConst1))
            BxByEndpoints.append((root2[0], root2[1], zConst1))
    # Bx=BZ=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + dx * zConst1,
        bx + fx * zConst1,
        cx + gx * zConst1,
        ex + hx * zConst1,
        az + dz * zConst1,
        bz + fz * zConst1,
        cz + gz * zConst1,
        ez + hz * zConst1,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxBzEndpoints.append((root1[0], root1[1], zConst1))
        else:
            BxBzEndpoints.append((root1[0], root1[1], zConst1))
            BxBzEndpoints.append((root2[0], root2[1], zConst1))
    # By=Bz=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ay + dy * zConst1,
        by + fy * zConst1,
        cy + gy * zConst1,
        ey + hy * zConst1,
        az + dz * zConst1,
        bz + fz * zConst1,
        cz + gz * zConst1,
        ez + hz * zConst1,
    )
    if root1 is not None:
        if is_close(root1, root2):
            ByBzEndpoints.append((root1[0], root1[1], zConst1))
        else:
            ByBzEndpoints.append((root1[0], root1[1], zConst1))
            ByBzEndpoints.append((root2[0], root2[1], zConst1))

    # Down Surface
    zConst2 = vspace[0][2][f000[0]][f000[1]][f000[2]]
    # Bx=By=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + dx * zConst2,
        bx + fx * zConst2,
        cx + gx * zConst2,
        ex + hx * zConst2,
        ay + dy * zConst2,
        by + fy * zConst2,
        cy + gy * zConst2,
        ey + hy * zConst2,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxByEndpoints.append((root1[0], root1[1], zConst2))
        else:
            BxByEndpoints.append((root1[0], root1[1], zConst2))
            BxByEndpoints.append((root2[0], root2[1], zConst2))
    # Bx=BZ=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ax + dx * zConst2,
        bx + fx * zConst2,
        cx + gx * zConst2,
        ex + hx * zConst2,
        az + dz * zConst2,
        bz + fz * zConst2,
        cz + gz * zConst2,
        ez + hz * zConst2,
    )
    if root1 is not None:
        if is_close(root1, root2):
            BxBzEndpoints.append((root1[0], root1[1], zConst2))
        else:
            BxBzEndpoints.append((root1[0], root1[1], zConst2))
            BxBzEndpoints.append((root2[0], root2[1], zConst2))
    # By=Bz=0 Curve Endpoint
    root1, root2 = bilinear_root(
        ay + dy * zConst2,
        by + fy * zConst2,
        cy + gy * zConst2,
        ey + hy * zConst2,
        az + dz * zConst2,
        bz + fz * zConst2,
        cz + gz * zConst2,
        ez + hz * zConst2,
    )
    if root1 is not None:
        if is_close(root1, root2):
            ByBzEndpoints.append((root1[0], root1[1], zConst2))
        else:
            ByBzEndpoints.append((root1[0], root1[1], zConst2))
            ByBzEndpoints.append((root2[0], root2[1], zConst2))

    xbound = vspace[0][0][f111[0]][f111[1]][f111[2]]
    ybound = vspace[0][1][f111[0]][f111[1]][f111[2]]
    zbound = vspace[0][2][f111[0]][f111[1]][f111[2]]

    def bound(epoint):
        a = (initial[0] < epoint[0] or is_close(initial[0], epoint[0])) and (
            epoint[0] < xbound or is_close(epoint[0], xbound)
        )
        b = (initial[1] < epoint[1] or is_close(initial[1], epoint[1])) and (
            epoint[1] < ybound or is_close(epoint[1], ybound)
        )
        c = (initial[2] < epoint[2] or is_close(initial[2], epoint[2])) and (
            epoint[2] < zbound or is_close(epoint[2], zbound)
        )
        return a and b and c

    BxByEndpoints = list(filter(bound, BxByEndpoints))
    BxBzEndpoints = list(filter(bound, BxBzEndpoints))
    ByBzEndpoints = list(filter(bound, ByBzEndpoints))

    tlApprox = trilinear_approx(vspace, cell)

    # Check on the Surfaces
    for p in BxByEndpoints:
        if np.linalg.norm(tlApprox(p[0], p[1], p[2])) < ATOL:
            return True
    for p in BxBzEndpoints:
        if np.linalg.norm(tlApprox(p[0], p[1], p[2])) < ATOL:
            return True
    for p in ByBzEndpoints:
        if np.linalg.norm(tlApprox(p[0], p[1], p[2])) < ATOL:
            return True

    # Check Grid Resolution
    if len(BxByEndpoints) == 0 and len(BxBzEndpoints) == 0 and len(ByBzEndpoints) == 0:
        warnings.warn("Possible Lack of Grid Resolution")
        return False

    if len(BxByEndpoints) != 2 or len(BxBzEndpoints) != 2 or len(ByBzEndpoints) != 2:
        return False

    isNullPoint = True

    # Checking sign of third component for Bx=By=0
    if (
        np.sign(
            tlApprox(BxByEndpoints[0][0], BxByEndpoints[0][1], BxByEndpoints[0][2])[2]
        )
        * np.sign(
            tlApprox(BxByEndpoints[1][0], BxByEndpoints[1][1], BxByEndpoints[1][2])[2]
        )
        > 0
    ):
        isNullPoint = False

    # Checking sign of third component for Bx=Bz=0
    if (
        np.sign(
            tlApprox(BxBzEndpoints[0][0], BxBzEndpoints[0][1], BxBzEndpoints[0][2])[1]
        )
        * np.sign(
            tlApprox(BxBzEndpoints[1][0], BxBzEndpoints[1][1], BxBzEndpoints[1][2])[1]
        )
        > 0
    ):
        isNullPoint = False

    # Checking sign of third component for By=Bz=0
    if (
        np.sign(
            tlApprox(ByBzEndpoints[0][0], ByBzEndpoints[0][1], ByBzEndpoints[0][2])[0]
        )
        * np.sign(
            tlApprox(ByBzEndpoints[1][0], ByBzEndpoints[1][1], ByBzEndpoints[1][2])[0]
        )
        > 0
    ):
        isNullPoint = False

    return isNullPoint


def locate_null_point(vspace, cell, n, err):
    r"""
    Return the coordinates of a nullpoint within
    a given grid cell in a vector space using the
    Newton-Rapshon method.
    Multiple initial position are tried until either
    one converges inside a the grid cell, or the maximum
    iteration is reached.
    If neither occurs, more starting positions are tried,
    by breaking up the cell into 8 smaller sub-grid cells,
    until one starting position does converge or stop inside
    the grid cell.
    This process is repeated a finite amount of times, after which
    the function returns None.


    Parameters
    ----------

    vspace: array_like
        The vector space as constructed by the vector_space function which is
        A 1 by 3 array with the first element containing the coordinates,
        the second element containing the vector values,
        and the third element containing the delta values for each dimension.

    cell: array_like of integers
        A grid cell, represented by a 1 by 3 array
        of integers, which correspond to a grid cell
        in the vector space.

    n: int
    The maximum number of times the iterative step
    of the Newton-Raphson method is repeated.

    err: float
    The threshold/error that determines if convergence has occured
    using the Newton-Raphson method.

    Returns
    -------
    array_like of floats
        A 1 by 3 array containing the converged coordinates of the
        null point.
    NoneType
        None if the coordinates of the null point could not be converged
        at a point inside the grid cell.

    Raises
    ------
    This function does not raise any exceptions.

    Warns
    -----
    :'UserWarning'
        If the maximum number of iteration has been
        reached, but convergence has not occurred.

    Notes
    -----
    N/A
    """

    global divide
    # Calculating the Jacobian and trillinear approximation functions for the cell
    tlApprox = trilinear_approx(vspace, cell)
    jcb = jacobian(vspace, cell)
    # Calculatiung the deltas
    deltax, deltay, deltaz = vspace[2]
    deltax = deltax[cell[0]]
    deltay = deltay[cell[1]]
    deltaz = deltaz[cell[2]]

    f000 = cell
    f001 = [cell[0], cell[1], cell[2] + 1]
    f010 = [cell[0], cell[1] + 1, cell[2]]
    f011 = [cell[0], cell[1] + 1, cell[2] + 1]
    f100 = [cell[0] + 1, cell[1], cell[2]]
    f101 = [cell[0] + 1, cell[1], cell[2] + 1]
    f110 = [cell[0] + 1, cell[1] + 1, cell[2]]
    f111 = [cell[0] + 1, cell[1] + 1, cell[2] + 1]
    corners = [f000, f001, f010, f011, f100, f101, f110, f111]

    # Critical Coordinates
    pos_000 = np.array(
        [
            vspace[0][0][f000[0]][f000[1]][f000[2]],
            vspace[0][1][f000[0]][f000[1]][f000[2]],
            vspace[0][2][f000[0]][f000[1]][f000[2]],
        ]
    )
    pos_111 = np.array(
        [
            vspace[0][0][f111[0]][f111[1]][f111[2]],
            vspace[0][1][f111[0]][f111[1]][f111[2]],
            vspace[0][2][f111[0]][f111[1]][f111[2]],
        ]
    )

    def inbound(pos):
        pos = pos.reshape(1, 3)[0]
        A = (np.isclose(pos_000[0], pos[0], atol=ATOL) or pos_000[0] < pos[0]) and (
            np.isclose(pos[0], pos_111[0], atol=ATOL) or pos[0] < pos_111[0]
        )
        B = (np.isclose(pos_000[1], pos[1], atol=ATOL) or pos_000[1] < pos[1]) and (
            np.isclose(pos[1], pos_111[1], atol=ATOL) or pos[1] < pos_111[1]
        )
        C = (np.isclose(pos_000[2], pos[2], atol=ATOL) or pos_000[2] < pos[2]) and (
            np.isclose(pos[2], pos_111[2], atol=ATOL) or pos[2] < pos_111[2]
        )
        return A and B and C

    starting_pos = []
    # Adding the Corners
    for point in corners:
        starting_pos.append(
            [
                vspace[0][0][point[0]][point[1]][point[2]],
                vspace[0][1][point[0]][point[1]][point[2]],
                vspace[0][2][point[0]][point[1]][point[2]],
            ]
        )
    # Adding the Mid Point
    starting_pos.append(
        [
            vspace[0][0][f000[0]][f000[1]][f000[2]] + deltax / 2.0,
            vspace[0][1][f000[0]][f000[1]][f000[2]] + deltay / 2.0,
            vspace[0][2][f000[0]][f000[1]][f000[2]] + deltaz / 2.0,
        ]
    )
    # Newton Iteration
    for x0 in starting_pos:
        x0 = np.array(x0)
        x0 = x0.reshape(3, 1)
        for i in range(n):
            locx = tlApprox(x0[0], x0[1], x0[2])[0]
            locy = tlApprox(x0[0], x0[1], x0[2])[1]
            locz = tlApprox(x0[0], x0[1], x0[2])[2]
            Bx0 = np.array([locx, locy, locz])
            Bx0 = Bx0.reshape(3, 1)
            prev_norm = np.linalg.norm(x0)
            x0 = np.subtract(
                x0, np.matmul(np.linalg.inv(jcb(x0[0], x0[1], x0[2])), Bx0)
            )
            norm = np.linalg.norm(x0)
            if np.abs(norm - prev_norm) < err and inbound(x0):
                return x0
        if inbound(x0):
            warnings.warn("Max Iterations Reached")
            return x0

    # Break Up the Cell into 8 smaller cells and try again
    divide = divide + 1
    if divide > MAX_DIVIDE:
        warnings.warn("Could Not Locate a possible Nullpoint")
        return None
    null_point_args = {
        "func": tlApprox,
        "x_range": [pos_000[0], pos_111[0]],
        "y_range": [pos_000[1], pos_111[1]],
        "z_range": [pos_000[2], pos_111[2]],
        "precision": [deltax / 2, deltay / 2, deltaz / 2],
    }

    return nullpoint_find(**null_point_args)


def nullpoint_find(
    x_arr=None,
    y_arr=None,
    z_arr=None,
    x_range=[0, 1],
    y_range=[0, 1],
    z_range=[0, 1],
    u_arr=None,
    v_arr=None,
    w_arr=None,
    func=(lambda x, y, z: [x, y, z]),
    precision=[0.05, 0.05, 0.05],
    MAX_ITERATIONS=500,
    err=10 ** (-10),
):
    r"""
    Returns an array of nullpoint object, representing
    the nullpoints of the given vector space.


    Parameters
    ----------
        x_arr: array_like
            The array representing the coordinates in the x-dimension.
            If not given, then range values are used to construct a
            uniform array on that interval.

        y_arr: array_like
            The array representing the coordinates in the y-dimension.
            If not given, then range values are used to construct a
            uniform array on that interval.

        z_arr: array_like
            The array representing the coordinates in the z-dimension.
            If not given, then range values are used to construct a
            uniform array on that interval.

        x_range: array_like
            A 1 by 2 array containing the range of x-vlaues for the vector spaces.
            If not given, the default interval [0,1] is assumed.

        y_range: array_like
            A 1 by 2 array containing the range of y-vlaues for the vector spaces.
            If not given, the default interval [0,1] is assumed.

        z_range: array_like
            A 1 by 2 array containing the range of z-vlaues for the vector spaces.
            If not given, the default interval [0,1] is assumed.

        u_arr: array_like
            A 3D array containing the x-component of the vector values for the vector
            space. If not given, the vector values are generated over the vector space
            using the function func.

        v_arr: array_like
            A 3D array containing the y-component of the vector values for the vector
            space. If not given, the vector values are generated over the vector space
            using the function func.

        w_arr: array_like
            A 3D array containing the z-component of the vector values for the vector
            space. If not given, the vector values are generated over the vector space
            using the function func.

        func: <class 'function'>
            A function that takes in 3 arguments, respectively representing a x, y, and z
            coordinate of a point and returns the vector value for that point in the form
            of a 1 by 3 array.

        precision: array_like
            A 1 by 3 array containing the approximate precision values for each dimension,
            in the case where uniform arrays are being used.
            The default value is [0.05, 0.05, 0.05].

        MAX_ITERATIONS: int
        The maximum iterations of the Newton-Raphson method.
        The default value is 500.

        err: float
        The threshold/error that determines if convergence has occured
        using the Newton-Raphson method.
        The default value is 10**(-10).


    Returns
    -------
        array_like of `~plasmapy.analysis.nullpoint.NullPoint`
            An array of NullPoint objects representing the nullpoints
            of the given vector space.

    Notes
    -----
    This method is described by :cite:t:`haynes:2007`
    """
    # Constructing the vspace
    vspace = vector_space(
        x_arr,
        y_arr,
        z_arr,
        x_range,
        y_range,
        z_range,
        u_arr,
        v_arr,
        w_arr,
        func,
        precision,
    )

    # Helper Function
    def in_null_list(elem, lst):
        for p in lst:
            if p.isEqual(elem):
                return True
        return False

    nullpoints = []
    for i in range(len(vspace[0][0]) - 1):
        for j in range(len(vspace[0][0][0]) - 1):
            for k in range(len(vspace[0][0][0][0]) - 1):
                if reduction(vspace, [i, j, k]):
                    if trilinear_analysis(vspace, [i, j, k]):
                        loc = locate_null_point(vspace, [i, j, k], MAX_ITERATIONS, err)
                        if not isinstance(loc, type(None)):
                            p = NullPoint(loc, "N/A")
                            if not in_null_list(p, nullpoints):
                                nullpoints.append(p)
    return nullpoints
