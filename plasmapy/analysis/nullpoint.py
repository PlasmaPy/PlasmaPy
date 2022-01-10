import matplotlib.pyplot as plt
import numpy as np
import warnings

# Declare Constants
ATOL = 10 ** (-10)
MAX_DIVIDE = 10


class Point:
    def __init__(self, loc, field):
        self.loc = loc
        self.type = type
        self.field = field

    def getLoc(self):
        return self.loc

    def getField(self):
        return self.field


class NullPoint(Point):
    def __init__(self, null_loc, type):
        super().__init__(null_loc, [0, 0, 0])
        self.type = type

    def getLoc(self):
        return self.loc

    def getType(self):
        return self.type

    def isEqual(self, point):
        return (
            np.isclose(self.getLoc()[0], point.getLoc()[0])
            and np.isclose(self.getLoc()[1], point.getLoc()[1])
            and np.isclose(self.getLoc()[2], point.getLoc()[2])
        )


def vector_space(func, x_range, y_range, z_range, precision=[0.01, 0.01, 0.01]):
    x_den = int(np.around((x_range[1] - x_range[0]) / precision[0]) + 1)
    y_den = int(np.around((y_range[1] - y_range[0]) / precision[1]) + 1)
    z_den = int(np.around((z_range[1] - z_range[0]) / precision[2]) + 1)
    x, y, z = np.meshgrid(
        np.linspace(x_range[0], x_range[1], x_den),
        np.linspace(y_range[0], y_range[1], y_den),
        np.linspace(z_range[0], z_range[1], z_den),
        indexing="ij",
    )
    u, v, w = func(x, y, z)
    # ax = plt.figure().add_subplot(projection='3d')
    # ax.quiver(x, y, z, u, v, w, length=0.2, normalize=True)
    # plt.show()
    dx = np.float128((x_range[1] - x_range[0]) / (x_den - 1))
    dy = np.float128((y_range[1] - y_range[0]) / (y_den - 1))
    dz = np.float128((z_range[1] - z_range[0]) / (z_den - 1))

    return np.array([x, y, z]), np.array([u, v, w]), np.array([dx, dy, dz])


def trilinear_coeff_cal(vspace, cell):
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
    x1 = float(x0 + deltax)
    y0 = float(vspace[0][1][f000[0]][f000[1]][f000[2]])
    y1 = float(y0 + deltay)
    z0 = float(vspace[0][2][f000[0]][f000[1]][f000[2]])
    z1 = float(z0 + deltaz)
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

    ax, bx, cx, dx, ex, fx, gx, hx = np.linalg.solve(A, sx)
    ay, by, cy, dy, ey, fy, gy, hy = np.linalg.solve(A, sy)
    az, bz, cz, dz, ez, fz, gz, hz = np.linalg.solve(A, sz)

    return np.array(
        [
            [ax, bx, cx, dx, ex, fx, gx, hx],
            [ay, by, cy, dy, ey, fy, gy, hy],
            [az, bz, cz, dz, ez, fz, gz, hz],
        ]
    )


def trilinear_approx(vspace, cell):
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


"""
    Return a true or false based on weather
    a grid cell passes the reduction phase.

    Parameters
    ----------

    vspace :

    cell:

    Returns
    -------
    doesPassReduction : bool
        True if a grid cell passes the reduction phase.
        False, otherwise.

    Raises
    ------


    Warns
    -----


    Notes
    -----


    Examples
    --------
"""


def reduction(vspace, cell):
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


"""
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
    roots : float[]
        A 1 by 2 array holding the two roots

    Raises
    ------


    Warns
    -----


    Notes
    -----


    Examples
    --------
"""


def bilinear_root(a1, b1, c1, d1, a2, b2, c2, d2):
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

    # Normalizing coefficients
    coeffs = [
        [ax, bx, cx, dx, ex, fx, gx, hx],
        [ay, by, cy, dy, ey, fy, gy, hy],
        [az, bz, cz, dz, ez, fz, gz, hz],
    ]
    normalized_coeffs = []
    for elem in coeffs:
        normalized_coeffs.append(list(map(lambda arr: float(arr[0]), elem)))

    ax, bx, cx, dx, ex, fx, gx, hx = normalized_coeffs[0]
    ay, by, cy, dy, ey, fy, gy, hy = normalized_coeffs[1]
    az, bz, cz, dz, ez, fz, gz, hz = normalized_coeffs[2]
    # print(ax, bx, cx, dx, ex, fx, gx, hx)

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


global divide
divide = 0


def locate_null_point(vspace, cell, n, err):
    global divide
    # def get_starting_pos():

    # Calculating the Jacobian and trillinear approximation functions for the cell
    tlApprox = trilinear_approx(vspace, cell)
    jcb = jacobian(vspace, cell)

    deltax, deltay, deltaz = vspace[2]
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
    new_vspace = vector_space(
        tlApprox,
        [pos_000[0], pos_111[0]],
        [pos_000[1], pos_111[1]],
        [pos_000[2], pos_111[2]],
        [deltax / 2, deltay / 2, deltaz / 2],
    )
    return nullpoint(new_vspace)


def inNullList(elem, lst):
    for p in lst:
        if p.isEqual(elem):
            return True
    return False


def nullpoint(vspace, MAX_ITERATIONS=500, err=10 ** (-10)):
    nullpoints = []
    for i in range(len(vspace[0][0]) - 1):
        for j in range(len(vspace[0][0][0]) - 1):
            for k in range(len(vspace[0][0][0][0]) - 1):
                if reduction(vspace, [i, j, k]):
                    if trillinear_analysis(vspace, [i, j, k]):
                        loc = locate_null_point(vspace, [i, j, k], MAX_ITERATIONS, err)
                        if not isinstance(loc, type(None)):
                            p = NullPoint(loc, "V")
                            if not inNullList(p, nullpoints):
                                nullpoints.append(p)
    return nullpoints


"""
Testing and Examples
"""

# def vector_Space_Func(x,y,z):
#     return [2*y-z-5.5, 3*x+z-22, x**2-11*x+y+24.75]
# def vector_Space_Func(x,y,z):
#     return [(y-2.5)**2-0.01, (z-2.5)**2-0.01, (x-2.5)**2-0.01]
# def vector_Space_Func(x, y, z):
#     return [-1.80 - 2.57 * x + 6.92 * y + 0.44 * z + 0.02 * y * z + 0.46 * x * z - 1.40 * x * y-0.10472,
#             0.44 - 3.05 * x + 2.09 * y + 0.20 * z - 0.46 * y * z - 0.34 * x * z + 1.46 * x * y - 0.235314,
#             -0.67 - 2.30 * x + 8.69 * y + 0.48 * z + 1.40 * y * z - 1.46 * x * z - 8.29 * x * y - 0.510349]
def vector_Space_Func(x, y, z):
    a = (
        -1.80
        - 2.57 * x
        + 6.92 * y
        + 0.44 * z
        + 0.02 * y * z
        + 0.46 * x * z
        - 1.40 * x * y
    )
    b = (
        0.44
        - 3.05 * x
        + 2.09 * y
        + 0.20 * z
        - 0.46 * y * z
        - 0.34 * x * z
        + 1.46 * x * y
    )
    c = (
        -0.67
        - 2.30 * x
        + 8.69 * y
        + 0.48 * z
        + 1.40 * y * z
        - 1.46 * x * z
        - 8.29 * x * y
    )
    return [a, b, c]


vspace1 = vector_space(vector_Space_Func, [0, 4], [0, 4], [0, 4], [1, 1, 1])
npoints = nullpoint(vspace1)
if len(npoints) == 0:
    print("No Nullpoints Found.")
# print(vspace1[2])
for points in npoints:
    print(points.getLoc())
