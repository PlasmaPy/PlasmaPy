import copy
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scipy
import warnings

from plasmapy.analysis.nullpoint import (
    _bilinear_root,
    _locate_null_point,
    _reduction,
    _trilinear_analysis,
    _trilinear_coeff_cal,
    _trilinear_jacobian,
    _vector_space,
    _vspace_iterator,
    MultipleNullPointWarning,
    NonZeroDivergence,
    null_point_find,
    NullPointError,
    NullPointWarning,
    trilinear_approx,
    uniform_null_point_find,
)


class FanPoint:
    """
    Abstract class for defining a Fan Point in 3D space.
    """

    def __init__(self, loc, index):
        self._loc = loc
        self._index = index

    def get_loc(self):
        r"""
        Returns the coordinates of the FanPoint object.
        """
        return self._loc

    def get_index(self):
        r"""
        Returns the index of the FanPoint object.
        """
        return self._index

    loc = property(get_loc)
    index = property(get_index)


class Separator:
    """
    Abstract class for defining a Separator in 3D space.
    """

    def __init__(self, points, starting_null, end_null):
        self._points = points
        self._starting_null = starting_null
        self._end_null = end_null

    def get_points(self):
        r"""
        Returns the points of the Separator object.
        """
        return self._points

    def get_starting_null(self):
        r"""
        Returns the starting null of the Separator object.
        """
        return self._starting_null

    def get_end_null(self):
        r"""
        Returns the ending null of the Separator object.
        """
        return self._end_null

    points = property(get_points)
    starting_null = property(get_starting_null)
    end_null = property(get_end_null)


class Fan:
    """
    Abstract class for defining a Fan in 3D space.
    """

    def __init__(self, ring_list, generating_null):
        self._ring_list = ring_list
        self._generating_null = generating_null

    def get_ring_list(self):
        r"""
        Returns the ring list of the Fan object.
        """
        return self._ring_list

    def get_generating_null(self):
        r"""
        Returns the generating null of the Fan object.
        """
        return self._generating_null

    ring_list = property(get_ring_list)
    generating_null = property(get_generating_null)


_EQUALITY_ATOL = 1e-10
_INDEX_CONST = 1048576
# Finds the eigen vectors that determine the fan plane
def eigen_handler(eigen_vals, eigen_vectors):
    lambda1_vector = None
    lambda2_vector = None
    lambda3_vector = None
    if (
        np.isreal(eigen_vals[0])
        and np.sign(np.real(eigen_vals[0])) * np.sign(np.real(eigen_vals[1])) < 0
        and np.sign(np.real(eigen_vals[0])) * np.sign(np.real(eigen_vals[2])) < 0
    ):
        lambda1_vector = eigen_vectors[:, 1]
        lambda2_vector = eigen_vectors[:, 2]
        lambda3_vector = eigen_vectors[:, 0]
    elif (
        np.isreal(eigen_vals[1])
        and np.sign(np.real(eigen_vals[1])) * np.sign(np.real(eigen_vals[0])) < 0
        and np.sign(np.real(eigen_vals[1])) * np.sign(np.real(eigen_vals[2])) < 0
    ):
        lambda1_vector = eigen_vectors[:, 0]
        lambda2_vector = eigen_vectors[:, 2]
        lambda3_vector = eigen_vectors[:, 1]
    elif (
        np.isreal(eigen_vals[2])
        and np.sign(np.real(eigen_vals[2])) * np.sign(np.real(eigen_vals[0])) < 0
        and np.sign(np.real(eigen_vals[2])) * np.sign(np.real(eigen_vals[1])) < 0
    ):
        lambda1_vector = eigen_vectors[:, 0]
        lambda2_vector = eigen_vectors[:, 1]
        lambda3_vector = eigen_vectors[:, 2]

    if np.iscomplex(lambda1_vector).any():
        lambda2_vector = np.imag(lambda1_vector)
        lambda1_vector = np.real(lambda1_vector)
    return np.array([[lambda1_vector, lambda2_vector], lambda3_vector])


# Approximates the field at any point using tl approximation
def B_approx(vspace, loc):
    # return np.array([np.NaN,np.NaN,np.NaN ]).reshape(3,)
    # Find the cell
    # Maybe use numpy.argwhere or numpy.where
    x, y, z = vspace[0]
    cell = [0, 0, 0]
    for i in range(len(x) - 1):
        if (x[i][0][0] <= loc[0]) and (loc[0] <= x[i + 1][0][0]):
            cell[0] = i
            break
    for j in range(len(y) - 1):
        if (y[0][j][0] <= loc[1]) and (loc[1] <= y[0][j + 1][0]):
            cell[1] = j
            break
    for k in range(len(z) - 1):
        if (z[0][0][k] <= loc[2]) and (loc[2] <= z[0][0][k + 1]):
            cell[2] = k
            break
    # Estimate
    tlApprox = trilinear_approx(vspace, cell)
    return tlApprox(loc[0], loc[1], loc[2])


# Returns list of spines for each nullpoint
def spine_find(
    vspace,
    nullpoint_list,
    alpha=1,
):
    spine_list = []
    # Maybe add a check so that you would not add the same spine twice
    for nullp in nullpoint_list:
        jcb = _trilinear_jacobian(vspace, nullp.cell)
        M = jcb(nullp.loc[0], nullp.loc[1], nullp.loc[2])
        # if not np.isclose(np.trace(M), 0, atol=_EQUALITY_ATOL):
        #     raise NonZeroDivergence()
        eigen_vals, eigen_vectors = np.linalg.eig(M)
        lambda3_vector = eigen_handler(eigen_vals, eigen_vectors)[1]
        # Generate seeds
        seed1 = nullp.loc + alpha * lambda3_vector.reshape(3, 1)
        seed2 = nullp.loc - alpha * lambda3_vector.reshape(3, 1)
        seed1 = seed1.reshape(
            3,
        )
        seed2 = seed2.reshape(
            3,
        )
        t1 = np.linalg.norm(seed1 - nullp.loc)
        t2 = np.linalg.norm(seed2 - nullp.loc)
        x, y, z = vspace[0]
        x_min = x[0][0][0]
        x_max = x[len(x) - 1][0][0]
        y_min = y[0][0][0]
        y_max = y[0][len(y) - 1][0]
        z_min = z[0][0][0]
        z_max = z[0][0][len(z) - 1]

        def event0(s, p):
            return p[0] - x_min

        def event1(s, p):
            return x_max - p[0]

        def event2(s, p):
            return p[1] - y_min

        def event3(s, p):
            return y_max - p[1]

        def event4(s, p):
            return p[2] - z_min

        def event5(s, p):
            return z_max - p[2]

        event0.terminal = True
        event1.terminal = True
        event2.terminal = True
        event3.terminal = True
        event4.terminal = True
        event5.terminal = True

        def f(s, p):
            approx = B_approx(vspace, p)
            result = np.array(
                [
                    approx[0] / np.linalg.norm(approx),
                    approx[1] / np.linalg.norm(approx),
                    approx[2] / np.linalg.norm(approx),
                ]
            ).reshape(
                3,
            )
            return result

        sol1 = scipy.integrate.solve_ivp(
            fun=f,
            y0=seed1,
            t_span=(t1, 10),
            events=[event0, event1, event2, event3, event4, event5],
        )
        sol2 = scipy.integrate.solve_ivp(
            fun=f,
            t_span=(t2, 10),
            y0=seed2,
            events=[event0, event1, event2, event3, event4, event5],
        )
        s1 = sol1.y.T
        s2 = sol2.y.T
        spine_list.append(s1)
        spine_list.append(s2)
    return np.array(spine_list)


def null_sign(vspace, nullp):
    jcb = _trilinear_jacobian(vspace, nullp.cell)
    M = jcb(nullp.loc[0], nullp.loc[1], nullp.loc[2])
    if np.isclose(np.linalg.det(M), 0, atol=_EQUALITY_ATOL):
        return 0
    elif np.linalg.det(M) < 0:
        return 1
    else:
        return -1


def seed_from_null(vspace, nullp, alpha):
    jcb = _trilinear_jacobian(vspace, nullp.cell)
    M = jcb(nullp.loc[0], nullp.loc[1], nullp.loc[2])
    eigen_vals, eigen_vectors = np.linalg.eig(M)
    lambda3_vector = eigen_handler(eigen_vals, eigen_vectors)[1]
    # Generate seeds
    seed1 = nullp.loc + alpha * lambda3_vector.reshape(3, 1)
    seed2 = nullp.loc - alpha * lambda3_vector.reshape(3, 1)
    return seed1, seed2


def can_null_break(starting_null, nullp, alpha, vspace, phi):
    # Must have opposite sign to the starting null point
    if null_sign(vspace, starting_null) * null_sign(vspace, nullp) >= 0:
        return False
    # Generate the seed points
    seed1, seed2 = seed_from_null(vspace, nullp, alpha)
    xk = nullp.loc.reshape(3, 1)
    # Getting eigen vectors
    jcb = _trilinear_jacobian(vspace, nullp.cell)
    M = jcb(nullp.loc[0], nullp.loc[1], nullp.loc[2])
    eigen_vals, eigen_vectors = np.linalg.eig(M)
    v3 = eigen_handler(eigen_vals, eigen_vectors)[1]
    # Check that seeds will be traced in opposite directions
    if (
        np.dot(np.squeeze(B_approx(vspace, seed1)), np.squeeze((seed1 - xk)), out=None)
        <= 0
    ):
        return False
    if (
        np.dot(np.squeeze(B_approx(vspace, seed2)), np.squeeze((seed2 - xk)), out=None)
        <= 0
    ):
        return False
    # Check the angle condition
    expression = np.dot(
        np.squeeze((seed2 - xk)), np.squeeze(v3), out=None
    ) / np.linalg.norm(np.squeeze((seed1 - xk)))
    print(expression)
    # if expression<= np.cos(phi):
    #     a=5/0
    #     return False
    # if expression >= (-1*np.cos(phi)):
    #     a = 5 / 0
    #     return False
    return True


def mid_point(p1, p2):
    temp = ((p1 + p2) / 2).reshape(
        3,
    )
    return temp


# May need to call multiple times in a for loop to ensure proper density
# Maybe seperate adding/removing points and then call both in while loops
def intrp_ring(ring_layer, min_density_tol, max_density_tol):
    # Populate if not enough points
    pop_ring_layer = add_intrp_ring(ring_layer, max_density_tol)
    # Remove if too many points
    rm_ring_layer = rm_intrp_ring(pop_ring_layer, min_density_tol)
    return rm_ring_layer


def add_intrp_ring(ring_layer, max_density_tol):
    pop_ring_layer = copy.deepcopy(ring_layer)
    # Populate if not enough points
    done = False
    while not done:
        done = True
        index = 0
        for i in range(len(ring_layer) - 1):
            if (
                np.linalg.norm(
                    ring_layer[i].loc - ring_layer[(i + 1) % len(ring_layer)].loc
                )
                > max_density_tol
            ):
                done = False
                index = index + 1
                point = FanPoint(
                    mid_point(
                        ring_layer[i].loc, ring_layer[(i + 1) % len(ring_layer)].loc
                    ),
                    (ring_layer[i].index + ring_layer[(i + 1) % len(ring_layer)].index)
                    / 2,
                )
                pop_ring_layer.insert(index, point)
            index = index + 1
        if (
            np.linalg.norm(ring_layer[len(ring_layer) - 1].loc - ring_layer[0].loc)
            > max_density_tol
        ):
            done = False
            index = index + 1
            point = FanPoint(
                mid_point(ring_layer[len(ring_layer) - 1].loc, ring_layer[0].loc),
                ring_layer[len(ring_layer) - 1].index + _INDEX_CONST,
            )
            pop_ring_layer.insert(index, point)
        ring_layer = copy.deepcopy(pop_ring_layer)
    return pop_ring_layer


def rm_intrp_ring(pop_ring_layer, min_density_tol):
    rm_ring_layer = copy.deepcopy(pop_ring_layer)
    # Remove if too many points
    done = False
    while not done:
        done = True
        index = 0
        for i in range(len(pop_ring_layer)):
            if (
                np.linalg.norm(
                    pop_ring_layer[i].loc
                    - pop_ring_layer[(i + 1) % len(pop_ring_layer)].loc
                )
                < min_density_tol
            ):
                done = False
                del rm_ring_layer[index]
                index = index - 1
            index = index + 1
        pop_ring_layer = copy.deepcopy(rm_ring_layer)
    return rm_ring_layer


def trace_one_step(point, step_size, vspace, backward=False):
    flag = 1
    if backward:
        flag = -1

    def f(s, p):
        approx = B_approx(vspace, p)
        result = np.array(
            [
                approx[0] / np.linalg.norm(approx),
                approx[1] / np.linalg.norm(approx),
                approx[2] / np.linalg.norm(approx),
            ]
        ).reshape(
            3,
        )
        return result

    # Make min_step a fraction of max_step so that solve_ivp can pick optimal step size
    sol = scipy.integrate.solve_ivp(
        fun=f,
        t_span=(0, flag * step_size),
        y0=point,
        min_step=(1 / 12) * step_size,
        max_step=step_size,
    ).y.T[-1]
    # Don't forget to change the index to -1 to get last point
    return sol


def init_ring(ring_list, vspace, alpha, seed_nums, nullp):
    # Initialization
    # Adding the first null point as the zeroth layer
    ring_list.append([FanPoint(nullp.loc, 0)])
    # Constructing the first ring layer
    jcb = _trilinear_jacobian(vspace, nullp.cell)
    M = jcb(nullp.loc[0], nullp.loc[1], nullp.loc[2])
    # if not np.isclose(np.trace(M), 0, atol=_EQUALITY_ATOL):
    #     raise NonZeroDivergence()
    eigen_vals, eigen_vectors = np.linalg.eig(M)
    lambda1_vector, lambda2_vector, = eigen_handler(
        eigen_vals, eigen_vectors
    )[0]

    normal_vector = np.cross(lambda1_vector, lambda2_vector).reshape(
        3,
    )
    normal_vector = (normal_vector / np.linalg.norm(normal_vector)).reshape(
        3,
    )
    v = alpha * lambda1_vector.reshape(
        3,
    )
    ring1 = []
    theta = 2 * np.pi / seed_nums

    for i in range(seed_nums):
        p = (
            nullp.loc.reshape(
                3,
            )
            + v
        )
        ring1.append(FanPoint(p, i * _INDEX_CONST))
        v = np.cos(theta) * v + np.sin(theta) * np.cross(normal_vector, v)
    ring_list.append(ring1)
    return None


# Returns list of rings
def fan_find(
    vspace,
    nullpoint_original_list,
    alpha,
    seed_nums=100,
    ring_nums=38,
    min_density_tol=10 ** -10,
    max_density_tol=10 ** -(2),
    phi=np.pi,
):
    # Make a copy of null point list
    nullpoint_list = copy.deepcopy(nullpoint_original_list)
    # Initilization of the ring list
    ring_list = []
    seperators = []
    starting_ring_layer_for_seperators = []
    nullp = nullpoint_list.pop(0)
    init_ring(ring_list, vspace, alpha, seed_nums, nullp)
    # Expansion and Null Breaking
    for i in range(1, ring_nums):
        # Expansion
        next_ring_layer = list(
            map(
                lambda elem: FanPoint(
                    trace_one_step(elem.loc, alpha, vspace), elem.index
                ),
                ring_list[i],
            )
        )
        # Density Check and Interpolation
        next_ring_layer = intrp_ring(next_ring_layer, min_density_tol, max_density_tol)
        # Null Breaking
        null_breaking_candidates = []
        prev = None
        next = None
        for k in range(len(next_ring_layer)):
            elem = next_ring_layer[k]
            for candidate_null in nullpoint_list:
                p = elem.loc
                norm = np.linalg.norm(
                    candidate_null.loc.reshape(
                        3,
                    )
                    - p.reshape(
                        3,
                    )
                )
                if np.isclose(norm, 0, atol=10 ** -2):
                    if can_null_break(nullp, candidate_null, alpha, vspace, phi):
                        prev = next_ring_layer[k - 1]
                        next = next_ring_layer[k + 1]
                        null_breaking_candidates.append([elem, candidate_null])
                        nullpoint_list.remove(candidate_null)
        for elem in null_breaking_candidates:
            key = elem[0].loc
            value = elem[1]
            # Removing p from next_ring layer
            for q in range(len(next_ring_layer)):
                if (
                    next_ring_layer[q].loc[0] == key[0]
                    and next_ring_layer[q].loc[1] == key[1]
                    and next_ring_layer[q].loc[2] == key[2]
                ):
                    index = q
                    next_ring_layer.pop(q)
                    break
            seperators.append([FanPoint(value.loc, elem[0].index)])
            starting_ring_layer_for_seperators.append(i)
            # next_ring_layer.append(value.loc)
            seed1, seed2 = seed_from_null(vspace, value, alpha)
            seed1 = seed1.reshape(
                3,
            )
            seed2 = seed2.reshape(
                3,
            )
            next_ring_layer.insert(
                index, FanPoint(seed1, prev.index + (next.index - prev.index) / 3)
            )
            next_ring_layer.insert(
                index + 1,
                FanPoint(seed2, prev.index + 2 * (next.index - prev.index) / 3),
            )
        ring_list.append(np.array(next_ring_layer))
    DATA = []
    color_arr = []
    for layer in ring_list:
        # print("###############")
        # print(len(layer))
        for p in layer:
            # print(p.index)
            DATA.append(p.loc)
            color_arr.append(p.index)
    # DATA = DATA[3:]
    Xs = []
    Ys = []
    Zs = []
    for p in DATA:
        Xs.append(p[0])
        Ys.append(p[1])
        Zs.append(p[2])
    trace = go.Scatter3d(
        x=Xs,
        y=Ys,
        z=Zs,
        mode="markers",
        marker=dict(
            size=12,
            color=color_arr,  # set color to an array/list of desired values
            colorscale="Turbo",
        ),
    )
    layout = go.Layout(title="3D Scatter plot")
    fig = go.Figure(data=[trace], layout=layout)
    fig.show()
    # Trace Back
    for n in range(len(seperators)):
        seperator_find(seperators[n], ring_list, starting_ring_layer_for_seperators[n])
    return ring_list, seperators


def seperator_find(seperator, ring_list, starting_ring_layer):
    def intrp_sep_point(ring_layer, fanpoint):
        for i in range(len(ring_layer)):
            if ring_layer[i].index == fanpoint.index:
                return ring_layer[i]
            elif ring_layer[i].index > fanpoint.index:
                return FanPoint(
                    ring_layer[i - 1].loc
                    + (ring_layer[i].loc - ring_layer[i - 1].loc)
                    * (
                        fanpoint.index / (ring_layer[i].index + ring_layer[i - 1].index)
                    ),
                    fanpoint.index,
                )
        return FanPoint(
            ring_layer[-1].loc
            + (ring_layer[0].loc - ring_layer[-1].loc)
            * (fanpoint.index / (_INDEX_CONST + 2 * ring_layer[-1].index)),
            fanpoint.index,
        )

    start = seperator[0]
    # Traverse the ring backwards
    for i in range(starting_ring_layer - 1, 0, -1):
        start = intrp_sep_point(ring_list[i], start)
        seperator.append(start)
    seperator.append(ring_list[0][0])
    return None


# Returns list of null points, spines, fans, and separators (The wrapper function provided to the user)
def magnetic_skeleton_find(
    x_arr=None,
    y_arr=None,
    z_arr=None,
    u_arr=None,
    v_arr=None,
    w_arr=None,
):
    vspace = _vector_space(
        x_arr,
        y_arr,
        z_arr,
        None,
        None,
        None,
        u_arr,
        v_arr,
        w_arr,
        None,
        None,
    )
    nullpoints = null_point_find(x_arr, y_arr, z_arr, u_arr, v_arr, w_arr)
    for p in nullpoints:
        print(p.loc)
    spines = spine_find(vspace, nullpoints, 0.1)
    fan, seperators = fan_find(vspace, nullpoints, 0.05)
    nullpoints = list(map(lambda elem: elem.loc, nullpoints))
    return nullpoints, spines, fan, seperators
