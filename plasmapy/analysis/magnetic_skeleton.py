import numpy as np
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

_EQUALITY_ATOL = 1e-10

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
        if not np.isclose(np.trace(M), 0, atol=_EQUALITY_ATOL):
            raise NonZeroDivergence()
        eigen_vals, eigen_vectors = np.linalg.eig(M)

        # Generate seeds
        lambda3 = None
        lambda3_vector = None
        if (
            np.isreal(eigen_vals[0])
            and np.sign(np.real(eigen_vals[0])) * np.sign(np.real(eigen_vals[1])) < 0
            and np.sign(np.real(eigen_vals[0])) * np.sign(np.real(eigen_vals[2])) < 0
        ):
            lambda3 = eigen_vals[0]
            lambda3_vector = eigen_vectors[:, 0]
        elif (
            np.isreal(eigen_vals[1])
            and np.sign(np.real(eigen_vals[1])) * np.sign(np.real(eigen_vals[0])) < 0
            and np.sign(np.real(eigen_vals[1])) * np.sign(np.real(eigen_vals[2])) < 0
        ):
            lambda3 = eigen_vals[1]
            lambda3_vector = eigen_vectors[:, 1]
        elif (
            np.isreal(eigen_vals[2])
            and np.sign(np.real(eigen_vals[2])) * np.sign(np.real(eigen_vals[0])) < 0
            and np.sign(np.real(eigen_vals[2])) * np.sign(np.real(eigen_vals[1])) < 0
        ):
            lambda3 = eigen_vals[2]
            lambda3_vector = eigen_vectors[:, 2]
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
            # print("#")
            # print(result)
            # print(result.shape)
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


# Returns list of rings
def fan_find(
    vspace,
    nullpoint_list,
    alpha,
    seed_nums=10,
    ring_nums=1000,
):
    ring_list = []
    # Initialization
    nullp = nullpoint_list.pop(0)
    # Adding the first null point as the zeroth layer
    ring_list.append(np.array(nullp.loc))
    # Constructing the first ring layer
    jcb = _trilinear_jacobian(vspace, nullp.cell)
    M = jcb(nullp.loc[0], nullp.loc[1], nullp.loc[2])
    if not np.isclose(np.trace(M), 0, atol=_EQUALITY_ATOL):
        raise NonZeroDivergence()
    eigen_vals, eigen_vectors = np.linalg.eig(M)
    # Generate seeds
    lambda1 = None
    lambda1_vector = None
    lambda2 = None
    lambda2_vector = None
    if (
        np.isreal(eigen_vals[0])
        and np.sign(np.real(eigen_vals[0])) * np.sign(np.real(eigen_vals[1])) < 0
        and np.sign(np.real(eigen_vals[0])) * np.sign(np.real(eigen_vals[2])) < 0
    ):
        lambda1 = eigen_vals[1]
        lambda1_vector = eigen_vectors[:, 1]
        lambda2 = eigen_vals[2]
        lambda2_vector = eigen_vectors[:, 2]
    elif (
        np.isreal(eigen_vals[1])
        and np.sign(np.real(eigen_vals[1])) * np.sign(np.real(eigen_vals[0])) < 0
        and np.sign(np.real(eigen_vals[1])) * np.sign(np.real(eigen_vals[2])) < 0
    ):
        lambda1 = eigen_vals[0]
        lambda1_vector = eigen_vectors[:, 0]
        lambda2 = eigen_vals[2]
        lambda2_vector = eigen_vectors[:, 2]
    elif (
        np.isreal(eigen_vals[2])
        and np.sign(np.real(eigen_vals[2])) * np.sign(np.real(eigen_vals[0])) < 0
        and np.sign(np.real(eigen_vals[2])) * np.sign(np.real(eigen_vals[1])) < 0
    ):
        lambda1 = eigen_vals[0]
        lambda1_vector = eigen_vectors[:, 0]
        lambda2 = eigen_vals[1]
        lambda2_vector = eigen_vectors[:, 1]
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
        ring1.append(p)
        v = np.cos(theta) * v + np.sin(theta) * np.cross(normal_vector, v)
    ring_list.append(np.array(ring1))
    # Expansion and Null Breaking
    for i in range(1, ring_nums):
        next_ring_layer = []
        for point in ring_list[i]:
            # Null Breaking
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

            sol = scipy.integrate.solve_ivp(
                fun=f,
                t_span=(0, 0.01),
                y0=point,
                min_step=0.01,
                max_step=0.01,
            ).y.T[1]
            next_ring_layer.append(sol)
        ring_list.append(np.array(next_ring_layer))
    for layer in ring_list:
        print("###################")
        print(layer)
    return np.array(ring_list)
    # Trace Back
    return None


# Returns list of nullpoints, spines, fans, and seperators (The wrapper function provided to the user)
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
    spines = spine_find(vspace, nullpoints, 0.1)
    fan = fan_find(vspace, nullpoints, 0.1)
    return nullpoints, spines, fan
