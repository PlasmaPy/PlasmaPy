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
    # Find the cell
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
):
    spine_list = np.array([])
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
        seed1 = nullp.loc + lambda3_vector.reshape(3, 1)
        seed2 = nullp.loc - lambda3_vector.reshape(3, 1)

        def f(s, p):
            approx = B_approx(vspace, p)
            result = np.concatenate(
                (approx[0] / approx, approx[1] / approx, approx[2] / approx), axis=1
            )
            print("#")
            print(result)
            print(result.shape)
            return result

        sol1 = scipy.integrate.solve_ivp(
            fun=f,
            t_span=[0, 10],
            y0=seed1.reshape(
                3,
            ),
            vectorized=True,
        )
        sol2 = scipy.integrate.solve_ivp(
            fun=f,
            t_span=[0, 10],
            y0=seed2.reshape(
                3,
            ),
            vectorized=True,
        )
        spine_list.append(np.concatenate((sol1, sol2)))
    return spine_list


# Returns list of rings
def fan_find(
    vspace,
    nullpoint_list,
):
    ring_list = np.array([])
    # Initializon
    # Expansion
    # Null Breaking
    # Trace Back
    return None


# Returns list of nullpoints, spines, fans, and seperaators (The wrapper function provided to the user)
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
    spines = spine_find(vspace, nullpoints)
    fan = fan_find(vspace, nullpoints)
    return nullpoints, spines, fan
