import numpy as np


def construct_integration_matrix(num_points, point_spacing):
    r"""Construct the matrix that is used for integrating a function, :math:`y(s)`, equally sampled in space.

    Parameters
    ----------
    num_points : `int`
        Number of sample points.
    point_spacing : `float`
        Distance between each sample point.

    Returns
    -------
    integration_matrix : `numpy.ndarray`
        2D-array of shape `(num_points, num_points)` of constants for each point, :math:`y(s_i)`.

    Notes
    -----
    This follows along the equations in appendix D from (D.20) onwards. We've
    changed the indexing to go from `0` to `n - 1` instead of `1` to `n`.

    Let :math:`Y(s)` be the integrated function where
    :math:`Y(s) = int_{s_0}^{s_{n - 1}} y(s) ds, n = num_points`. Then we can
    get :math:`Y(s_i)` for all :math:`i \in {0, ..., n - 1}` by doing
    `Y_values = integration_matrix @ y_values + integration_constant`,
    where :math:`integration_constant = Y(s_0)`.
    """
    if num_points < 3:
        raise ValueError("Can't create integration matrix for fewer than 3 points.")

    integration_matrix = np.zeros((num_points, num_points), dtype=float)
    for Y_index in range(1, num_points):
        if Y_index == 1:  # Y_1 is a parabolic arc through y_0, y_1, and y_2.
            integration_matrix[Y_index, : Y_index + 2] = np.array([5, 8, -1]) / 12
        elif (
            Y_index < num_points - 1
        ):  # Y_i is a cubic arc through y_{i - 2}, y_{i - 1}, y_{i}, and y_{i + 1}.
            integration_matrix[Y_index, : Y_index + 1] = integration_matrix[
                Y_index - 1, : Y_index + 1
            ]
            integration_matrix[Y_index, Y_index - 2 : Y_index + 2] += (
                np.array([-1, 13, 13, -1]) / 24
            )
        elif (
            Y_index == num_points - 1
        ):  # Y_{n - 1} is a parabolic arc through y_{n - 3}, y_{n - 2}, and y_{n - 1}.
            integration_matrix[Y_index, : Y_index + 1] = integration_matrix[
                Y_index - 1, : Y_index + 1
            ]
            integration_matrix[Y_index, Y_index - 2 : Y_index + 1] += (
                np.array([-1, 8, 5]) / 12
            )

    integration_matrix *= point_spacing
    return integration_matrix
