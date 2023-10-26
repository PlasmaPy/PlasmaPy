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

    Examples
    --------
    >>> import numpy as np
    >>> num_points = 20
    >>> x_points, dx = np.linspace(0, 10, num_points, retstep=True)
    >>> integration_matrix = construct_integration_matrix(num_points, dx)
    >>> y_points = x_points**2
    >>> integrated_y_points = integration_matrix @ y_points
    """
    # TODO: Change this so that it takes in an array of points and then integrates instead of passing around this integration object.
    # TODO: This should be a slightly different matrix as the last term is -1 but should actually follow the last equation of (D.22) for all lines.
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


def integrate(integrand, start, end, point_spacing):
    f"""Integrate an array where the start and end indeces can be floats.
    
    This will only calculate the integrand where it is real.
    The integration is defined in the `CAL` function on 
    page 37 - 39 of the computer program listing, appendix I.

    Parameters
    ----------
    integrand : `numpy.ndarray`
        2D-array to integrate. Of shape `(num_x_points, num_y_points)`.
    starts, ends : `numpy.ndarray`
        1D-array of indeces of the `integrand` array to integrate between. Of 
        shape `(num_x_points,)`. Additionally, `0 <= starts <= num_y_points - 1` 
        and `0 <= ends <= num_y_points - 1`.
    point_spacing : `float`
        Distance between each sample point.

    Notes
    -----
    This operates by using the integration described in appendix D for 
    integration between the indeces closest to the `starts` and `ends`. If the 
    start or end point is not an integer then a Taylor expansion is done around 
    either the closest integer index or one away from that depending on if the 
    integrand is a complex number.
    """
    start_float = np.where(start <= end, start, end)
    end_float = np.where(start <= end, end, start)
    reverse_integration = np.where(start <= end, 1, -1)

    start_index = np.ceil(start_float)
    end_index = np.floor(end_float)

    integral = np.zeros(integrand.shape[0], dtype=float)

    # Integrate up to the start index.
    interpolated_value = np.zeros_like(integral)
    # Line 231.
    interpolated_value[start_float == start_index] = integrand[start_float == start_index]
    np.logical_and(integrand[])
    integral += np.where(start_float != start_index, )
