"""Calculate the initial potential and net charge density profiles to start iterations from."""
import numpy as np


def evaluate_polynomial(coefficients, normalized_probe_potential, normalized_probe_radius, x_points, dx_ds_points, spherical=True):
    r"""
    Parameters
    ----------
    coefficients : `numpy.ndarray`
        An array of coefficients for a polynomial in increasing order.
    normalized_probe_potential : `float`
        The normalized probe potential as defined in `~plasmapy.diagnostics.brl.normalizations.get_normalized_potential`.
    normalized_probe_radius : `float`
        The radius of the probe normalized to the attracted particle Debye length as defined in `~plasmapy.diagnostics.brl.normalizations.get_normalized_probe_radius`.
    x_points : `numpy.ndarray`
        The normalized inverse distance from the origin.
    dx_ds_points : `numpy.ndarray`
        The derivative of `x_points` with respect to `s`.
    spherical : `bool`, optional
        Whether the polynomial is in spherical coordinates. Default is `True`.

    Returns
    -------
    z : `numpy.ndarray`
        This is used to calculate the potential.
    dxi_ds, eta : `numpy.ndarray`
        The derivative of the potential and the net charge density at each point.

    Notes
    -----
    `coefficients` is the same as `A` in the thesis. This follows the code on 
    page 25 of the code in the thesis.

    This operates by computing a polynomial for :math:`z(s)` of order from 2 to 
    8 (inclusive) if `spherical` or from 1 to 7 (inclusive) in 
    :math:`x = 1 / r`. Let the minimum and maximum order be given by :math:`p_i` 
    and :math:`p_f`. The `coefficients`, :math:`c_i,\, i \in \{1, 2, \dots, 8\}`, 
    are normalized by their sum divided by the normalized probe potential,

    .. math::

        N &= \sum_{i=p_i}^{p_f} \frac{c_i}{\chi_p} \\
        \overline{c}_i &= \frac{c_i}{N}.

    Then we can calculate :math:`z(s)`, :math:`\frac{d\xi}{ds}`, and :math:`\eta` as

    .. math::

        z(s) &= \sum_{i=p_i}^{p_f} \overline{c}_i x(s)^i, \\
        \frac{d\xi}{ds} &= \sum_{i=p_i}^{p_f} i \overline{c}_i x(s)^{i-1} \frac{dx}{ds}, \\
        \eta &= \frac{-1}{R_p^2} \sum_{i=p_i}^{p_f} a_i i \overline{c}_i x(s)^{i + 2}.

    :math:`a_i` is :math:`i - 1` if `spherical` and :math:`i` if not.
    """
    min_order = 2 if spherical else 1 # `JAA` in the thesis.
    max_order = 8 if spherical else 7 # `JBB` in the thesis.
    # Line 24. `normalization_inverse` is `ERROR` in the thesis.
    normalization_inverse = normalized_probe_potential / np.sum(coefficients[min_order - 1:max_order], dtype=float)
    # Line 27.
    new_coefficients = np.array(coefficients, dtype=float, copy=True)
    new_coefficients[min_order - 1:max_order] *= normalization_inverse

    z = np.zeros_like(x_points)
    dxi_ds = np.zeros_like(x_points)
    eta = np.zeros_like(x_points)
    for power in range(min_order, max_order + 1):
        z += new_coefficients[power - 1] * x_points**power
        dxi_ds += power * new_coefficients[power - 1] * x_points**(power - 1) * dx_ds_points
        eta -= power * (power - spherical) * new_coefficients[power - 1] * x_points**(power + 2) / normalized_probe_radius**2

    return z, dxi_ds, eta


def evaluate_debye_and_power_law(
    power_law_coefficient, normalized_probe_potential, 
):
    pass


if __name__ == "__main__":
    coefficients = np.ones(9)
    # coefficients = np.arange(9, dtype=float)
    normalized_probe_potential = 1
    gamma = 1
    num_points = 10
    spherical = True

    from plasmapy.diagnostics.brl.net_spacing import get_s_points, get_x_and_dx_ds
    s_points = get_s_points(num_points, 0.4)
    x_points, dx_ds_points = get_x_and_dx_ds(s_points, gamma**0.5)

    import matplotlib.pyplot as plt
    z, dxi_ds, eta = evaluate_polynomial(
        coefficients, normalized_probe_potential, gamma, num_points, 
        x_points, dx_ds_points, spherical=spherical
    )

    fig, ax = plt.subplots()
    ax.plot(1 / x_points, z, label="Z", marker="o")
    ax.plot(1 / x_points, dxi_ds, label="Potential derivative", marker="o")
    ax.plot(1 / x_points, eta, label="Net charge density", marker="o")
    ax.set_xlabel("r")
    ax.legend()

    print(f"{coefficients = }")
    print(f"{normalized_probe_potential = }")
    print(f"{gamma = }")
    print(f"{num_points = }")
    print(f"{x_points = }")
    print(f"{dx_ds_points = }")
    print(f"{spherical = }")
    print(f"{z = }")
    print(f"{dxi_ds = }")
    print(f"{eta = }")

    plt.show()

