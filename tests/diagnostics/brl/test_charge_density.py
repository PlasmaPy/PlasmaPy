"""Tests for calculating charge density."""

import numpy as np
import pytest
from hypothesis import example, given
from hypothesis.extra.numpy import arrays
from hypothesis.strategies import (
    booleans,
    composite,
    floats,
    integers,
    just,
    shared,
    tuples,
)
from scipy.integrate import quad
from scipy.special import erfc

from plasmapy.diagnostics.brl import charge_density

num_array_elements_including_empty = shared(integers(min_value=0, max_value=100))
num_array_elements_non_empty = shared(integers(min_value=1, max_value=100))
num_array_elements_small = shared(integers(min_value=1, max_value=10))


class Test_g:
    """Test the function `g`."""

    @pytest.fixture(autouse=True)
    def _caplog_fixture(self, caplog):
        self._caplog = caplog

    @given(
        x=arrays(float, num_array_elements_non_empty, elements=floats()).filter(
            lambda x: np.any(x < 0)
        )
    )
    def test_g_negative_x(self, x):
        """Test that the function raises a warning when x is negative."""
        charge_density._g(x)
        assert "WARNING" in [rec.levelname for rec in self._caplog.records]

    @given(
        x=arrays(
            float, num_array_elements_including_empty, elements=floats(min_value=0)
        )
    )
    def test_g_non_negative_x(self, x):
        """Test that the function returns the expected value for non-negative x."""
        result = charge_density._g(x)
        assert result.size == x.size

        if x.size > 0:
            assert np.all(charge_density._g(x) >= 0)

    @given(
        x=arrays(
            float,
            num_array_elements_non_empty,
            elements=floats(min_value=0, max_value=20),
        )
    )
    def test_g_bounded_x(self, x):
        """Test that the function returns the expected value for bounded x."""
        assert np.allclose(
            charge_density._g(x), np.pi**0.5 / 2 * np.exp(x**2) * erfc(x)
        )


class TestEta1:
    """Test the function `eta_1`."""

    @given(
        A=arrays(
            float,
            num_array_elements_non_empty,
            elements=floats(max_value=0, exclude_max=True),
        ),
        chi=arrays(float, num_array_elements_non_empty, elements=floats()),
        spherical=booleans(),
    )
    def test_eta_1_negative_A(self, A, chi, spherical):
        """Test that the function raises a ValueError when A is negative."""
        with pytest.raises(ValueError):
            charge_density.eta_1(A, chi, spherical)

    @given(
        A_chi_array=arrays(float, tuples(just(2), num_array_elements_non_empty))
        .map(lambda arr: np.sort(arr, axis=0))
        .filter(lambda arr: np.any(arr[0] < arr[1])),
        spherical=booleans(),
    )
    def test_eta_1_A_less_than_chi(self, A_chi_array, spherical):
        """Test that the function raises a ValueError when A is less than chi."""
        A, chi = A_chi_array
        with pytest.raises(ValueError):
            charge_density.eta_1(A, chi, spherical)

    @given(
        chi_A_array=arrays(
            float,
            tuples(just(2), num_array_elements_including_empty),
            elements=floats(min_value=0),
        ).map(lambda arr: np.sort(arr, axis=0)),
    )
    def test_eta_1_cylindrical(self, chi_A_array):
        """Test that the function returns the correct value for cylindrical coordinates."""
        chi, A = chi_A_array

        is_spherical = False
        if A.size == 0:
            assert charge_density.eta_1(A, chi, is_spherical).size == 0
        else:
            assert np.all(charge_density.eta_1(A, chi, is_spherical) == 0)

    @given(
        chi_A_array=arrays(
            float,
            tuples(just(2), num_array_elements_including_empty),
            elements=floats(min_value=0),
        ).map(lambda arr: np.sort(arr, axis=0)),
    )
    def test_eta_1_spherical(self, chi_A_array):
        """Test that the function returns a value within the expected bounds."""
        chi, A = chi_A_array

        is_spherical = True
        eta_1_result = charge_density.eta_1(A, chi, is_spherical)

        assert eta_1_result.size == A.size
        if A.size > 0:
            assert np.all(eta_1_result <= 0)
            theoretical_value = np.where(
                A == np.inf,
                0,
                -np.exp(-A)
                / np.pi**0.5
                * ((A - chi) ** 0.5 + charge_density._g((A - chi) ** 0.5)),
            )
            assert np.allclose(eta_1_result, theoretical_value)


x_some_not_one = shared(
    arrays(
        float,
        num_array_elements_non_empty,
        elements=floats(min_value=0, exclude_min=True, max_value=1, exclude_max=False),
    ).filter(lambda x: np.any(x < 1))
)


def kappa(chi, chi_p, x):
    return np.where(1 - x**2 == 0, 0, (chi - x**2 * chi_p) / (1 - x**2))


@composite
def eta_2_input_values(
    draw,
    num_array_elements,
    x_elements=None,
    A_minus_kappa_elements=None,
):
    """
    Draw an input value for calculating eta 2.

    Parameters
    ----------
    draw : func
        Function passed to this `composite` to allow drawing of examples.
    num_array_elements : int
        Number of array elements to use in each array.
    x_elements : `hypothesis.strategy`
        Elements to put in the `x` array. Should be between 0 and 1 inclusive.
    A_minus_kappa_elements : `hypothesis.strategy`
        Elements to put in the `A` array after subtracting by `TestEta2.kappa`.

    Returns
    -------
    x : np.array[float]
    chi_p : float
    chi : np.array[float]
    A : np.array[float]
    spherical : bool
    """
    if x_elements is None:
        x_elements = floats(min_value=0, max_value=1)
    if A_minus_kappa_elements is None:
        A_minus_kappa_elements = floats(min_value=0, max_value=10**6)

    x = draw(arrays(float, num_array_elements, elements=x_elements))
    chi_p = draw(floats(min_value=-(10**6), max_value=10**6).filter(lambda x: x != 0))
    chi = draw(
        arrays(
            float,
            num_array_elements,
            elements=floats(min_value=-(10**6), max_value=10**6),
        )
    )
    chi[np.isclose(x, 1)] = chi_p
    A_minus_kappa = draw(
        arrays(float, num_array_elements, elements=A_minus_kappa_elements)
    )
    spherical = draw(booleans())
    return (
        x,
        chi_p,
        chi,
        A_minus_kappa + np.max([kappa(chi, chi_p, x), np.full_like(x, chi_p)], axis=0),
        spherical,
    )


class TestEta2:
    """Test the `eta_2` function."""

    @given(
        eta_2_input_values(
            num_array_elements_non_empty,
            A_minus_kappa_elements=floats(allow_subnormal=False),
        ).filter(
            lambda result: np.any((result[3] < kappa(*result[2::-1]))[result[0] != 1])
        )
    )
    def test_eta_2_A_less_than_kappa(self, result):
        """Test that the function raises a ValueError when A is less than kappa."""
        x, chi_p, chi, A, spherical = result
        with pytest.raises(ValueError):
            charge_density.eta_2(A, chi, chi_p, x, spherical)

    @given(eta_2_input_values(0))
    def test_eta_2_empty_array(self, result):
        """Test that the function returns an empty array when the input arrays are empty."""
        x, chi_p, chi, A, spherical = result
        assert charge_density.eta_2(A, chi, chi_p, x, spherical).size == 0

    @given(eta_2_input_values(num_array_elements_non_empty, x_elements=just(1)))
    def test_eta_2_all_x_are_one(self, result):
        """Test that the function returns the expected value when all x are 1."""
        x, chi_p, chi, A, spherical = result
        assert np.all(
            charge_density.eta_2(A, chi, chi_p, x, spherical)
            == charge_density.eta_5(A, spherical=spherical)
        )

    @staticmethod
    def manually_integrate_eta_2_spherical(A, chi, chi_p, x):
        """Integrate equation (E.2) at a single point."""

        def integrand(beta):
            return np.exp(-beta) * (beta - chi - x**2 * (beta - chi_p)) ** 0.5

        result = quad(integrand, A, np.inf)
        return -1 / np.pi**0.5 * result[0]

    @staticmethod
    def manually_integrate_eta_2_cylindrical(A, chi, chi_p, x):
        """Integrate equation (E.10) at a single point."""

        def integrand(beta):
            return np.exp(-beta) * np.arcsin(
                ((x**2 * (beta - chi_p)) / (beta - chi)) ** 0.5
            )

        if A == np.inf:
            return 0.0
        result = quad(integrand, A, np.inf)
        return 1 / np.pi * result[0]

    @given(
        eta_2_input_values(
            num_array_elements_small,
            x_elements=floats(
                min_value=0, max_value=1, allow_subnormal=False, exclude_min=True
            ),
        )
    )
    # @example(result=(np.array([0.5**0.5]), 0.0, np.array([0.]), np.array([0.]), False))
    # @example(result=(np.array([0.]), 1.0, np.array([0.]), np.array([0.]), False))
    # @example(result=(np.array([1., 0.5]), 1.0, np.array([1., 0.]), np.array([1., 1.]), False))
    @example(
        result=(
            np.array([0.99999]),
            0.0625,
            np.array([0.0625]),
            np.array([0.0625]),
            False,
        )
    )
    @example(result=(np.array([0.5]), 1.0, np.array([-710.0]), np.array([1.0]), False))
    @example(
        result=(
            np.array([0.5, 1.0, 1.0]),
            -2130.0,
            np.array([0.0, -2130.0, -2130.0]),
            np.array([710.0, 0.0, 0.0]),
            False,
        )
    )
    @example(result=(np.array([1.0]), 1.0, np.array([1.0]), np.array([1.0]), False))
    @example(
        result=(
            np.array([1.0, 1.0]),
            -1.5,
            np.array([-1.5, -1.5]),
            np.array([-1.5, 0.0]),
            False,
        )
    )
    def test_eta_2_same_as_manual_integration(self, result):
        """
        Test that the function gives the same result as manually integrating
        each term using the initial equations.
        """
        x, chi_p, chi, A, spherical = result

        eta_2_actual_values = []
        for A_value, chi_value, x_value in zip(A, chi, x, strict=False):
            if spherical:
                eta_2_actual_value = self.manually_integrate_eta_2_spherical(
                    A_value, chi_value, chi_p, x_value
                )
            else:
                eta_2_actual_value = self.manually_integrate_eta_2_cylindrical(
                    A_value, chi_value, chi_p, x_value
                )

            eta_2_actual_values.append(eta_2_actual_value)
        eta_2_actual_values = np.array(eta_2_actual_values)

        eta_2_calculated = charge_density.eta_2(A, chi, chi_p, x, spherical=spherical)
        assert np.allclose(eta_2_calculated, eta_2_actual_values, atol=1e-4, rtol=1e-3)


@pytest.mark.parametrize(
    ("beta_M", "correct_index"),
    [(0, 8), (6, 0), (4.5, 3), (3.5, 3), (2.5, 4), (0.5, 6)],
)
def test_find_beta_M_index(beta_M, correct_index):
    """Test that the beta_M index is correct."""
    beta_G = np.array([3, 4, 5, 3, 2, 1, 0.5, 0.1])

    assert charge_density.find_beta_M_index(beta_M, beta_G) == correct_index


@pytest.mark.parametrize(
    ("beta_M", "beta_M_index", "beta_G_func", "omega_M", "omega_G_func"),
    [
        (1, 1, lambda x: 2 - x, 1, lambda x: x),
        (1.8, 0, lambda x: 2 - x, 0.2, lambda x: x),
        (1.8, 1, lambda x: 2 - x, 0.2, lambda x: x),
        (0.5, 8, lambda x: (x - 9) ** 2, 9 - 0.5**0.5, lambda x: x),
        (0.5, 9, lambda x: (x - 9) ** 2, 9 - 0.5**0.5, lambda x: x),
        (-5, 2, lambda x: -(x**2), 5**0.5, lambda x: x),
        (1, 1, lambda x: 2 - x, 2, lambda x: 2 * x),
        (1.8, 1, lambda x: 2 - x, 1.2, lambda x: 1 + x),
        (0.5, 9, lambda x: (x - 9) ** 2, (9 - 0.5**0.5) ** 2, lambda x: x**2),
        (-5, 2, lambda x: -(x**2), 1 - 5**0.5, lambda x: 1 - x),
    ],
)
def test_estimate_omega_M(beta_M, beta_M_index, beta_G_func, omega_M, omega_G_func):
    """Test that the returned omega_M is close."""
    sample_points = np.arange(0, 10)
    # beta_M_index = charge_density.find_beta_M_index(beta_M, beta_G_func(sample_points))
    assert np.isclose(
        charge_density.estimate_omega_M(
            beta_M,
            beta_M_index,
            beta_G_func(sample_points),
            omega_G_func(sample_points),
        ),
        omega_M,
    )


@pytest.mark.parametrize(
    ("chi", "x", "spherical", "expected_charge_density"),
    [
        (10 * np.ones(10), np.linspace(1, 0.1, num=10), False, np.zeros(10)),
        (10 * np.ones(10), np.linspace(1, 0.1, num=10), True, np.zeros(10)),
    ],
)
def test_delta_function_charge_density_high_potential(
    chi, x, spherical, expected_charge_density
):
    """Test that the delta function charge density is correct."""
    omega_G = -1 / (2 * x) * np.gradient(chi, x)
    beta_G = chi - x / 2 * np.gradient(chi, x)

    eta = charge_density.delta_function_charge_density(
        chi, x, omega_G, beta_G, spherical
    )
    assert eta.size == chi.size
    assert np.all(np.isclose(eta, expected_charge_density))


@pytest.mark.parametrize(
    ("chi", "x", "spherical"),
    [
        (np.zeros(10), np.linspace(1, 0.1, num=10), False),
        (np.zeros(10), np.linspace(1, 0.1, num=10), True),
    ],
)
def test_delta_function_charge_density_low_potential(chi, x, spherical):
    """Test that the delta function charge density is greater than 0 and increasing."""
    omega_G = -1 / (2 * x) * np.gradient(chi, x)
    beta_G = chi - x / 2 * np.gradient(chi, x)

    eta = charge_density.delta_function_charge_density(
        chi, x, omega_G, beta_G, spherical
    )
    assert eta.size == chi.size
    assert np.all(eta >= 0)
    assert np.all(np.diff(eta) >= 0)
