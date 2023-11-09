"""Tests for doing integration."""

import numpy as np
import pytest

from plasmapy.diagnostics.brl import integration


class Test__integrate:
    r"""Test the integrate function in integration.py"""

    num_points = 10
    sample_points = np.arange(num_points)
    sample_spacing = 1
    constant_integrand = np.ones(num_points)[np.newaxis, :]

    def constant_integral_func(x):
        return x

    real_integrand = (sample_points**2 + 3)[np.newaxis, :]
    real_integral = (1 / 3 * sample_points**3 + 3 * sample_points)[np.newaxis, :]

    def real_integral_func(x):
        return 1 / 3 * x**3 + 3 * x

    def test_invalid_start_after_end(self):
        r"""Test error integral if start before end."""
        start = np.array([3.3])
        end = start - 1

        with pytest.raises(ValueError):
            integration.integrate(
                self.constant_integrand,
                start,
                end,
                np.zeros_like(start),
                np.zeros_like(start),
                self.sample_spacing,
            )

    def test_multiple_integrals(self):
        r"""Test multiple integrals are computed if multiple start points given."""
        start = np.array([0, 2, 5])
        end = np.array([0, 3, 8])
        integrand = np.tile(self.constant_integrand, [start.size, 1])

        assert np.allclose(
            integration.integrate(
                integrand,
                start,
                end,
                np.zeros_like(start),
                np.zeros_like(start),
                self.sample_spacing,
            ),
            end - start,
        )

    @pytest.mark.parametrize(
        ("start_index", "end_index", "error"),
        [
            # Test many combinations of spacing between points as different formulas are used for each.
            (2, 2, 10**-6),
            (2, 4, 10**-6),
            (2, 5, 10**-6),
            (2, 6, 10**-6),
            (2, 7, 10**-6),
            (2, 8, 10**-6),
            (1.3, 2, 10**-6),
            (1.3, 3, 10**-6),
            (1.3, 4, 10**-6),
            (2, 2.3, 10**-6),
            (2, 3.3, 10**-6),
            (2, 4.3, 10**-6),
            (1.3, 1.7, 10**-6),
            (1.3, 2.7, 10**-6),
            (1.3, 3.7, 10**-6),
            (1.3, 4.7, 10**-6),
            # The large error for the starting and ending points being one apart is due to using the trapezoidal rule for the integration.
            (2, 3, 2 * 10**-1),
        ],
    )
    def test_non_complex_start_and_end_combinations(
        self, start_index, end_index, error
    ):
        r"""Test integration for start and end points that aren't complex."""
        start = np.array([start_index])
        end = np.array([end_index])

        assert np.isclose(
            integration.integrate(
                self.constant_integrand,
                start,
                end,
                np.zeros_like(start),
                np.zeros_like(start),
                self.sample_spacing,
            ),
            (
                Test__integrate.constant_integral_func(end)
                - Test__integrate.constant_integral_func(start)
            ),
            atol=error,
        )
        assert np.isclose(
            integration.integrate(
                self.real_integrand,
                start,
                end,
                np.zeros_like(start),
                np.zeros_like(start),
                self.sample_spacing,
            ),
            (
                Test__integrate.real_integral_func(end)
                - Test__integrate.real_integral_func(start)
            ),
            atol=error,
        )
