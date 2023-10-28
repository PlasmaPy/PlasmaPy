"""Tests for doing integration."""

import numpy as np

from plasmapy.diagnostics.brl import integration


class Test__integrate:
    r"""Test the integrate function in integration.py"""

    @staticmethod
    def test_equal_start_and_end_point():
        r"""Test zero integral if start and end point are equal."""
        integrand = np.ones(10)[np.newaxis, :]
        start = np.array([3.3])

        assert (
            integration.integrate(
                integrand, start, start, np.zeros_like(start), np.zeros_like(start), 1
            )
            == 0
        )
