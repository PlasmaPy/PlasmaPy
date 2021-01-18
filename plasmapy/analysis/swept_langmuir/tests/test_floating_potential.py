"""
Tests for functionality contained in
`plasmapy.analysis.swept_langmuir.floating_potential`.
"""

import numpy as np
import pytest

from unittest import mock

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis import swept_langmuir as _sl
from plasmapy.analysis.swept_langmuir.floating_potential import (
    find_floating_potential,
    FloatingPotentialResults,
)


class TestFindFloatingPotential:
    """
    Tests for function
    `~plasmapy.analysis.swept_langmuir.floating_potential.find_floating_potential`.
    """

    def test_returned_namedtuple(self):
        """Test structure of the namedtuple used to return computed data."""

        assert issubclass(FloatingPotentialResults, tuple)
        assert hasattr(FloatingPotentialResults, "_fields")
        assert FloatingPotentialResults._fields == (
            "vf",
            "vf_err",
            "rsq",
            "func",
            "islands",
            "indices",
        )
        assert hasattr(FloatingPotentialResults, "_fields_defaults")
        assert FloatingPotentialResults._field_defaults == {}

    def test_call_of_check_sweep(self):
        """
        Test `find_floating_potential` appropriately calls
        `plasmapy.analysis.swept_langmuir.helpers.check_sweep` so we can relay on
        the `check_sweep` tests.
        """
        varr = np.linspace(-20.0, 20.0, 100)
        carr = np.linspace(-20.0, 20.0, 100)

        assert _sl.helpers.check_sweep is _sl.floating_potential.check_sweep

        with mock.patch(_sl.floating_potential.__name__ + ".check_sweep") as mock_cs:
            find_floating_potential(voltage=varr, current=carr, fit_type="linear")

            assert mock_cs.call_count == 1
            assert mock_cs.call_args.kwargs == {}
            assert len(mock_cs.call_args.args) == 2
            assert np.array_equal(mock_cs.call_args.args[0], varr)
            assert np.array_equal(mock_cs.call_args.args[1], carr)
