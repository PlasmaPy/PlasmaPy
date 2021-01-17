"""
Tests for functionality contained in
`plasmapy.analysis.swept_langmuir.floating_potential`.
"""

import numpy as np
import pytest

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis.swept_langmuir.floating_potential import (
    FloatingPotentialResults,
    find_floating_potential,
)


class TestFindFloatingPotential:
    """
    Tests for function
    `~plasmapy.analysis.swept_langmuir.floating_potential.find_floating_potential`.
    """

    def test_return_namedtuple(self):
        """Test structure of the namedtuple used to return computed data."""

        assert issubclass(FloatingPotentialResults, tuple)
        assert hasattr(FloatingPotentialResults, "_fields")
        assert (
            FloatingPotentialResults._fields
            == ("vf", "vf_err", "rsq", "func", "islands", "indices")
        )
        assert hasattr(FloatingPotentialResults, "_fields_defaults")
        assert FloatingPotentialResults._field_defaults == {}
