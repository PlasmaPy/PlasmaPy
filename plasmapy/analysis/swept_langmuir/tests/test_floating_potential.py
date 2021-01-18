"""
Tests for functionality contained in
`plasmapy.analysis.swept_langmuir.floating_potential`.
"""

import numpy as np
import pytest
import sys

from unittest import mock

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis import swept_langmuir as _sl
from plasmapy.analysis.swept_langmuir.floating_potential import (
    find_floating_potential,
    FloatingPotentialResults,
)


def test_floating_potential_namedtuple():
    """
    Test structure of the namedtuple used to return computed floating potential
    data.
    """

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
    if sys.version_info >= (3, 7):
        # TODO: remove if clause when python 3.6 support is dropped
        assert hasattr(FloatingPotentialResults, "_field_defaults")
        assert FloatingPotentialResults._field_defaults == {}


class TestFindFloatingPotential:
    """
    Tests for function
    `~plasmapy.analysis.swept_langmuir.floating_potential.find_floating_potential`.
    """

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

            # passed args
            assert len(mock_cs.call_args[0]) == 2
            assert np.array_equal(mock_cs.call_args[0][0], varr)
            assert np.array_equal(mock_cs.call_args[0][1], carr)

            # passed kwargs
            assert mock_cs.call_args[1] == {}

    @pytest.mark.parametrize(
        "kwargs, _error",
        [
            # errors on kwarg fit_type
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "fit_type": "wrong",
                },
                ValueError,
            ),
            #
            # errors on kwarg min_points
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "min_points": "wrong",
                },
                TypeError,
            ),
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "min_points": -1,
                },
                ValueError,
            ),
        ],
    )
    def test_raises(self, kwargs, _error):
        """Test scenarios that raise `Exception`s."""
        with pytest.raises(_error):
            find_floating_potential(**kwargs)
