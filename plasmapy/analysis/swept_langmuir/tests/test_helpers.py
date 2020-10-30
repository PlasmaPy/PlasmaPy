"""Tests for `plasmapy.analysis.swept_langmuir.helpers`."""

import numpy as np
import pytest

# ExitStack can be replaced with nullcontext when we require >= python 3.7
from contextlib import ExitStack as does_not_raise

from plasmapy.analysis.swept_langmuir.helpers import check_sweep


@pytest.mark.parametrize(
    "voltage, current, with_context",
    [
        # the one that works
        (np.linspace(-40.0, 40, 100), np.linspace(-10.0, 30, 100), does_not_raise()),
        # -- voltage cases --
        # not the right type
        ("not a numpy array", np.linspace(-10.0, 30, 100), pytest.raises(TypeError)),
        # not 1D
        (
            np.empty((2, 2), dtype=np.float64),
            np.linspace(-10.0, 30, 100),
            pytest.raises(ValueError),
        ),
        # not linearly increasing
        (
            np.linspace(40.0, -40, 100),
            np.linspace(-10.0, 30, 100),
            pytest.raises(ValueError),
        ),
        # -- current cases --
        # not the right type
        (np.linspace(-40.0, 40, 100), "not a numpy array", pytest.raises(TypeError)),
        # not 1D
        (
            np.linspace(-40.0, 40, 100),
            np.empty((2, 2), dtype=np.float64),
            pytest.raises(ValueError),
        ),
        # no floating potential (i.e. current never crosses zero)
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(10.0, 30, 100),
            pytest.raises(ValueError),
        ),
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(-30.0, -5, 100),
            pytest.raises(ValueError),
        ),
        # current needs to start from negative and go positive
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(30.0, -5, 100),
            pytest.raises(ValueError),
        ),
        # -- mixed cases --
        # voltage and current must have the same size
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(-5.0, 30, 150),
            pytest.raises(ValueError),
        ),
    ],
)
def test_check_sweep(voltage, current, with_context):
    """Test functionality of `plasmapy.analysis.swept_langmuir.helpers.check_sweep`."""
    with with_context:
        check_sweep(voltage=voltage, current=current)
