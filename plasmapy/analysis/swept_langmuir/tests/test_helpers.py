"""Tests for `plasmapy.analysis.swept_langmuir.helpers`."""

import numpy as np
import pytest

from contextlib import ExitStack as does_not_raise

from plasmapy.analysis.swept_langmuir.helpers import check_sweep


@pytest.mark.parametrize(
    "voltage, current, with_context",
    [
        (np.linspace(-40., 40., 100), np.linspace(-10., 30, 100), does_not_raise()),
        ("not a numpy array", np.linspace(-10., 30, 100), pytest.raises(TypeError)),
        # voltage not 1D
        (
            np.empty((2, 2), dtype=np.float64),
            np.linspace(-10., 30, 100),
            pytest.raises(ValueError),
        ),
        # voltage not linearly increasing
        (
            np.linspace(40., -40., 100),
            np.linspace(-10., 30, 100),
            pytest.raises(ValueError),
        ),
    ],
)
def test_check_sweep(voltage, current, with_context):
    """Test functionality of `plasmapy.analysis.swept_langmuir.helpers.check_sweep`."""
    with with_context:
        check_sweep(voltage=voltage, current=current)
