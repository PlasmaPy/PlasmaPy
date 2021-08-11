"""Tests for `plasmapy.analysis.swept_langmuir.helpers`."""

import astropy.units as u
import numpy as np
import pytest

from contextlib import nullcontext as does_not_raise

from plasmapy.analysis.swept_langmuir.helpers import check_sweep


@pytest.mark.parametrize(
    "voltage, current, kwargs, with_context, expected",
    [
        # the one that works
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(-10.0, 30, 100),
            {},
            does_not_raise(),
            "expected same as inputs",
        ),
        # -- voltage cases --
        # not the right type
        (
            "not a numpy array",
            np.linspace(-10.0, 30, 100),
            {},
            pytest.raises(TypeError),
            None,
        ),
        # not 1D
        (
            np.empty((2, 2), dtype=np.float64),
            np.linspace(-10.0, 30, 100),
            {},
            pytest.raises(ValueError),
            None,
        ),
        # not linearly increasing
        (
            np.linspace(40.0, -40, 100),
            np.linspace(-10.0, 30, 100),
            {},
            pytest.raises(ValueError),
            None,
        ),
        # list to array
        (
            [-5.0, -4, -3, -2, -1, 0, 1, 2, 3, 4],
            np.linspace(-10.0, 30, 10),
            {},
            does_not_raise(),
            (
                np.array([-5.0, -4, -3, -2, -1, 0, 1, 2, 3, 4]),
                np.linspace(-10.0, 30, 10),
            ),
        ),
        # wrong dtype
        (
            np.array(["one", "two", "three", "four", "five"]),
            np.linspace(-10.0, 30, 5),
            {},
            pytest.raises(ValueError),
            None,
        ),
        (
            ["one", "two", "three", "four", "five"],
            np.linspace(-10.0, 30, 5),
            {},
            pytest.raises(ValueError),
            None,
        ),
        # -- current cases --
        # not the right type
        (
            np.linspace(-40.0, 40, 100),
            "not a numpy array",
            {},
            pytest.raises(TypeError),
            None,
        ),
        # not 1D
        (
            np.linspace(-40.0, 40, 100),
            np.empty((2, 2), dtype=np.float64),
            {},
            pytest.raises(ValueError),
            None,
        ),
        # no floating potential (i.e. current never crosses zero)
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(10.0, 30, 100),
            {},
            pytest.raises(ValueError),
            None,
        ),
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(-30.0, -5, 100),
            {},
            pytest.raises(ValueError),
            None,
        ),
        # current needs to start from negative and go positive
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(30.0, -5, 100),
            {},
            pytest.raises(ValueError),
            None,
        ),
        # list to array
        (
            np.linspace(-40.0, 40, 10),
            [-5.0, -4, -3, -2, -1, 0, 1, 2, 3, 4],
            {},
            does_not_raise(),
            (
                np.linspace(-40.0, 40, 10),
                np.array([-5.0, -4, -3, -2, -1, 0, 1, 2, 3, 4]),
            ),
        ),
        # wrong dtype
        (
            np.linspace(-40.0, 40, 5),
            np.array(["one", "two", "three", "four", "five"]),
            {},
            pytest.raises(ValueError),
            None,
        ),
        (
            np.linspace(-40.0, 40, 5),
            ["one", "two", "three", "four", "five"],
            {},
            pytest.raises(ValueError),
            None,
        ),
        # -- mixed cases --
        # voltage and current must have the same size
        (
            np.linspace(-40.0, 40, 100),
            np.linspace(-5.0, 30, 150),
            {},
            pytest.raises(ValueError),
            None,
        ),
        # -- cases with units --
        # Note: u.Quantity is a subclass of np.ndarray, so all the tests above
        #       still work for an input argument that is a u.Quantity
        (
            np.linspace(-40.0, 40, 100) * u.volt,
            np.linspace(-10.0, 30, 100) * u.amp,
            {"strip_units": False},
            does_not_raise(),
            "expected same as inputs",
        ),
        (
            np.linspace(-40.0, 40, 100) * u.volt,
            np.linspace(-10.0, 30, 100) * u.amp,
            {},
            does_not_raise(),
            (np.linspace(-40.0, 40, 100), np.linspace(-10.0, 30, 100)),
        ),
    ],
)
def test_check_sweep(voltage, current, kwargs, with_context, expected):
    """Test functionality of `plasmapy.analysis.swept_langmuir.helpers.check_sweep`."""
    with with_context:
        rtn_voltage, rtn_current = check_sweep(
            voltage=voltage,
            current=current,
            **kwargs,
        )

        if expected == "expected same as inputs":
            assert np.allclose(rtn_voltage, voltage)
            assert np.allclose(rtn_current, current)
        else:
            assert np.allclose(rtn_voltage, expected[0])
            assert np.allclose(rtn_current, expected[1])
