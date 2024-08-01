"""Tests for functionality contained in `plasmapy.simulation.test_resolution_constraints`."""

import warnings

import astropy.units as u
import numpy as np
import pytest

from plasmapy.simulation.resolution_constraints import CFL_limit_electromagnetic_yee


class TestCFL_limit_electromagnetic_yee:
    """Tests for `plasmapy.simulation.resolution_constraints.CFL_limit_electromagnetic_yee`."""

    @pytest.mark.parametrize(
        ("args", "kwargs", "_error"),
        [
            ((1 * u.s,), {}, u.UnitTypeError),
            ((-1 * u.m,), {}, ValueError),
        ],
    )
    def test_raises(self, args, kwargs, _error) -> None:
        """Test scenarios that raise an exception."""

        with warnings.catch_warnings(), pytest.raises(_error):
            # we don't care about warnings for these tests
            warnings.simplefilter("ignore")
            CFL_limit_electromagnetic_yee(*args, **kwargs)

    @pytest.mark.parametrize(
        ("args", "kwargs", "expected", "atol"),
        [
            (
                (10 * u.nm,),
                {},
                [3.335640951981521e-17] * u.s,
                None,
            ),
            (
                (np.array([5, 10, 15]) * u.nm,),
                {},
                [1.4295604079920803e-17] * u.s,
                None,
            ),
        ],
    )
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_values(self, args, kwargs, expected, atol) -> None:
        if atol is None:
            atol = 1e-8

        rc = CFL_limit_electromagnetic_yee(*args, **kwargs)
        assert np.allclose(rc, expected, atol=atol)
