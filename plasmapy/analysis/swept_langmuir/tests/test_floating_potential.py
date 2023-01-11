"""
Tests for functionality contained in
`plasmapy.analysis.swept_langmuir.floating_potential`.
"""

import numpy as np
import pytest

from unittest import mock

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis import swept_langmuir as sla
from plasmapy.analysis.swept_langmuir.floating_potential import (
    find_floating_potential,
    find_vf_,
    VFExtras,
)
from plasmapy.utils.exceptions import PlasmaPyWarning


def test_floating_potential_namedtuple():
    """
    Test structure of the namedtuple used to return computed floating potential
    data.
    """

    assert issubclass(VFExtras, tuple)
    assert hasattr(VFExtras, "_fields")
    assert VFExtras._fields == (
        "vf_err",
        "rsq",
        "fitted_func",
        "islands",
        "fitted_indices",
    )
    assert hasattr(VFExtras, "_field_defaults")
    assert VFExtras._field_defaults == {}


class TestFindFloatingPotential:
    """
    Tests for function
    `~plasmapy.analysis.swept_langmuir.floating_potential.find_floating_potential`.
    """

    _null_result = {
        **VFExtras(
            vf_err=np.nan, rsq=None, fitted_func=None, islands=None, fitted_indices=None
        )._asdict(),
        "vf": np.nan,
    }

    _voltage = np.linspace(-10.0, 15, 70)
    _linear_current = np.linspace(-3.1, 4.1, 70)
    _linear_p_sine_current = _linear_current + 1.2 * np.sin(1.2 * _voltage)
    _exp_current = -1.3 + 2.2 * np.exp(_voltage)

    def test_alias(self):
        """Test the associated alias(es) is(are) defined correctly."""
        assert find_vf_ is find_floating_potential

    def test_call_of_check_sweep(self):
        """
        Test `find_floating_potential` appropriately calls
        `plasmapy.analysis.swept_langmuir.helpers.check_sweep` so we can relay on
        the `check_sweep` tests.
        """
        varr = np.linspace(-20.0, 20.0, 100)
        carr = np.linspace(-20.0, 20.0, 100)

        assert sla.helpers.check_sweep is sla.floating_potential.check_sweep

        with mock.patch(f"{sla.floating_potential.__name__}.check_sweep") as mock_cs:
            mock_cs.return_value = varr, carr
            find_floating_potential(voltage=varr, current=carr, fit_type="linear")

            assert mock_cs.call_count == 1

            # passed args
            assert len(mock_cs.call_args[0]) == 2
            assert np.array_equal(mock_cs.call_args[0][0], varr)
            assert np.array_equal(mock_cs.call_args[0][1], carr)

            # passed kwargs
            assert mock_cs.call_args[1] == {"strip_units": True}

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
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "min_points": 0,
                },
                ValueError,
            ),
            #
            # errors on kwarg threshold
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "threshold": -1,
                },
                ValueError,
            ),
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "threshold": "wrong type",
                },
                TypeError,
            ),
            #
            # TypeError on voltage/current arrays from check_sweep
            (
                {
                    "voltage": "not an array",
                    "current": np.array([-1.0, 0, 1, 2]),
                    "fit_type": "linear",
                },
                TypeError,
            ),
            #
            # ValueError on voltage/current arrays from check_sweep
            #   (not linearly increasing)
            (
                {
                    "voltage": np.array([2.0, 1, 0, -1]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "fit_type": "linear",
                },
                ValueError,
            ),
        ],
    )
    def test_raises(self, kwargs, _error):
        """Test scenarios that raise `Exception`s."""
        with pytest.raises(_error):
            find_floating_potential(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, expected, _warning",
        [
            # too many crossing islands
            (
                {
                    "voltage": _voltage,
                    "current": _linear_p_sine_current,
                    "fit_type": "linear",
                },
                {
                    **_null_result,
                    "fitted_func": ffuncs.Linear(),
                    "islands": [slice(27, 29), slice(36, 38), slice(39, 41)],
                },
                PlasmaPyWarning,
            ),
            #
            # min_points is larger than array size
            (
                {
                    "voltage": _voltage,
                    "current": _linear_p_sine_current,
                    "fit_type": "linear",
                    "threshold": 8,
                    "min_points": 80,
                },
                {
                    **_null_result,
                    "vf": 0.6355491,
                    "vf_err": 0.03306472,
                    "rsq": 0.8446441,
                    "fitted_func": ffuncs.Linear(),
                    "islands": [slice(27, 41)],
                    "fitted_indices": slice(0, 70),
                },
                PlasmaPyWarning,
            ),
        ],
    )
    def test_warnings(self, kwargs, expected, _warning):
        """Test scenarios that issue warnings."""
        with pytest.warns(_warning):
            vf, extras = find_floating_potential(**kwargs)
            assert isinstance(extras, VFExtras)

        for key, val in expected.items():
            rtn_val = vf if key == "vf" else getattr(extras, key)
            if val is None:
                assert rtn_val is None
            elif key == "fitted_func":
                assert isinstance(rtn_val, val.__class__)
            elif np.isscalar(val):
                if np.isnan(val):
                    assert np.isnan(rtn_val)
                else:
                    assert np.isclose(rtn_val, val)
            else:
                assert rtn_val == val

    @pytest.mark.parametrize(
        "min_points, fit_type, islands, indices",
        [
            (np.inf, "linear", [slice(29, 31)], slice(0, 70)),
            (1, "linear", [slice(29, 31)], slice(29, 31)),
            (15, "linear", [slice(29, 31)], slice(22, 38)),
            (16, "linear", [slice(29, 31)], slice(22, 38)),
            (0.14, "linear", [slice(29, 31)], slice(25, 35)),
            #
            # rely on default min_points
            # - linear -> 0.1 * array size & ceiling to the nearest even
            # - exponential -> 0.2 * array size & ceiling to the nearest even
            (None, "linear", [slice(29, 31)], slice(26, 34)),
            (None, "exponential", [slice(26, 28)], slice(20, 34)),
        ],
    )
    def test_kwarg_min_points(self, min_points, fit_type, islands, indices):
        """
        Test functionality of keyword `min_points` and how it affects the
        size of the crossing-point island.
        """
        voltage = self._voltage
        current = self._linear_current if fit_type == "linear" else self._exp_current
        vf, extras = find_floating_potential(
            voltage,
            current,
            min_points=min_points,
            fit_type=fit_type,
        )
        assert isinstance(extras, VFExtras)

        assert extras.islands == islands
        assert extras.fitted_indices == indices

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            # simple linear
            (
                {
                    "voltage": _voltage,
                    "current": _linear_current,
                    "fit_type": "linear",
                    "min_points": 16,
                },
                {
                    **_null_result,
                    "vf": 0.7638889,
                    "vf_err": 0.0,
                    "rsq": 1.0,
                    "fitted_func": ffuncs.Linear(),
                    "islands": [slice(29, 31)],
                    "fitted_indices": slice(22, 38),
                },
            ),
            #
            # multiple islands merged with min_points
            (
                {
                    "voltage": _voltage,
                    "current": _linear_p_sine_current,
                    "fit_type": "linear",
                    "min_points": 16,
                },
                {
                    **_null_result,
                    "vf": -8.8243208,
                    "vf_err": 032.9961,
                    "rsq": 0.005084178,
                    "fitted_func": ffuncs.Linear(),
                    "islands": [slice(27, 29), slice(36, 38), slice(39, 41)],
                    "fitted_indices": slice(26, 42),
                },
            ),
            #
            # crossing-point near front of the array
            (
                {
                    "voltage": _voltage,
                    "current": _linear_current + 2.5,
                    "fit_type": "linear",
                    "min_points": 16,
                },
                {
                    **_null_result,
                    "vf": -7.91666667,
                    "vf_err": 3.153378e-8,
                    "rsq": 1.0,
                    "fitted_func": ffuncs.Linear(),
                    "islands": [slice(5, 7)],
                    "fitted_indices": slice(0, 16),
                },
            ),
            #
            # crossing-point near end of the array
            (
                {
                    "voltage": _voltage,
                    "current": _linear_current - 4.0,
                    "fit_type": "linear",
                    "min_points": 16,
                },
                {
                    **_null_result,
                    "vf": 14.6527778,
                    "vf_err": 0.0,
                    "rsq": 1.0,
                    "fitted_func": ffuncs.Linear(),
                    "islands": [slice(68, 70)],
                    "fitted_indices": slice(54, 70),
                },
            ),
        ],
    )
    def test_island_finding(self, kwargs, expected):
        """
        Test scenarios related to the identification of crossing-point islands.
        """
        vf, extras = find_floating_potential(**kwargs)
        assert isinstance(extras, VFExtras)

        for key, val in expected.items():
            rtn_val = vf if key == "vf" else getattr(extras, key)
            if val is None:
                assert rtn_val is None
            elif key == "fitted_func":
                assert isinstance(rtn_val, val.__class__)
            elif np.isscalar(val):
                if np.isnan(val):
                    assert np.isnan(rtn_val)
                else:
                    assert np.isclose(rtn_val, val, atol=1e-7)
            else:
                assert rtn_val == val

    @pytest.mark.parametrize("m, b", [(2.0, 0.0), (1.33, -0.1), (0.5, -0.1)])
    def test_perfect_linear(self, m, b):
        """Test calculated fit parameters on a few perfectly linear cases."""
        voltage = self._voltage
        current = m * voltage + b

        vf, extras = find_floating_potential(
            voltage=voltage,
            current=current,
            fit_type="linear",
            min_points=0.8,
        )

        assert isinstance(extras, VFExtras)
        assert np.isclose(vf, -b / m)
        assert np.isclose(extras.vf_err, 0.0)
        assert np.isclose(extras.rsq, 1.0)
        assert isinstance(extras.fitted_func, ffuncs.Linear)
        assert np.allclose(extras.fitted_func.params, (m, b))
        assert np.allclose(extras.fitted_func.param_errors, (0.0, 0.0), atol=2e-8)

    @pytest.mark.parametrize(
        "a, alpha, b",
        [(1.0, 0.2, -0.2), (2.7, 0.2, -10.0), (6.0, 0.6, -10.0)],
    )
    def test_perfect_exponential(self, a, alpha, b):
        """Test calculated fit parameters on a few perfectly exponential cases."""
        voltage = self._voltage
        current = a * np.exp(alpha * voltage) + b

        vf, extras = find_floating_potential(
            voltage=voltage,
            current=current,
            fit_type="exponential",
            min_points=0.8,
        )

        assert isinstance(extras, VFExtras)
        assert np.isclose(vf, np.log(-b / a) / alpha)
        assert np.isclose(extras.vf_err, 0.0, 1e-7)
        assert np.isclose(extras.rsq, 1.0)
        assert isinstance(extras.fitted_func, ffuncs.ExponentialPlusOffset)
        assert np.allclose(extras.fitted_func.params, (a, alpha, b))
        assert np.allclose(extras.fitted_func.param_errors, (0.0, 0.0, 0.0), atol=2e-8)
