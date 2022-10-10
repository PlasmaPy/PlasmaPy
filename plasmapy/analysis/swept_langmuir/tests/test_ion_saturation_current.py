"""
Tests for functionality contained in
`plasmapy.analysis.swept_langmuir.ion_saturation_current`.
"""

import numpy as np
import pytest

from pathlib import Path

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis.swept_langmuir.ion_saturation_current import (
    find_ion_saturation_current,
    find_isat_,
    ISatExtras,
)


def test_ion_saturation_current_namedtuple():
    """
    Test structure of the namedtuple used to return the computed ion saturation
    current data.
    """

    assert issubclass(ISatExtras, tuple)
    assert hasattr(ISatExtras, "_fields")
    assert ISatExtras._fields == (
        "rsq",
        "fitted_func",
        "fitted_indices",
    )
    assert hasattr(ISatExtras, "_field_defaults")
    assert ISatExtras._field_defaults == {}


class TestFindIonSaturationCurrent:
    """
    Tests for function
    `~plasmapy.analysis.swept_langmuir.ion_saturation_current.find_ion_saturation_current`.
    """

    analytical_funcs = {
        "linear": ffuncs.Linear(params=(0.0004, -0.014)),
        "exp_offset": ffuncs.ExponentialPlusOffset(params=(0.001, 0.1, -0.01)),
        "exp_linear": ffuncs.ExponentialPlusLinear(params=(0.001, 0.1, 0.00005, -0.01)),
    }
    analytical_data = {"voltage": np.linspace(-30.0, 35, 100)}
    analytical_data.update(
        {
            "current_linear": analytical_funcs["linear"](analytical_data["voltage"]),
            "current_exp_offset": analytical_funcs["exp_offset"](
                analytical_data["voltage"]
            ),
            "current_exp_linear": analytical_funcs["exp_linear"](
                analytical_data["voltage"]
            ),
        },
    )

    _null_result = (
        None,
        ISatExtras(rsq=None, fitted_func=None, fitted_indices=None)._asdict(),
    )
    _voltage = np.linspace(-10.0, 15, 70)
    _linear_current = np.linspace(-3.1, 4.1, 70)
    _linear_p_sine_current = _linear_current + 1.2 * np.sin(1.2 * _voltage)
    _exp_current = -1.3 + 2.2 * np.exp(_voltage)

    def test_alias(self):
        """Test the associated alias(es) is(are) defined correctly."""
        assert find_isat_ is find_ion_saturation_current

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
            # current_bound and voltage_bound specified
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "current_bound": 0.5,
                    "voltage_bound": 0.0,
                },
                ValueError,
            ),
            #
            # current_bound/voltage_bound not the right type
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "current_bound": "not a number",
                },
                TypeError,
            ),
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "voltage_bound": "not a number",
                },
                TypeError,
            ),
            #
            # arguments result in no collected indices to fit
            (
                {
                    "voltage": np.array([1.0, 2, 3, 4]),
                    "current": np.array([-1.0, 0, 1, 2]),
                    "voltage_bound": 0.0,
                },
                ValueError,
            ),
        ],
    )
    def test_raises(self, kwargs, _error):
        """Test scenarios that raise exceptions."""
        with pytest.raises(_error):
            find_ion_saturation_current(**kwargs)

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            # linear fit to linear analytical data
            (
                {
                    "voltage": analytical_data["voltage"],
                    "current": analytical_data["current_linear"],
                    "fit_type": "linear",
                },
                (
                    analytical_funcs["linear"],
                    ISatExtras(
                        fitted_func=analytical_funcs["linear"],
                        rsq=None,
                        fitted_indices=slice(0, 40),
                    ),
                ),
            ),
            (
                {
                    "voltage": analytical_data["voltage"],
                    "current": analytical_data["current_linear"],
                    "fit_type": "linear",
                    "current_bound": 0.2,
                },
                (
                    analytical_funcs["linear"],
                    ISatExtras(
                        fitted_func=analytical_funcs["linear"],
                        rsq=None,
                        fitted_indices=slice(0, 20),
                    ),
                ),
            ),
            (
                {
                    "voltage": analytical_data["voltage"],
                    "current": analytical_data["current_linear"],
                    "fit_type": "linear",
                    "voltage_bound": 0.0,
                },
                (
                    analytical_funcs["linear"],
                    ISatExtras(
                        fitted_func=analytical_funcs["linear"],
                        rsq=None,
                        fitted_indices=slice(0, 46),
                    ),
                ),
            ),
            # exponential plus offset fit to exponential plus offset analytical data
            (
                {
                    "voltage": analytical_data["voltage"],
                    "current": analytical_data["current_exp_offset"],
                    "fit_type": "exp_plus_offset",
                },
                (
                    ffuncs.Linear(
                        params=(0.0, analytical_funcs["exp_offset"].params.b)
                    ),
                    ISatExtras(
                        fitted_func=analytical_funcs["exp_offset"],
                        rsq=None,
                        fitted_indices=slice(0, 81),
                    ),
                ),
            ),
            (
                {
                    "voltage": analytical_data["voltage"],
                    "current": analytical_data["current_exp_offset"],
                    "fit_type": "exp_plus_offset",
                    "voltage_bound": 30,
                },
                (
                    ffuncs.Linear(
                        params=(0.0, analytical_funcs["exp_offset"].params.b)
                    ),
                    ISatExtras(
                        fitted_func=analytical_funcs["exp_offset"],
                        rsq=None,
                        fitted_indices=slice(0, 92),
                    ),
                ),
            ),
            # exponential plus linear fit to exponential plus linear analytical data
            (
                {
                    "voltage": analytical_data["voltage"],
                    "current": analytical_data["current_exp_linear"],
                    "fit_type": "exp_plus_linear",
                },
                (
                    ffuncs.Linear(
                        params=(
                            analytical_funcs["exp_linear"].params.m,
                            analytical_funcs["exp_linear"].params.b,
                        )
                    ),
                    ISatExtras(
                        fitted_func=analytical_funcs["exp_linear"],
                        rsq=None,
                        fitted_indices=slice(0, 79),
                    ),
                ),
            ),
            (
                {
                    "voltage": analytical_data["voltage"],
                    "current": analytical_data["current_exp_linear"],
                    "fit_type": "exp_plus_linear",
                    "voltage_bound": 30,
                },
                (
                    ffuncs.Linear(
                        params=(
                            analytical_funcs["exp_linear"].params.m,
                            analytical_funcs["exp_linear"].params.b,
                        )
                    ),
                    ISatExtras(
                        fitted_func=analytical_funcs["exp_linear"],
                        rsq=None,
                        fitted_indices=slice(0, 92),
                    ),
                ),
            ),
        ],
    )
    def test_analytical_fits(self, kwargs, expected):
        """Test functionality on analytical traces."""
        isat, extras = find_ion_saturation_current(**kwargs)

        # assertions on isat
        assert isinstance(isat, type(expected[0]))
        assert np.allclose(isat.params, expected[0].params)

        # assertions on extras
        assert isinstance(extras, ISatExtras)
        assert isinstance(extras.fitted_func, type(expected[1].fitted_func))
        assert np.allclose(extras.fitted_func.params, expected[1].fitted_func.params)
        assert np.isclose(extras.rsq, 1.0)
        assert extras.fitted_indices == expected[1].fitted_indices

    def test_on_pace_data(self):
        """
        Test functionality on D. Pace data.

        Data was obtained from: https://davidpace.com/example-of-langmuir-probe-analysis/
        """
        filepath = (Path(__file__).parent / "Pace2015.npy").resolve()
        voltage, current = np.load(filepath)

        isort = np.argsort(voltage)
        voltage = voltage[isort]
        current = current[isort]

        isat, extras = find_ion_saturation_current(
            voltage, current, fit_type="exp_plus_linear", current_bound=3.6
        )

        assert np.isclose(isat.params.m, 3.81079e-6)
        assert np.isclose(isat.params.b, 0.000110284)
        assert np.isclose(extras.rsq, 0.982, atol=0.001)
        assert np.isclose(np.min(isat(voltage)), -0.00014275)
