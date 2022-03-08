"""
Tests for functionality contained in
`plasmapy.analysis.swept_langmuir.ion_saturation_current`.
"""

import numpy as np
import pytest

from unittest import mock

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis import swept_langmuir as _sl
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
            # #
            # # errors on kwarg upper_bound
            # (
            #     {
            #         "voltage": np.array([1.0, 2, 3, 4]),
            #         "current": np.array([-1.0, 0, 1, 2]),
            #         "min_points": "wrong",
            #     },
            #     TypeError,
            # ),
            # (
            #     {
            #         "voltage": np.array([1.0, 2, 3, 4]),
            #         "current": np.array([-1.0, 0, 1, 2]),
            #         "min_points": -1,
            #     },
            #     ValueError,
            # ),
            # (
            #     {
            #         "voltage": np.array([1.0, 2, 3, 4]),
            #         "current": np.array([-1.0, 0, 1, 2]),
            #         "min_points": 0,
            #     },
            #     ValueError,
            # ),
            # #
            # # errors on kwarg threshold
            # (
            #     {
            #         "voltage": np.array([1.0, 2, 3, 4]),
            #         "current": np.array([-1.0, 0, 1, 2]),
            #         "threshold": -1,
            #     },
            #     ValueError,
            # ),
            # (
            #     {
            #         "voltage": np.array([1.0, 2, 3, 4]),
            #         "current": np.array([-1.0, 0, 1, 2]),
            #         "threshold": "wrong type",
            #     },
            #     TypeError,
            # ),
            # #
            # # TypeError on voltage/current arrays from check_sweep
            # (
            #     {
            #         "voltage": "not an array",
            #         "current": np.array([-1.0, 0, 1, 2]),
            #         "fit_type": "linear",
            #     },
            #     TypeError,
            # ),
            # #
            # # ValueError on voltage/current arrays from check_sweep
            # #   (not linearly increasing)
            # (
            #     {
            #         "voltage": np.array([2.0, 1, 0, -1]),
            #         "current": np.array([-1.0, 0, 1, 2]),
            #         "fit_type": "linear",
            #     },
            #     ValueError,
            # ),
        ],
    )
    def test_raises(self, kwargs, _error):
        """Test scenarios that raise `Exception`s."""
        with pytest.raises(_error):
            find_ion_saturation_current(**kwargs)
