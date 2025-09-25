"""
Tests for functionality contained in
`plasmapy.analysis.swept_langmuir.plasma_potential` that calculates the
peak slope of the langmuir drives (i.e. peak dI/dV).
"""

from unittest import mock

import operator
import numpy as np
import pytest

from plasmapy.analysis import fit_functions as ffuncs
from plasmapy.analysis import swept_langmuir as sla
from plasmapy.analysis.swept_langmuir.plasma_potential import (
    dIdVExtras,
    find_didv_peak_location,
)
from plasmapy.utils.exceptions import PlasmaPyWarning


@pytest.mark.parametrize(
    ("assert_op", "arg1", "expected"),
    [
        (issubclass, dIdVExtras, tuple),
        (hasattr, dIdVExtras, "_fields"),
        (hasattr, dIdVExtras, "_field_defaults"),
        (
            operator.eq,
            set(dIdVExtras._fields),
            {"std", "data_slice", "savgol_windows", "savgol_peaks"},
        ),
        (operator.eq, dIdVExtras._field_defaults, {}),
    ],
)
def test_didv_namedtuple(assert_op, arg1, expected):
    assert assert_op(arg1, expected)


@pytest.mark.parametrize(
    ("index", "field_name"),
    [(0, "std"), (1, "data_slice"), (2, "savgol_windows"), (3, "savgol_peaks")],
)
def test_didv_namedtuple_index_field_mapping(index, field_name):
    extras = dIdVExtras(
        std=0.1,
        data_slice=slice(4, 20, 2),
        savgol_windows=[4, 6, 8],
        savgol_peaks=[1.1, 1.2, 1.3],
    )

    assert extras[index] == getattr(extras, field_name)


class TestFinddIdVPeakLocation:

    @pytest.mark.parametrize(
        ("helper", "used_callable"),
        [
            (
                sla.helpers._condition_voltage_window,
                sla.plasma_potential._condition_voltage_window,
            ),
            (sla.helpers.check_sweep, sla.plasma_potential.check_sweep),
            (
                sla.helpers.merge_voltage_clusters,
                sla.plasma_potential.merge_voltage_clusters,
            ),
        ],
    )
    def test_helper_callable_consistency(self, helper, used_callable):
        assert helper is used_callable

    def test_call_of_check_sweep(self) -> None:
        varr = np.linspace(-20.0, 20.0, 100)
        carr = np.linspace(-20.0, 20.0, 100)

        with mock.patch(f"{sla.plasma_potential.__name__}.check_sweep") as mock_cs:
            mock_cs.return_value = varr, carr
            find_didv_peak_location(voltage=varr, current=carr)

            assert mock_cs.call_count == 1

            # passed args
            assert len(mock_cs.call_args[0]) == 2
            assert np.array_equal(mock_cs.call_args[0][0], varr)
            assert np.array_equal(mock_cs.call_args[0][1], carr)

            # passed kwargs
            assert mock_cs.call_args[1] == {"strip_units": True}

    def test_call_of_merge_voltage_clusters(self) -> None:
        varr = np.linspace(-20.0, 20.0, 100)
        carr = np.linspace(-20.0, 20.0, 100)

        with (
            mock.patch(f"{sla.plasma_potential.__name__}.merge_voltage_clusters")
            as mock_mvc
        ):
            mock_mvc.return_value = varr, carr
            find_didv_peak_location(voltage=varr, current=carr)

            assert mock_mvc.call_count == 1

            # passed args
            assert len(mock_mvc.call_args[0]) == 2
            assert np.array_equal(mock_mvc.call_args[0][0], varr)
            assert np.array_equal(mock_mvc.call_args[0][1], carr)

            # passed kwargs
            assert mock_mvc.call_args[1] == {"voltage_step_size": 0}

    def test_call_of_condition_voltage_window(self) -> None:
        varr = np.linspace(-20.0, 20.0, 100)
        carr = np.linspace(-20.0, 20.0, 100)

        with (
            mock.patch(
                f"{sla.plasma_potential.__name__}._condition_voltage_window",
                side_effect=sla.plasma_potential._condition_voltage_window,
            )
            as mock_cvw
        ):
            # mock_cvw.return_value = varr, carr
            find_didv_peak_location(voltage=varr, current=carr)

            assert mock_cvw.call_count == 1

            # passed args
            assert len(mock_cvw.call_args[0]) == 2
            assert np.array_equal(mock_cvw.call_args[0][0], varr)
            assert mock_cvw.call_args[0][1] is None  # voltage_window arg

            # passed kwargs
            assert mock_cvw.call_args[1] == {}
