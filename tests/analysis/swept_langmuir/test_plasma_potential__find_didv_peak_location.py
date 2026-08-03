"""
Tests for functionality contained in
`plasmapy.analysis.swept_langmuir.plasma_potential` that calculates the
peak slope of the langmuir drives (i.e. peak dI/dV).
"""

import operator
from unittest import mock

import numpy as np
import pytest

from plasmapy.analysis import swept_langmuir as sla
from plasmapy.analysis.swept_langmuir.plasma_potential import (
    dIdVExtras,
    find_didv_peak_location,
)


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

    def test_call_of_check_sweep(self, simple_voltage, simple_current) -> None:
        with mock.patch(f"{sla.plasma_potential.__name__}.check_sweep") as mock_cs:
            mock_cs.return_value = simple_voltage, simple_current
            find_didv_peak_location(voltage=simple_voltage, current=simple_current)

            assert mock_cs.call_count == 1

            # passed args
            assert len(mock_cs.call_args[0]) == 2
            assert np.array_equal(mock_cs.call_args[0][0], simple_voltage)
            assert np.array_equal(mock_cs.call_args[0][1], simple_current)

            # passed kwargs
            assert mock_cs.call_args[1] == {"strip_units": True}

    def test_call_of_merge_voltage_clusters(
        self, simple_voltage, simple_current
    ) -> None:
        with mock.patch(
            f"{sla.plasma_potential.__name__}.merge_voltage_clusters"
        ) as mock_mvc:
            mock_mvc.return_value = simple_voltage, simple_current
            find_didv_peak_location(voltage=simple_voltage, current=simple_current)

            assert mock_mvc.call_count == 1

            # passed args
            assert len(mock_mvc.call_args[0]) == 2
            assert np.array_equal(mock_mvc.call_args[0][0], simple_voltage)
            assert np.array_equal(mock_mvc.call_args[0][1], simple_current)

            # passed kwargs
            assert mock_mvc.call_args[1] == {"voltage_step_size": 0}

    def test_call_of_condition_voltage_window(
        self, simple_voltage, simple_current
    ) -> None:
        with mock.patch(
            f"{sla.plasma_potential.__name__}._condition_voltage_window",
            side_effect=sla.plasma_potential._condition_voltage_window,
        ) as mock_cvw:
            # mock_cvw.return_value = varr, carr
            find_didv_peak_location(voltage=simple_voltage, current=simple_current)

            assert mock_cvw.call_count == 1

            # passed args
            assert len(mock_cvw.call_args[0]) == 2
            assert np.array_equal(mock_cvw.call_args[0][0], simple_voltage)
            assert mock_cvw.call_args[0][1] is None  # voltage_window arg

            # passed kwargs
            assert mock_cvw.call_args[1] == {}

    def test_call_of_condition_smooth_fractions(
        self, simple_voltage, simple_current
    ) -> None:
        with mock.patch(
            f"{sla.plasma_potential.__name__}._condition_smooth_fractions",
            side_effect=sla.plasma_potential._condition_smooth_fractions,
        ) as mock_csf:
            find_didv_peak_location(voltage=simple_voltage, current=simple_current)

            assert mock_csf.call_count == 1

            # passed args
            assert len(mock_csf.call_args[0]) == 2
            assert mock_csf.call_args[0][0] is None  # smooth_fractions arg
            assert mock_csf.call_args[0][1] == simple_voltage.size  # data_size arg

            # passed kwargs
            assert mock_csf.call_args[1] == {}

    @pytest.mark.usefixtures("request")
    @pytest.mark.parametrize(
        ("_raises", "voltage", "current", "voltage_window", "smooth_fractions"),
        [
            # voltage_window is too small (<3 points)
            (
                pytest.raises(ValueError),
                "simple_voltage",
                "simple_current",
                [4, 4.5],
                None,
            ),
            # no valid computer plasma potentials
            (
                pytest.raises(RuntimeError),
                "simple_voltage",
                np.append(
                    np.append(np.linspace(0, 10.0, 50), np.linspace(0, 10.0, 100)),
                    np.linspace(0, 10, 50),
                ),
                None,
                [0.0],
            ),
        ],
    )
    def test_raises(
        self, _raises, voltage, current, voltage_window, smooth_fractions, request
    ):
        if isinstance(voltage, str):
            # assume fixture name
            voltage = request.getfixturevalue(voltage)

        if isinstance(current, str):
            # assume fixture name
            current = request.getfixturevalue(current)

        with (
            mock.patch(
                f"{sla.plasma_potential.__name__}.check_sweep",
                return_value=(voltage, current),
            ),
            mock.patch(
                f"{sla.plasma_potential.__name__}.merge_voltage_clusters",
                return_value=(voltage, current),
            ),
            _raises,
        ):
            find_didv_peak_location(
                voltage=voltage,
                current=current,
                voltage_window=voltage_window,
                smooth_fractions=smooth_fractions,
            )

    @pytest.mark.usefixtures("request")
    @pytest.mark.parametrize(
        ("voltage", "current", "kwargs", "expected"),
        [
            (
                "simple_voltage",
                "simple_current_two_arches",
                {"smooth_fractions": [0.0]},  # no smoothing
                {"vp": 0.0},
            ),
        ],
    )
    def test_designed_didv(self, voltage, current, kwargs, expected, request):
        if isinstance(voltage, str):
            # assume fixture name
            voltage = request.getfixturevalue(voltage)

        if isinstance(current, str):
            # assume fixture name
            current = request.getfixturevalue(current)

        vp, extras = find_didv_peak_location(voltage, current, **kwargs)

        assert isinstance(vp, float)
        assert isinstance(extras, dIdVExtras)

        dv = float(np.average(np.diff(voltage)))
        vp_expected = expected["vp"]

        assert -0.5 * dv < (vp - vp_expected) < 0.5 * dv

    def test_on_pace_data(self, pace_voltage, pace_current):
        vp, extras = find_didv_peak_location(
            pace_voltage,
            pace_current,
            voltage_window=[-33.814, None],
        )

        # assertions on computed plasma potential
        assert isinstance(vp, float)
        assert np.isclose(vp, -16.6526)

        # assertions on extras from computation
        assert isinstance(extras, dIdVExtras)
        assert np.isclose(extras.std, 0.348483)
        assert extras.data_slice == slice(78, None, 1)
        assert extras.savgol_windows == [
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
        ]
        assert extras.savgol_peaks == [
            -16.4,
            -17.2,
            -16.8,
            -17.2,
            -17.2,
            -16.4,
            -16.4,
            -16.4,
            -16.0,
            -16.8,
            -16.8,
            -16.8,
            -16.8,
            -16.8,
            -16.8,
            -16.8,
            -16.4,
            -16.4,
            -16.0,
        ]

    @pytest.mark.parametrize(
        ("argmax_return", "expected"),
        [
            (np.array([10, 11, 12], dtype=np.int64), -17.78894472),
            (np.array([45], dtype=np.int64), -10.95477386),
            (np.int64(45), -10.95477386),
        ],
    )
    def test_peak_index_identification(
        self, argmax_return, expected, simple_voltage, simple_current
    ):
        with mock.patch("numpy.where") as mock_argmax:
            mock_argmax.return_value = (argmax_return,)

            vp, extras = find_didv_peak_location(
                simple_voltage,
                simple_current,
                voltage_window=None,
                smooth_fractions=[0.1],
            )

            assert np.isclose(vp, expected)

            mock_argmax.reset_mock()
