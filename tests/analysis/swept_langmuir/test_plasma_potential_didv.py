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
        (operator.eq, dIdVExtras._field_defaults, dict()),
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

