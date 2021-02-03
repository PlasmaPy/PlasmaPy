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
    IonSaturationCurrentResults,
)


def test_ion_saturation_current_namedtuple():
    """
    Test structure of the namedtuple used to return the computed ion saturation
    current data.
    """

    assert issubclass(IonSaturationCurrentResults, tuple)
    assert hasattr(IonSaturationCurrentResults, "_fields")
    assert IonSaturationCurrentResults._fields == (
        "isat_func",
        "rsq",
        "func",
        "indices",
    )
    assert hasattr(IonSaturationCurrentResults, "_field_defaults")
    assert IonSaturationCurrentResults._field_defaults == {}

