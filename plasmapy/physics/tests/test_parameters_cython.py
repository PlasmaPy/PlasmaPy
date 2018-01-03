"""Tests for functions that calculate plasma parameters using cython."""

import numpy as np
import pytest
from astropy import units as u
from warnings import simplefilter

from ...utils.exceptions import RelativityWarning, RelativityError
from ...utils.exceptions import PhysicsError
from ...constants import c, m_p, m_e, e, mu0

from parameters_cython import (fib,
                               )


def test_fib():
    r"""Test fibonacci dummy function to check if Cython works."""
    fibMeth = fib(10)
    fibTrue = 55
    testTrue = fibMeth == fibTrue
    exceptStr = f'fib() value is {fibMeth}, should be {fibTrue}.'
    assert testTrue, exceptStr