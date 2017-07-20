"""Tests for functions that calculate transport coefficients."""

import numpy as np
import pytest
from astropy import units as u

from ...constants import c, m_p, m_e, e, mu0

from ..transport import (Coulomb_logarithm)


def test_Coulomb_logarithm():
    Coulomb_logarithm()
