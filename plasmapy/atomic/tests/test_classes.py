import pytest
from ..classes import Particle


def test_Particle():
    r"""Original, temporary tests of Particle class."""
    p = Particle('H', mass_numb=1, Z=1)
    assert p._generation is None
