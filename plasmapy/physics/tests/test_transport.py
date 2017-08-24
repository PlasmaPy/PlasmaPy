"""Tests for functions that calculate transport coefficients."""


import numpy as np
import pytest
from astropy import units as u

from ...constants import c, m_p, m_e, e, mu0

from ..transport import (Coulomb_logarithm)


def test_Coulomb_logarithm():

    n_e = np.array([1e9, 1e9, 1e24])*u.cm**-3
    T = np.array([1e2, 1e7, 1e8])*u.K
    Lambda = np.array([5.97, 21.66, 6.69])
    particles = ('e', 'p')

    Coulomb_logarithm(5*u.m**-3, 1*u.eV, ('e', 'e'))

    for i in range(3):
        assert np.isclose(Coulomb_logarithm(n_e[i], T[i], particles),
                          Lambda[i], atol=0.01)

    assert np.isclose(Coulomb_logarithm(1e9*u.cm**-3, 1e2*u.K, ('e', 'p')),
                      5.97, atol=0.01)

    assert np.isclose(Coulomb_logarithm(1e9*u.cm**-3, 1e7*u.K, ('e', 'p')),
                      21.6, atol=0.1)

    assert np.isclose(Coulomb_logarithm(1e24*u.cm**-3, 1e8*u.K, ('e', 'p')),
                      6.69, atol=0.01)

    assert np.allclose(Coulomb_logarithm(n_e, T, particles), Lambda, atol=0.01)

    with pytest.raises(UserWarning):
        Coulomb_logarithm(1*u.m**-3, 1e5*u.K, ('e', 'p'), 299792458*u.m/u.s)

    with pytest.raises(u.UnitConversionError):
        Coulomb_logarithm(1*u.m**-3, 1e5*u.g, ('e', 'p'), 29979245*u.m/u.s)

    with pytest.raises(ValueError):
        Coulomb_logarithm(5*u.m**-3, 1*u.K, ('e'))

    with pytest.raises(ValueError):
        Coulomb_logarithm(5*u.m**-3, 1*u.K, ('e', 'g'))

    with pytest.raises(ValueError):
        Coulomb_logarithm(5*u.m**-3, 1*u.K, ('e', 'D'))
