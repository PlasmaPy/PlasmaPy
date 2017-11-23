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

    for i in range(3):
        assert np.isclose(Coulomb_logarithm(T[i], n_e[i], particles),
                          Lambda[i], atol=0.01)

    assert np.isclose(Coulomb_logarithm(1*u.eV, 5*u.m**-3, ('e', 'e')),
                      Coulomb_logarithm(11604.5220*u.K, 5*u.m**-3, ('e', 'e')))

    assert np.isclose(Coulomb_logarithm(1e2*u.K, 1e9*u.cm**-3, ('e', 'p')),
                      5.97, atol=0.01)

    assert np.isclose(Coulomb_logarithm(1e7*u.K, 1e9*u.cm**-3, ('e', 'p')),
                      21.6, atol=0.1)

    assert np.isclose(Coulomb_logarithm(1e8*u.K, 1e24*u.cm**-3, ('e', 'p')),
                      6.69, atol=0.01)

    assert np.allclose(Coulomb_logarithm(T, n_e, particles), Lambda, atol=0.01)

    assert np.isclose(Coulomb_logarithm(1e5*u.K, 5*u.m**-3, ('e', 'e'),
                                        V=1e4*u.m/u.s), 21.379082011)

    with pytest.warns(UserWarning):
        Coulomb_logarithm(1e5*u.K, 1*u.m**-3, ('e', 'p'), 299792458*u.m/u.s)

    with pytest.raises(u.UnitConversionError):
        Coulomb_logarithm(1e5*u.g, 1*u.m**-3, ('e', 'p'), 29979245*u.m/u.s)

    with pytest.raises(ValueError):
        Coulomb_logarithm(1*u.K, 5*u.m**-3, ('e'))

    with pytest.raises(ValueError):
        Coulomb_logarithm(1*u.K, 5*u.m**-3, ('e', 'g'))

    with pytest.raises(ValueError):
        Coulomb_logarithm(1*u.K, 5*u.m**-3, ('e', 'D'))
