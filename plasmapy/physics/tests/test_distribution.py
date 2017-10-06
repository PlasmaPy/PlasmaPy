"""Tests for functions that uses Distribution functions."""

import numpy as np
import pytest
from astropy import units as u


from ...constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, e)
from ...atomic import (ion_mass, charge_state)
from ..distribution import (Maxwellian_1D)

T_e = 30000*u.K
V_drift = 1000000*u.m/u.s

start = -5000
stop = - start
dv = 10000 * u.m/u.s

v_vect = np.arange(start, stop, dtype='float64') * dv


def test_Maxwellian_1D():
    """test the 1D maxwellian distribution function"""

    max_index = Maxwellian_1D(v_vect, T=T_e, particle='e', V_drift=0*u.m/u.s
                              ).argmax()
    assert np.isclose(v_vect[max_index].value, 0.0)

    max_index = Maxwellian_1D(v_vect, T=T_e, particle='e', V_drift=V_drift
                              ).argmax()
    assert np.isclose(v_vect[max_index].value, V_drift.value)

    # integral of the distribution over v_vect
    integral = (Maxwellian_1D(v_vect, T=30000*u.K, particle='e')).sum()*dv
    assert np.isclose(integral, 1.0)

    std = (Maxwellian_1D(v_vect, T=T_e, particle='e')*v_vect**2*dv).sum()
    std = np.sqrt(std)
    T_distri = (std**2/k_B*m_e).to(u.K)

    assert np.isclose(T_distri.value, T_e.value)

    with pytest.raises(ValueError):
        Maxwellian_1D(1*u.m/u.s, T=1*u.K, particle='XXX')
