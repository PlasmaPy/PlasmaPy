"""Tests for functions that uses Distribution functions."""

import numpy as np
import pytest
from astropy import units as u
import scipy.integrate as spint


from ...constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, e)
from ...atomic import (ion_mass, charge_state)
from ..distribution import (Maxwellian_1D,
                            Maxwellian_speed_1D,
                            Maxwellian_velocity_3D,
                            Maxwellian_speed_3D)
from ..parameters import thermal_speed

T_e = 30000*u.K
V_drift = 1000000*u.m/u.s

start = -5000
stop = - start
dv = 10000 * u.m/u.s

v_vect = np.arange(start, stop, dtype='float64') * dv


def test_Maxwellian_1D():
    r"""test the 1D maxwellian distribution function"""

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

def test_Maxwellian_speed_1D():
    r"""test the 1D Maxwellian speed distribution function"""
    T = 1.0 * u.eV
    particle = 'H'
    # get thermal velocity and thermal velocity squared
    vTh = thermal_speed(T, particle=particle, method="most_probable")
    # setting up integration from 0 to 10*vTh
    xData1D = np.arange(0, 10.01, 0.01) * vTh
    yData1D = Maxwellian_speed_1D(v=xData1D,T=T,particle=particle)
    # integrating, this should be close to 1
    integ = spint.trapz(y=yData1D, x=xData1D)
    exceptStr = "Integral of distribution function should be 1."
    assert np.isclose(integ.value, 1), exceptStr
    
    
def test_Maxwellian_velocity_3D():
    r"""test the 3D Maxwellian velocity distribution function"""
    T = 1.0 * u.eV
    particle = 'H'
    # get thermal velocity and thermal velocity squared
    vTh = thermal_speed(T, particle=particle, method="most_probable")
    vTh = vTh.si.value
    # setting up integration from -10*vTh to 10*vTh, which is close to Inf
    infApprox = (10 * vTh)
    # integrating, this should be close to 1
    integ = spint.tplquad(Maxwellian_velocity_3D,
                          -infApprox,
                          infApprox,
                          lambda z: -infApprox,
                          lambda z: infApprox,
                          lambda z, y: -infApprox,
                          lambda z, y: infApprox,
                          args=(T,particle,0,0,0, vTh, "unitless"),
                          epsabs=1e0,
                          epsrel=1e0,
                          )
    # value returned from tplquad is (integral, error), we just need the 1st
    integVal = integ[0]
    exceptStr = ("Integral of distribution function should be 1 "
                 f"and not {integVal}.")
    assert np.isclose(integVal,
                      1,
                      rtol=1e-3,
                      atol=0.0), exceptStr
    
def test_Maxwellian_speed_3D():
    r"""test the 3D Maxwellian speed distribution function"""
    T = 1.0 * u.eV
    particle = 'H'
    # get thermal velocity and thermal velocity squared
    vTh = thermal_speed(T, particle=particle, method="most_probable")
    vTh = vTh.si.value
    # setting up integration from 0 to 10*vTh, which is close to Inf
    infApprox = (10 * vTh)
    # integral should be close to 1
    integ = spint.tplquad(Maxwellian_speed_3D,
                          0,
                          infApprox,
                          lambda z: 0,
                          lambda z: infApprox,
                          lambda z, y: 0,
                          lambda z, y: infApprox,
                          args=(T,particle,0,0,0, vTh, "unitless"),
                          epsabs=1e0,
                          epsrel=1e0,
                          )
    # value returned from tplquad is (integral, error), we just need the 1st
    integVal = integ[0]
    exceptStr = ("Integral of distribution function should be 1 "
                 f"and not {integVal}")
    assert np.isclose(integVal,
                      1,
                      rtol=1e-3,
                      atol=0.0), exceptStr