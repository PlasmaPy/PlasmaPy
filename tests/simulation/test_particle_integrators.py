"""
Tests for the radiation-reaction methods of RelativisticBorisIntegrator in
particle_integrators.py (issue #3306).

Two parametrized tests, each verifying a known limiting case:

1. E = 0 with v perpendicular to B ==> rrf_full reduces to the
synchrotron-drag force proposed in PlasmaPy discussion #3306, checked across
several Lorentz factors.
2. E = 0 with v parallel to B ==> No change in acceleration ==>
nothing radiates ==> RelativisticBorisIntegratorRRF.push() reduces to 
the regular RelativisticBorisIntegrator.push() for different particles and speeds.
"""

import astropy.constants as const
import numpy as np
import pytest

from plasmapy.simulation.particle_integrators import (
    RelativisticBorisIntegrator,
    RelativisticBorisIntegratorRRF,
)

# Physical constants as raw SI floats (the integrator methods take SI floats,
# not astropy Quantity objects).
_c = const.c.si.value
_mu0 = const.mu0.si.value
_e = const.e.si.value
_m_e = const.m_e.si.value
_m_p = const.m_p.si.value


@pytest.mark.parametrize(
    "gamma",
    [
        2.0,  # mildly relativistic
        10.0,  # relativistic
        100.0,  # ultra-relativistic 
    ],
)
def test_rrf_full_reduces_to_discussion_synchrotron_drag(gamma) -> None:
    """
    With E = 0 and v perpendicular to B, the Landau-Lifshitz force must reduce
    to the synchrotron-drag model proposed in discussion #3306:
    """
    # establish an electron moving perpendicular to B at some Lorentz factor
    q, m = -_e, _m_e
    speed = _c * np.sqrt(1 - 1 / gamma**2)  # speed corresponding to gamma
    B_magnitude = 1.0
    v = np.array([[speed, 0.0, 0.0]])  # velocity along x
    B = np.array([[0.0, 0.0, B_magnitude]])  # B along z  =>  v perpendicular to B
    E = np.zeros((1, 3))  # no electric field

    # evaluate the radiation-reaction force from established conditions
    f_R = RelativisticBorisIntegratorRRF.rrf_full(v, B, E, q, m)

    # establish theoretically predicted eqn for RRF force proposed in discussion #3306
    f0 = _mu0 * q**4 * B_magnitude**2 / (6 * np.pi * m**3 * _c)
    p_perp = gamma * m * v
    expected = -f0 * gamma * p_perp
    assert np.allclose(f_R, expected, rtol=1e-12)


@pytest.mark.parametrize(
    "gamma,q,m",
    [
        (2.0, -_e, _m_e),  # relativistic electron
        (10.0, -_e, _m_e),  # ultra-relativistic electron
        (5.0, _e, _m_p),  # relativistic proton (different charge and mass)
    ],
)
def test_push_reduces_to_relativistic_boris_when_v_parallel_B(gamma, q, m) -> None:
    """
    With E = 0 and v parallel to B, v x B = 0 so the particle is unaccelerated.
    Should reproduce RelativisticBorisIntegrator.push() for range of particles 
    and velocities
    """
    # establish a particle moving ALONG B (so v x B = 0), and no electric field
    speed = _c * np.sqrt(1 - 1 / gamma**2)
    x = np.array([[0.0, 0.0, 0.0]])
    v = np.array([[0.0, 0.0, speed]])  # velocity along z
    B = np.array([[0.0, 0.0, 2.0]])  # B along z  =>  v parallel to B
    E = np.zeros((1, 3))  # no electric field
    dt = 1.0e-11

    # finding velocity and position push with rrf and without
    x_rr, v_rr = RelativisticBorisIntegratorRRF.push(x, v, B, E, q, m, dt)
    x_boris, v_boris = RelativisticBorisIntegrator.push(x, v, B, E, q, m, dt)

    # the radiation-reaction kick is zero, each particle pusher should equate
    assert np.allclose(x_rr, x_boris, rtol=1e-12)
    assert np.allclose(v_rr, v_boris, rtol=1e-12)
