"""
Tests for the radiation-reaction methods of ``RelativisticBorisIntegrator`` in
``particle_integrators.py`` (issue #3310).

Two parametrized tests, each verifying a known limiting case:

1. ``E = 0`` with ``v`` perpendicular to ``B`` -- ``rrf_full`` reduces to the
   synchrotron-drag force proposed in PlasmaPy discussion #3306, checked across
   several Lorentz factors.
2. ``E = 0`` with ``v`` parallel to ``B`` -- nothing radiates, so
   ``push_including_rrf`` reduces to the plain relativistic ``push``, checked
   for several particles and speeds.

Following Nick Murphy's pytest tutorial (github.com/namurphy/pytest-tutorial),
cases are supplied with ``@pytest.mark.parametrize`` and floating-point results
are compared using ``numpy.allclose`` wrapped in ``assert``.
"""

import astropy.constants as const
import numpy as np
import pytest

from plasmapy.simulation.particle_integrators import RelativisticBorisIntegrator

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
        100.0,  # ultra-relativistic (checks the gamma**2 scaling at high gamma)
    ],
)
def test_rrf_full_reduces_to_discussion_synchrotron_drag(gamma) -> None:
    """
    With E = 0 and v perpendicular to B, the Landau-Lifshitz force must reduce
    to the synchrotron-drag model proposed in discussion #3306:
    ``f_R = -f0 * gamma * p_perp`` with ``f0 = mu0 q**4 B**2 / (6 pi m**3 c)``
    and ``p_perp = gamma m v``.
    """
    # --- Arrange: an electron moving perpendicular to B at this Lorentz factor
    q, m = -_e, _m_e
    speed = _c * np.sqrt(1 - 1 / gamma**2)  # speed corresponding to gamma
    B_magnitude = 1.0
    v = np.array([[speed, 0.0, 0.0]])  # velocity along x
    B = np.array([[0.0, 0.0, B_magnitude]])  # B along z  =>  v perpendicular to B
    E = np.zeros((1, 3))  # no electric field

    # --- Act: evaluate the radiation-reaction force
    f_R = RelativisticBorisIntegrator.rrf_full(v, B, E, q, m)

    # --- Assert: it equals the force proposed in discussion #3306
    f0 = _mu0 * q**4 * B_magnitude**2 / (6 * np.pi * m**3 * _c)
    p_perp = gamma * m * v
    expected = -f0 * gamma * p_perp
    # rtol is 1e-9 (not tighter): recomputing gamma from a near-light speed
    # loses ~1e-12 to floating-point cancellation at gamma = 100.
    assert np.allclose(f_R, expected, rtol=1e-9)


@pytest.mark.parametrize(
    "gamma, q, m",
    [
        (2.0, -_e, _m_e),  # relativistic electron
        (10.0, -_e, _m_e),  # ultra-relativistic electron
        (5.0, _e, _m_p),  # relativistic proton (different charge and mass)
    ],
)
def test_push_reduces_to_relativistic_boris_when_v_parallel_B(gamma, q, m) -> None:
    """
    With E = 0 and v parallel to B, v x B = 0 so the particle is unaccelerated
    and nothing radiates. ``push_including_rrf`` must therefore reproduce the
    plain relativistic ``push`` exactly, for any particle moving along B.
    """
    # --- Arrange: a particle moving ALONG B (so v x B = 0), no electric field
    speed = _c * np.sqrt(1 - 1 / gamma**2)
    x = np.array([[0.0, 0.0, 0.0]])
    v = np.array([[0.0, 0.0, speed]])  # velocity along z
    B = np.array([[0.0, 0.0, 2.0]])  # B along z  =>  v parallel to B
    E = np.zeros((1, 3))  # no electric field
    dt = 1.0e-11

    # --- Act: one step with radiation reaction, one plain relativistic step
    x_rr, v_rr = RelativisticBorisIntegrator.push_including_rrf(x, v, B, E, q, m, dt)
    x_boris, v_boris = RelativisticBorisIntegrator.push(x, v, B, E, q, m, dt)

    # --- Assert: the radiation-reaction kick is zero, so the two agree exactly
    assert np.allclose(x_rr, x_boris, rtol=1e-12)
    assert np.allclose(v_rr, v_boris, rtol=1e-12)
