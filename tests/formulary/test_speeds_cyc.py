import inspect
from collections.abc import Callable
from typing import cast

import astropy.units as u
import numpy as np

from plasmapy.formulary.speeds import ion_sound_speed
from plasmapy.particles.particle_class import Particle


def _raw_ion_sound_speed() -> Callable[..., u.Quantity]:
    """Get undecorated ion_sound_speed to avoid particle-input decorators."""
    return cast("Callable[..., u.Quantity]", inspect.unwrap(ion_sound_speed))


def test_ion_sound_speed_uses_Z_when_provided_unwrapped() -> None:
    """
    Regression test for the LOGIC bug: Z was ignored in the formula.
    We bypass decorators so Z isn't interpreted as particle charge override.
    """
    f = _raw_ion_sound_speed()

    T_e = 1e5 * u.K
    T_i = 1e5 * u.K
    ion = Particle("p+")  # already fully specified

    v_Z1 = f(T_e=T_e, T_i=T_i, ion=ion, Z=1.0)
    v_Z15 = f(T_e=T_e, T_i=T_i, ion=ion, Z=1.5)

    # Before fix: v_Z15 == v_Z1 (since Z ignored)
    # After fix: v_Z15 > v_Z1
    assert v_Z15 > v_Z1
    assert not u.isclose(v_Z15, v_Z1)

    # Optional: analytic ratio check when k=None and n_e=None (=> klD2 = 0)
    gamma_e = 1
    gamma_i = 3
    expected_ratio = np.sqrt(
        ((gamma_e * 1.5 * T_e + gamma_i * T_i) / (gamma_e * 1.0 * T_e + gamma_i * T_i))
        .decompose()
        .value
    )
    got_ratio = (v_Z15 / v_Z1).decompose().value
    assert np.isclose(got_ratio, expected_ratio, rtol=1e-12, atol=0)


def test_ion_sound_speed_Z_default_equals_charge_number_unwrapped() -> None:
    """
    Z=None should default to ion.charge_number (docstring behavior).
    """
    f = _raw_ion_sound_speed()

    T_e = 2e5 * u.K
    T_i = 5e4 * u.K
    ion = Particle("He 2+")

    v_default = f(T_e=T_e, T_i=T_i, ion=ion)  # Z=None implicit
    v_Z2 = f(T_e=T_e, T_i=T_i, ion=ion, Z=2)  # explicit

    assert u.isclose(v_default, v_Z2)

    v_Z1 = f(T_e=T_e, T_i=T_i, ion=ion, Z=1)
    assert not u.isclose(v_default, v_Z1)
