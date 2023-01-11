"""Tests for functionality contained in `plasmapy.formulary.frequencies`."""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.frequencies import (
    gyrofrequency,
    lower_hybrid_frequency,
    oc_,
    plasma_frequency,
    upper_hybrid_frequency,
    wc_,
    wlh_,
    wuh_,
)
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

Z = 1
ion = "p"
n_i = 5e19 * u.m**-3
n_e = Z * 5e19 * u.m**-3

B = 1.0 * u.T
B_nanarr = np.array([0.001, np.nan]) * u.T


@pytest.mark.parametrize(
    "alias, parent",
    [
        (oc_, gyrofrequency),
        (wc_, gyrofrequency),
        (wlh_, lower_hybrid_frequency),
        (wuh_, upper_hybrid_frequency),
    ],
)
def test_aliases(alias, parent):
    """Test all aliases defined in frequencies.py"""
    assert alias is parent


def test_gyrofrequency():
    r"""Test the gyrofrequency function in frequencies.py."""

    assert gyrofrequency(B, "e-").unit.is_equivalent(u.rad / u.s)

    assert gyrofrequency(B, "e-", to_hz=True).unit.is_equivalent(u.Hz)

    assert np.isclose(gyrofrequency(1 * u.T, "e-").value, 175882008784.72018)

    assert np.isclose(gyrofrequency(2.4 * u.T, "e-").value, 422116821083.3284)

    assert np.isclose(
        gyrofrequency(1 * u.T, "e-", to_hz=True).value, 27992490076.528206
    )

    assert np.isclose(
        gyrofrequency(2.4 * u.T, "e-", signed=True).value, -422116821083.3284
    )

    assert np.isclose(gyrofrequency(1 * u.G, "e-").cgs.value, 1.76e7, rtol=1e-3)

    with pytest.raises(TypeError):
        with pytest.warns(u.UnitsWarning):
            gyrofrequency(u.m, "e-")

    with pytest.raises(u.UnitTypeError):
        gyrofrequency(u.m * 1, "e-")

    assert np.isnan(gyrofrequency(B_nanarr, "e-")[-1])

    # The following is a test to check that equivalencies from astropy
    # are working.
    omega_ce = gyrofrequency(2.2 * u.T, "e-")
    f_ce = (omega_ce / (2 * np.pi)) / u.rad
    f_ce_use_equiv = omega_ce.to(u.Hz, equivalencies=[(u.cy / u.s, u.Hz)])
    assert np.isclose(f_ce.value, f_ce_use_equiv.value)

    with pytest.warns(u.UnitsWarning):
        assert gyrofrequency(5.0, "e-") == gyrofrequency(5.0 * u.T, "e-")

    assert gyrofrequency(B, particle=ion).unit.is_equivalent(u.rad / u.s)

    assert np.isclose(gyrofrequency(1 * u.T, particle="p").value, 95788335.834874)

    assert np.isclose(gyrofrequency(2.4 * u.T, particle="p").value, 229892006.00369796)

    assert np.isclose(gyrofrequency(1 * u.G, particle="p").cgs.value, 9.58e3, rtol=2e-3)

    assert gyrofrequency(-5 * u.T, "p") == gyrofrequency(5 * u.T, "p")

    # Case when Z=1 is assumed
    # assert gyrofrequency(B, particle='p+') == gyrofrequency(B, particle='H-1')

    assert gyrofrequency(B, particle="e+") == gyrofrequency(B, "e-")

    with pytest.warns(u.UnitsWarning):
        gyrofrequency(8, "p")

    with pytest.raises(u.UnitTypeError):
        gyrofrequency(5 * u.m, "p")

    with pytest.raises(InvalidParticleError):
        gyrofrequency(8 * u.T, particle="asdfasd")

    with pytest.warns(u.UnitsWarning):
        # TODO this should be WARNS, not RAISES. and it's probably still raised
        assert gyrofrequency(5.0, "p") == gyrofrequency(5.0 * u.T, "p")

    gyrofrequency(1 * u.T, particle="p")
    # testing for user input Z
    testMeth1 = gyrofrequency(1 * u.T, particle="p", Z=0.8).si.value
    testTrue1 = 76630665.79318453
    errStr = f"gyrofrequency() gave {testMeth1}, should be {testTrue1}."
    assert np.isclose(testMeth1, testTrue1, atol=0.0, rtol=1e-5), errStr

    assert_can_handle_nparray(gyrofrequency, kwargs={"signed": True})

    assert_can_handle_nparray(gyrofrequency, kwargs={"signed": False})


def test_lower_hybrid_frequency():
    r"""Test the lower_hybrid_frequency function in frequencies.py."""

    ion = "He-4 1+"
    omega_ci = gyrofrequency(B, particle=ion)
    omega_pi = plasma_frequency(n=n_i, particle=ion)
    omega_ce = gyrofrequency(B, "e-")
    omega_lh = lower_hybrid_frequency(B, n_i=n_i, ion=ion)
    omega_lh_hz = lower_hybrid_frequency(B, n_i=n_i, ion=ion, to_hz=True)
    assert omega_ci.unit.is_equivalent(u.rad / u.s)
    assert omega_pi.unit.is_equivalent(u.rad / u.s)
    assert omega_ce.unit.is_equivalent(u.rad / u.s)
    assert omega_lh.unit.is_equivalent(u.rad / u.s)
    left_hand_side = omega_lh**-2
    right_hand_side = (
        1 / (omega_ci**2 + omega_pi**2) + omega_ci**-1 * omega_ce**-1
    )
    assert np.isclose(left_hand_side.value, right_hand_side.value)

    assert np.isclose(omega_lh_hz.value, 299878691.3223296)

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=5e19 * u.m**-3, ion="asdfasd")

    with pytest.raises(ValueError):
        lower_hybrid_frequency(0.2 * u.T, n_i=-5e19 * u.m**-3, ion="asdfasd")

    with pytest.raises(ValueError):
        lower_hybrid_frequency(np.nan * u.T, n_i=-5e19 * u.m**-3, ion="asdfasd")

    with pytest.warns(u.UnitsWarning):
        assert lower_hybrid_frequency(1.3, 1e19, "p+") == lower_hybrid_frequency(
            1.3 * u.T, 1e19 * u.m**-3, "p+"
        )
    assert_can_handle_nparray(lower_hybrid_frequency)


def test_upper_hybrid_frequency():
    r"""Test the upper_hybrid_frequency function in frequencies.py."""

    omega_uh = upper_hybrid_frequency(B, n_e=n_e)
    omega_uh_hz = upper_hybrid_frequency(B, n_e=n_e, to_hz=True)
    omega_ce = gyrofrequency(B, "e-")
    omega_pe = plasma_frequency(n=n_e, particle="e-")
    assert omega_ce.unit.is_equivalent(u.rad / u.s)
    assert omega_pe.unit.is_equivalent(u.rad / u.s)
    assert omega_uh.unit.is_equivalent(u.rad / u.s)
    assert omega_uh_hz.unit.is_equivalent(u.Hz)
    left_hand_side = omega_uh**2
    right_hand_side = omega_ce**2 + omega_pe**2
    assert np.isclose(left_hand_side.value, right_hand_side.value)

    assert np.isclose(omega_uh_hz.value, 69385868857.90918)

    with pytest.raises(ValueError):
        upper_hybrid_frequency(5 * u.T, n_e=-1 * u.m**-3)

    with pytest.warns(u.UnitsWarning):
        assert upper_hybrid_frequency(1.2, 1.3) == upper_hybrid_frequency(
            1.2 * u.T, 1.3 * u.m**-3
        )

    with pytest.warns(u.UnitsWarning):
        assert upper_hybrid_frequency(1.4 * u.T, 1.3) == upper_hybrid_frequency(
            1.4, 1.3 * u.m**-3
        )

    assert_can_handle_nparray(upper_hybrid_frequency)
