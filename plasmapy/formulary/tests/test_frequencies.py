"""Tests for functionality contained in `plasmapy.formulary.frequencies`."""
import astropy.units as u
import numpy as np
import pytest

from plasmapy.formulary.frequencies import gyrofrequency, oc_, wc_
from plasmapy.particles.exceptions import InvalidParticleError
from plasmapy.utils.pytest_helpers import assert_can_handle_nparray

ion = "p"

B = 1.0 * u.T
B_nanarr = np.array([0.001, np.nan]) * u.T


@pytest.mark.parametrize(
    "alias, parent",
    [
        (oc_, gyrofrequency),
        (wc_, gyrofrequency),
    ],
)
def test_parameters_aliases(alias, parent):
    """Test all aliases defined in parameters.py"""
    assert alias is parent


def test_gyrofrequency():
    r"""Test the gyrofrequency function in parameters.py."""

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
