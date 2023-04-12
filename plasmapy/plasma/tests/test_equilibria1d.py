from astropy import units as u

from plasmapy.plasma.equilibria1d import HarrisSheet


def test_HarrisSheet():
    B0 = 1 * u.T
    delta = 1 * u.m
    P0 = 0 * u.Pa
    hs = HarrisSheet(B0, delta, P0)
    B = hs.magnetic_field(0 * u.m)
    assert u.isclose(
        B, 0 * u.T, atol=1e-9 * u.T
    ), "Magnetic field is supposed to be zero at Y=0"


def test_fails():
    assert False
