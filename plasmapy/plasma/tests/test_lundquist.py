import astropy.units as u

from plasmapy.plasma import Lundquist


def test_B_z():
    B0 = 2 * u.T
    a = 3 * (1 / u.m)
    r = 0 * u.m
    flux_rope = Lundquist(B0=B0, a=a)
    B_axis = flux_rope.B_z(r=r)
    assert u.allclose(B_axis, B0, atol=1e-9 * u.T)
