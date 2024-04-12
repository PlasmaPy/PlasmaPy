import astropy.units as u

from plasmapy.plasma import ForceFreeFluxRope


def test_B_z() -> None:
    B0 = 2 * u.T
    alpha = 3 * (1 / u.m)
    r = [0, 5] * u.m
    flux_rope = ForceFreeFluxRope(B0=B0, alpha=alpha)
    B_z = flux_rope.B_z(r=r)
    assert u.allclose(B_z, [B0, -0.028448945653561195 * u.T], atol=1e-9 * u.T)


def test_B_theta() -> None:
    B0 = 2 * u.T
    alpha = 3 * (1 / u.m)
    r = [0, 5] * u.m
    flux_rope = ForceFreeFluxRope(B0=B0, alpha=alpha)
    B_theta = flux_rope.B_theta(r=r)
    assert u.allclose(
        B_theta,
        [
            0 * u.T,
            0.41020807722704555 * u.T,
        ],
        atol=1e-9 * u.T,
    )


def test_B_magnitude() -> None:
    B0 = 2 * u.T
    alpha = 3 * (1 / u.m)
    r = [0, 5] * u.m
    flux_rope = ForceFreeFluxRope(B0=B0, alpha=alpha)
    B_magnitude = flux_rope.B_magnitude(r=r)
    assert u.allclose(B_magnitude, [B0, 0.4111933962639831 * u.T], atol=1e-9 * u.T)
