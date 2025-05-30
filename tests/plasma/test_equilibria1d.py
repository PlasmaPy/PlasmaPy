import astropy.units as u
import numpy as np

from plasmapy.formulary import magnetic_pressure
from plasmapy.plasma.equilibria1d import HarrisSheet


def test_HarrisSheet() -> None:
    B0 = 1 * u.T
    delta = 1 * u.m
    P0 = 0 * u.Pa
    hs = HarrisSheet(B0, delta, P0)
    B = hs.magnetic_field(0 * u.m)
    assert u.isclose(B, 0 * u.T, atol=1e-9 * u.T), (
        "Magnetic field is supposed to be zero at y=0"
    )


def test_HarrisSheet_pressure_balance() -> None:
    B0 = 1 * u.T
    delta = 1 * u.m
    P0 = 0 * u.Pa
    hs = HarrisSheet(B0, delta, P0)
    y = [-7, -3, 0, 2, 47] * u.m
    B = hs.magnetic_field(y)
    P = hs.plasma_pressure(y)
    p_b = magnetic_pressure(B)
    total_pressure = P + p_b
    assert u.allclose(total_pressure, total_pressure[0], atol=1e-9 * u.Pa)


def test_HarrisSheet_current_density() -> None:
    B0 = 1 * u.T
    delta = 1 * u.m
    P0 = 0 * u.Pa
    hs = HarrisSheet(B0, delta, P0)
    y = [-2, 0, 2] * u.m
    J = hs.current_density(y)
    correct_J = [-56222.1400445, -795774.715459, -56222.1400445] * u.A / u.m**2
    assert u.allclose(J, correct_J, atol=1e-8 * u.A / u.m**2)


def test_HarrisSheet_magnetic_field() -> None:
    B0 = 1 * u.T
    delta = 1 * u.m
    P0 = 0 * u.Pa
    hs = HarrisSheet(B0, delta, P0)
    y = [-2, 0, 2] * u.m
    B = hs.magnetic_field(y)
    correct_B = [-0.96402758007, 0, 0.96402758007] * u.T
    assert u.allclose(B, correct_B, atol=1e-9 * u.T)


def test_HarrisSheet_limits() -> None:
    y = [-np.inf, np.inf] * u.m
    B0 = 1 * u.T
    delta = 1 * u.m
    P0 = 0 * u.Pa
    hs = HarrisSheet(B0, delta, P0)
    B = hs.magnetic_field(y)
    P = hs.plasma_pressure(y)
    J = hs.current_density(y)
    assert u.allclose(B, [-B0, B0], atol=1e-9 * u.T)
    assert u.allclose(P, [P0, P0], atol=1e-9 * u.Pa)
    assert u.allclose(J, [0, 0] * u.amp / u.m**2, atol=1e-9 * u.amp / u.m**2)
