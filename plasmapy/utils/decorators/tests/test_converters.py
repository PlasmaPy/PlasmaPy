import astropy.units as u
import numpy as np

from plasmapy.utils.decorators.converter import angular_freq_to_hz


def test_to_hz():
    @angular_freq_to_hz
    def func():
        return 2 * np.pi * u.rad / u.s

    assert func().unit == (
        u.rad / u.s
    ), f"Unit expected is {(u.rad / u.s)} instead of {func().unit}"
    assert (
        func().value == 2 * np.pi
    ), f"Value expected is {2 * np.pi} instead of {func().value}"
    assert (
        func(to_hz=True).unit == u.Hz
    ), f"Unit expected is {u.Hz} instead of {func(to_hz=True).unit}"
    assert (
        func(to_hz=True).value == 1
    ), f"Value expected is 1 instead of {func(to_hz=True).value}"


def test_to_hz_complicated_signature():
    """
    Test that `angular_freq_to_hz` can decorate a function with
    positional-only, postional, var-positional,
    """

    @angular_freq_to_hz
    def func2(a, /, b, *args, c, d=2, **kwargs):  # noqa: ARG001
        return 2 * np.pi * u.rad / u.s

    result_rad_per_s = func2(1, 2, 3, 4, c=5, d=6, e=7)
    result_hz = func2(1, 2, 3, 4, c=5, d=6, e=7, to_hz=True)

    assert result_rad_per_s.unit == (
        u.rad / u.s
    ), f"Unit expected is {(u.rad / u.s)} instead of {result_rad_per_s.unit}"
    assert (
        result_rad_per_s.value == 2 * np.pi
    ), f"Value expected is {2 * np.pi} instead of {result_rad_per_s.value}"

    assert (
        result_hz.unit == u.Hz
    ), f"Unit expected is {u.Hz} instead of {result_hz.unit}"
    assert result_hz.value == 1, f"Value expected is 1 instead of {result_hz.value}"
