import inspect

import astropy.units as u
import numpy as np

from plasmapy.particles import ParticleLike, particle_input
from plasmapy.utils.decorators import validate_quantities
from plasmapy.utils.decorators.converter import angular_freq_to_hz


def test_to_hz() -> None:
    @angular_freq_to_hz
    def func():
        return 2 * np.pi * u.rad / u.s

    assert func().unit == (u.rad / u.s), (
        f"Unit expected is {u.rad / u.s} instead of {func().unit}"
    )
    assert np.isclose(func().value, 2 * np.pi), (
        f"Value expected is {2 * np.pi} instead of {func().value}"
    )
    assert func(to_hz=True).unit == u.Hz, (
        f"Unit expected is {u.Hz} instead of {func(to_hz=True).unit}"
    )
    assert func(to_hz=True).value == 1, (
        f"Value expected is 1 instead of {func(to_hz=True).value}"
    )


def test_to_hz_complicated_signature() -> None:
    """
    Test that `angular_freq_to_hz` can decorate a function with
    positional-only, positional, var-positional, keyword, keyword-only,
    or var-keyword.
    """

    @angular_freq_to_hz
    def func2(a, /, b, *args, c, d: int = 2, **kwargs):
        return 2 * np.pi * u.rad / u.s

    result_rad_per_s = func2(1, 2, 3, 4, c=5, d=6, e=7)
    result_hz = func2(1, 2, 3, 4, c=5, d=6, e=7, to_hz=True)

    assert result_rad_per_s.unit == u.rad / u.s, (
        f"Unit expected is {(u.rad / u.s)} instead of {result_rad_per_s.unit}"
    )
    assert np.isclose(result_rad_per_s.value, 2 * np.pi), (
        f"Value expected is {2 * np.pi} instead of {result_rad_per_s.value}"
    )

    assert result_hz.unit == u.Hz, (
        f"Unit expected is {u.Hz} instead of {result_hz.unit}"
    )
    assert result_hz.value == 1, f"Value expected is 1 instead of {result_hz.value}"


def test_to_hz_stacked_decorators() -> None:
    """Test that @angular_freq_to_hz can be stacked with multiple decorators."""

    @particle_input
    @validate_quantities
    @angular_freq_to_hz
    def func(particle: ParticleLike | None = None):
        return 2 * np.pi * u.rad / u.s

    assert u.isclose(func(), 2 * np.pi * u.rad / u.s)
    assert u.isclose(func(to_hz=True), 1 * u.Hz)


def test_angular_freq_to_hz_preserves_signature() -> None:
    """
    Tests if the angular_freq_to_hz decorator preserves the signature of
    the wrapped function.
    """

    @angular_freq_to_hz
    def test_func(
        pos_only, /, arg, *args, required_kwarg, optional_kwarg: int = 2, **kwargs
    ):
        return 2 * u.rad / u.s

    def func_with_expected_signature(
        pos_only,
        /,
        arg,
        *args,
        required_kwarg,
        optional_kwarg: int = 2,
        to_hz: bool = False,
        **kwargs,
    ):
        return 2 * u.rad / u.s

    original_signature = inspect.signature(test_func)
    expected_signature = inspect.signature(func_with_expected_signature)

    assert original_signature == expected_signature, (
        f"Expected signature: {expected_signature}, but got: {original_signature}"
    )
