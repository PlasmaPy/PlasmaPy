import numpy as np
import pytest

from hypothesis.extra.numpy import arrays
from hypothesis.strategies import floats, integers


def test_fs_Bmag(num_regression, flux_surface):
    num_regression.check(
        {"fs-averaged mod B": flux_surface.flux_surface_average(flux_surface.Bmag)}
    )


import hypothesis

from hypothesis import given, infer, settings
from hypothesis.strategies import composite


@composite
def contour_shaped_array(draw, elements):
    def c_shaped_array(flux_surface):
        array = draw(arrays(float, flux_surface.R.shape, elements=elements))
        return array

    return c_shaped_array


# @given(
#     build=contour_shaped_array(
#         floats(-1e10, 1e10, allow_nan=False, allow_infinity=False)
#     )
# )
# def test_fs_flux_surface_average(build, flux_surface):
#     array = build(flux_surface)
#     hypothesis.note(f"array = {repr(array)}")  # extra reporting

#     average = flux_surface.flux_surface_average(array)
#     assert array.min() <= average
#     assert average <= array.max()


@given(m=integers(-100, 100).filter(lambda x: x != 0))
def test_fs_flux_surface_sine_modes(m: int, flux_surface):
    mode = np.sin(flux_surface.Theta * m)
    mode /= np.linalg.norm(mode, ord=2)
    assert abs(flux_surface.flux_surface_average(mode)) < 0.01


@given(m=integers(1, 100))
def test_fs_flux_surface_cosine_modes(m: int, flux_surface):
    mode = np.cos(flux_surface.Theta * m)
    mode /= np.linalg.norm(mode, ord=2)
    assert abs(flux_surface.flux_surface_average(mode)) < 0.01


# @pytest.mark.xfail(raises=AssertionError,
#                    reason="This is only really true for periodic functions IIRC")
# @given(
#     array=arrays(
#         float, (flux_surface.R.size, 2), elements=floats(allow_nan=False, allow_infinity=False)
#     )
# )
# def test_fs_flux_surface_average_annihilates_Bdot(array, flux_surface):
#     under_average = np.array([B @ ai for B, ai in zip(flux_surface.Bvectors.T, array)])
#     assert abs(flux_surface.flux_surface_average(under_average)) < 0.01


