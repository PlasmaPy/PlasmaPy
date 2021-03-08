import numpy as np
import pytest

from hypothesis.extra.numpy import arrays
from hypothesis.strategies import floats, integers

from plasmapy.plasma.fluxsurface import FluxSurface
from plasmapy.plasma.symbolicequilibrium import SymbolicEquilibrium


@pytest.fixture(scope="module")
def equilibrium():
    import plasmaboundaries

    params = plasmaboundaries.ITER.copy()
    equilibrium = SymbolicEquilibrium(**params, B0=5.2, config="single-null")
    return equilibrium


@pytest.fixture(scope="module")
def flux_surface(equilibrium, psi_value=-0.01):
    return equilibrium.get_flux_surface(psi_value)


def test_fs_Bmag(num_regression, flux_surface):
    num_regression.check(
        {"fs-averaged mod B": flux_surface.flux_surface_average(flux_surface.Bmag)}
    )


from hypothesis import given, infer, settings

# @given(
#     array=arrays(
#         float,
#         flux_surface.R.size,
#         elements=floats(-1e60, 1e60, allow_nan=False, allow_infinity=False),
#     )
# )
# def test_fs_flux_surface_average(array):
#     average = fs.flux_surface_average(array)
#     assert average.min() <= average <= average.max()


@given(m=integers(-100, 100).filter(lambda x: x != 0))
def test_fs_flux_surface_sine_modes(m: int, flux_surface):
    mode = np.sin(flux_surface.Theta * m)
    assert abs(flux_surface.flux_surface_average(mode)) < 0.1


@given(m=integers(1, 100))
def test_fs_flux_surface_cosine_modes(m: int, flux_surface):
    mode = np.cos(flux_surface.Theta * m)
    assert abs(flux_surface.flux_surface_average(mode)) < 0.1


# @pytest.mark.xfail(reason="This is only really true for periodic functions IIRC")
# @given(
#     array=arrays(
#         float, (flux_surface.R.size, 2), elements=floats(allow_nan=False, allow_infinity=False)
#     )
# )
# def test_fs_flux_surface_average_annihilates_Bdot(array, flux_surface):
#     under_average = np.array([B @ ai for B, ai in zip(flux_surface.Bvectors.T, array)])
#     assert abs(flux_surface.flux_surface_average(under_average)) < 0.01


@pytest.mark.xfail(reason="Currently negative")
def test_fs_trapped_fraction(num_regression, flux_surface):
    f_t = flux_surface.trapped_fraction()
    num_regression.check({"f_t": f_t})
    assert 0 < f_t < 0.5
