import numpy as np
import pytest

from hypothesis.extra.numpy import arrays
from hypothesis.strategies import floats, integers

from plasmapy.plasma.fluxsurface import FluxSurface


def Bt_func(R, Z):
    return 5.3


def FS(psi_value=-0.01, Bt_func=Bt_func):
    import plasmaboundaries

    from skimage import measure

    # plasma parameters
    params = plasmaboundaries.ITER

    # compute magnetic flux psi(R, z)
    psi = plasmaboundaries.compute_psi(params, config="non-null")

    # plot the results
    import matplotlib.pyplot as plt
    import numpy as np

    rmin, rmax = 0.6, 1.4
    zmin, zmax = -0.6, 0.6
    r = np.arange(rmin, rmax, step=0.01)
    z = np.arange(zmin, zmax, step=0.01)
    R, Z = np.meshgrid(r, z)
    PSI = psi(R, Z)  # compute magnetic flux
    import sympy

    Rsym, Zsym = sympy.symbols("R Z")
    psisym = psi(Rsym, Zsym, pkg="sp")
    Br = -psisym.diff(Zsym) / Rsym
    Bz = psisym.diff(Rsym) / Rsym
    B = sympy.sqrt(Br ** 2 + Bz ** 2)
    Bdiff_r = B.diff(Rsym)
    Bdiff_z = B.diff(Zsym)
    Brfunc = sympy.lambdify((Rsym, Zsym), Br)
    Bzfunc = sympy.lambdify((Rsym, Zsym), Bz)
    Brdifffunc = sympy.lambdify((Rsym, Zsym), Bdiff_r)
    Bzdifffunc = sympy.lambdify((Rsym, Zsym), Bdiff_z)
    contours = measure.find_contours(PSI, psi_value, positive_orientation="high")
    if len(contours) != 1:
        print(f"Could not find contour for psi = {psi_value} ({len(contours)=})")
    assert len(contours) == 1, (len(contours), psi_value, i)
    contour = contours[0]
    RcontourArrayUnits, ZcontourArrayUnits = contour[:, 1], contour[:, 0]

    Zcontour = ZcontourArrayUnits / PSI.shape[0] * (zmax - zmin) + zmin
    Rcontour = RcontourArrayUnits / PSI.shape[1] * (rmax - rmin) + rmin
    Brvals = Brfunc(Rcontour, Zcontour)
    Bzvals = Bzfunc(Rcontour, Zcontour)
    Bprimervals = Brdifffunc(Rcontour, Zcontour)
    Bprimezvals = Bzdifffunc(Rcontour, Zcontour)

    dZ = np.gradient(Zcontour)
    dR = np.gradient(Rcontour)
    dL = np.sqrt(dZ ** 2 + dR ** 2)
    plt.plot(Rcontour, Zcontour)
    Brvals = Brfunc(Rcontour, Zcontour)
    Bzvals = Bzfunc(Rcontour, Zcontour)
    Bprimervals = Brdifffunc(Rcontour, Zcontour)
    Bprimezvals = Bzdifffunc(Rcontour, Zcontour)
    fs = FluxSurface(
        Rcontour, Zcontour, psi_value, Brvals, Bzvals, Bt_func, Bprimervals, Bprimezvals
    )
    #     fs.flux_surface_average(fs.Bmag)
    return fs


def test_fs_Bmag(num_regression):
    num_regression.check({"fs-averaged mod B": fs.flux_surface_average(fs.Bmag)})


fs = FS()
from hypothesis import given, infer, settings


@given(
    array=arrays(
        float,
        fs.R.size,
        elements=floats(-1e60, 1e60, allow_nan=False, allow_infinity=False),
    )
)
def test_fs_flux_surface_average(array):
    average = fs.flux_surface_average(array)
    assert average.min() <= average <= average.max()


@given(m=integers(-100, 100))
def test_fs_flux_surface_sine_modes(m: int):
    mode = np.sin(fs.Theta * m)
    assert abs(fs.flux_surface_average(mode)) < 0.1


@given(m=integers(-100, 100).filter(lambda x: x != 0))
def test_fs_flux_surface_cosine_modes(m: int):
    mode = np.cos(fs.Theta * m)
    assert abs(fs.flux_surface_average(mode)) < 0.1


@pytest.mark.xfail(reason="This is only really true for periodic functions IIRC")
@given(
    array=arrays(
        float, (fs.R.size, 2), elements=floats(allow_nan=False, allow_infinity=False)
    )
)
def test_fs_flux_surface_average_annihilates_Bdot(array):
    under_average = np.array([B @ ai for B, ai in zip(fs.Bvectors.T, array)])
    assert abs(fs.flux_surface_average(under_average)) < 0.01


@pytest.mark.xfail(reason="Currently negative")
def test_fs_trapped_fraction(num_regression):
    f_t = fs.trapped_fraction()
    num_regression.check({"f_t": f_t})
    assert 0 < f_t < 0.5
