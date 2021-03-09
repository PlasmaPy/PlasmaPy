"""
Tests for Thomson scattering analysis functions
"""

import astropy.units as u
import copy
import numpy as np
import pytest

from lmfit import Parameters

from plasmapy.diagnostics import thomson
from plasmapy.particles import Particle


def width_at_value(x, y, val):
    """
    Calculates the width of a curve at a given value.
    """
    above = np.where(y > val, x, np.nan)
    return np.abs(np.nanmax(above) - np.nanmin(above))


def gen_collective_spectrum():
    """
    Generates an example Thomson scattering spectrum in the collective regime
    """
    wavelengths = np.arange(520, 545, 0.01) * u.nm
    probe_wavelength = 532 * u.nm
    n = 5e17 * u.cm ** -3
    probe_vec = np.array([1, 0, 0])
    scatter_vec = np.array([0, 1, 0])
    Te = 10 * u.eV
    Ti = np.array([10]) * u.eV
    ion_species = ["C-12 5+"]

    alpha, Skw = thomson.spectral_density(
        wavelengths,
        probe_wavelength,
        n,
        Te,
        Ti,
        ion_species=ion_species,
        probe_vec=probe_vec,
        scatter_vec=scatter_vec,
    )

    return alpha, wavelengths, Skw


def gen_multiple_ion_species_spectrum():
    """
    Generates an example Thomson scattering spectrum for multiple ion species
    that also have drift velocities. Parameters are set to be in the
    collective regime where ion species are important.
    """
    wavelengths = np.arange(520, 545, 0.01) * u.nm
    probe_wavelength = 532 * u.nm
    n = 5e17 * u.cm ** -3
    probe_vec = np.array([1, 0, 0])
    scatter_vec = np.array([0, 1, 0])
    ifract = np.array([0.7, 0.3])
    Te = 10 * u.eV
    Ti = np.array([5, 5]) * u.eV
    electron_vel = np.array([[300, 0, 0]]) * u.km / u.s
    ion_vel = np.array([[-500, 0, 0], [0, 500, 0]]) * u.km / u.s

    # Use this to also test passing in ion species as Particle objects
    ion_species = [Particle("p+"), Particle("C-12 5+")]

    alpha, Skw = thomson.spectral_density(
        wavelengths,
        probe_wavelength,
        n,
        Te,
        Ti,
        ifract=ifract,
        ion_species=ion_species,
        probe_vec=probe_vec,
        scatter_vec=scatter_vec,
        electron_vel=electron_vel,
        ion_vel=ion_vel,
    )

    return alpha, wavelengths, Skw


def gen_non_collective_spectrum():
    """
    Generates an example Thomson scattering spectrum in the non-collective
    regime
    """
    wavelengths = np.arange(500, 570, 0.01) * u.nm
    probe_wavelength = 532 * u.nm
    n = 5e15 * u.cm ** -3
    probe_vec = np.array([1, 0, 0])
    scatter_vec = np.array([0, 1, 0])
    Te = 100 * u.eV
    Ti = np.array([10]) * u.eV
    ion_species = ["H+"]

    alpha, Skw = thomson.spectral_density(
        wavelengths,
        probe_wavelength,
        n,
        Te,
        Ti,
        ion_species=ion_species,
        probe_vec=probe_vec,
        scatter_vec=scatter_vec,
    )

    return alpha, wavelengths, Skw


def gen_collective_epw_spectrum():
    """
    Generates an example EPW-only spectrum
    """
    wavelengths = np.arange(520, 545, 0.01) * u.nm
    probe_wavelength = 532 * u.nm
    n = 5e17 * u.cm ** -3
    probe_vec = np.array([1, 0, 0])
    scatter_vec = np.array([0, 1, 0])
    Te = 10 * u.eV

    alpha, Skw = thomson.spectral_density_epw(
        wavelengths,
        probe_wavelength,
        n,
        Te,
        probe_vec=probe_vec,
        scatter_vec=scatter_vec,
    )

    return alpha, wavelengths, Skw


def test_different_input_types():

    # Define some constants
    wavelengths = np.arange(520, 545, 0.01) * u.nm
    probe_wavelength = 532 * u.nm
    n = 5e17 * u.cm ** -3
    probe_vec = np.array([1, 0, 0])
    scatter_vec = np.array([0, 1, 0])
    ifract = np.array([1.0])
    Te = np.array([10]) * u.eV
    Ti = np.array([10]) * u.eV
    ion_species = "C-12 5+"

    # Raise a ValueError with inconsistent ion array lengths
    with pytest.raises(ValueError):
        alpha, Skw = thomson.spectral_density(
            wavelengths,
            probe_wavelength,
            n,
            Te,
            Ti,
            ifract=np.array([0.5, 0.5]),
            ion_species=ion_species,
            probe_vec=probe_vec,
            scatter_vec=scatter_vec,
        )

    # Raise a ValueError with inconsistent ion temperature array
    with pytest.raises(ValueError):
        alpha, Skw = thomson.spectral_density(
            wavelengths,
            probe_wavelength,
            n,
            Te,
            np.array([5, 5]) * u.eV,
            ifract=ifract,
            ion_species=ion_species,
            probe_vec=probe_vec,
            scatter_vec=scatter_vec,
        )

    # Raise a ValueError with empty ion_species
    with pytest.raises(ValueError):
        alpha, Skw = thomson.spectral_density(
            wavelengths,
            probe_wavelength,
            n,
            Te,
            Ti,
            ifract=ifract,
            ion_species=[],
            probe_vec=probe_vec,
            scatter_vec=scatter_vec,
        )

    # Raise a Value Error with inconsistent electron array lengths
    # Te.size != efract.size
    with pytest.raises(ValueError):
        alpha, Skw = thomson.spectral_density(
            wavelengths,
            probe_wavelength,
            n,
            np.array([1, 10]) * u.eV,
            Ti,
            efract=np.array([0.5, 0.2, 0.3]),
            probe_vec=probe_vec,
            scatter_vec=scatter_vec,
        )

    # Electron vel shape not compatible with efract.size
    with pytest.raises(ValueError):
        alpha, Skw = thomson.spectral_density(
            wavelengths,
            probe_wavelength,
            n,
            Te,
            Ti,
            efract=np.array([0.5, 0.5]),
            electron_vel=np.array([[100, 0, 0]]) * u.km / u.s,
            probe_vec=probe_vec,
            scatter_vec=scatter_vec,
        )


def test_collective_spectrum():
    """
    Compares the generated spectrum to previously determined values
    """
    alpha, wavelength, Skw = gen_collective_spectrum()

    # Check that alpha is correct
    assert np.isclose(
        alpha.value, 1.801, atol=0.01
    ), "Collective case alpha returns {alpha} instead of expected 1.801"

    i_width = width_at_value(wavelength.value, Skw.value, 2e-13)
    e_width = width_at_value(wavelength.value, Skw.value, 0.2e-13)

    # Check that the widths of the ion and electron features match expectations
    assert np.isclose(i_width, 0.1599, 1e-3), (
        "Collective case ion feature "
        f"width is {i_width}"
        "instead of expected 0.1599"
    )

    assert np.isclose(e_width, 17.7899, 1e-3), (
        "Collective case electron "
        f"feature width is {e_width} "
        "instead of expected 17.7899"
    )


def test_non_collective_spectrum():
    """
    Compares the generated spectrum to previously determined values
    """
    alpha, wavelength, Skw = gen_non_collective_spectrum()

    # Check that alpha is correct
    assert np.isclose(
        alpha.value, 0.05707, atol=0.01
    ), "Non-collective case alpha returns {alpha} instead of expected 0.05707"

    e_width = width_at_value(wavelength.value, Skw.value, 0.2e-13)

    # Check that the widts of the electron feature matchs expectations
    assert np.isclose(e_width, 22.6699, 1e-3), (
        "Non-collective case electron "
        f"feature width is {e_width} "
        "instead of expected 22.6699"
    )


def test_multiple_ion_species_spectrum():
    """
    Compares the generated spectrum to previously determined values
    """

    alpha, wavelength, Skw = gen_multiple_ion_species_spectrum()

    # Compute the width and max of the spectrum, and the wavelength
    # of the max (sensitive to ion vel)
    width = width_at_value(wavelength.value, Skw.value, 0.2e-11)
    max_skw = np.max(Skw.value)
    max_wavelength = wavelength.value[np.argmax(Skw.value)]

    # Check width
    assert np.isclose(width, 0.049999, 1e-2), (
        f"Multiple ion species case spectrum width is {width} instead of "
        "expected 0.04999"
    )

    # Check max value
    assert np.isclose(max_skw, 2.4e-11, 1e-11), (
        f"Multiple ion species case spectrum max is {max_skw} instead of "
        "expected 2.4e-11"
    )

    # Check max peak location
    assert np.isclose(max_wavelength, 530.799, 1e-2), (
        "Multiple ion species case spectrum peak wavelength is "
        f"{max_wavelength} instead of expected 530.79"
    )


def test_split_populations():
    """
    This test makes sure that splitting a single population of ions or electrons
    into two identical halves returns the same result.
    """

    wavelengths = np.arange(520, 545, 0.01) * u.nm
    probe_wavelength = 532 * u.nm
    n = 5e17 * u.cm ** -3
    probe_vec = np.array([1, 0, 0])
    scatter_vec = np.array([0, 1, 0])

    # Combined
    Te = np.array([10]) * u.eV
    Ti = np.array([10]) * u.eV
    ion_species = ["H+"]
    ifract = np.array([1.0])
    efract = np.array([1.0])

    alpha, Skw0 = thomson.spectral_density(
        wavelengths,
        probe_wavelength,
        n,
        Te,
        Ti,
        ifract=ifract,
        efract=efract,
        ion_species=ion_species,
        probe_vec=probe_vec,
        scatter_vec=scatter_vec,
    )

    # Split e and i populations into two parts
    # this should not change the results since the parts are identical
    Te = np.array([10, 10]) * u.eV
    Ti = np.array([10, 10]) * u.eV
    ion_species = ["H+", "H+"]
    ifract = np.array([0.2, 0.8])
    efract = np.array([0.8, 0.2])

    alpha, Skw1 = thomson.spectral_density(
        wavelengths,
        probe_wavelength,
        n,
        Te,
        Ti,
        ifract=ifract,
        efract=efract,
        ion_species=ion_species,
        probe_vec=probe_vec,
        scatter_vec=scatter_vec,
    )

    # Calculate the deviation between the two spectra
    # (any differences should be in the noise)
    deviation = (Skw0 - Skw1) / Skw0 * 100

    assert np.all(deviation < 1e-6), "Failed split populations test"


def test_param_to_array_fcns():
    """
    Tests a few low-level routines used to convert lmfit scalar parameters
    into array input for `spectral_density` based on a naming convention
    """
    params = Parameters()

    # Create two groups of test variabels, one scalars and one vectors
    prefix = "Te"
    for i in range(3):
        params.add(prefix + f"_{i}", value=2)

    prefix = "ion_vel"
    for i in range(2):
        for j in ["x", "y", "z"]:
            params.add(prefix + f"_{j}_{i}", value=2)

    arr = thomson._params_to_array(params, "Te", vector=False)
    assert arr.shape == (3,)
    assert np.mean(arr) == 2

    arr = thomson._params_to_array(params, "ion_vel", vector=True)
    assert arr.shape == (2, 3)
    assert np.mean(arr) == 2


def run_fit(
    wavelengths,
    params,
    settings,
    noise_amp=0.05,
    notch=None,
    fit_method="differential_evolution",
    fit_kws={},
    max_iter=None,
    check_errors=True,
    require_redchi=1,
):
    """
    This function takes a Parameters object, generates some synthetic data near it,
    perturbs the initial values, then tries a fit

    """

    true_params = copy.deepcopy(params)

    skeys = list(settings.keys())
    pkeys = list(params.keys())

    # Fill any missing settings
    if "electron_vdir" not in skeys:
        settings["electron_vdir"] = None

    if "ion_vdir" not in skeys:
        settings["ion_vdir"] = None

    # Fill any missing required parameters
    if "efract_0" not in pkeys:
        params.add("efract_0", value=1.0, vary=False)

    if "ifract_0" not in pkeys:
        params.add("ifract_0", value=1.0, vary=False)

    if "electron_speed" not in pkeys:
        params.add("electron_speed_0", value=0.0, vary=False)

    if "ion_speed" not in pkeys:
        params.add("ion_speed_0", value=0.0, vary=False)

    # LOAD FROM SETTINGS
    ion_species = settings["ion_species"]
    probe_vec = settings["probe_vec"]
    scatter_vec = settings["scatter_vec"]
    electron_vdir = settings["electron_vdir"]
    ion_vdir = settings["ion_vdir"]
    probe_wavelength = settings["probe_wavelength"]

    # LOAD FROM PARAMS
    n = params["n"] * u.cm ** -3
    Te = thomson._params_to_array(params, "Te") * u.eV
    Ti = thomson._params_to_array(params, "Ti") * u.eV
    efract = thomson._params_to_array(params, "efract")
    ifract = thomson._params_to_array(params, "ifract")
    electron_speed = thomson._params_to_array(params, "electron_speed") * u.m / u.s
    ion_speed = thomson._params_to_array(params, "ion_speed") * u.m / u.s

    # Create the synthetic data
    alpha, Skw = thomson.spectral_density(
        wavelengths,
        probe_wavelength,
        n,
        Te,
        Ti,
        ifract=ifract,
        efract=efract,
        ion_species=ion_species,
        probe_vec=probe_vec,
        scatter_vec=scatter_vec,
        ion_vdir=ion_vdir,
        ion_speed=ion_speed,
        electron_vdir=electron_vdir,
        electron_speed=electron_speed,
    )

    data = Skw.value

    if notch is not None:
        x0 = np.argmin(np.abs(wavelengths.value - notch[0]))
        x1 = np.argmin(np.abs(wavelengths.value - notch[1]))
        data[x0:x1] = np.nan

    data *= 1 + np.random.normal(loc=0, scale=noise_amp, size=wavelengths.size)
    data *= 1 / np.nanmax(data)

    # Randomly choose the starting values of the parameters within the
    # search space (to make the algorithm do some work!)
    for p in list(params.keys()):
        if params[p].vary:
            params[p].value = np.random.uniform(
                low=params[p].min, high=params[p].max, size=1
            )

    # Make the model, then perform the fit
    model = thomson.thomson_model(wavelengths, settings, params)

    result = model.fit(
        data,
        params,
        wavelengths=wavelengths,
        method=fit_method,
        max_nfev=max_iter,
        fit_kws=fit_kws,
    )

    # Assert that the fit reduced chi2 is under the requirement specified
    assert result.redchi < require_redchi

    # Compute confidence intervals using the Monte-Carlo function and
    # check that the true values fall within the confidence intervals
    if check_errors:

        sigmas = thomson.mc_error(
            wavelengths,
            data,
            model,
            result,
            verbose=False,
            chisqr_cutoff=0.02,
            nsamples=300,
        )

        for k in list(sigmas.keys()):
            true = true_params[k].value
            floor = sigmas[k][0]
            ceiling = sigmas[k][1]

            assert true > floor
            assert true < ceiling


@pytest.mark.slow
def test_fitting():
    """
    This function sets up a bunch of fit Parameters then runs the run_fit
    function on them
    """

    # ***************************************************
    # Setup a number of different settings
    # ***************************************************

    probe_wavelength = 532 * u.nm
    probe_vec = np.array([1, 0, 0])
    scattering_angle = np.deg2rad(63)
    scatter_vec = np.array([np.cos(scattering_angle), np.sin(scattering_angle), 0])

    settings = {}
    settings["probe_wavelength"] = probe_wavelength
    settings["probe_vec"] = probe_vec
    settings["scatter_vec"] = scatter_vec
    settings["ion_species"] = ["H+"]
    settings["ion_vdir"] = np.array([[1, 0, 0]])
    settings["electron_vdir"] = np.array([[1, 0, 0]])
    single_species = copy.deepcopy(settings)

    settings = {}
    settings["probe_wavelength"] = probe_wavelength
    settings["probe_vec"] = probe_vec
    settings["scatter_vec"] = scatter_vec
    settings["ion_species"] = ["H+", "H+", "C-12 +4"]
    settings["ion_vdir"] = np.array([[0, 1, 0]])
    settings["electron_vdir"] = np.array([[0, 1, 0]])
    multi_species = copy.deepcopy(settings)

    # ***************************************************
    # Setup a number of different parameters
    # ***************************************************

    params = Parameters()
    params.add("n", value=2e17, vary=True, min=8e16, max=6e17)
    params.add("Te_0", value=10, vary=True, min=5, max=20)
    params.add("Ti_0", value=20, vary=False, min=5, max=70)
    epw_params = copy.deepcopy(params)

    params = Parameters()
    params.add("n", value=2e17, vary=False)
    params.add("Te_0", value=10, vary=False, min=5, max=20)
    params.add("Ti_0", value=20, vary=True, min=5, max=70)
    params.add("ifract_0", value=1.0, vary=False, min=0.2, max=0.8)
    params.add("ion_speed_0", value=0, vary=False)
    params.add("electron_speed_0", value=0, vary=False)
    iaw_params = copy.deepcopy(params)

    params = Parameters()
    params.add("n", value=1e19, vary=False)
    params.add("Te_0", value=500, vary=False, min=5, max=1000)
    params.add("Ti_0", value=200, vary=True, min=5, max=1000)
    params.add("Ti_1", value=500, vary=True, min=5, max=1000)
    params.add("Ti_2", value=400, vary=False, min=5, max=1000)
    params.add("ifract_0", value=0.4, vary=False, min=0.2, max=0.8)
    params.add("ifract_1", value=0.3, vary=False, min=0.2, max=0.8)
    params.add("ifract_2", value=0.3, vary=False, min=0.2, max=0.8)
    params.add("ion_speed_0", value=0, vary=False)
    params.add("ion_speed_1", value=1e5, vary=True, min=0, max=5e5)
    params.add("ion_speed_2", value=2e5, vary=False, min=0, max=5e5)
    params.add("electron_speed_0", value=0, vary=False)
    iaw_multispecies_params = copy.deepcopy(params)

    params = Parameters()
    params.add("n", value=2e17, vary=True, min=8e16, max=6e17)
    params.add("Te_0", value=10, vary=True, min=5, max=20)
    params.add("Ti_0", value=120, vary=False, min=5, max=70)
    params.add("efract_0", value=1.0, vary=False)
    params.add("electron_speed_0", value=0, vary=False)
    nc_params = copy.deepcopy(params)

    # ***************************************************
    # Setup a number of different wavelength ranges
    # ***************************************************

    w0 = probe_wavelength.value
    iaw_wavelengths = np.linspace(w0 - 5, w0 + 5, num=512) * u.nm
    epw_wavelengths = np.linspace(w0 - 40, w0 + 40, num=512) * u.nm
    nc_wavelengths = np.linspace(w0 - 60, w0 + 60, num=512) * u.nm

    # ***************************************************
    # Run a bunch of test fits
    # ***************************************************

    print("Running EPW test fit")
    run_fit(epw_wavelengths, epw_params, single_species, notch=(531, 533))

    print("Running IAW test fit")
    run_fit(iaw_wavelengths, iaw_params, single_species)

    print("Running IAW multi-species test fit")
    run_fit(iaw_wavelengths, iaw_multispecies_params, multi_species)

    print("Running Non-Collective test fit")
    run_fit(nc_wavelengths, nc_params, single_species)


if __name__ == "__main__":
    test_fitting()
    # test_collective_spectrum()
    # test_non_collective_spectrum()
    # test_different_input_types()
    # test_split_populations()
    # test_param_to_array_fcns()
    pass
