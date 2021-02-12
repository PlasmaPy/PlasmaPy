"""
Defines the Thomson scattering analysis module as
part of the diagnostics package.
"""

__all__ = [
    "spectral_density",
    "thomson_model",
    "mc_error",
]

import astropy.constants as const
import astropy.units as u
import numpy as np
import warnings

from collections import namedtuple
from lmfit import Model
from typing import List, Tuple, Union

from plasmapy.formulary.dielectric import permittivity_1D_Maxwellian
from plasmapy.formulary.parameters import plasma_frequency, thermal_speed
from plasmapy.particles import Particle
from plasmapy.utils.decorators import validate_quantities

# Define some constants
C = const.c.si  # speed of light


# TODO: interface for inputting a multi-species configuration could be
# simplified using the plasmapy.classes.plasma_base class if that class
# included ion and electron drift velocities and information about the ion
# atomic species.


# TODO: If we can make this object pickle-able then we can set
# workers=-1 as a kw to differential_evolution to parallelize execution for fitting!
# The probem is a lambda function used in the Particle class...


@validate_quantities(
    wavelengths={"can_be_negative": False},
    probe_wavelength={"can_be_negative": False},
    n={"can_be_negative": False},
    Te={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    Ti={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def spectral_density(
    wavelengths: u.nm,
    probe_wavelength: u.nm,
    n: u.m ** -3,
    Te: u.K,
    Ti: u.K,
    efract: np.ndarray = None,
    ifract: np.ndarray = None,
    ion_species: Union[str, List[str], Particle, List[Particle]] = "H+",
    electron_vel: u.m / u.s = None,
    electron_vdir: np.ndarray = None,
    electron_speed: u.m / u.s = None,
    ion_vel: u.m / u.s = None,
    ion_vdir: np.ndarray = None,
    ion_speed: np.ndarray = None,
    probe_vec=np.array([1, 0, 0]),
    scatter_vec=np.array([0, 1, 0]),
) -> Tuple[Union[np.floating, np.ndarray], np.ndarray]:
    r"""
    Calculate the spectral density function for Thomson scattering of a
    probe laser beam by a multi-species Maxwellian plasma.

    This function calculates the spectral density function for Thomson
    scattering of a probe laser beam by a plasma consisting of one or more ion
    species and a one or more thermal electron populations (the entire plasma
    is assumed to be quasi-neutral)

    .. math::
        S(k,\omega) = \sum_e \frac{2\pi}{k}
        \bigg |1 - \frac{\chi_e}{\epsilon} \bigg |^2
        f_{e0,e} \bigg (\frac{\omega}{k} \bigg ) +
        \sum_i \frac{2\pi Z_i}{k}
        \bigg |\frac{\chi_e}{\epsilon} \bigg |^2 f_{i0,i}
        \bigg ( \frac{\omega}{k} \bigg )

    where :math:`\chi_e` is the electron component susceptibility of the
    plasma and :math:`\epsilon = 1 + \sum_e \chi_e + \sum_i \chi_i` is the total
    plasma dielectric  function (with :math:`\chi_i` being the ion component
    of the susceptibility), :math:`Z_i` is the charge of each ion, :math:`k`
    is the scattering wavenumber, :math:`\omega` is the scattering frequency,
    and :math:`f_{e0,e}` and :math:`f_{i0,i}` are the electron and ion velocity
    distribution functions respectively. In this function the electron and ion
    velocity distribution functions are assumed to be Maxwellian, making this
    function equivalent to Eq. 3.4.6 in `Sheffield`_.

    Parameters
    ----------

    wavelengths : `~astropy.units.Quantity`
        Array of wavelengths over which the spectral density function
        will be calculated. (convertible to nm)

    probe_wavelength : `~astropy.units.Quantity`
        Wavelength of the probe laser. (convertible to nm)

    n : `~astropy.units.Quantity`
        Mean (0th order) density of all plasma components combined.
        (convertible to cm^-3.)

    Te : `~astropy.units.Quantity`, shape (Ne, )
        Temperature of each electron component. Shape (Ne, ) must be equal to the
        number of electron components Ne. (in K or convertible to eV)

    Ti : `~astropy.units.Quantity`, shape (Ni, )
        Temperature of each ion component. Shape (Ni, ) must be equal to the
        number of ion components Ni. (in K or convertible to eV)

    efract : array_like, shape (Ne, ), optional
        An array-like object where each element represents the fraction (or ratio)
        of the electron component number density to the total electron number density.
        Must sum to 1.0. Default is a single electron component.

    ifract : array_like, shape (Ni, ), optional
        An array-like object where each element represents the fraction (or ratio)
        of the ion component number density to the total ion number density.
        Must sum to 1.0. Default is a single ion species.

    ion_species : str or `~plasmapy.particles.Particle`, shape (Ni, ), optional
        A list or single instance of `~plasmapy.particles.Particle`, or strings
        convertible to `~plasmapy.particles.Particle`. Default is `'H+'`
        corresponding to a single species of hydrogen ions.

    electron_vel : `~astropy.units.Quantity`, shape (Ne, 3), optional
        Velocity of each electron component in the rest frame. (convertible to m/s)
        If set, overrides electron_vdir and electron_speed.
        Defaults to a stationary plasma [0, 0, 0] m/s.

    electron_vdir : np.ndarray, shape (Ne,3), optional
        Unit vectors describing the velocity of each electron population.
        Setting electron_vel overrides this keyword.

    electron_speed : `~astropy.units.Quantity`, shape (Ne), optional
        A scalar speed for each electron population. Must be used along with
        electron_vdir. Setting electron_vel overrides this keyword.

    ion_vel : `~astropy.units.Quantity`, shape (Ni, 3), optional
        Velocity vectors for each electron population in the rest frame
        (convertible to m/s). If set, overrides ion_vdir and ion_speed.
        Defaults zero drift for all specified ion species.

    ion_vdir : np.ndarray, shape (Ne,3), optional
        Unit vectors describing the velocity of each ion population.
        Setting ion_vel overrides this keyword.

    ion_speed : `~astropy.units.Quantity`, shape (Ne), optional
        A scalar speed for each ion population. Must be used along with
        ion_vdir. Setting ion_vel overrides this keyword.

    probe_vec : float `~numpy.ndarray`, shape (3, )
        Unit vector in the direction of the probe laser. Defaults to
        [1, 0, 0].

    scatter_vec : float `~numpy.ndarray`, shape (3, )
        Unit vector pointing from the scattering volume to the detector.
        Defaults to [0, 1, 0] which, along with the default `probe_vec`,
        corresponds to a 90 degree scattering angle geometry.

    Returns
    -------
    alpha : float
        Mean scattering parameter, where `alpha` > 1 corresponds to collective
        scattering and `alpha` < 1 indicates non-collective scattering. The
        scattering parameter is calculated based on the total plasma density n.

    Skw : `~astropy.units.Quantity`
        Computed spectral density function over the input `wavelengths` array
        with units of s/rad.

    Notes
    -----

    For details, see "Plasma Scattering of Electromagnetic Radiation" by
    Sheffield et al. `ISBN 978\\-0123748775`_. This code is a modified version
    of the program described therein.

    For a concise summary of the relevant physics, see Chapter 5 of Derek
    Schaeffer's thesis, DOI: `10.5281/zenodo.3766933`_.

    .. _`ISBN 978\\-0123748775`: https://www.sciencedirect.com/book/9780123748775/plasma-scattering-of-electromagnetic-radiation
    .. _`10.5281/zenodo.3766933`: https://doi.org/10.5281/zenodo.3766933
    .. _`Sheffield`: https://doi.org/10.1016/B978-0-12-374877-5.00003-8
    """

    # Validate efract
    if efract is None:
        efract = np.ones(1)
    else:
        efract = np.asarray(efract, dtype=np.float64)

    # Validate ifract
    if ifract is None:
        ifract = np.ones(1)
    else:
        ifract = np.asarray(ifract, dtype=np.float64)

    # TODO: Write tests and update docstring for these different ways
    # of specifying velocities

    # Condition the electron velocity keywords
    if electron_vel is not None:
        pass
    elif (electron_speed is not None) and (electron_vdir is not None):
        norm = np.linalg.norm(electron_vdir, axis=1, keepdims=True)
        if all(norm == 0):
            electron_vel = np.zeros(3) * u.m / u.s
        else:
            electron_vel = np.outer(electron_speed, np.ones(3)) * electron_vdir
    else:
        electron_vel = np.zeros([efract.size, 3]) * u.m / u.s

    # Condition the electron velocity keywords
    if ion_vel is not None:
        pass
    elif (ion_speed is not None) and (ion_vdir is not None):
        norm = np.linalg.norm(ion_vdir, axis=1, keepdims=True)
        if all(norm == 0):
            ion_vel = np.zeros(3) * u.m / u.s
        else:
            ion_vel = np.outer(ion_speed, np.ones(3)) * ion_vdir
    else:
        ion_vel = np.zeros([ifract.size, 3]) * u.m / u.s

    # Condition ion_species
    if isinstance(ion_species, (str, Particle)):
        ion_species = [ion_species]
    if len(ion_species) == 0:
        raise ValueError("At least one ion species needs to be defined.")
    for ii, ion in enumerate(ion_species):
        if isinstance(ion, Particle):
            continue
        ion_species[ii] = Particle(ion)

    # Condition Ti
    if Ti.size == 1:
        # If a single quantity is given, put it in an array so it's iterable
        # If Ti.size != len(ion_species), assume same temp. for all species
        Ti = [Ti.value] * len(ion_species) * Ti.unit
    elif Ti.size != len(ion_species):
        raise ValueError(
            f"Got {Ti.size} ion temperatures and expected {len(ion_species)}."
        )

    # Make sure the sizes of ion_species, ifract, ion_vel, and Ti all match
    if (
        (len(ion_species) != ifract.size)
        or (ion_vel.shape[0] != ifract.size)
        or (Ti.size != ifract.size)
    ):
        raise ValueError(
            f"Inconsistent number of species in ifract ({ifract}), "
            f"ion_species ({len(ion_species)}), Ti ({Ti.size}), "
            f"and/or ion_vel ({ion_vel.shape[0]})."
        )

    # Condition Te
    if Te.size == 1:
        # If a single quantity is given, put it in an array so it's iterable
        # If Te.size != len(efract), assume same temp. for all species
        Te = [Te.value] * len(efract) * Te.unit
    elif Te.size != len(efract):
        raise ValueError(
            f"Got {Te.size} electron temperatures and expected {len(efract)}."
        )

    # Make sure the sizes of efract, electron_vel, and Te all match
    if (electron_vel.shape[0] != efract.size) or (Te.size != efract.size):
        raise ValueError(
            f"Inconsistent number of electron populations in efract ({efract.size}), "
            f"Te ({Te.size}), or electron velocity ({electron_vel.shape[0]})."
        )

    probe_vec = probe_vec / np.linalg.norm(probe_vec)
    scatter_vec = scatter_vec / np.linalg.norm(scatter_vec)
    scattering_angle = np.arccos(np.dot(probe_vec, scatter_vec))

    # Calculate plasma parameters
    vTe = thermal_speed(Te, particle="e-")
    vTi, ion_z = [], []
    for T, ion in zip(Ti, ion_species):
        vTi.append(thermal_speed(T, particle=ion).value)
        ion_z.append(ion.integer_charge * u.dimensionless_unscaled)
    vTi = vTi * vTe.unit
    zbar = np.sum(ifract * ion_z)

    # Compute electron and ion densities
    # np.nan_to_num and warning filter catch nan's generated when
    # efract or ifract == 0 (when calculating EPW or IAW spectra separately)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ne = np.nan_to_num(efract * n, nan=0)
        ni = np.nan_to_num(ifract * n / zbar, nan=0)  # ne/zbar = sum(ni)

    # wpe is calculated for the entire plasma (all electron populations combined)
    wpe = plasma_frequency(n=n, particle="e-")

    # Convert wavelengths to angular frequencies (electromagnetic waves, so
    # phase speed is c)
    ws = (2 * np.pi * u.rad * C / wavelengths).to(u.rad / u.s)
    wl = (2 * np.pi * u.rad * C / probe_wavelength).to(u.rad / u.s)

    # Compute the frequency shift (required by energy conservation)
    w = ws - wl

    # Compute the wavenumbers in the plasma
    # See Sheffield Sec. 1.8.1 and Eqs. 5.4.1 and 5.4.2
    ks = np.sqrt(ws ** 2 - wpe ** 2) / C
    kl = np.sqrt(wl ** 2 - wpe ** 2) / C

    # Compute the wavenumber shift (required by momentum conservation)\
    # Eq. 1.7.10 in Sheffield
    k = np.sqrt(ks ** 2 + kl ** 2 - 2 * ks * kl * np.cos(scattering_angle))
    # Normal vector along k
    k_vec = (scatter_vec - probe_vec) * u.dimensionless_unscaled

    # Compute Doppler-shifted frequencies for both the ions and electrons
    # Matmul is simultaneously conducting dot product over all wavelengths
    # and ion components
    w_e = w - np.matmul(electron_vel, np.outer(k, k_vec).T)
    w_i = w - np.matmul(ion_vel, np.outer(k, k_vec).T)

    # Compute the scattering parameter alpha
    # expressed here using the fact that v_th/w_p = root(2) * Debye length
    alpha = np.sqrt(2) * wpe / np.outer(k, vTe)

    # Calculate the normalized phase velocities (Sec. 3.4.2 in Sheffield)
    xe = (np.outer(1 / vTe, 1 / k) * w_e).to(u.dimensionless_unscaled)
    xi = (np.outer(1 / vTi, 1 / k) * w_i).to(u.dimensionless_unscaled)

    # Calculate the susceptibilities
    chiE = np.zeros([efract.size, w.size], dtype=np.complex128)
    for i, fract in enumerate(efract):
        chiE[i, :] = permittivity_1D_Maxwellian(w_e[i, :], k, Te[i], ne[i], "e-")

    # Treatment of multiple species is an extension of the discussion in
    # Sheffield Sec. 5.1
    chiI = np.zeros([ifract.size, w.size], dtype=np.complex128)
    for i, ion in enumerate(ion_species):
        chiI[i, :] = permittivity_1D_Maxwellian(
            w_i[i, :], k, Ti[i], ni[i], ion, z_mean=ion_z[i]
        )

    # Calculate the longitudinal dielectric function
    epsilon = 1 + np.sum(chiE, axis=0) + np.sum(chiI, axis=0)

    econtr = np.zeros([efract.size, w.size], dtype=np.complex128) * u.s / u.rad
    for m in range(efract.size):
        econtr[m, :] = efract[m] * (
            2
            * np.sqrt(np.pi)
            / k
            / vTe[m]
            * np.power(np.abs(1 - np.sum(chiE, axis=0) / epsilon), 2)
            * np.exp(-xe[m, :] ** 2)
        )

    icontr = np.zeros([ifract.size, w.size], dtype=np.complex128) * u.s / u.rad
    for m in range(ifract.size):
        icontr[m, :] = ifract[m] * (
            2
            * np.sqrt(np.pi)
            * ion_z[m]
            / k
            / vTi[m]
            * np.power(np.abs(np.sum(chiE, axis=0) / epsilon), 2)
            * np.exp(-xi[m, :] ** 2)
        )

    # Recast as real: imaginary part is already zero
    Skw = np.real(np.sum(econtr, axis=0) + np.sum(icontr, axis=0))

    return np.mean(alpha), Skw


# ***************************************************************************
# These functions are necessary to interface scalar Parameter objects with
# the array inputs of spectral_density
# ***************************************************************************


def _count_populations_in_params(keys, prefix):
    """
    Counts the number of entries matching the pattern prefix_i in a
    list of keys
    """

    i = 0
    while True:
        if prefix + f"_{i}" in keys:
            i += 1
        else:
            break
    return i


def _params_to_array(params, prefix, vector=False):
    """
    Takes a list of parameters and returns an array of the values corresponding
    to a key, based on the following naming convention:

    Each parameter should be named prefix_i
    Where i is an integer (starting at 0)

    This function allows lmfit.Parameter inputs to be converted into the
    array-type inputs required by the spectral density function

    """
    keys = list(params.keys())

    if vector:
        npop = _count_populations_in_params(keys, prefix + "_x")
        output = np.zeros([npop, 3])
        for i in range(npop):
            for j, ax in enumerate(["x", "y", "z"]):
                output[i, j] = params[prefix + f"_{ax}_{i}"].value

    else:
        npop = _count_populations_in_params(keys, prefix)
        output = np.zeros([npop])
        for i in range(npop):
            output[i] = params[prefix + f"_{i}"]

    return output


# ***************************************************************************
# Fitting functions
# ***************************************************************************


def _thomson_model(wavelengths, settings=None, **params):
    """
    lmfit Model function for fitting Thomson spectra

    For descriptions of arguments, see the `thomson_model` function.

    """

    # LOAD FROM SETTINGS
    ion_species = settings["ion_species"]
    probe_vec = settings["probe_vec"]
    scatter_vec = settings["scatter_vec"]
    electron_vdir = settings["electron_vdir"]
    ion_vdir = settings["ion_vdir"]
    probe_wavelength = settings["probe_wavelength"]

    # LOAD FROM PARAMS
    n = params["n"] * u.cm ** -3
    Te = _params_to_array(params, "Te") * u.eV
    Ti = _params_to_array(params, "Ti") * u.eV
    efract = _params_to_array(params, "efract")
    ifract = _params_to_array(params, "ifract")
    electron_speed = _params_to_array(params, "electron_speed") * u.m / u.s
    ion_speed = _params_to_array(params, "ion_speed") * u.m / u.s

    alpha, model_Skw = spectral_density(
        wavelengths,
        probe_wavelength,
        n,
        Te,
        Ti,
        ion_species=ion_species,
        efract=efract,
        ifract=ifract,
        electron_vdir=electron_vdir,
        electron_speed=electron_speed,
        ion_vdir=ion_vdir,
        ion_speed=ion_speed,
        probe_vec=probe_vec,
        scatter_vec=scatter_vec,
    )

    # Strip units and normalize
    model_Skw = model_Skw.to(u.s / u.rad).value
    model_Skw *= 1 / np.max(model_Skw)

    return model_Skw


def thomson_model(wavelengths, settings, params):
    """
    Returns a `lmfit.Model` function for Thomson spectrum


    Parameters
    ----------

    wavelength : 'astropy.units.Quantity' array
        Array of wavelengths over which to to evaluate the model


    settings : dict
        A dictionary of non-variable inputs to the spectral density function
        which must include the following:

            - probe_wavelength: Probe wavelength in nm
            - probe_vec : (3,) unit vector in the probe direction
            - scatter_vec: (3,) unit vector in the scattering direction
            - ion_species : list of Particle strings describing each ion species

        and may contain the following optional variables
            - electron_vdir : (e#, 3) array of electron velocity unit vectors
            - ion_vdir : (e#, 3) array of ion velocity unit vectors

        These quantities cannot be varied.


    params : `lmfit.Parameters` object
        A Parameters object that must contains the following variables
            - n: 0th order density in cm^-3
            - Te_e#
            - Ti_i#

        and may contain the following optional variables
            - efract_e# : Fraction of each electron population (must sum to 1) (optional)
            - ifract_i# : Fraction of each ion population (must sum to 1) (optional)
            - electron_speed_e# : Electron speed in m/s (optional)
            - ion_speed_i# : Ion speed in m/s (optional)

        where i# and e# are the number of electron and ion populations,
        zero-indexed, respectively (eg. 0,1,2...).

        These quantities can be either fixed or varying.


    Returns
    -------

    Spectral density (optimization function)


    """

    skeys = list(settings.keys())
    pkeys = list(params.keys())

    req_settings = ["probe_wavelength", "probe_vec", "scatter_vec", "ion_species"]
    for k in req_settings:
        if k not in skeys:
            raise KeyError(f"{k} was not provided in settings, but is required.")

    req_params = ["n", "Te_0", "Ti_0"]
    for k in req_params:
        if k not in pkeys:
            raise KeyError(f"{k} was not provided in parameters, but is required.")

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

    # Automatically add an expression to the last efract parameter to
    # indicate that it depends on the others (so they sum to 1.0)
    # The resulting expression for the last of three will look like
    # efract_2.expr = "1.0 - efract_0 - efract_1"
    num_e = _count_populations_in_params(pkeys, "efract")
    if num_e > 1:
        nums = ["efract_" + str(i) for i in range(num_e - 1)]
        nums.insert(0, "1.0")
        params["efract_" + str(num_e - 1)].expr = " - ".join(nums)

    num_i = _count_populations_in_params(pkeys, "ifract")
    if num_i > 1:
        nums = ["ifract_" + str(i) for i in range(num_i - 1)]
        nums.insert(0, "1.0")
        params["ifract_" + str(num_i - 1)].expr = " - ".join(nums)

    # Create a lmfit.Model
    # nan_policy='omit' automatically ignores NaN values in data, allowing those
    # to be used to represnt regions of missing data
    # the "settings" dict is an additional kwarg that will be passed to the model function on every call
    model = Model(
        _thomson_model,
        independent_vars=["wavelengths"],
        nan_policy="omit",
        settings=settings,
    )

    return model


def mc_error(
    wavelengths,
    data,
    model,
    result,
    nsamples=1000,
    chisqr_cutoff=0.05,
    max_iter=15,
    verbose=False,
):

    # Create a mask that sets NaN => 0
    not_nan = np.argwhere(np.logical_not(np.isnan(data)))
    data = data[not_nan]

    # Estimate the degrees of freedom (non-trivial datapoints)
    # This allows calculation of reduced chi2, without which the chisq_cutoff
    # appropriate will change with the length of the data sample.
    DOF = np.sum(np.where(data > 0.01, 1, 0))

    # Get the best-fit parameters
    params = result.init_params
    best_chisqr = result.redchi

    if best_chisqr > chisqr_cutoff:
        raise ValueError("Best fit chisqr > chisqr_cutoff. That won't work!")

    # Find all of the free parameters in params, add to a dict
    free_param_keys = []
    for p in list(params.keys()):
        if params[p].vary:
            free_param_keys.append(p)

    # Create a dictionary of ranges to explore
    ranges = {}
    Range = namedtuple("Range", ["min", "max", "center", "init_range"])
    for p in free_param_keys:
        min_val = params[p].min
        max_val = params[p].max
        center_val = result.best_values[p]
        init_range = 0.5 * (max_val - min_val)
        ranges[p] = Range(
            min=min_val, max=max_val, center=center_val, init_range=init_range
        )

    iteration = 0
    range_mult = 1
    free_params = {}

    final_nsamples = nsamples
    quick_samples = np.min([nsamples, 100])
    nsamples = quick_samples

    while iteration < max_iter:
        iteration += 1

        if verbose:
            print(f"Iteration: {iteration}, Nsamples: {nsamples}")

        # Create Monte-Carlo distributions of values for the free params
        for p in free_param_keys:
            r = ranges[p]
            lower = np.clip(r.center - range_mult * r.init_range, r.min, r.max)
            upper = np.clip(r.center + range_mult * r.init_range, r.min, r.max)

            # print(f"{p}-bounds: [{lower},{upper}]")

            free_params[p] = np.random.uniform(low=lower, high=upper, size=nsamples)

        # Create a list of dicts that can be used as keywords to override params
        changed_params = []
        for i in range(nsamples):
            changed_params.append({k: v[i] for (k, v) in free_params.items()})

        chisq = np.zeros(nsamples)
        for i in range(nsamples):
            # Evaluate the model at the parameters, normalize the result then
            # calculate the Chi2 against the data
            evaluated = model.eval(params, **changed_params[i], wavelengths=wavelengths)
            evaluated *= 1 / np.nanmax(evaluated[not_nan])

            chisq[i] = np.sum((evaluated[not_nan] - data) ** 2 / data) / DOF

        # Find all samples that are within the chisq cutoff
        ind = np.argwhere(chisq < chisqr_cutoff)

        fract_in = len(ind) / nsamples

        if verbose:
            print(f"Min chi2, Max chi2: {np.min(chisq)}, {np.max(chisq)}")
            print(f"Fract in: {fract_in}")

        # If all the chi2 are too big, reduce the range multiplier
        if fract_in < 0.3:
            # print("Decreasing range")
            range_mult *= 0.6
            nsamples = quick_samples
        # If all the chi2 are too small, increase the range multiplier
        elif fract_in > 0.7:
            # print("Increasing range")
            range_mult *= 1.5
            nsamples = quick_samples
        else:
            if nsamples == final_nsamples:
                # print(f"Terminating with {len(ind)} inside bounds")
                break
            else:
                # Repeat one more time, now with the total number of samples
                nsamples = final_nsamples

    if iteration >= max_iter:
        print("WARNING: did not finish properly!")

    # Calculate the min and max values for each parameter that satisfy the chi2 test
    for p in list(free_params.keys()):
        v = free_params[p][ind]
        free_params[p] = (np.min(v), np.max(v))

    return free_params
