"""
Defines the Thomson scattering analysis module as
part of the diagnostics package.
"""

__all__ = [
    "spectral_density",
]

import astropy.constants as const
import astropy.units as u
import numpy as np

from plasmapy.formulary.parameters import plasma_frequency, thermal_speed
from plasmapy.formulary.dielectric import permittivity_1D_Maxwellian
from plasmapy.particles import Particle
from plasmapy.utils.decorators import validate_quantities
from typing import List, Tuple, Union


# TODO: interface for inputting a multi-species configuration could be
# simplified using the plasmapy.classes.plasma_base class if that class
# included ion and electron drift velocities and information about the ion
# atomic species.


@validate_quantities(
    wavelengths={"can_be_negative": False},
    probe_wavelength={"can_be_negative": False},
    ne={"can_be_negative": False},
    Te={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    Ti={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def spectral_density(
    wavelengths: u.nm,
    probe_wavelength: u.nm,
    ne: u.m ** -3,
    Te: u.K,
    Ti: u.K,
    fract: np.ndarray = np.ones(1),
    ion_species: Union[str, List[str], Particle, List[Particle]] = "H+",
    fluid_vel: u.m / u.s = np.zeros(3) * u.m / u.s,
    ion_vel: u.m / u.s = None,
    probe_vec=np.array([1, 0, 0]),
    scatter_vec=np.array([0, 1, 0]),
) -> Tuple[Union[np.floating, np.ndarray], np.ndarray]:
    r"""
    Calculate the spectral density function for Thomson scattering of a
    probe laser beam by a multi-species Maxwellian plasma.

    This function calculates the spectral density function for Thomson
    scattering of a probe laser beam by a plasma consisting of one or more ion
    species and a neutralizing electron fluid:

    .. math::
        S(k,\omega) = 1 - \frac{2\pi}{k}
        \bigg |1 - \frac{\chi_e}{\epsilon} \bigg |^2
        f_{e0} \bigg (\frac{\omega}{k} \bigg ) +
        \sum_i \frac{2\pi Z_i}{k}
        \bigg |\frac{\chi_e}{\epsilon} \bigg |^2 f_{i0,i}
        \bigg ( \frac{\omega}{k} \bigg )

    where :math:`\chi_e` is the electron component susceptibility of the
    plasma and :math:`\epsilon = 1 + \chi_e + \sum_i \chi_i` is the total
    plasma dielectric  function (with :math:`\chi_i` being the ion component
    of the susceptibility), :math:`Z_i` is the charge of each ion, :math:`k`
    is the scattering wavenumber, :math:`\omega` is the scattering frequency,
    and :math:`f_{e0}` and :math:`f_{i0,i}` are the electron and ion velocity
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

    ne : `~astropy.units.Quantity`
        Mean (0th order) electron density of all plasma components combined.
        (convertible to cm^-3.)

    Te : `~astropy.units.Quantity`
        Temperature of the electron component. (in K or convertible to eV)

    Ti : `~astropy.units.Quantity`, shape (N, )
        Temperature of each ion component. Shape (N, ) must be equal to the
        number of ion components N. (in K or convertible to eV)

    fract : array_like, shape (N, ), optional
        An array-like object where each element represents the fraction (or ratio)
        of the ion component number density to the total ion number density.
        Must sum to 1.0. Default is a single ion species.

    ion_species : str or `~plasmapy.particles.Particle`, shape (N, ), optional
        A list or single instance of `~plasmapy.particles.Particle`, or strings
        convertible to `~plasmapy.particles.Particle`. Default is `'H+'`
        corresponding to a single species of hydrogen ions.

    fluid_vel : `~astropy.units.Quantity`, shape (3, ), optional
        Electron fluid velocity in the rest frame. (convertible to m/s)
        Defaults to a stationary plasma [0, 0, 0] m/s.

    ion_vel : `~astropy.units.Quantity`, shape (N, 3), optional
        Velocity vectors for each ion population relative to the
        fluid velocity. (convertible to m/s) Defaults zero drift
        for all specified ion species.

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
        scattering and `alpha` < 1 indicates non-collective scattering.

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
    fract = np.asarray(fract, dtype=np.float64)

    # If ion drift velocity is not specified, create an array corresponding
    # to zero drift
    if ion_vel is None:
        ion_vel = np.zeros([fract.size, 3]) * u.m / u.s

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
            f"Got {Ti.size} ion temperatures and expected " f"{len(ion_species)}."
        )

    if (
        (len(ion_species) != fract.size)
        or (ion_vel.shape[0] != fract.size)
        or (Ti.size != fract.size)
    ):
        raise ValueError(
            "Inconsistent number of species in fract, "
            "ion_species, Ti, and/or ion_vel."
        )

    # Ensure unit vectors are normalized
    probe_vec = probe_vec / np.linalg.norm(probe_vec)
    scatter_vec = scatter_vec / np.linalg.norm(scatter_vec)

    # Define some constants
    C = const.c.si  # speed of light

    # Calculate plasma parameters
    vTe = thermal_speed(Te, particle="e-")
    vTi, ion_z = [], []
    for T, ion in zip(Ti, ion_species):
        vTi.append(thermal_speed(T, particle=ion).value)
        ion_z.append(ion.integer_charge * u.dimensionless_unscaled)
    vTi = vTi * vTe.unit
    zbar = np.sum(fract * ion_z)
    ni = fract * ne / zbar  # ne/zbar = sum(ni)
    wpe = plasma_frequency(n=ne, particle="e-")

    # Compute the ion velocity in the rest frame
    ion_vel = fluid_vel + ion_vel

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

    # Compute the wavenumber shift (required by momentum conservation)
    scattering_angle = np.arccos(np.dot(probe_vec, scatter_vec))
    # Eq. 1.7.10 in Sheffield
    k = np.sqrt(ks ** 2 + kl ** 2 - 2 * ks * kl * np.cos(scattering_angle))
    # Normal vector along k
    k_vec = (scatter_vec - probe_vec) * u.dimensionless_unscaled

    # Compute Doppler-shifted frequencies for both the ions and electrons
    # Matmul is simultaneously conducting dot product over all wavelengths
    # and ion components
    w_e = w - k * np.dot(fluid_vel, k_vec)
    w_i = w - np.matmul(ion_vel, np.outer(k, k_vec).T)

    # Compute the scattering parameter alpha
    # expressed here using the fact that v_th/w_p = root(2) * Debye length
    alpha = np.sqrt(2) * wpe / (k * vTe)

    # Calculate the normalized phase velocities (Sec. 3.4.2 in Sheffield)
    xe = (w_e / (k * vTe)).to(u.dimensionless_unscaled)
    xi = (np.outer(1 / vTi, 1 / k) * w_i).to(u.dimensionless_unscaled)

    # Calculate the susceptibilities
    chiE = permittivity_1D_Maxwellian(w_e, k, Te, ne, "e-")

    # Treatment of multiple species is an extension of the discussion in
    # Sheffield Sec. 5.1
    chiI = np.zeros([fract.size, w.size], dtype=np.complex128)
    for i, ion in enumerate(ion_species):
        chiI[i, :] = permittivity_1D_Maxwellian(
            w_i[i, :], k, Ti[i], ni[i], ion, z_mean=ion_z[i]
        )

    # Calculate the longitudinal dielectric function
    epsilon = 1 + chiE + np.sum(chiI, axis=0)

    # Calculate the contributions to the spectral density function
    econtr = (
        2
        * np.sqrt(np.pi)
        / k
        / vTe
        * np.power(np.abs(1 - chiE / epsilon), 2)
        * np.exp(-(xe ** 2))
    )

    icontr = np.zeros([fract.size, w.size], dtype=np.complex128) * u.s / u.rad
    for m in range(fract.size):
        icontr[m, :] = (
            2
            * np.sqrt(np.pi)
            * ion_z[m]
            / k
            / vTi[m]
            * np.power(np.abs(chiE / epsilon), 2)
            * np.exp(-xi[m, :] ** 2)
        )

    # Recast as real: imaginary part is already zero
    Skw = np.real(econtr + np.sum(icontr, axis=0))

    return np.mean(alpha), Skw
