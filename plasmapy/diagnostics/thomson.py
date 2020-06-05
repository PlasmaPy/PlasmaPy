"""
Defines the Thomson scattering analysis module as
part of the diagnostics package.
"""

__all__ = [
    "spectral_density",
]

import numpy as np
from astropy import units as u
import astropy.constants as const

from plasmapy.utils.decorators import validate_quantities
from plasmapy.formulary.parameters import plasma_frequency as wp

from plasmapy.formulary.dielectric import permittivity_1D_Maxwellian

from plasmapy.particles import Particle


# TODO: interface for inputting a multi-species configuration could be
# simplified using the plasmapy.classes.plasma_base class if that class
# included ion and electron drift velocities and information about the ion
# atomic species.

@validate_quantities(wavelengths=u.nm, probe_wavelength=u.nm, ne=u.cm**-3,
                  kTe=u.eV, kTi=u.eV, ion_vel=u.cm/u.s, fluid_vel=u.cm/u.s)
def spectral_density(wavelengths, probe_wavelength, ne, kTe, kTi,
                     fract=np.ones(1), ion_species=['H+'],
                     fluid_vel=np.zeros([3])*u.cm/u.s,
                     ion_vel=None,
                     probe_vec=np.array([1, 0, 0]),
                     scatter_vec=np.array([0, 1, 0])):
    r"""
    Calculate the spectral dispersion function for Thomson scattering of a
    probe laser beam by a multi-species Maxwellian plasma.

    Parameters
    ---------

    wavelengths : astropy.units.Quantity ndarray (required)
        Array of wavelengths over which the spectral density function
        will be calculated, convertable to nm.

    probe_wavelength : astropy.units.Quantity (required)
        Wavelength of the probe laser convertable to nm

    ne : astropy.units.Quantity, ndarray (required)
        mean (0th order) electron density of all plasma components combined
        convertable to cm^-3.

    kTe : astropy.units.Quantity (required)
        Temperature of the electron component convertable to eV.

    kTi : astropy.units.Quantity ndarray, shape [N] (required)
        Temperature of each ion component convertable to eV.

    fract : float ndarray, shape [N]
        Fraction (by number) of the total number of ions made up by each ion
        species. Must sum to 1.0. Default is a single ion species.

    ion_species : str ndarray, shape [N]
        Strings representing each ion species in the format interpretable by
        the plasmapy Particle class. Default is ['H+'] corresponding to a
        single species of hydrogen ions.

    fluid_vel : astropy.units.Quantity ndarray shape [3]
        Velocity of the fluid or electron component in the rest frame in
        units convertable to cm/s. Defaults to [0, 0, 0], representing
        a stationary plasma.

    ion_vel : astropy.units.Quantity ndarray, shape [N,3]
        Velocity vectors for each ion population relative to the
        fluid velocit in units convertable to cm/s. Defaults zero drift
        for all specified ion species.

    probe_vec : float ndarray, shape [3]
        Unit vector in the direction of the probe laser. Defaults to
        [1, 0, 0].

    scatter_vec : float ndarray, shape [3]
        Unit vector pointing from the scattering volume to the detector.
        Defaults to [0, 1, 0] which, along with the default probe_vec,
        corresponds to a 90 degree scattering angle geometry.

    Returns
    -------
    alpha : float
        mean scattering parameter. alpha > 1 corresponds to collective
        scattering, while alpha < 1 indicates non-collective scattering.

    Skw : astropy.units.Quantity ndarray
        Computed spectral density function over the input wavelength range
        with units of nm.

    Notes
    -----
    For details, see "Plasma Scattering of Electromagnetic Radiation" by
    Sheffield et al. ISBN 978-0123748775. This code is a modified version of
    the program described therein.

    For a concise summary of the relevant physics, see Chapter 5 of Derek
    Schaeffer's thesis, DOI: 10.5281/zenodo.3766933

    """

    if np.sum(fract) != 1.0:
        print("WARNING: Sum(fract) != 1.0. Normalizing array.")
        fract = fract/np.sum(fract)

    # If ion drift velocity is not specified, create an array corresponding
    # to zero drift
    if ion_vel is None:
        ion_vel = np.zeros([fract.size, 3])*u.cm/u.s

    #Check to make sure the ion arrays are compatible
    if (ion_species.size != fract.size) or \
        (ion_vel.shape[0] != fract.size) or (kTi.size != fract.size):
        raise ValueError("WARNING: Inconsistent number of species in fract, "
              "ion_species, kTi, and/or ion_vel.")

    # Ensure unit vectors are normalized
    probe_vec = probe_vec/np.linalg.norm(probe_vec)
    scatter_vec = scatter_vec/np.linalg.norm(scatter_vec)

    # Define some constants
    C = const.c.to(u.cm/u.s)  # speed of light
    m_e = const.m_e.to(u.eV/(u.cm/u.s)**2)
    Te = (kTe/const.k_B).to(u.K)
    Ti = (kTi/const.k_B).to(u.K)

    #Calculate the mass and charge of each ion species
    particles = [Particle(s) for s in ion_species]
    Mi = u.Quantity([p.mass.to(u.eV/(u.cm/u.s)**2) for p in particles])
    ion_z = u.Quantity([p.integer_charge for p in particles])
    zbar = np.sum(fract*ion_z)
    ni = fract*ne/zbar # ne/zbar = sum(ni)

    vTe = np.sqrt(2*kTe/m_e)  # Electron thermal velocity
    vTi = np.sqrt(2*kTi/Mi)  # Ion thermal velocity
    # Electron plasma frequency
    wpe = wp(n=ne)

    # Compute the ion velocity in the rest frame
    ion_vel = fluid_vel + ion_vel

    # Convert wavelengths to angular frequencies (electromagnetic waves, so
    # phase speed is c)
    ws = (2*np.pi*u.rad*C/wavelengths).to(u.rad/u.s)
    wl = (2*np.pi*u.rad*C/probe_wavelength).to(u.rad/u.s)

    # Compute the frequency shift (required by energy conservation)
    w = ws - wl

    # Compute the wavenumbers in the plasma
    # See Sheffield Sec. 1.8.1 and Eqs. 5.4.1 and 5.4.2
    ks = np.sqrt(ws**2 - wpe**2)/C
    kl = np.sqrt(wl**2 - wpe**2)/C

    # Compute the wavenumber shift (required by momentum conservation)
    scattering_angle = np.arccos(np.dot(probe_vec, scatter_vec))
    # Eq. 1.7.10 in Sheffield
    k = np.sqrt(ks**2 + kl**2 - 2*ks*kl*np.cos(scattering_angle))
    # Normal vector along k
    k_vec = (scatter_vec - probe_vec)* u.dimensionless_unscaled

    # Compute Doppler-shifted frequencies for both the ions and electrons
    # Matmul is simultaneously conducting dot product over all wavelengths
    # and ion components
    w_e = w - k*np.dot(fluid_vel, k_vec)
    w_i = w - np.matmul(ion_vel, np.outer(k, k_vec).T)

    # Compute the scattering parameter alpha
    # expressed here using the fact that v_th/w_p = root(2) * Debye length
    alpha = np.sqrt(2)*wpe/(k*vTe)

    # Calculate the normalized phase velocities (Sec. 3.4.2 in Sheffield)
    xe = (w_e/(k*vTe)).to(u.dimensionless_unscaled)
    xi = (np.outer(1/vTi, 1/k)*w_i).to(u.dimensionless_unscaled)

    # Calculate the succeptabilities
    chiE = permittivity_1D_Maxwellian(w_e, k, Te, ne, 'e-')

    # Treatment of multiple species is an extension of the discussion in
    # Sheffield Sec. 5.1
    chiI = np.zeros([fract.size, w.size], dtype=np.complex128)
    for i,s in enumerate(particles):
        chiI[i, :] = permittivity_1D_Maxwellian(w_i[i, :], k, Ti[i], ni[i],
                                                ion_species[i],
                                                z_mean=ion_z[i])

    # Calculate the logitudinal dielectric function
    epsilon = 1 + chiE + np.sum(chiI, axis=0)

    # Calculate the contributions to the spectral density function
    econtr = 2*np.sqrt(np.pi)/k/vTe* \
        np.power(np.abs(1 - chiE/epsilon), 2)*np.exp(-xe**2)

    icontr = np.zeros([fract.size, w.size], dtype=np.complex128)*u.s/u.rad
    for m in range(fract.size):
        icontr[m, :] = 2*np.sqrt(np.pi)*ion_z[m]/k/vTi[m]* \
            np.power(np.abs(chiE/epsilon), 2)*np.exp(-xi[m, :]**2)

    # Recast as real: imaginary part is already zero
    Skw = np.real(econtr + np.sum(icontr, axis=0))

    return np.mean(alpha), Skw
