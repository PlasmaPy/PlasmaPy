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

from plasmapy.formulary.dispersionfunction \
    import plasma_dispersion_func_deriv as ZPrime

from plasmapy.utils.decorators import validate_quantities


# TODO: interface for inputting a multi-species configuration could be
# simplified using the plasmapy.classes.plasma_base class if that class
# included ion and electron drift velocities and information about the ion
# atomic species.

@validate_quantities(wavelength=u.nm, probe_wavelength=u.nm, ne=u.cm**-3,
                  Te=u.eV, Ti=u.eV, ion_vel=u.cm/u.s, fluid_vel=u.cm/u.s)
def spectral_density(wavelength, probe_wavelength=532*u.nm, ne=1e15*u.cm**-3,
                     fract=np.ones(1), Te=1*u.eV, Ti=np.ones(1)*u.eV,
                     ion_z=np.ones(1), ion_mu=np.ones(1),
                     fluid_vel=np.zeros([3])*u.cm/u.s,
                     ion_vel=np.zeros([1, 3])*u.cm/u.s,
                     probe_n=np.array([1, 0, 0]),
                     scatter_n=np.array([0, 1, 0])):
    r"""
    Calculate the spectral dispersion function for Thomson scattering of a
    probe laser beam by a multi-species Maxwellian plasma.

    Positional Arguments
    --------------------
    wavelength : astropy.units.Quantity ndarray
        Array of wavelengths over which the spectral density function
        will be calculated, convertable to nm.

    Keyword Arguments
    -----------------
    probe_wavelength : astropy.units.Quantity
        Wavelength of the probe laser convertable to nm

    ne : astropy.units.Quantity, ndarray
        mean (0th order) electron density of all plasma components combined
        convertable to cm^-3.

    fract : float ndarray, shape [N]
        Relative fractions of the overall electron charge density
        corresponding to each ion population such that ne = -SUM(fract*Z*ni)

    Te : astropy.units.Quantity
        Temperature of the electron component convertable to eV.

    Ti : astropy.units.Quantity ndarray, shape [N]
        Temperature of each ion component convertable to eV.

    ion_z : float ndarray, shape [N]
        Charge of each ion species normalized to the proton charge.

    ion_mu : float ndarray, shape [N]
        Mass of each ion species normalized to the proton mass.

    fluid_vel : astropy.units.Quantity ndarray shape [3]
        Velocity of the fluid or electron component in the rest frame in
        units convertable to cm/s

    ion_vel : astropy.units.Quantity ndarray, shape [N,3]
        Velocity vectors for each ion population relative to the
        fluid velocit in units convertable to cm/s.

    probe_n : float ndarray, shape [3]
        Unit vector in the direction of the probe laser.

    scatter_n : float ndarray, shape [3]
        Unit vector pointing from the scattering volume to the detector.

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

    # Ensure unit vectors are normalized
    probe_n = probe_n/np.linalg.norm(probe_n)
    scatter_n = scatter_n/np.linalg.norm(scatter_n)

    # Define some constants
    C = const.c.to(u.cm/u.s)  # speed of light
    m_e = const.m_e.to(u.eV/(u.cm/u.s)**2)  # Electron mass in eV per (cm/s)^2
    m_p = const.m_p.to(u.eV/(u.cm/u.s)**2)  # Proton mass in same units
    Mi = ion_mu*m_p  # mass of each ion species in same units
    vTe = np.sqrt(Te/m_e)  # Electron thermal velocity
    vTi = np.sqrt(Ti/Mi)  # Ion thermal velocity
    # Electron plasma frequency per the NRL formulary
    wpe = 5.64e4*(u.cm**1.5*u.rad/u.s)*np.sqrt(ne)

    # Compute the ion velocity in the rest frame
    ion_vel = fluid_vel + ion_vel

    # Convert wavelengths to angular frequencies (electromagnetic waves, so
    # phase speed is c)
    ws = (2*np.pi*u.rad*C/wavelength).to(u.rad/u.s)
    wl = (2*np.pi*u.rad*C/probe_wavelength).to(u.rad/u.s)

    # Compute the frequency shift (required by energy conservation)
    w = ws - wl

    # Compute the wavenumbers in the plasma
    # See Sheffield Sec. 1.8.1 and Eqs. 5.4.1 and 5.4.2
    ks = np.sqrt(ws**2 - wpe**2)/C
    kl = np.sqrt(wl**2 - wpe**2)/C

    # Compute the wavenumber shift (required by momentum conservation)
    scattering_angle = np.arccos(np.dot(probe_n, scatter_n))
    k = np.sqrt(ks**2 + kl**2 - 2*ks*kl*np.cos(scattering_angle))  # Eq. 1.7.10 in Sheffield
    k_n = scatter_n - probe_n  # Normal vector along k

    # Compute Doppler-shifted frequencies for both the ions and electrons
    # Matmul is simultaneously conducting dot product over all wavelengths
    # and ion components
    w_e = w - k*np.dot(fluid_vel, k_n)
    w_i = w - np.matmul(ion_vel, np.outer(k, k_n).T)

    # Compute the scattering parameter alpha
    # expressed here using the fact that v_th/w_p = root(2) * Debye length
    alpha = wpe/(np.sqrt(2)*k*vTe)

    # Calculate the normalized phase velocities (Sec. 3.4.2 in Sheffield)
    xe = w_e/(k*np.sqrt(2)*vTe)
    xi = 1/np.sqrt(2)*np.outer(1/vTi, 1/k)*w_i

    # Calculate the succeptabilities
    # Treatment of multiple species is an extension of the discussion in
    # Sheffield Sec. 5.1
    chiE = -0.5*np.power(alpha, 2)*ZPrime(xe)

    chiI = np.zeros([fract.size, w.size], dtype=np.complex128)
    for m in range(fract.size):
        chiI[m, :] = -0.5*fract[m]*alpha**2*ion_z[m]*(Te/Ti[m])* \
            ZPrime(xi[m, :])

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
