# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 15:01:23 2020

@author: Peter
"""

import matplotlib.pyplot as plt

import astropy.constants as const
import astropy.units as u
import numpy as np

from typing import List, Tuple, Union

from plasmapy.formulary.parameters import plasma_frequency
from plasmapy.particles import Particle
from plasmapy.utils.decorators import validate_quantities

from scipy.special import exp1

@validate_quantities(
    frequencies={"can_be_negative": False},
    ne={"can_be_negative": False},
    ni={"can_be_negative": False},
    Te={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def maxwellian_bremsstrahlung(
    frequencies: u.Hz,
    ne: u.m ** -3,
    ni: u.m ** -3,
    Te: u.K,
    ion_species: Union[str, Particle] = "H+",
    kmax : u.m = None,
) -> np.ndarray:
    r"""
    Calculate the Bremsstrahlung emission spectrum for a Maxwellian plasma
    in the Rayleigh-Jeans limit hbar*w << k_B*Te

    Parameters
    ----------

    frequencies : `~astropy.units.Quantity`
        Array of frequencies over which the bremsstrahlung spectrum will be
        calculated (convertable to Hz).

    ne : `~astropy.units.Quantity`
        Electron number density in the plasma (convertable to m^-3)

    ni : `~astropy.units.Quantity`
        Ion number density in the plasma (convertable to m^-3)

    Te : `~astropy.units.Quantity`
        Temperature of the electrons (in K or convertible to eV)

    ion_species : str or `~plasmapy.particles.Particle`, optional
        An instance of `~plasmapy.particles.Particle`, or a string
        convertible to `~plasmapy.particles.Particle`. Default is `'H+'`
        corresponding to hydrogen ions.

    kmax :  `~astropy.units.Quantity`
        Cutoff wavenumber (convertable to u.rad/u.m). Defaults to the inverse
        of the electron de Broglie wavelength.

    Returns
    -------


    spectrum : `~astropy.units.Quantity`
        Computed bremsstrahlung spectrum over the frequencies provided.

    Notes
    -----



    """

    # TODO: Error check that all frequencies are GREATER than wpe

    # TODO: Error check that hw << kTe (Rayleigh-Jeans limit)

    # Condition ion_species
    if isinstance(ion_species, str):
        ion_species = Particle(ion_species)

    # Default value of kmax is the electrom thermal de Broglie wavelength
    if kmax is None:
        kmax = (np.sqrt(const.m_e.si*const.k_B.si*Te)/const.hbar.si).to(1/u.m)

    # Convert frequencies to angular frequencies
    w = (frequencies*2*np.pi*u.rad).to(u.rad / u.s)

    wpe = plasma_frequency(n=ne, particle="e-")

    c1 = (8/3)*np.sqrt(2/np.pi) \
         *(const.e.si**2/(4*np.pi*const.eps0.si))**3 \
         *1/(const.m_e.si*const.c.si**2)**1.5

    Zi = ion_species.integer_charge
    c2 = np.sqrt(1 - wpe**2/w**2)*Zi**2*ni*ne/np.sqrt(const.k_B.si*Te)

    print(c1.unit)
    print(c2.unit)

    # Dimensionless argument for exponential integral
    arg = 0.5*w**2*const.m_e.si/(kmax**2*const.k_B.si*Te)/u.rad**2
    # Remove units, get ndarray of values
    arg = (arg.to(u.dimensionless_unscaled)).value

    return c1*c2*exp1(arg)



frequencies = np.arange(14, 17, .05)
frequencies = (10**frequencies)/u.s



energy = (frequencies*const.h.si).to(u.eV)



ne = 1e22*u.cm**-3
ni=ne
Te = 100*u.eV

spectrum = maxwellian_bremsstrahlung(frequencies, ne, ni, Te)



plt.plot(energy, spectrum)
plt.xlabel('Energy (eV)')
plt.show()