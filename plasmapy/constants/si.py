"""
Physical constants in SI units.
"""

import numpy as np
import warnings
from astropy.constants import Constant, EMConstant

# Constant gives a warning when a variable is already defined in
# Astropy, even though these are not defined in Astropy

from astropy.utils.exceptions import AstropyWarning
warnings.filterwarnings('ignore', category=AstropyWarning, append=True)

# MATHEMATICAL CONSTANTS

pi = 3.141592653589793238462643383279

# PHYSICAL CONSTANTS

h = Constant('h', "Planck constant", 6.62606957e-34, 'J s', 0.00000029e-34,
             'CODATA 2010', system='si')

hbar = Constant('hbar', "Reduced Planck constant", h.value * 0.5 / np.pi,
                'J s', h.uncertainty * 0.5 / np.pi, h.reference, system='si')

k_B = Constant('k_B', "Boltzmann constant", 1.3806488e-23, 'J / (K)',
               0.0000013e-23, 'CODATA 2010', system='si')

c = Constant('c', "Speed of light in vacuum", 2.99792458e8, 'm / (s)', 0.,
             'CODATA 2010', system='si')

G = Constant('G', "Gravitational constant", 6.67384e-11, 'm3 / (kg s2)',
             0.00080e-11, 'CODATA 2010', system='si')

g0 = Constant('g0', "Standard acceleration of gravity", 9.80665, 'm / s2', 0.0,
              'CODATA 2010', system='si')

m_p = Constant('m_p', "Proton mass", 1.672621777e-27, 'kg', 0.000000074e-27,
               'CODATA 2010', system='si')

m_n = Constant('m_n', "Neutron mass", 1.674927351e-27, 'kg', 0.000000074e-27,
               'CODATA 2010', system='si')

m_e = Constant('m_e', "Electron mass", 9.10938291e-31, 'kg', 0.00000040e-31,
               'CODATA 2010', system='si')

u = Constant('u', "Atomic mass unit", 1.660538921e-27, 'kg', 0.000000073e-27,
             'CODATA 2010', system='si') 

sigma_sb = Constant('sigma_sb', "Stefan-Boltzmann constant", 5.670373e-8,
                    'W / (K4 m2)', 0.000021e-8, 'CODATA 2010', system='si')

# EM constants require a system to be specified
e = EMConstant('e', 'Electron charge', 1.602176565e-19, 'C', 0.000000035e-19,
               'CODATA 2010', system='si')

eps0 = EMConstant('eps0', 'Electric constant', 8.854187817e-12, 'F/m', 0.0,
                  'CODATA 2010', system='si')

N_A = Constant('N_A', "Avogadro's number", 6.02214129e23, '1 / (mol)',
               0.00000027e23, 'CODATA 2010', system='si')

R = Constant('R', "Gas constant", 8.3144621, 'J / (K mol)', 0.0000075,
             'CODATA 2010', system='si')

Ryd = Constant('Ryd', 'Rydberg constant', 10973731.568539, '1 / (m)', 0.000055,
               'CODATA 2010', system='si')

a0 = Constant('a0', "Bohr radius", 0.52917721092e-10, 'm', 0.00000000017e-10,
              'CODATA 2010', system='si')

muB = Constant('muB', "Bohr magneton", 927.400968e-26, 'J/T', 0.00002e-26,
               'CODATA 2010', system='si')

atm = Constant('atmosphere', "Atmosphere", 101325, 'Pa', 0.0,
               'CODATA 2010', system='si')

mu0 = Constant('mu0', "Magnetic constant", 4.0e-7 * np.pi, 'N/A2', 0.0,
               'CODATA 2010', system='si')

sigma_T = Constant('sigma_T', "Thomson scattering cross-section", 
                   0.6652458734e-28, 'm2', 0.0000000013e-28, 
                   'CODATA 2010', system='si')

# DISTANCE

au = Constant('au', "Astronomical Unit", 1.49597870700e11, 'm', 0.0,
              "IAU 2012 Resolution B2", system='si')

pc = Constant('pc', "Parsec", au.value / np.tan(np.radians(1. / 3600.)), 'm',
              au.uncertainty / np.tan(np.radians(1. / 3600.)),
              "Derived from au", system='si')

kpc = Constant('kpc', "Kiloparsec",
               1000. * au.value / np.tan(np.radians(1. / 3600.)), 'm',
               1000. * au.uncertainty / np.tan(np.radians(1. / 3600.)),
               "Derived from au", system='si')

# SOLAR QUANTITIES

L_sun = Constant('L_sun', "Solar luminosity", 3.846e26, 'W', 0.0005e26,
                 "Allen's Astrophysical Quantities 4th Ed.", system='si')

M_sun = Constant('M_sun', "Solar mass", 1.9891e30, 'kg', 0.00005e30,
                 "Allen's Astrophysical Quantities 4th Ed.", system='si')

R_sun = Constant('R_sun', "Solar radius", 6.95508e8, 'm', 0.00026e8,
                 "Allen's Astrophysical Quantities 4th Ed.", system='si')

# OTHER SOLAR SYSTEM QUANTITIES

M_earth = Constant('M_earth', "Earth mass", 5.9742e24, 'kg', 0.00005e24,
                   "Allen's Astrophysical Quantities 4th Ed.", system='si')

R_earth = Constant('R_earth', "Earth equatorial radius", 6.378136e6, 'm',
                   0.0000005e6, "Allen's Astrophysical Quantities 4th Ed.",
                   system='si')

warnings.filters.pop() # should be same as before 
