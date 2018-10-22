# coding: utf-8
"""
Langmuir probe data analysis
============================

Let's analyze a few Langmuir probe characteristics using the
`diagnostics.langmuir` subpackage. First we need to import the module and some
basics.
"""

from plasmapy.diagnostics.langmuir import Characteristic, swept_probe_analysis
import astropy.units as u
import numpy as np
import os

######################################################
# The first characteristic we analyze is a simple single-probe measurement in
# a low (ion) temperature, low density plasma with a cylindrical probe. This
# allows us to utilize OML theory implemented in `swept_probe_analysis`.
# The data has been preprocessed with some smoothing, which allows us to obtain
# a Electron Energy Distribution Function (EEDF) as well.

# Load the bias and current values stored in the .p pickle file.
path = os.path.join("langmuir_samples", "Beckers2017.npy")
bias, current = np.load(path)

# Create the Characteristic object, taking into account the correct units
characteristic = Characteristic(np.array(bias) * u.V,
                                np.array(current)*1e3 * u.mA)

# Calculate the cylindrical probe surface area
probe_length = 1.145 * u.mm
probe_diameter = 1.57 * u.mm
probe_area = (probe_length * np.pi * probe_diameter +
              np.pi * 0.25 * probe_diameter**2)

######################################################
# Now we can actually perform the analysis. Since the plasma is in Helium an
# ion mass number of 4 is entered. The results are visualized and the obtained
# EEDF is also shown.
print(swept_probe_analysis(characteristic,
                           probe_area, 'He-4+',
                           visualize=True,
                           plot_EEDF=True))

######################################################
# The cyan and yellow lines indicate the fitted electron and ion currents,
# respectively. The green line is the sum of these and agrees nicely with the
# data. This indicates a successful analysis.

######################################################
# The next sample probe data is provided by David Pace. It is also obtained
# from a low relatively ion temperature and density plasma, in Argon.

# Load the data from a file and create the Characteristic object
path = os.path.join("langmuir_samples", "Pace2015.npy")
bias, current = np.load(path)
characteristic = Characteristic(np.array(bias) * u.V,
                                np.array(current) * 1e3 * u.mA)

######################################################
# Initially the electrons are assumed to be Maxwellian. To check this the fit
# of the electron growth region will be plotted.
swept_probe_analysis(characteristic,
                     0.738 * u.cm**2,
                     'Ar-40 1+',
                     bimaxwellian=False,
                     plot_electron_fit=True)

######################################################
# It can be seen that this plasma is slightly bi-Maxwellian, as there are two
# distinct slopes in the exponential section. The analysis is now performed
# with bimaxwellian set to True, which yields improved results.
print(swept_probe_analysis(characteristic,
                           0.738 * u.cm**2,
                           'Ar-40 1+',
                           bimaxwellian=True,
                           visualize=True,
                           plot_electron_fit=True))

######################################################
# The probe current resolution of the raw data is relatively poor, but the
# analysis still performs well in the ion current region. The bi-Maxwellian
# properties are not significant but do make a difference. Check this analysis
# without setting `bimaxwellian` to True!
# This is reflected in the results, which indicate that the temperatures of
# the cold and hot electron population are indeed different, but relatively
# close.

######################################################
# This Helium plasma is fully bi-Maxwellian.

# Import probe data and calculate probe surface area.
path = os.path.join("langmuir_samples", "Beckers2017b.npy")
bias, current = np.load(path)
characteristic = Characteristic(np.array(bias) * u.V,
                                np.array(current) * 1e3 * u.mA)
probe_length = 1.145 * u.mm
probe_diameter = 1.57 * u.mm
probe_area = (probe_length * np.pi * probe_diameter +
              np.pi * 0.25 * probe_diameter**2)

######################################################
# `plot_electron_fit` is set to True to check the bi-Maxwellian properties.
# The fit converges nicely to the two slopes of the electron growth region.
print(swept_probe_analysis(characteristic,
                           probe_area,
                           'He-4+',
                           bimaxwellian=True,
                           plot_electron_fit=True,
                           visualize=True))
