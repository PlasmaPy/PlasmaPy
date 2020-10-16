"""Tests for the two fluid dispersion solver"""

import numpy as np
import pytest

from astropy import units as u

from plasmapy.dispersion.two_fluid_dispersion_solver import
two_fluid_dispersion_solution

k = 0.01 * u.m ** -1
theta = 30 * u.deg
B = 8.3E-9 * u.T
T_e = 1.6e6 * u.K
T_i = 4.e5 * u.K
z = 1

k_arr = np.linspace(10**-7, 10**-2, 1E4) * u.m ** -1
theta_arr = np.linspace(5, 85, 100) * u.deg
n = 5 * u.cm ** -3
B = 8.3E-9 * u.T
T_e = 1.6e6 * u.K
T_i = 4.e5 * u.K
z = 1
c = 3.e8 * u.m/u.s
