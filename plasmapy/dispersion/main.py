from astropy import units as u
from numerical.cold_plasma_function_solver import cold_plasma_function_solver

inputs = {
   "B": 8.3e-9 * u.T,
   "k": [0.001] * u.rad / u.m,
   "omega_p": [1.6e5]* u.rad / u.s,
   "omega_e": [4.0e5] * u.rad / u.s,
   "omega_alpha": [3.0e5] * u.rad / u.s,
   "theta": 30 * u.deg,
}

w = cold_plasma_function_solver(**inputs)
for i in range(len(w[0.001])):
    print(w[0.001][i])
