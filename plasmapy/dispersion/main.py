from astropy import units as u
from numerical.cold_plasma_function_solver import cold_plasma_function_solver

k = [0.001,0.005] 

inputs = {
   "B": 8.3e-9 * u.T,
   "k": k* u.rad / u.m,
   "ions": ['e-','H+'],
   "omega_ions": [4.0e5,2.0e5] * u.rad / u.s,
   "theta": 30 * u.deg,
}

w = cold_plasma_function_solver(**inputs)
for i in range(len(k)):
    print(w[k[i]])
    print(w[k[i]][0])

