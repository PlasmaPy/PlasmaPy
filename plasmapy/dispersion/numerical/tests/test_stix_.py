from astropy import units as u

from ..stix_ import stix

inputs = {
    "B": 8.3e-9 * u.T,
    "k": 0.001* u.rad / u.m,
    "ions": ['e-','H+'],
    "omega_ions": [4.0e5,2.0e5] * u.rad / u.s,
    "theta": 30 * u.deg,
}
w =stix(**inputs)
print(w[0.001])
