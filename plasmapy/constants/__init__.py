"""Physical and mathematical constants for use within PlasmaPy"""

from numpy import pi

from astropy.constants.si import (h, hbar, k_B, c, G, g0, m_p, m_n, m_e,
                                  u, sigma_sb, e, eps0, N_A, R, Ryd, a0,
                                  muB, mu0, sigma_T, au, pc, kpc, L_sun,
                                  M_sun, R_sun, M_earth, R_earth)

try:
    from astropy.constants import atm  # astropy 2.0 and later
except:
    from astropy.constants import atmosphere as atm  # astropy 1.x and before

from .atomic import (element_symbol, isotope_symbol, atomic_number,
                     mass_number, element_name, standard_atomic_weight,
                     isotope_mass, ion_mass, nuclear_binding_energy, half_life,
                     energy_from_nuclear_reaction, is_isotope_stable,
                     known_isotopes, common_isotopes, stable_isotopes,
                     isotopic_abundance, charge_state)

from .elements import Elements

from .isotopes import Isotopes

from . import tests
