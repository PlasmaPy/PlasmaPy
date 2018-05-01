"""
Contains physical and mathematical constants commonly used in plasma
physics.

"""

from numpy import pi

from astropy.constants.si import (
    e,
    mu0,
    eps0,
    k_B,
    c,
    G,
    h,
    hbar,
    m_p,
    m_n,
    m_e,
    u,
    sigma_sb,
    N_A,
    R,
    Ryd,
    a0,
    muB,
    sigma_T,
    au,
    pc,
    kpc,
    g0,
    L_sun,
    M_sun,
    R_sun,
    M_earth,
    R_earth,
)

from astropy.constants import atm

# The following code is modified from astropy.constants to produce a
# table containing information on the constants contained with PlasmaPy.
# Mathematical constants can be just entered.

_lines = [
    'The following constants are available:\n',
    '========== ================= ================ ============================================',
    'Name       Value             Units            Description',
    '========== ================= ================ ============================================',
    "   pi      3.141592653589793                  Ratio of circumference to diameter of circle",
]

_constants = [eval(item) for item in dir() if item[0] != '_' and item != 'pi']
for _const in _constants:
    _lines.append('{0:^10} {1:^17.12g} {2:^16} {3}'
                  .format(_const.abbrev, _const.value, _const._unit_string, _const.name))

_lines.append(_lines[1])

__doc__ += '\n'.join(_lines)

del _lines, _const, _constants
