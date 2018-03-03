"""
Create a dictionary containing basic information for isotopes and
neutrons.
"""

import astropy.units as u

_Isotopes = {
    'n': {
        'name': 'neutron',
        'atomic number': 0,
        'mass number': 1,
        'mass': 1.00866491588 * u.u,
        'stable': False,
        'half-life': 881.5 * u.s,
    },

    'H-1': {
        'name': 'protium',
        'atomic number': 1,
        'mass number': 1,
        'mass': 1.00782503223 * u.u,
        'stable': True,
        'abundance': 0.999885,
    },

    'D': {
        'name': 'deuterium',
        'atomic number': 1,
        'mass number': 2,
        'mass': 2.01410177812 * u.u,
        'stable': True,
        'abundance': 0.000115,
    },

    'T': {
        'name': 'tritium',
        'atomic number': 1,
        'mass number': 3,
        'mass': 3.0160492779 * u.u,
        'stable': False,
        'half-life': 388800000.0 * u.s,
    },

    'H-4': {
        'atomic number': 1,
        'mass number': 4,
        'mass': 4.02643 * u.u,
        'stable': False,
        'half-life': 1.39e-22 * u.s,
    },

    'H-5': {
        'atomic number': 1,
        'mass number': 5,
        'mass': 5.035311 * u.u,
        'stable': False,
        'half-life': '>910 ys',
    },

    'H-6': {
        'atomic number': 1,
        'mass number': 6,
        'mass': 6.04496 * u.u,
        'stable': False,
        'half-life': 2.9e-22 * u.s,
    },

    'H-7': {
        'atomic number': 1,
        'mass number': 7,
        'mass': 7.0527 * u.u,
        'stable': False,
        'half-life': '500# ys',
    },

    'He-3': {
        'atomic number': 2,
        'mass number': 3,
        'mass': 3.0160293201 * u.u,
        'stable': True,
        'abundance': 1.34e-06,
    },

    'He-4': {
        'atomic number': 2,
        'mass number': 4,
        'mass': 4.00260325413 * u.u,
        'stable': True,
        'abundance': 0.99999866,
    },

    'He-5': {
        'atomic number': 2,
        'mass number': 5,
        'mass': 5.012057 * u.u,
        'stable': False,
        'half-life': 7e-22 * u.s,
    },

    'He-6': {
        'atomic number': 2,
        'mass number': 6,
        'mass': 6.018885891 * u.u,
        'stable': False,
        'half-life': 0.80692 * u.s,
    },

    'He-7': {
        'atomic number': 2,
        'mass number': 7,
        'mass': 7.0279907 * u.u,
        'stable': False,
        'half-life': 2.51e-21 * u.s,
    },

    'He-8': {
        'atomic number': 2,
        'mass number': 8,
        'mass': 8.03393439 * u.u,
        'stable': False,
        'half-life': 0.1191 * u.s,
    },

    'He-9': {
        'atomic number': 2,
        'mass number': 9,
        'mass': 9.043946 * u.u,
        'stable': False,
        'half-life': 2.5e-21 * u.s,
    },

    'He-10': {
        'atomic number': 2,
        'mass number': 10,
        'mass': 10.05279 * u.u,
        'stable': False,
        'half-life': 3.1e-21 * u.s,
    },

    'Li-3': {
        'atomic number': 3,
        'mass number': 3,
        'mass': 3.0308 * u.u,
        'stable': False,
    },

    'Li-4': {
        'atomic number': 3,
        'mass number': 4,
        'mass': 4.02719 * u.u,
        'stable': False,
        'half-life': 9.1e-23 * u.s,
    },

    'Li-5': {
        'atomic number': 3,
        'mass number': 5,
        'mass': 5.012538 * u.u,
        'stable': False,
        'half-life': 3.7e-22 * u.s,
    },

    'Li-6': {
        'atomic number': 3,
        'mass number': 6,
        'mass': 6.0151228874 * u.u,
        'stable': True,
        'abundance': 0.0759,
    },

    'Li-7': {
        'atomic number': 3,
        'mass number': 7,
        'mass': 7.0160034366 * u.u,
        'stable': True,
        'abundance': 0.9241,
    },

    'Li-8': {
        'atomic number': 3,
        'mass number': 8,
        'mass': 8.022486246 * u.u,
        'stable': False,
        'half-life': 0.8394 * u.s,
    },

    'Li-9': {
        'atomic number': 3,
        'mass number': 9,
        'mass': 9.02679019 * u.u,
        'stable': False,
        'half-life': 0.1783 * u.s,
    },

    'Li-10': {
        'atomic number': 3,
        'mass number': 10,
        'mass': 10.035483 * u.u,
        'stable': False,
        'half-life': 2e-21 * u.s,
    },

    'Li-11': {
        'atomic number': 3,
        'mass number': 11,
        'mass': 11.04372358 * u.u,
        'stable': False,
        'half-life': 0.00875 * u.s,
    },

    'Li-12': {
        'atomic number': 3,
        'mass number': 12,
        'mass': 12.052517 * u.u,
        'stable': False,
        'half-life': '<10 ns',
    },

    'Li-13': {
        'atomic number': 3,
        'mass number': 13,
        'mass': 13.06263 * u.u,
        'stable': False,
        'half-life': 3.3e-21 * u.s,
    },

    'Be-5': {
        'atomic number': 4,
        'mass number': 5,
        'mass': 5.0399 * u.u,
        'stable': False,
    },

    'Be-6': {
        'atomic number': 4,
        'mass number': 6,
        'mass': 6.0197264 * u.u,
        'stable': False,
        'half-life': 5e-21 * u.s,
    },

    'Be-7': {
        'atomic number': 4,
        'mass number': 7,
        'mass': 7.016928717 * u.u,
        'stable': False,
        'half-life': 4598208.0 * u.s,
    },

    'Be-8': {
        'atomic number': 4,
        'mass number': 8,
        'mass': 8.005305102 * u.u,
        'stable': False,
        'half-life': 6.7e-17 * u.s,
    },

    'Be-9': {
        'atomic number': 4,
        'mass number': 9,
        'mass': 9.012183065 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Be-10': {
        'atomic number': 4,
        'mass number': 10,
        'mass': 10.013534695 * u.u,
        'stable': False,
        'half-life': 47650958260000.0 * u.s,
    },

    'Be-11': {
        'atomic number': 4,
        'mass number': 11,
        'mass': 11.02166108 * u.u,
        'stable': False,
        'half-life': 13.76 * u.s,
    },

    'Be-12': {
        'atomic number': 4,
        'mass number': 12,
        'mass': 12.0269221 * u.u,
        'stable': False,
        'half-life': 0.0215 * u.s,
    },

    'Be-13': {
        'atomic number': 4,
        'mass number': 13,
        'mass': 13.036135 * u.u,
        'stable': False,
        'half-life': 1e-21 * u.s,
    },

    'Be-14': {
        'atomic number': 4,
        'mass number': 14,
        'mass': 14.04289 * u.u,
        'stable': False,
        'half-life': 0.00435 * u.s,
    },

    'Be-15': {
        'atomic number': 4,
        'mass number': 15,
        'mass': 15.05342 * u.u,
        'stable': False,
        'half-life': 7.9e-22 * u.s,
    },

    'Be-16': {
        'atomic number': 4,
        'mass number': 16,
        'mass': 16.06167 * u.u,
        'stable': False,
        'half-life': 6.5e-22 * u.s,
    },

    'B-6': {
        'atomic number': 5,
        'mass number': 6,
        'mass': 6.0508 * u.u,
        'stable': False,
    },

    'B-7': {
        'atomic number': 5,
        'mass number': 7,
        'mass': 7.029712 * u.u,
        'stable': False,
        'half-life': 5.7e-22 * u.s,
    },

    'B-8': {
        'atomic number': 5,
        'mass number': 8,
        'mass': 8.0246073 * u.u,
        'stable': False,
        'half-life': 0.77 * u.s,
    },

    'B-9': {
        'atomic number': 5,
        'mass number': 9,
        'mass': 9.01332965 * u.u,
        'stable': False,
        'half-life': 8e-19 * u.s,
    },

    'B-10': {
        'atomic number': 5,
        'mass number': 10,
        'mass': 10.01293695 * u.u,
        'stable': True,
        'abundance': 0.199,
    },

    'B-11': {
        'atomic number': 5,
        'mass number': 11,
        'mass': 11.00930536 * u.u,
        'stable': True,
        'abundance': 0.801,
    },

    'B-12': {
        'atomic number': 5,
        'mass number': 12,
        'mass': 12.0143527 * u.u,
        'stable': False,
        'half-life': 0.0202 * u.s,
    },

    'B-13': {
        'atomic number': 5,
        'mass number': 13,
        'mass': 13.0177802 * u.u,
        'stable': False,
        'half-life': 0.01733 * u.s,
    },

    'B-14': {
        'atomic number': 5,
        'mass number': 14,
        'mass': 14.025404 * u.u,
        'stable': False,
        'half-life': 0.0125 * u.s,
    },

    'B-15': {
        'atomic number': 5,
        'mass number': 15,
        'mass': 15.031088 * u.u,
        'stable': False,
        'half-life': 0.00993 * u.s,
    },

    'B-16': {
        'atomic number': 5,
        'mass number': 16,
        'mass': 16.039842 * u.u,
        'stable': False,
        'half-life': '>4.6 zs',
    },

    'B-17': {
        'atomic number': 5,
        'mass number': 17,
        'mass': 17.04699 * u.u,
        'stable': False,
        'half-life': 0.00508 * u.s,
    },

    'B-18': {
        'atomic number': 5,
        'mass number': 18,
        'mass': 18.05566 * u.u,
        'stable': False,
        'half-life': '<26 ns',
    },

    'B-19': {
        'atomic number': 5,
        'mass number': 19,
        'mass': 19.0631 * u.u,
        'stable': False,
        'half-life': 0.00292 * u.s,
    },

    'B-20': {
        'atomic number': 5,
        'mass number': 20,
        'mass': 20.07207 * u.u,
        'stable': False,
    },

    'B-21': {
        'atomic number': 5,
        'mass number': 21,
        'mass': 21.08129 * u.u,
        'stable': False,
    },

    'C-8': {
        'atomic number': 6,
        'mass number': 8,
        'mass': 8.037643 * u.u,
        'stable': False,
        'half-life': 3.5e-21 * u.s,
    },

    'C-9': {
        'atomic number': 6,
        'mass number': 9,
        'mass': 9.0310372 * u.u,
        'stable': False,
        'half-life': 0.1265 * u.s,
    },

    'C-10': {
        'atomic number': 6,
        'mass number': 10,
        'mass': 10.01685331 * u.u,
        'stable': False,
        'half-life': 19.3009 * u.s,
    },

    'C-11': {
        'atomic number': 6,
        'mass number': 11,
        'mass': 11.0114336 * u.u,
        'stable': False,
        'half-life': 1221.84 * u.s,
    },

    'C-12': {
        'atomic number': 6,
        'mass number': 12,
        'mass': 12.0 * u.u,
        'stable': True,
        'abundance': 0.9893,
    },

    'C-13': {
        'atomic number': 6,
        'mass number': 13,
        'mass': 13.00335483507 * u.u,
        'stable': True,
        'abundance': 0.0107,
    },

    'C-14': {
        'atomic number': 6,
        'mass number': 14,
        'mass': 14.0032419884 * u.u,
        'stable': False,
        'half-life': 180825048000.0 * u.s,
    },

    'C-15': {
        'atomic number': 6,
        'mass number': 15,
        'mass': 15.01059926 * u.u,
        'stable': False,
        'half-life': 2.449 * u.s,
    },

    'C-16': {
        'atomic number': 6,
        'mass number': 16,
        'mass': 16.0147013 * u.u,
        'stable': False,
        'half-life': 0.747 * u.s,
    },

    'C-17': {
        'atomic number': 6,
        'mass number': 17,
        'mass': 17.022577 * u.u,
        'stable': False,
        'half-life': 0.193 * u.s,
    },

    'C-18': {
        'atomic number': 6,
        'mass number': 18,
        'mass': 18.026751 * u.u,
        'stable': False,
        'half-life': 0.092 * u.s,
    },

    'C-19': {
        'atomic number': 6,
        'mass number': 19,
        'mass': 19.0348 * u.u,
        'stable': False,
        'half-life': 0.0462 * u.s,
    },

    'C-20': {
        'atomic number': 6,
        'mass number': 20,
        'mass': 20.04032 * u.u,
        'stable': False,
        'half-life': 0.016 * u.s,
    },

    'C-21': {
        'atomic number': 6,
        'mass number': 21,
        'mass': 21.049 * u.u,
        'stable': False,
    },

    'C-22': {
        'atomic number': 6,
        'mass number': 22,
        'mass': 22.05753 * u.u,
        'stable': False,
        'half-life': 0.0062 * u.s,
    },

    'C-23': {
        'atomic number': 6,
        'mass number': 23,
        'mass': 23.0689 * u.u,
        'stable': False,
    },

    'N-10': {
        'atomic number': 7,
        'mass number': 10,
        'mass': 10.04165 * u.u,
        'stable': False,
        'half-life': 2e-22 * u.s,
    },

    'N-11': {
        'atomic number': 7,
        'mass number': 11,
        'mass': 11.026091 * u.u,
        'stable': False,
        'half-life': 5.5e-22 * u.s,
    },

    'N-12': {
        'atomic number': 7,
        'mass number': 12,
        'mass': 12.0186132 * u.u,
        'stable': False,
        'half-life': 0.011 * u.s,
    },

    'N-13': {
        'atomic number': 7,
        'mass number': 13,
        'mass': 13.00573861 * u.u,
        'stable': False,
        'half-life': 597.9 * u.s,
    },

    'N-14': {
        'atomic number': 7,
        'mass number': 14,
        'mass': 14.00307400443 * u.u,
        'stable': True,
        'abundance': 0.99636,
    },

    'N-15': {
        'atomic number': 7,
        'mass number': 15,
        'mass': 15.00010889888 * u.u,
        'stable': True,
        'abundance': 0.00364,
    },

    'N-16': {
        'atomic number': 7,
        'mass number': 16,
        'mass': 16.0061019 * u.u,
        'stable': False,
        'half-life': 7.13 * u.s,
    },

    'N-17': {
        'atomic number': 7,
        'mass number': 17,
        'mass': 17.008449 * u.u,
        'stable': False,
        'half-life': 4.173 * u.s,
    },

    'N-18': {
        'atomic number': 7,
        'mass number': 18,
        'mass': 18.014078 * u.u,
        'stable': False,
        'half-life': 0.6192 * u.s,
    },

    'N-19': {
        'atomic number': 7,
        'mass number': 19,
        'mass': 19.017022 * u.u,
        'stable': False,
        'half-life': 0.336 * u.s,
    },

    'N-20': {
        'atomic number': 7,
        'mass number': 20,
        'mass': 20.023366 * u.u,
        'stable': False,
        'half-life': 0.136 * u.s,
    },

    'N-21': {
        'atomic number': 7,
        'mass number': 21,
        'mass': 21.02711 * u.u,
        'stable': False,
        'half-life': 0.084 * u.s,
    },

    'N-22': {
        'atomic number': 7,
        'mass number': 22,
        'mass': 22.03439 * u.u,
        'stable': False,
        'half-life': 0.023 * u.s,
    },

    'N-23': {
        'atomic number': 7,
        'mass number': 23,
        'mass': 23.04114 * u.u,
        'stable': False,
        'half-life': 0.0139 * u.s,
    },

    'N-24': {
        'atomic number': 7,
        'mass number': 24,
        'mass': 24.05039 * u.u,
        'stable': False,
    },

    'N-25': {
        'atomic number': 7,
        'mass number': 25,
        'mass': 25.0601 * u.u,
        'stable': False,
    },

    'O-12': {
        'atomic number': 8,
        'mass number': 12,
        'mass': 12.034262 * u.u,
        'stable': False,
        'half-life': '>6.3 zs',
    },

    'O-13': {
        'atomic number': 8,
        'mass number': 13,
        'mass': 13.024815 * u.u,
        'stable': False,
        'half-life': 0.00858 * u.s,
    },

    'O-14': {
        'atomic number': 8,
        'mass number': 14,
        'mass': 14.00859636 * u.u,
        'stable': False,
        'half-life': 70.62 * u.s,
    },

    'O-15': {
        'atomic number': 8,
        'mass number': 15,
        'mass': 15.00306562 * u.u,
        'stable': False,
        'half-life': 122.24 * u.s,
    },

    'O-16': {
        'atomic number': 8,
        'mass number': 16,
        'mass': 15.99491461957 * u.u,
        'stable': True,
        'abundance': 0.99757,
    },

    'O-17': {
        'atomic number': 8,
        'mass number': 17,
        'mass': 16.9991317565 * u.u,
        'stable': True,
        'abundance': 0.00038,
    },

    'O-18': {
        'atomic number': 8,
        'mass number': 18,
        'mass': 17.99915961286 * u.u,
        'stable': True,
        'abundance': 0.00205,
    },

    'O-19': {
        'atomic number': 8,
        'mass number': 19,
        'mass': 19.003578 * u.u,
        'stable': False,
        'half-life': 26.47 * u.s,
    },

    'O-20': {
        'atomic number': 8,
        'mass number': 20,
        'mass': 20.00407535 * u.u,
        'stable': False,
        'half-life': 13.51 * u.s,
    },

    'O-21': {
        'atomic number': 8,
        'mass number': 21,
        'mass': 21.008655 * u.u,
        'stable': False,
        'half-life': 3.42 * u.s,
    },

    'O-22': {
        'atomic number': 8,
        'mass number': 22,
        'mass': 22.009966 * u.u,
        'stable': False,
        'half-life': 2.25 * u.s,
    },

    'O-23': {
        'atomic number': 8,
        'mass number': 23,
        'mass': 23.015696 * u.u,
        'stable': False,
        'half-life': 0.097 * u.s,
    },

    'O-24': {
        'atomic number': 8,
        'mass number': 24,
        'mass': 24.01986 * u.u,
        'stable': False,
        'half-life': 0.0774 * u.s,
    },

    'O-25': {
        'atomic number': 8,
        'mass number': 25,
        'mass': 25.02936 * u.u,
        'stable': False,
        'half-life': 5.18e-21 * u.s,
    },

    'O-26': {
        'atomic number': 8,
        'mass number': 26,
        'mass': 26.03729 * u.u,
        'stable': False,
        'half-life': 4.2e-12 * u.s,
    },

    'O-27': {
        'atomic number': 8,
        'mass number': 27,
        'mass': 27.04772 * u.u,
        'stable': False,
    },

    'O-28': {
        'atomic number': 8,
        'mass number': 28,
        'mass': 28.05591 * u.u,
        'stable': False,
    },

    'F-14': {
        'atomic number': 9,
        'mass number': 14,
        'mass': 14.034315 * u.u,
        'stable': False,
        'half-life': 5e-22 * u.s,
    },

    'F-15': {
        'atomic number': 9,
        'mass number': 15,
        'mass': 15.018043 * u.u,
        'stable': False,
        'half-life': 1.1e-21 * u.s,
    },

    'F-16': {
        'atomic number': 9,
        'mass number': 16,
        'mass': 16.0114657 * u.u,
        'stable': False,
        'half-life': 1.1e-20 * u.s,
    },

    'F-17': {
        'atomic number': 9,
        'mass number': 17,
        'mass': 17.00209524 * u.u,
        'stable': False,
        'half-life': 64.37 * u.s,
    },

    'F-18': {
        'atomic number': 9,
        'mass number': 18,
        'mass': 18.00093733 * u.u,
        'stable': False,
        'half-life': 6586.236 * u.s,
    },

    'F-19': {
        'atomic number': 9,
        'mass number': 19,
        'mass': 18.99840316273 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'F-20': {
        'atomic number': 9,
        'mass number': 20,
        'mass': 19.999981252 * u.u,
        'stable': False,
        'half-life': 11.163 * u.s,
    },

    'F-21': {
        'atomic number': 9,
        'mass number': 21,
        'mass': 20.9999489 * u.u,
        'stable': False,
        'half-life': 4.158 * u.s,
    },

    'F-22': {
        'atomic number': 9,
        'mass number': 22,
        'mass': 22.002999 * u.u,
        'stable': False,
        'half-life': 4.23 * u.s,
    },

    'F-23': {
        'atomic number': 9,
        'mass number': 23,
        'mass': 23.003557 * u.u,
        'stable': False,
        'half-life': 2.23 * u.s,
    },

    'F-24': {
        'atomic number': 9,
        'mass number': 24,
        'mass': 24.008115 * u.u,
        'stable': False,
        'half-life': 0.384 * u.s,
    },

    'F-25': {
        'atomic number': 9,
        'mass number': 25,
        'mass': 25.012199 * u.u,
        'stable': False,
        'half-life': 0.08 * u.s,
    },

    'F-26': {
        'atomic number': 9,
        'mass number': 26,
        'mass': 26.020038 * u.u,
        'stable': False,
        'half-life': 0.0082 * u.s,
    },

    'F-27': {
        'atomic number': 9,
        'mass number': 27,
        'mass': 27.02644 * u.u,
        'stable': False,
        'half-life': 0.0049 * u.s,
    },

    'F-28': {
        'atomic number': 9,
        'mass number': 28,
        'mass': 28.03534 * u.u,
        'stable': False,
        'half-life': 4.6e-20 * u.s,
    },

    'F-29': {
        'atomic number': 9,
        'mass number': 29,
        'mass': 29.04254 * u.u,
        'stable': False,
        'half-life': 0.0025 * u.s,
    },

    'F-30': {
        'atomic number': 9,
        'mass number': 30,
        'mass': 30.05165 * u.u,
        'stable': False,
    },

    'F-31': {
        'atomic number': 9,
        'mass number': 31,
        'mass': 31.05971 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'Ne-16': {
        'atomic number': 10,
        'mass number': 16,
        'mass': 16.02575 * u.u,
        'stable': False,
        'half-life': '>5.7 zs',
    },

    'Ne-17': {
        'atomic number': 10,
        'mass number': 17,
        'mass': 17.01771396 * u.u,
        'stable': False,
        'half-life': 0.1092 * u.s,
    },

    'Ne-18': {
        'atomic number': 10,
        'mass number': 18,
        'mass': 18.0057087 * u.u,
        'stable': False,
        'half-life': 1.6642 * u.s,
    },

    'Ne-19': {
        'atomic number': 10,
        'mass number': 19,
        'mass': 19.00188091 * u.u,
        'stable': False,
        'half-life': 17.274 * u.s,
    },

    'Ne-20': {
        'atomic number': 10,
        'mass number': 20,
        'mass': 19.9924401762 * u.u,
        'stable': True,
        'abundance': 0.9048,
    },

    'Ne-21': {
        'atomic number': 10,
        'mass number': 21,
        'mass': 20.993846685 * u.u,
        'stable': True,
        'abundance': 0.0027,
    },

    'Ne-22': {
        'atomic number': 10,
        'mass number': 22,
        'mass': 21.991385114 * u.u,
        'stable': True,
        'abundance': 0.0925,
    },

    'Ne-23': {
        'atomic number': 10,
        'mass number': 23,
        'mass': 22.99446691 * u.u,
        'stable': False,
        'half-life': 37.14 * u.s,
    },

    'Ne-24': {
        'atomic number': 10,
        'mass number': 24,
        'mass': 23.99361065 * u.u,
        'stable': False,
        'half-life': 202.8 * u.s,
    },

    'Ne-25': {
        'atomic number': 10,
        'mass number': 25,
        'mass': 24.997789 * u.u,
        'stable': False,
        'half-life': 0.602 * u.s,
    },

    'Ne-26': {
        'atomic number': 10,
        'mass number': 26,
        'mass': 26.000515 * u.u,
        'stable': False,
        'half-life': 0.197 * u.s,
    },

    'Ne-27': {
        'atomic number': 10,
        'mass number': 27,
        'mass': 27.007553 * u.u,
        'stable': False,
        'half-life': 0.0315 * u.s,
    },

    'Ne-28': {
        'atomic number': 10,
        'mass number': 28,
        'mass': 28.01212 * u.u,
        'stable': False,
        'half-life': 0.02 * u.s,
    },

    'Ne-29': {
        'atomic number': 10,
        'mass number': 29,
        'mass': 29.01975 * u.u,
        'stable': False,
        'half-life': 0.0147 * u.s,
    },

    'Ne-30': {
        'atomic number': 10,
        'mass number': 30,
        'mass': 30.02473 * u.u,
        'stable': False,
        'half-life': 0.00722 * u.s,
    },

    'Ne-31': {
        'atomic number': 10,
        'mass number': 31,
        'mass': 31.0331 * u.u,
        'stable': False,
        'half-life': 0.0034 * u.s,
    },

    'Ne-32': {
        'atomic number': 10,
        'mass number': 32,
        'mass': 32.03972 * u.u,
        'stable': False,
        'half-life': 0.0035 * u.s,
    },

    'Ne-33': {
        'atomic number': 10,
        'mass number': 33,
        'mass': 33.04938 * u.u,
        'stable': False,
    },

    'Ne-34': {
        'atomic number': 10,
        'mass number': 34,
        'mass': 34.05673 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'Na-18': {
        'atomic number': 11,
        'mass number': 18,
        'mass': 18.02688 * u.u,
        'stable': False,
        'half-life': 1.3e-21 * u.s,
    },

    'Na-19': {
        'atomic number': 11,
        'mass number': 19,
        'mass': 19.01388 * u.u,
        'stable': False,
        'half-life': '>1 as',
    },

    'Na-20': {
        'atomic number': 11,
        'mass number': 20,
        'mass': 20.0073544 * u.u,
        'stable': False,
        'half-life': 0.4479 * u.s,
    },

    'Na-21': {
        'atomic number': 11,
        'mass number': 21,
        'mass': 20.99765469 * u.u,
        'stable': False,
        'half-life': 22.422 * u.s,
    },

    'Na-22': {
        'atomic number': 11,
        'mass number': 22,
        'mass': 21.99443741 * u.u,
        'stable': False,
        'half-life': 82163808.0 * u.s,
    },

    'Na-23': {
        'atomic number': 11,
        'mass number': 23,
        'mass': 22.989769282 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Na-24': {
        'atomic number': 11,
        'mass number': 24,
        'mass': 23.99096295 * u.u,
        'stable': False,
        'half-life': 53824.32 * u.s,
    },

    'Na-25': {
        'atomic number': 11,
        'mass number': 25,
        'mass': 24.989954 * u.u,
        'stable': False,
        'half-life': 59.1 * u.s,
    },

    'Na-26': {
        'atomic number': 11,
        'mass number': 26,
        'mass': 25.9926346 * u.u,
        'stable': False,
        'half-life': 1.07128 * u.s,
    },

    'Na-27': {
        'atomic number': 11,
        'mass number': 27,
        'mass': 26.9940765 * u.u,
        'stable': False,
        'half-life': 0.301 * u.s,
    },

    'Na-28': {
        'atomic number': 11,
        'mass number': 28,
        'mass': 27.998939 * u.u,
        'stable': False,
        'half-life': 0.0305 * u.s,
    },

    'Na-29': {
        'atomic number': 11,
        'mass number': 29,
        'mass': 29.0028771 * u.u,
        'stable': False,
        'half-life': 0.0441 * u.s,
    },

    'Na-30': {
        'atomic number': 11,
        'mass number': 30,
        'mass': 30.0090979 * u.u,
        'stable': False,
        'half-life': 0.0484 * u.s,
    },

    'Na-31': {
        'atomic number': 11,
        'mass number': 31,
        'mass': 31.013163 * u.u,
        'stable': False,
        'half-life': 0.01735 * u.s,
    },

    'Na-32': {
        'atomic number': 11,
        'mass number': 32,
        'mass': 32.02019 * u.u,
        'stable': False,
        'half-life': 0.0129 * u.s,
    },

    'Na-33': {
        'atomic number': 11,
        'mass number': 33,
        'mass': 33.02573 * u.u,
        'stable': False,
        'half-life': 0.0082 * u.s,
    },

    'Na-34': {
        'atomic number': 11,
        'mass number': 34,
        'mass': 34.03359 * u.u,
        'stable': False,
        'half-life': 0.0055 * u.s,
    },

    'Na-35': {
        'atomic number': 11,
        'mass number': 35,
        'mass': 35.04062 * u.u,
        'stable': False,
        'half-life': 0.0015 * u.s,
    },

    'Na-36': {
        'atomic number': 11,
        'mass number': 36,
        'mass': 36.04929 * u.u,
        'stable': False,
    },

    'Na-37': {
        'atomic number': 11,
        'mass number': 37,
        'mass': 37.05705 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'Mg-19': {
        'atomic number': 12,
        'mass number': 19,
        'mass': 19.034169 * u.u,
        'stable': False,
        'half-life': 5e-12 * u.s,
    },

    'Mg-20': {
        'atomic number': 12,
        'mass number': 20,
        'mass': 20.01885 * u.u,
        'stable': False,
        'half-life': 0.093 * u.s,
    },

    'Mg-21': {
        'atomic number': 12,
        'mass number': 21,
        'mass': 21.011716 * u.u,
        'stable': False,
        'half-life': 0.1186 * u.s,
    },

    'Mg-22': {
        'atomic number': 12,
        'mass number': 22,
        'mass': 21.99957065 * u.u,
        'stable': False,
        'half-life': 3.8755 * u.s,
    },

    'Mg-23': {
        'atomic number': 12,
        'mass number': 23,
        'mass': 22.99412421 * u.u,
        'stable': False,
        'half-life': 11.317 * u.s,
    },

    'Mg-24': {
        'atomic number': 12,
        'mass number': 24,
        'mass': 23.985041697 * u.u,
        'stable': True,
        'abundance': 0.7899,
    },

    'Mg-25': {
        'atomic number': 12,
        'mass number': 25,
        'mass': 24.985836976 * u.u,
        'stable': True,
        'abundance': 0.1,
    },

    'Mg-26': {
        'atomic number': 12,
        'mass number': 26,
        'mass': 25.982592968 * u.u,
        'stable': True,
        'abundance': 0.1101,
    },

    'Mg-27': {
        'atomic number': 12,
        'mass number': 27,
        'mass': 26.984340624 * u.u,
        'stable': False,
        'half-life': 566.1 * u.s,
    },

    'Mg-28': {
        'atomic number': 12,
        'mass number': 28,
        'mass': 27.9838767 * u.u,
        'stable': False,
        'half-life': 75294.0 * u.s,
    },

    'Mg-29': {
        'atomic number': 12,
        'mass number': 29,
        'mass': 28.988617 * u.u,
        'stable': False,
        'half-life': 1.3 * u.s,
    },

    'Mg-30': {
        'atomic number': 12,
        'mass number': 30,
        'mass': 29.9904629 * u.u,
        'stable': False,
        'half-life': 0.313 * u.s,
    },

    'Mg-31': {
        'atomic number': 12,
        'mass number': 31,
        'mass': 30.996648 * u.u,
        'stable': False,
        'half-life': 0.236 * u.s,
    },

    'Mg-32': {
        'atomic number': 12,
        'mass number': 32,
        'mass': 31.9991102 * u.u,
        'stable': False,
        'half-life': 0.086 * u.s,
    },

    'Mg-33': {
        'atomic number': 12,
        'mass number': 33,
        'mass': 33.0053271 * u.u,
        'stable': False,
        'half-life': 0.0905 * u.s,
    },

    'Mg-34': {
        'atomic number': 12,
        'mass number': 34,
        'mass': 34.008935 * u.u,
        'stable': False,
        'half-life': 0.02 * u.s,
    },

    'Mg-35': {
        'atomic number': 12,
        'mass number': 35,
        'mass': 35.01679 * u.u,
        'stable': False,
        'half-life': 0.07 * u.s,
    },

    'Mg-36': {
        'atomic number': 12,
        'mass number': 36,
        'mass': 36.02188 * u.u,
        'stable': False,
        'half-life': 0.0039 * u.s,
    },

    'Mg-37': {
        'atomic number': 12,
        'mass number': 37,
        'mass': 37.03037 * u.u,
        'stable': False,
        'half-life': 0.008 * u.s,
    },

    'Mg-38': {
        'atomic number': 12,
        'mass number': 38,
        'mass': 38.03658 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'Mg-39': {
        'atomic number': 12,
        'mass number': 39,
        'mass': 39.04538 * u.u,
        'stable': False,
    },

    'Mg-40': {
        'atomic number': 12,
        'mass number': 40,
        'mass': 40.05218 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'Al-21': {
        'atomic number': 13,
        'mass number': 21,
        'mass': 21.02897 * u.u,
        'stable': False,
    },

    'Al-22': {
        'atomic number': 13,
        'mass number': 22,
        'mass': 22.01954 * u.u,
        'stable': False,
        'half-life': 0.0911 * u.s,
    },

    'Al-23': {
        'atomic number': 13,
        'mass number': 23,
        'mass': 23.00724435 * u.u,
        'stable': False,
        'half-life': 0.47 * u.s,
    },

    'Al-24': {
        'atomic number': 13,
        'mass number': 24,
        'mass': 23.9999489 * u.u,
        'stable': False,
        'half-life': 2.053 * u.s,
    },

    'Al-25': {
        'atomic number': 13,
        'mass number': 25,
        'mass': 24.9904281 * u.u,
        'stable': False,
        'half-life': 7.183 * u.s,
    },

    'Al-26': {
        'atomic number': 13,
        'mass number': 26,
        'mass': 25.986891904 * u.u,
        'stable': False,
        'half-life': 22626315942000.0 * u.s,
    },

    'Al-27': {
        'atomic number': 13,
        'mass number': 27,
        'mass': 26.98153853 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Al-28': {
        'atomic number': 13,
        'mass number': 28,
        'mass': 27.98191021 * u.u,
        'stable': False,
        'half-life': 134.7 * u.s,
    },

    'Al-29': {
        'atomic number': 13,
        'mass number': 29,
        'mass': 28.9804565 * u.u,
        'stable': False,
        'half-life': 393.6 * u.s,
    },

    'Al-30': {
        'atomic number': 13,
        'mass number': 30,
        'mass': 29.98296 * u.u,
        'stable': False,
        'half-life': 3.62 * u.s,
    },

    'Al-31': {
        'atomic number': 13,
        'mass number': 31,
        'mass': 30.983945 * u.u,
        'stable': False,
        'half-life': 0.644 * u.s,
    },

    'Al-32': {
        'atomic number': 13,
        'mass number': 32,
        'mass': 31.988085 * u.u,
        'stable': False,
        'half-life': 0.033 * u.s,
    },

    'Al-33': {
        'atomic number': 13,
        'mass number': 33,
        'mass': 32.990909 * u.u,
        'stable': False,
        'half-life': 0.0417 * u.s,
    },

    'Al-34': {
        'atomic number': 13,
        'mass number': 34,
        'mass': 33.996705 * u.u,
        'stable': False,
        'half-life': 0.0563 * u.s,
    },

    'Al-35': {
        'atomic number': 13,
        'mass number': 35,
        'mass': 34.999764 * u.u,
        'stable': False,
        'half-life': 0.0372 * u.s,
    },

    'Al-36': {
        'atomic number': 13,
        'mass number': 36,
        'mass': 36.00639 * u.u,
        'stable': False,
        'half-life': 0.09 * u.s,
    },

    'Al-37': {
        'atomic number': 13,
        'mass number': 37,
        'mass': 37.01053 * u.u,
        'stable': False,
        'half-life': 0.0115 * u.s,
    },

    'Al-38': {
        'atomic number': 13,
        'mass number': 38,
        'mass': 38.0174 * u.u,
        'stable': False,
        'half-life': 0.009 * u.s,
    },

    'Al-39': {
        'atomic number': 13,
        'mass number': 39,
        'mass': 39.02254 * u.u,
        'stable': False,
        'half-life': 0.0076 * u.s,
    },

    'Al-40': {
        'atomic number': 13,
        'mass number': 40,
        'mass': 40.03003 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Al-41': {
        'atomic number': 13,
        'mass number': 41,
        'mass': 41.03638 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'Al-42': {
        'atomic number': 13,
        'mass number': 42,
        'mass': 42.04384 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'Al-43': {
        'atomic number': 13,
        'mass number': 43,
        'mass': 43.05147 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'Si-22': {
        'atomic number': 14,
        'mass number': 22,
        'mass': 22.03579 * u.u,
        'stable': False,
        'half-life': 0.029 * u.s,
    },

    'Si-23': {
        'atomic number': 14,
        'mass number': 23,
        'mass': 23.02544 * u.u,
        'stable': False,
        'half-life': 0.0423 * u.s,
    },

    'Si-24': {
        'atomic number': 14,
        'mass number': 24,
        'mass': 24.011535 * u.u,
        'stable': False,
        'half-life': 0.14 * u.s,
    },

    'Si-25': {
        'atomic number': 14,
        'mass number': 25,
        'mass': 25.004109 * u.u,
        'stable': False,
        'half-life': 0.22 * u.s,
    },

    'Si-26': {
        'atomic number': 14,
        'mass number': 26,
        'mass': 25.99233384 * u.u,
        'stable': False,
        'half-life': 2.2453 * u.s,
    },

    'Si-27': {
        'atomic number': 14,
        'mass number': 27,
        'mass': 26.98670481 * u.u,
        'stable': False,
        'half-life': 4.15 * u.s,
    },

    'Si-28': {
        'atomic number': 14,
        'mass number': 28,
        'mass': 27.97692653465 * u.u,
        'stable': True,
        'abundance': 0.92223,
    },

    'Si-29': {
        'atomic number': 14,
        'mass number': 29,
        'mass': 28.9764946649 * u.u,
        'stable': True,
        'abundance': 0.04685,
    },

    'Si-30': {
        'atomic number': 14,
        'mass number': 30,
        'mass': 29.973770136 * u.u,
        'stable': True,
        'abundance': 0.03092,
    },

    'Si-31': {
        'atomic number': 14,
        'mass number': 31,
        'mass': 30.975363194 * u.u,
        'stable': False,
        'half-life': 9441.6 * u.s,
    },

    'Si-32': {
        'atomic number': 14,
        'mass number': 32,
        'mass': 31.97415154 * u.u,
        'stable': False,
        'half-life': 4828209678.0 * u.s,
    },

    'Si-33': {
        'atomic number': 14,
        'mass number': 33,
        'mass': 32.97797696 * u.u,
        'stable': False,
        'half-life': 6.18 * u.s,
    },

    'Si-34': {
        'atomic number': 14,
        'mass number': 34,
        'mass': 33.978576 * u.u,
        'stable': False,
        'half-life': 2.77 * u.s,
    },

    'Si-35': {
        'atomic number': 14,
        'mass number': 35,
        'mass': 34.984583 * u.u,
        'stable': False,
        'half-life': 0.78 * u.s,
    },

    'Si-36': {
        'atomic number': 14,
        'mass number': 36,
        'mass': 35.986695 * u.u,
        'stable': False,
        'half-life': 0.45 * u.s,
    },

    'Si-37': {
        'atomic number': 14,
        'mass number': 37,
        'mass': 36.992921 * u.u,
        'stable': False,
        'half-life': 0.09 * u.s,
    },

    'Si-38': {
        'atomic number': 14,
        'mass number': 38,
        'mass': 37.995523 * u.u,
        'stable': False,
        'half-life': '90# ms',
    },

    'Si-39': {
        'atomic number': 14,
        'mass number': 39,
        'mass': 39.002491 * u.u,
        'stable': False,
        'half-life': 0.0475 * u.s,
    },

    'Si-40': {
        'atomic number': 14,
        'mass number': 40,
        'mass': 40.00583 * u.u,
        'stable': False,
        'half-life': 0.033 * u.s,
    },

    'Si-41': {
        'atomic number': 14,
        'mass number': 41,
        'mass': 41.01301 * u.u,
        'stable': False,
        'half-life': 0.02 * u.s,
    },

    'Si-42': {
        'atomic number': 14,
        'mass number': 42,
        'mass': 42.01778 * u.u,
        'stable': False,
        'half-life': 0.0125 * u.s,
    },

    'Si-43': {
        'atomic number': 14,
        'mass number': 43,
        'mass': 43.0248 * u.u,
        'stable': False,
        'half-life': '15# ms',
    },

    'Si-44': {
        'atomic number': 14,
        'mass number': 44,
        'mass': 44.03061 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Si-45': {
        'atomic number': 14,
        'mass number': 45,
        'mass': 45.03995 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'P-24': {
        'atomic number': 15,
        'mass number': 24,
        'mass': 24.03577 * u.u,
        'stable': False,
    },

    'P-25': {
        'atomic number': 15,
        'mass number': 25,
        'mass': 25.02119 * u.u,
        'stable': False,
    },

    'P-26': {
        'atomic number': 15,
        'mass number': 26,
        'mass': 26.01178 * u.u,
        'stable': False,
        'half-life': 0.0437 * u.s,
    },

    'P-27': {
        'atomic number': 15,
        'mass number': 27,
        'mass': 26.999224 * u.u,
        'stable': False,
        'half-life': 0.26 * u.s,
    },

    'P-28': {
        'atomic number': 15,
        'mass number': 28,
        'mass': 27.9923266 * u.u,
        'stable': False,
        'half-life': 0.2703 * u.s,
    },

    'P-29': {
        'atomic number': 15,
        'mass number': 29,
        'mass': 28.98180079 * u.u,
        'stable': False,
        'half-life': 4.142 * u.s,
    },

    'P-30': {
        'atomic number': 15,
        'mass number': 30,
        'mass': 29.97831375 * u.u,
        'stable': False,
        'half-life': 149.88 * u.s,
    },

    'P-31': {
        'atomic number': 15,
        'mass number': 31,
        'mass': 30.97376199842 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'P-32': {
        'atomic number': 15,
        'mass number': 32,
        'mass': 31.973907643 * u.u,
        'stable': False,
        'half-life': 1232323.2 * u.s,
    },

    'P-33': {
        'atomic number': 15,
        'mass number': 33,
        'mass': 32.9717257 * u.u,
        'stable': False,
        'half-life': 2190240.0 * u.s,
    },

    'P-34': {
        'atomic number': 15,
        'mass number': 34,
        'mass': 33.97364589 * u.u,
        'stable': False,
        'half-life': 12.43 * u.s,
    },

    'P-35': {
        'atomic number': 15,
        'mass number': 35,
        'mass': 34.9733141 * u.u,
        'stable': False,
        'half-life': 47.3 * u.s,
    },

    'P-36': {
        'atomic number': 15,
        'mass number': 36,
        'mass': 35.97826 * u.u,
        'stable': False,
        'half-life': 5.6 * u.s,
    },

    'P-37': {
        'atomic number': 15,
        'mass number': 37,
        'mass': 36.979607 * u.u,
        'stable': False,
        'half-life': 2.31 * u.s,
    },

    'P-38': {
        'atomic number': 15,
        'mass number': 38,
        'mass': 37.984252 * u.u,
        'stable': False,
        'half-life': 0.64 * u.s,
    },

    'P-39': {
        'atomic number': 15,
        'mass number': 39,
        'mass': 38.986227 * u.u,
        'stable': False,
        'half-life': 0.282 * u.s,
    },

    'P-40': {
        'atomic number': 15,
        'mass number': 40,
        'mass': 39.99133 * u.u,
        'stable': False,
        'half-life': 0.15 * u.s,
    },

    'P-41': {
        'atomic number': 15,
        'mass number': 41,
        'mass': 40.994654 * u.u,
        'stable': False,
        'half-life': 0.101 * u.s,
    },

    'P-42': {
        'atomic number': 15,
        'mass number': 42,
        'mass': 42.00108 * u.u,
        'stable': False,
        'half-life': 0.0485 * u.s,
    },

    'P-43': {
        'atomic number': 15,
        'mass number': 43,
        'mass': 43.00502 * u.u,
        'stable': False,
        'half-life': 0.0358 * u.s,
    },

    'P-44': {
        'atomic number': 15,
        'mass number': 44,
        'mass': 44.01121 * u.u,
        'stable': False,
        'half-life': 0.0185 * u.s,
    },

    'P-45': {
        'atomic number': 15,
        'mass number': 45,
        'mass': 45.01645 * u.u,
        'stable': False,
        'half-life': '8# ms',
    },

    'P-46': {
        'atomic number': 15,
        'mass number': 46,
        'mass': 46.02446 * u.u,
        'stable': False,
        'half-life': '4# ms',
    },

    'P-47': {
        'atomic number': 15,
        'mass number': 47,
        'mass': 47.03139 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'S-26': {
        'atomic number': 16,
        'mass number': 26,
        'mass': 26.02907 * u.u,
        'stable': False,
    },

    'S-27': {
        'atomic number': 16,
        'mass number': 27,
        'mass': 27.01828 * u.u,
        'stable': False,
        'half-life': 0.0155 * u.s,
    },

    'S-28': {
        'atomic number': 16,
        'mass number': 28,
        'mass': 28.00437 * u.u,
        'stable': False,
        'half-life': 0.125 * u.s,
    },

    'S-29': {
        'atomic number': 16,
        'mass number': 29,
        'mass': 28.996611 * u.u,
        'stable': False,
        'half-life': 0.188 * u.s,
    },

    'S-30': {
        'atomic number': 16,
        'mass number': 30,
        'mass': 29.98490703 * u.u,
        'stable': False,
        'half-life': 1.1759 * u.s,
    },

    'S-31': {
        'atomic number': 16,
        'mass number': 31,
        'mass': 30.97955701 * u.u,
        'stable': False,
        'half-life': 2.5534 * u.s,
    },

    'S-32': {
        'atomic number': 16,
        'mass number': 32,
        'mass': 31.9720711744 * u.u,
        'stable': True,
        'abundance': 0.9499,
    },

    'S-33': {
        'atomic number': 16,
        'mass number': 33,
        'mass': 32.9714589098 * u.u,
        'stable': True,
        'abundance': 0.0075,
    },

    'S-34': {
        'atomic number': 16,
        'mass number': 34,
        'mass': 33.967867004 * u.u,
        'stable': True,
        'abundance': 0.0425,
    },

    'S-35': {
        'atomic number': 16,
        'mass number': 35,
        'mass': 34.96903231 * u.u,
        'stable': False,
        'half-life': 7548768.0 * u.s,
    },

    'S-36': {
        'atomic number': 16,
        'mass number': 36,
        'mass': 35.96708071 * u.u,
        'stable': True,
        'abundance': 0.0001,
    },

    'S-37': {
        'atomic number': 16,
        'mass number': 37,
        'mass': 36.97112551 * u.u,
        'stable': False,
        'half-life': 303.0 * u.s,
    },

    'S-38': {
        'atomic number': 16,
        'mass number': 38,
        'mass': 37.9711633 * u.u,
        'stable': False,
        'half-life': 10218.0 * u.s,
    },

    'S-39': {
        'atomic number': 16,
        'mass number': 39,
        'mass': 38.975134 * u.u,
        'stable': False,
        'half-life': 11.5 * u.s,
    },

    'S-40': {
        'atomic number': 16,
        'mass number': 40,
        'mass': 39.9754826 * u.u,
        'stable': False,
        'half-life': 8.8 * u.s,
    },

    'S-41': {
        'atomic number': 16,
        'mass number': 41,
        'mass': 40.9795935 * u.u,
        'stable': False,
        'half-life': 1.99 * u.s,
    },

    'S-42': {
        'atomic number': 16,
        'mass number': 42,
        'mass': 41.9810651 * u.u,
        'stable': False,
        'half-life': 1.016 * u.s,
    },

    'S-43': {
        'atomic number': 16,
        'mass number': 43,
        'mass': 42.9869076 * u.u,
        'stable': False,
        'half-life': 0.265 * u.s,
    },

    'S-44': {
        'atomic number': 16,
        'mass number': 44,
        'mass': 43.9901188 * u.u,
        'stable': False,
        'half-life': 0.1 * u.s,
    },

    'S-45': {
        'atomic number': 16,
        'mass number': 45,
        'mass': 44.99572 * u.u,
        'stable': False,
        'half-life': 0.068 * u.s,
    },

    'S-46': {
        'atomic number': 16,
        'mass number': 46,
        'mass': 46.00004 * u.u,
        'stable': False,
        'half-life': 0.05 * u.s,
    },

    'S-47': {
        'atomic number': 16,
        'mass number': 47,
        'mass': 47.00795 * u.u,
        'stable': False,
        'half-life': '20# ms',
    },

    'S-48': {
        'atomic number': 16,
        'mass number': 48,
        'mass': 48.0137 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'S-49': {
        'atomic number': 16,
        'mass number': 49,
        'mass': 49.02276 * u.u,
        'stable': False,
    },

    'Cl-28': {
        'atomic number': 17,
        'mass number': 28,
        'mass': 28.02954 * u.u,
        'stable': False,
    },

    'Cl-29': {
        'atomic number': 17,
        'mass number': 29,
        'mass': 29.01478 * u.u,
        'stable': False,
    },

    'Cl-30': {
        'atomic number': 17,
        'mass number': 30,
        'mass': 30.00477 * u.u,
        'stable': False,
    },

    'Cl-31': {
        'atomic number': 17,
        'mass number': 31,
        'mass': 30.992414 * u.u,
        'stable': False,
        'half-life': 0.19 * u.s,
    },

    'Cl-32': {
        'atomic number': 17,
        'mass number': 32,
        'mass': 31.98568464 * u.u,
        'stable': False,
        'half-life': 0.298 * u.s,
    },

    'Cl-33': {
        'atomic number': 17,
        'mass number': 33,
        'mass': 32.97745199 * u.u,
        'stable': False,
        'half-life': 2.5038 * u.s,
    },

    'Cl-34': {
        'atomic number': 17,
        'mass number': 34,
        'mass': 33.973762485 * u.u,
        'stable': False,
        'half-life': 1.5266 * u.s,
    },

    'Cl-35': {
        'atomic number': 17,
        'mass number': 35,
        'mass': 34.968852682 * u.u,
        'stable': True,
        'abundance': 0.7576,
    },

    'Cl-36': {
        'atomic number': 17,
        'mass number': 36,
        'mass': 35.968306809 * u.u,
        'stable': False,
        'half-life': 9508101803800.0 * u.s,
    },

    'Cl-37': {
        'atomic number': 17,
        'mass number': 37,
        'mass': 36.965902602 * u.u,
        'stable': True,
        'abundance': 0.2424,
    },

    'Cl-38': {
        'atomic number': 17,
        'mass number': 38,
        'mass': 37.96801044 * u.u,
        'stable': False,
        'half-life': 2234.4 * u.s,
    },

    'Cl-39': {
        'atomic number': 17,
        'mass number': 39,
        'mass': 38.9680082 * u.u,
        'stable': False,
        'half-life': 3372.0 * u.s,
    },

    'Cl-40': {
        'atomic number': 17,
        'mass number': 40,
        'mass': 39.970415 * u.u,
        'stable': False,
        'half-life': 81.0 * u.s,
    },

    'Cl-41': {
        'atomic number': 17,
        'mass number': 41,
        'mass': 40.970685 * u.u,
        'stable': False,
        'half-life': 38.4 * u.s,
    },

    'Cl-42': {
        'atomic number': 17,
        'mass number': 42,
        'mass': 41.97325 * u.u,
        'stable': False,
        'half-life': 6.8 * u.s,
    },

    'Cl-43': {
        'atomic number': 17,
        'mass number': 43,
        'mass': 42.97389 * u.u,
        'stable': False,
        'half-life': 3.13 * u.s,
    },

    'Cl-44': {
        'atomic number': 17,
        'mass number': 44,
        'mass': 43.97787 * u.u,
        'stable': False,
        'half-life': 0.56 * u.s,
    },

    'Cl-45': {
        'atomic number': 17,
        'mass number': 45,
        'mass': 44.98029 * u.u,
        'stable': False,
        'half-life': 0.413 * u.s,
    },

    'Cl-46': {
        'atomic number': 17,
        'mass number': 46,
        'mass': 45.98517 * u.u,
        'stable': False,
        'half-life': 0.232 * u.s,
    },

    'Cl-47': {
        'atomic number': 17,
        'mass number': 47,
        'mass': 46.98916 * u.u,
        'stable': False,
        'half-life': 0.101 * u.s,
    },

    'Cl-48': {
        'atomic number': 17,
        'mass number': 48,
        'mass': 47.99564 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Cl-49': {
        'atomic number': 17,
        'mass number': 49,
        'mass': 49.00123 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'Cl-50': {
        'atomic number': 17,
        'mass number': 50,
        'mass': 50.00905 * u.u,
        'stable': False,
        'half-life': '20# ms',
    },

    'Cl-51': {
        'atomic number': 17,
        'mass number': 51,
        'mass': 51.01554 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'Ar-30': {
        'atomic number': 18,
        'mass number': 30,
        'mass': 30.02307 * u.u,
        'stable': False,
    },

    'Ar-31': {
        'atomic number': 18,
        'mass number': 31,
        'mass': 31.01212 * u.u,
        'stable': False,
        'half-life': 0.0151 * u.s,
    },

    'Ar-32': {
        'atomic number': 18,
        'mass number': 32,
        'mass': 31.9976378 * u.u,
        'stable': False,
        'half-life': 0.098 * u.s,
    },

    'Ar-33': {
        'atomic number': 18,
        'mass number': 33,
        'mass': 32.98992555 * u.u,
        'stable': False,
        'half-life': 0.173 * u.s,
    },

    'Ar-34': {
        'atomic number': 18,
        'mass number': 34,
        'mass': 33.98027009 * u.u,
        'stable': False,
        'half-life': 0.8438 * u.s,
    },

    'Ar-35': {
        'atomic number': 18,
        'mass number': 35,
        'mass': 34.97525759 * u.u,
        'stable': False,
        'half-life': 1.7756 * u.s,
    },

    'Ar-36': {
        'atomic number': 18,
        'mass number': 36,
        'mass': 35.967545105 * u.u,
        'stable': True,
        'abundance': 0.003336,
    },

    'Ar-37': {
        'atomic number': 18,
        'mass number': 37,
        'mass': 36.96677633 * u.u,
        'stable': False,
        'half-life': 3024950.4 * u.s,
    },

    'Ar-38': {
        'atomic number': 18,
        'mass number': 38,
        'mass': 37.96273211 * u.u,
        'stable': True,
        'abundance': 0.000629,
    },

    'Ar-39': {
        'atomic number': 18,
        'mass number': 39,
        'mass': 38.964313 * u.u,
        'stable': False,
        'half-life': 8488813094.0 * u.s,
    },

    'Ar-40': {
        'atomic number': 18,
        'mass number': 40,
        'mass': 39.9623831237 * u.u,
        'stable': True,
        'abundance': 0.996035,
    },

    'Ar-41': {
        'atomic number': 18,
        'mass number': 41,
        'mass': 40.96450057 * u.u,
        'stable': False,
        'half-life': 6576.6 * u.s,
    },

    'Ar-42': {
        'atomic number': 18,
        'mass number': 42,
        'mass': 41.9630457 * u.u,
        'stable': False,
        'half-life': 1038222865.4 * u.s,
    },

    'Ar-43': {
        'atomic number': 18,
        'mass number': 43,
        'mass': 42.9656361 * u.u,
        'stable': False,
        'half-life': 322.2 * u.s,
    },

    'Ar-44': {
        'atomic number': 18,
        'mass number': 44,
        'mass': 43.9649238 * u.u,
        'stable': False,
        'half-life': 712.2 * u.s,
    },

    'Ar-45': {
        'atomic number': 18,
        'mass number': 45,
        'mass': 44.96803973 * u.u,
        'stable': False,
        'half-life': 21.48 * u.s,
    },

    'Ar-46': {
        'atomic number': 18,
        'mass number': 46,
        'mass': 45.968083 * u.u,
        'stable': False,
        'half-life': 8.4 * u.s,
    },

    'Ar-47': {
        'atomic number': 18,
        'mass number': 47,
        'mass': 46.972935 * u.u,
        'stable': False,
        'half-life': 1.23 * u.s,
    },

    'Ar-48': {
        'atomic number': 18,
        'mass number': 48,
        'mass': 47.97591 * u.u,
        'stable': False,
        'half-life': 0.415 * u.s,
    },

    'Ar-49': {
        'atomic number': 18,
        'mass number': 49,
        'mass': 48.9819 * u.u,
        'stable': False,
        'half-life': 0.236 * u.s,
    },

    'Ar-50': {
        'atomic number': 18,
        'mass number': 50,
        'mass': 49.98613 * u.u,
        'stable': False,
        'half-life': 0.106 * u.s,
    },

    'Ar-51': {
        'atomic number': 18,
        'mass number': 51,
        'mass': 50.9937 * u.u,
        'stable': False,
        'half-life': '60# ms',
    },

    'Ar-52': {
        'atomic number': 18,
        'mass number': 52,
        'mass': 51.99896 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Ar-53': {
        'atomic number': 18,
        'mass number': 53,
        'mass': 53.00729 * u.u,
        'stable': False,
        'half-life': '3# ms',
    },

    'K-32': {
        'atomic number': 19,
        'mass number': 32,
        'mass': 32.02265 * u.u,
        'stable': False,
    },

    'K-33': {
        'atomic number': 19,
        'mass number': 33,
        'mass': 33.00756 * u.u,
        'stable': False,
    },

    'K-34': {
        'atomic number': 19,
        'mass number': 34,
        'mass': 33.99869 * u.u,
        'stable': False,
    },

    'K-35': {
        'atomic number': 19,
        'mass number': 35,
        'mass': 34.98800541 * u.u,
        'stable': False,
        'half-life': 0.178 * u.s,
    },

    'K-36': {
        'atomic number': 19,
        'mass number': 36,
        'mass': 35.98130201 * u.u,
        'stable': False,
        'half-life': 0.341 * u.s,
    },

    'K-37': {
        'atomic number': 19,
        'mass number': 37,
        'mass': 36.97337589 * u.u,
        'stable': False,
        'half-life': 1.2365 * u.s,
    },

    'K-38': {
        'atomic number': 19,
        'mass number': 38,
        'mass': 37.96908112 * u.u,
        'stable': False,
        'half-life': 458.16 * u.s,
    },

    'K-39': {
        'atomic number': 19,
        'mass number': 39,
        'mass': 38.9637064864 * u.u,
        'stable': True,
        'abundance': 0.932581,
    },

    'K-40': {
        'atomic number': 19,
        'mass number': 40,
        'mass': 39.963998166 * u.u,
        'stable': False,
        'abundance': 0.000117,
    },

    'K-41': {
        'atomic number': 19,
        'mass number': 41,
        'mass': 40.9618252579 * u.u,
        'stable': True,
        'abundance': 0.067302,
    },

    'K-42': {
        'atomic number': 19,
        'mass number': 42,
        'mass': 41.96240231 * u.u,
        'stable': False,
        'half-life': 44478.0 * u.s,
    },

    'K-43': {
        'atomic number': 19,
        'mass number': 43,
        'mass': 42.9607347 * u.u,
        'stable': False,
        'half-life': 80280.0 * u.s,
    },

    'K-44': {
        'atomic number': 19,
        'mass number': 44,
        'mass': 43.96158699 * u.u,
        'stable': False,
        'half-life': 1327.8 * u.s,
    },

    'K-45': {
        'atomic number': 19,
        'mass number': 45,
        'mass': 44.96069149 * u.u,
        'stable': False,
        'half-life': 1068.0 * u.s,
    },

    'K-46': {
        'atomic number': 19,
        'mass number': 46,
        'mass': 45.96198159 * u.u,
        'stable': False,
        'half-life': 105.0 * u.s,
    },

    'K-47': {
        'atomic number': 19,
        'mass number': 47,
        'mass': 46.9616616 * u.u,
        'stable': False,
        'half-life': 17.5 * u.s,
    },

    'K-48': {
        'atomic number': 19,
        'mass number': 48,
        'mass': 47.96534119 * u.u,
        'stable': False,
        'half-life': 6.8 * u.s,
    },

    'K-49': {
        'atomic number': 19,
        'mass number': 49,
        'mass': 48.96821075 * u.u,
        'stable': False,
        'half-life': 1.26 * u.s,
    },

    'K-50': {
        'atomic number': 19,
        'mass number': 50,
        'mass': 49.97238 * u.u,
        'stable': False,
        'half-life': 0.472 * u.s,
    },

    'K-51': {
        'atomic number': 19,
        'mass number': 51,
        'mass': 50.975828 * u.u,
        'stable': False,
        'half-life': 0.365 * u.s,
    },

    'K-52': {
        'atomic number': 19,
        'mass number': 52,
        'mass': 51.98224 * u.u,
        'stable': False,
        'half-life': 0.11 * u.s,
    },

    'K-53': {
        'atomic number': 19,
        'mass number': 53,
        'mass': 52.98746 * u.u,
        'stable': False,
        'half-life': 0.03 * u.s,
    },

    'K-54': {
        'atomic number': 19,
        'mass number': 54,
        'mass': 53.99463 * u.u,
        'stable': False,
        'half-life': 0.01 * u.s,
    },

    'K-55': {
        'atomic number': 19,
        'mass number': 55,
        'mass': 55.00076 * u.u,
        'stable': False,
        'half-life': '3# ms',
    },

    'K-56': {
        'atomic number': 19,
        'mass number': 56,
        'mass': 56.00851 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'Ca-34': {
        'atomic number': 20,
        'mass number': 34,
        'mass': 34.01487 * u.u,
        'stable': False,
    },

    'Ca-35': {
        'atomic number': 20,
        'mass number': 35,
        'mass': 35.00514 * u.u,
        'stable': False,
        'half-life': 0.0257 * u.s,
    },

    'Ca-36': {
        'atomic number': 20,
        'mass number': 36,
        'mass': 35.993074 * u.u,
        'stable': False,
        'half-life': 0.1012 * u.s,
    },

    'Ca-37': {
        'atomic number': 20,
        'mass number': 37,
        'mass': 36.98589785 * u.u,
        'stable': False,
        'half-life': 0.1811 * u.s,
    },

    'Ca-38': {
        'atomic number': 20,
        'mass number': 38,
        'mass': 37.97631922 * u.u,
        'stable': False,
        'half-life': 0.4437 * u.s,
    },

    'Ca-39': {
        'atomic number': 20,
        'mass number': 39,
        'mass': 38.97071081 * u.u,
        'stable': False,
        'half-life': 0.8603 * u.s,
    },

    'Ca-40': {
        'atomic number': 20,
        'mass number': 40,
        'mass': 39.962590863 * u.u,
        'stable': True,
        'abundance': 0.96941,
    },

    'Ca-41': {
        'atomic number': 20,
        'mass number': 41,
        'mass': 40.96227792 * u.u,
        'stable': False,
        'half-life': 3136758444400.0 * u.s,
    },

    'Ca-42': {
        'atomic number': 20,
        'mass number': 42,
        'mass': 41.95861783 * u.u,
        'stable': True,
        'abundance': 0.00647,
    },

    'Ca-43': {
        'atomic number': 20,
        'mass number': 43,
        'mass': 42.95876644 * u.u,
        'stable': True,
        'abundance': 0.00135,
    },

    'Ca-44': {
        'atomic number': 20,
        'mass number': 44,
        'mass': 43.95548156 * u.u,
        'stable': True,
        'abundance': 0.02086,
    },

    'Ca-45': {
        'atomic number': 20,
        'mass number': 45,
        'mass': 44.95618635 * u.u,
        'stable': False,
        'half-life': 14049504.0 * u.s,
    },

    'Ca-46': {
        'atomic number': 20,
        'mass number': 46,
        'mass': 45.953689 * u.u,
        'stable': True,
        'abundance': 4e-05,
    },

    'Ca-47': {
        'atomic number': 20,
        'mass number': 47,
        'mass': 46.9545424 * u.u,
        'stable': False,
        'half-life': 391910.4 * u.s,
    },

    'Ca-48': {
        'atomic number': 20,
        'mass number': 48,
        'mass': 47.95252276 * u.u,
        'stable': False,
        'abundance': 0.00187,
    },

    'Ca-49': {
        'atomic number': 20,
        'mass number': 49,
        'mass': 48.95566274 * u.u,
        'stable': False,
        'half-life': 523.08 * u.s,
    },

    'Ca-50': {
        'atomic number': 20,
        'mass number': 50,
        'mass': 49.9574992 * u.u,
        'stable': False,
        'half-life': 13.9 * u.s,
    },

    'Ca-51': {
        'atomic number': 20,
        'mass number': 51,
        'mass': 50.960989 * u.u,
        'stable': False,
        'half-life': 10.0 * u.s,
    },

    'Ca-52': {
        'atomic number': 20,
        'mass number': 52,
        'mass': 51.963217 * u.u,
        'stable': False,
        'half-life': 4.6 * u.s,
    },

    'Ca-53': {
        'atomic number': 20,
        'mass number': 53,
        'mass': 52.96945 * u.u,
        'stable': False,
        'half-life': 0.461 * u.s,
    },

    'Ca-54': {
        'atomic number': 20,
        'mass number': 54,
        'mass': 53.9734 * u.u,
        'stable': False,
        'half-life': 0.09 * u.s,
    },

    'Ca-55': {
        'atomic number': 20,
        'mass number': 55,
        'mass': 54.9803 * u.u,
        'stable': False,
        'half-life': 0.022 * u.s,
    },

    'Ca-56': {
        'atomic number': 20,
        'mass number': 56,
        'mass': 55.98508 * u.u,
        'stable': False,
        'half-life': 0.011 * u.s,
    },

    'Ca-57': {
        'atomic number': 20,
        'mass number': 57,
        'mass': 56.99262 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Ca-58': {
        'atomic number': 20,
        'mass number': 58,
        'mass': 57.99794 * u.u,
        'stable': False,
        'half-life': '3# ms',
    },

    'Sc-36': {
        'atomic number': 21,
        'mass number': 36,
        'mass': 36.01648 * u.u,
        'stable': False,
    },

    'Sc-37': {
        'atomic number': 21,
        'mass number': 37,
        'mass': 37.00374 * u.u,
        'stable': False,
    },

    'Sc-38': {
        'atomic number': 21,
        'mass number': 38,
        'mass': 37.99512 * u.u,
        'stable': False,
    },

    'Sc-39': {
        'atomic number': 21,
        'mass number': 39,
        'mass': 38.984785 * u.u,
        'stable': False,
        'half-life': '<300 ns',
    },

    'Sc-40': {
        'atomic number': 21,
        'mass number': 40,
        'mass': 39.9779673 * u.u,
        'stable': False,
        'half-life': 0.1823 * u.s,
    },

    'Sc-41': {
        'atomic number': 21,
        'mass number': 41,
        'mass': 40.969251105 * u.u,
        'stable': False,
        'half-life': 0.5963 * u.s,
    },

    'Sc-42': {
        'atomic number': 21,
        'mass number': 42,
        'mass': 41.96551653 * u.u,
        'stable': False,
        'half-life': 0.68079 * u.s,
    },

    'Sc-43': {
        'atomic number': 21,
        'mass number': 43,
        'mass': 42.9611505 * u.u,
        'stable': False,
        'half-life': 14007.6 * u.s,
    },

    'Sc-44': {
        'atomic number': 21,
        'mass number': 44,
        'mass': 43.9594029 * u.u,
        'stable': False,
        'half-life': 14551.2 * u.s,
    },

    'Sc-45': {
        'atomic number': 21,
        'mass number': 45,
        'mass': 44.95590828 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Sc-46': {
        'atomic number': 21,
        'mass number': 46,
        'mass': 45.95516826 * u.u,
        'stable': False,
        'half-life': 7242998.4 * u.s,
    },

    'Sc-47': {
        'atomic number': 21,
        'mass number': 47,
        'mass': 46.9524037 * u.u,
        'stable': False,
        'half-life': 289370.88 * u.s,
    },

    'Sc-48': {
        'atomic number': 21,
        'mass number': 48,
        'mass': 47.9522236 * u.u,
        'stable': False,
        'half-life': 157212.0 * u.s,
    },

    'Sc-49': {
        'atomic number': 21,
        'mass number': 49,
        'mass': 48.9500146 * u.u,
        'stable': False,
        'half-life': 3430.8 * u.s,
    },

    'Sc-50': {
        'atomic number': 21,
        'mass number': 50,
        'mass': 49.952176 * u.u,
        'stable': False,
        'half-life': 102.5 * u.s,
    },

    'Sc-51': {
        'atomic number': 21,
        'mass number': 51,
        'mass': 50.953592 * u.u,
        'stable': False,
        'half-life': 12.4 * u.s,
    },

    'Sc-52': {
        'atomic number': 21,
        'mass number': 52,
        'mass': 51.95688 * u.u,
        'stable': False,
        'half-life': 8.2 * u.s,
    },

    'Sc-53': {
        'atomic number': 21,
        'mass number': 53,
        'mass': 52.95909 * u.u,
        'stable': False,
        'half-life': 2.4 * u.s,
    },

    'Sc-54': {
        'atomic number': 21,
        'mass number': 54,
        'mass': 53.96393 * u.u,
        'stable': False,
        'half-life': 0.526 * u.s,
    },

    'Sc-55': {
        'atomic number': 21,
        'mass number': 55,
        'mass': 54.96782 * u.u,
        'stable': False,
        'half-life': 0.096 * u.s,
    },

    'Sc-56': {
        'atomic number': 21,
        'mass number': 56,
        'mass': 55.97345 * u.u,
        'stable': False,
        'half-life': 0.026 * u.s,
    },

    'Sc-57': {
        'atomic number': 21,
        'mass number': 57,
        'mass': 56.97777 * u.u,
        'stable': False,
        'half-life': 0.022 * u.s,
    },

    'Sc-58': {
        'atomic number': 21,
        'mass number': 58,
        'mass': 57.98403 * u.u,
        'stable': False,
        'half-life': 0.012 * u.s,
    },

    'Sc-59': {
        'atomic number': 21,
        'mass number': 59,
        'mass': 58.98894 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Sc-60': {
        'atomic number': 21,
        'mass number': 60,
        'mass': 59.99565 * u.u,
        'stable': False,
        'half-life': '3# ms',
    },

    'Sc-61': {
        'atomic number': 21,
        'mass number': 61,
        'mass': 61.001 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'Ti-38': {
        'atomic number': 22,
        'mass number': 38,
        'mass': 38.01145 * u.u,
        'stable': False,
    },

    'Ti-39': {
        'atomic number': 22,
        'mass number': 39,
        'mass': 39.00236 * u.u,
        'stable': False,
        'half-life': 0.0285 * u.s,
    },

    'Ti-40': {
        'atomic number': 22,
        'mass number': 40,
        'mass': 39.9905 * u.u,
        'stable': False,
        'half-life': 0.0524 * u.s,
    },

    'Ti-41': {
        'atomic number': 22,
        'mass number': 41,
        'mass': 40.983148 * u.u,
        'stable': False,
        'half-life': 0.0819 * u.s,
    },

    'Ti-42': {
        'atomic number': 22,
        'mass number': 42,
        'mass': 41.97304903 * u.u,
        'stable': False,
        'half-life': 0.20865 * u.s,
    },

    'Ti-43': {
        'atomic number': 22,
        'mass number': 43,
        'mass': 42.9685225 * u.u,
        'stable': False,
        'half-life': 0.509 * u.s,
    },

    'Ti-44': {
        'atomic number': 22,
        'mass number': 44,
        'mass': 43.95968995 * u.u,
        'stable': False,
        'half-life': 1914105600.0 * u.s,
    },

    'Ti-45': {
        'atomic number': 22,
        'mass number': 45,
        'mass': 44.95812198 * u.u,
        'stable': False,
        'half-life': 11088.0 * u.s,
    },

    'Ti-46': {
        'atomic number': 22,
        'mass number': 46,
        'mass': 45.95262772 * u.u,
        'stable': True,
        'abundance': 0.0825,
    },

    'Ti-47': {
        'atomic number': 22,
        'mass number': 47,
        'mass': 46.95175879 * u.u,
        'stable': True,
        'abundance': 0.0744,
    },

    'Ti-48': {
        'atomic number': 22,
        'mass number': 48,
        'mass': 47.94794198 * u.u,
        'stable': True,
        'abundance': 0.7372,
    },

    'Ti-49': {
        'atomic number': 22,
        'mass number': 49,
        'mass': 48.94786568 * u.u,
        'stable': True,
        'abundance': 0.0541,
    },

    'Ti-50': {
        'atomic number': 22,
        'mass number': 50,
        'mass': 49.94478689 * u.u,
        'stable': True,
        'abundance': 0.0518,
    },

    'Ti-51': {
        'atomic number': 22,
        'mass number': 51,
        'mass': 50.94661065 * u.u,
        'stable': False,
        'half-life': 345.6 * u.s,
    },

    'Ti-52': {
        'atomic number': 22,
        'mass number': 52,
        'mass': 51.946893 * u.u,
        'stable': False,
        'half-life': 102.0 * u.s,
    },

    'Ti-53': {
        'atomic number': 22,
        'mass number': 53,
        'mass': 52.94973 * u.u,
        'stable': False,
        'half-life': 32.7 * u.s,
    },

    'Ti-54': {
        'atomic number': 22,
        'mass number': 54,
        'mass': 53.95105 * u.u,
        'stable': False,
        'half-life': 2.1 * u.s,
    },

    'Ti-55': {
        'atomic number': 22,
        'mass number': 55,
        'mass': 54.95527 * u.u,
        'stable': False,
        'half-life': 1.3 * u.s,
    },

    'Ti-56': {
        'atomic number': 22,
        'mass number': 56,
        'mass': 55.95791 * u.u,
        'stable': False,
        'half-life': 0.2 * u.s,
    },

    'Ti-57': {
        'atomic number': 22,
        'mass number': 57,
        'mass': 56.96364 * u.u,
        'stable': False,
        'half-life': 0.095 * u.s,
    },

    'Ti-58': {
        'atomic number': 22,
        'mass number': 58,
        'mass': 57.9666 * u.u,
        'stable': False,
        'half-life': 0.055 * u.s,
    },

    'Ti-59': {
        'atomic number': 22,
        'mass number': 59,
        'mass': 58.97247 * u.u,
        'stable': False,
        'half-life': 0.0285 * u.s,
    },

    'Ti-60': {
        'atomic number': 22,
        'mass number': 60,
        'mass': 59.97603 * u.u,
        'stable': False,
        'half-life': 0.0222 * u.s,
    },

    'Ti-61': {
        'atomic number': 22,
        'mass number': 61,
        'mass': 60.98245 * u.u,
        'stable': False,
        'half-life': 0.015 * u.s,
    },

    'Ti-62': {
        'atomic number': 22,
        'mass number': 62,
        'mass': 61.98651 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Ti-63': {
        'atomic number': 22,
        'mass number': 63,
        'mass': 62.99375 * u.u,
        'stable': False,
        'half-life': '3# ms',
    },

    'V-40': {
        'atomic number': 23,
        'mass number': 40,
        'mass': 40.01276 * u.u,
        'stable': False,
    },

    'V-41': {
        'atomic number': 23,
        'mass number': 41,
        'mass': 41.00021 * u.u,
        'stable': False,
    },

    'V-42': {
        'atomic number': 23,
        'mass number': 42,
        'mass': 41.99182 * u.u,
        'stable': False,
    },

    'V-43': {
        'atomic number': 23,
        'mass number': 43,
        'mass': 42.980766 * u.u,
        'stable': False,
        'half-life': 0.0793 * u.s,
    },

    'V-44': {
        'atomic number': 23,
        'mass number': 44,
        'mass': 43.97411 * u.u,
        'stable': False,
        'half-life': 0.111 * u.s,
    },

    'V-45': {
        'atomic number': 23,
        'mass number': 45,
        'mass': 44.9657748 * u.u,
        'stable': False,
        'half-life': 0.547 * u.s,
    },

    'V-46': {
        'atomic number': 23,
        'mass number': 46,
        'mass': 45.96019878 * u.u,
        'stable': False,
        'half-life': 0.42264 * u.s,
    },

    'V-47': {
        'atomic number': 23,
        'mass number': 47,
        'mass': 46.95490491 * u.u,
        'stable': False,
        'half-life': 1956.0 * u.s,
    },

    'V-48': {
        'atomic number': 23,
        'mass number': 48,
        'mass': 47.9522522 * u.u,
        'stable': False,
        'half-life': 1380110.4 * u.s,
    },

    'V-49': {
        'atomic number': 23,
        'mass number': 49,
        'mass': 48.9485118 * u.u,
        'stable': False,
        'half-life': 28512000.0 * u.s,
    },

    'V-50': {
        'atomic number': 23,
        'mass number': 50,
        'mass': 49.94715601 * u.u,
        'stable': False,
        'abundance': 0.0025,
    },

    'V-51': {
        'atomic number': 23,
        'mass number': 51,
        'mass': 50.94395704 * u.u,
        'stable': True,
        'abundance': 0.9975,
    },

    'V-52': {
        'atomic number': 23,
        'mass number': 52,
        'mass': 51.94477301 * u.u,
        'stable': False,
        'half-life': 224.58 * u.s,
    },

    'V-53': {
        'atomic number': 23,
        'mass number': 53,
        'mass': 52.9443367 * u.u,
        'stable': False,
        'half-life': 92.58 * u.s,
    },

    'V-54': {
        'atomic number': 23,
        'mass number': 54,
        'mass': 53.946439 * u.u,
        'stable': False,
        'half-life': 49.8 * u.s,
    },

    'V-55': {
        'atomic number': 23,
        'mass number': 55,
        'mass': 54.94724 * u.u,
        'stable': False,
        'half-life': 6.54 * u.s,
    },

    'V-56': {
        'atomic number': 23,
        'mass number': 56,
        'mass': 55.95048 * u.u,
        'stable': False,
        'half-life': 0.216 * u.s,
    },

    'V-57': {
        'atomic number': 23,
        'mass number': 57,
        'mass': 56.95252 * u.u,
        'stable': False,
        'half-life': 0.35 * u.s,
    },

    'V-58': {
        'atomic number': 23,
        'mass number': 58,
        'mass': 57.95672 * u.u,
        'stable': False,
        'half-life': 0.191 * u.s,
    },

    'V-59': {
        'atomic number': 23,
        'mass number': 59,
        'mass': 58.95939 * u.u,
        'stable': False,
        'half-life': 0.095 * u.s,
    },

    'V-60': {
        'atomic number': 23,
        'mass number': 60,
        'mass': 59.96431 * u.u,
        'stable': False,
        'half-life': 0.122 * u.s,
    },

    'V-61': {
        'atomic number': 23,
        'mass number': 61,
        'mass': 60.96725 * u.u,
        'stable': False,
        'half-life': 0.0482 * u.s,
    },

    'V-62': {
        'atomic number': 23,
        'mass number': 62,
        'mass': 61.97265 * u.u,
        'stable': False,
        'half-life': 0.0336 * u.s,
    },

    'V-63': {
        'atomic number': 23,
        'mass number': 63,
        'mass': 62.97639 * u.u,
        'stable': False,
        'half-life': 0.0196 * u.s,
    },

    'V-64': {
        'atomic number': 23,
        'mass number': 64,
        'mass': 63.98264 * u.u,
        'stable': False,
        'half-life': 0.015 * u.s,
    },

    'V-65': {
        'atomic number': 23,
        'mass number': 65,
        'mass': 64.9875 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'V-66': {
        'atomic number': 23,
        'mass number': 66,
        'mass': 65.99398 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Cr-42': {
        'atomic number': 24,
        'mass number': 42,
        'mass': 42.0067 * u.u,
        'stable': False,
        'half-life': 0.0133 * u.s,
    },

    'Cr-43': {
        'atomic number': 24,
        'mass number': 43,
        'mass': 42.99753 * u.u,
        'stable': False,
        'half-life': 0.0211 * u.s,
    },

    'Cr-44': {
        'atomic number': 24,
        'mass number': 44,
        'mass': 43.98536 * u.u,
        'stable': False,
        'half-life': 0.0428 * u.s,
    },

    'Cr-45': {
        'atomic number': 24,
        'mass number': 45,
        'mass': 44.97905 * u.u,
        'stable': False,
        'half-life': 0.0609 * u.s,
    },

    'Cr-46': {
        'atomic number': 24,
        'mass number': 46,
        'mass': 45.968359 * u.u,
        'stable': False,
        'half-life': 0.2243 * u.s,
    },

    'Cr-47': {
        'atomic number': 24,
        'mass number': 47,
        'mass': 46.9628974 * u.u,
        'stable': False,
        'half-life': 0.5 * u.s,
    },

    'Cr-48': {
        'atomic number': 24,
        'mass number': 48,
        'mass': 47.9540291 * u.u,
        'stable': False,
        'half-life': 77616.0 * u.s,
    },

    'Cr-49': {
        'atomic number': 24,
        'mass number': 49,
        'mass': 48.9513333 * u.u,
        'stable': False,
        'half-life': 2538.0 * u.s,
    },

    'Cr-50': {
        'atomic number': 24,
        'mass number': 50,
        'mass': 49.94604183 * u.u,
        'stable': True,
        'abundance': 0.04345,
    },

    'Cr-51': {
        'atomic number': 24,
        'mass number': 51,
        'mass': 50.94476502 * u.u,
        'stable': False,
        'half-life': 2393366.4 * u.s,
    },

    'Cr-52': {
        'atomic number': 24,
        'mass number': 52,
        'mass': 51.94050623 * u.u,
        'stable': True,
        'abundance': 0.83789,
    },

    'Cr-53': {
        'atomic number': 24,
        'mass number': 53,
        'mass': 52.94064815 * u.u,
        'stable': True,
        'abundance': 0.09501,
    },

    'Cr-54': {
        'atomic number': 24,
        'mass number': 54,
        'mass': 53.93887916 * u.u,
        'stable': True,
        'abundance': 0.02365,
    },

    'Cr-55': {
        'atomic number': 24,
        'mass number': 55,
        'mass': 54.94083843 * u.u,
        'stable': False,
        'half-life': 209.82 * u.s,
    },

    'Cr-56': {
        'atomic number': 24,
        'mass number': 56,
        'mass': 55.9406531 * u.u,
        'stable': False,
        'half-life': 356.4 * u.s,
    },

    'Cr-57': {
        'atomic number': 24,
        'mass number': 57,
        'mass': 56.943613 * u.u,
        'stable': False,
        'half-life': 21.1 * u.s,
    },

    'Cr-58': {
        'atomic number': 24,
        'mass number': 58,
        'mass': 57.94435 * u.u,
        'stable': False,
        'half-life': 7.0 * u.s,
    },

    'Cr-59': {
        'atomic number': 24,
        'mass number': 59,
        'mass': 58.94859 * u.u,
        'stable': False,
        'half-life': 1.05 * u.s,
    },

    'Cr-60': {
        'atomic number': 24,
        'mass number': 60,
        'mass': 59.95008 * u.u,
        'stable': False,
        'half-life': 0.49 * u.s,
    },

    'Cr-61': {
        'atomic number': 24,
        'mass number': 61,
        'mass': 60.95442 * u.u,
        'stable': False,
        'half-life': 0.243 * u.s,
    },

    'Cr-62': {
        'atomic number': 24,
        'mass number': 62,
        'mass': 61.9561 * u.u,
        'stable': False,
        'half-life': 0.206 * u.s,
    },

    'Cr-63': {
        'atomic number': 24,
        'mass number': 63,
        'mass': 62.96165 * u.u,
        'stable': False,
        'half-life': 0.129 * u.s,
    },

    'Cr-64': {
        'atomic number': 24,
        'mass number': 64,
        'mass': 63.96408 * u.u,
        'stable': False,
        'half-life': 0.043 * u.s,
    },

    'Cr-65': {
        'atomic number': 24,
        'mass number': 65,
        'mass': 64.96996 * u.u,
        'stable': False,
        'half-life': 0.0275 * u.s,
    },

    'Cr-66': {
        'atomic number': 24,
        'mass number': 66,
        'mass': 65.97366 * u.u,
        'stable': False,
        'half-life': 0.0238 * u.s,
    },

    'Cr-67': {
        'atomic number': 24,
        'mass number': 67,
        'mass': 66.98016 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Cr-68': {
        'atomic number': 24,
        'mass number': 68,
        'mass': 67.98403 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Mn-44': {
        'atomic number': 25,
        'mass number': 44,
        'mass': 44.00715 * u.u,
        'stable': False,
    },

    'Mn-45': {
        'atomic number': 25,
        'mass number': 45,
        'mass': 44.99449 * u.u,
        'stable': False,
    },

    'Mn-46': {
        'atomic number': 25,
        'mass number': 46,
        'mass': 45.98609 * u.u,
        'stable': False,
        'half-life': 0.0362 * u.s,
    },

    'Mn-47': {
        'atomic number': 25,
        'mass number': 47,
        'mass': 46.975775 * u.u,
        'stable': False,
        'half-life': 0.088 * u.s,
    },

    'Mn-48': {
        'atomic number': 25,
        'mass number': 48,
        'mass': 47.96852 * u.u,
        'stable': False,
        'half-life': 0.1581 * u.s,
    },

    'Mn-49': {
        'atomic number': 25,
        'mass number': 49,
        'mass': 48.959595 * u.u,
        'stable': False,
        'half-life': 0.382 * u.s,
    },

    'Mn-50': {
        'atomic number': 25,
        'mass number': 50,
        'mass': 49.95423778 * u.u,
        'stable': False,
        'half-life': 0.28319 * u.s,
    },

    'Mn-51': {
        'atomic number': 25,
        'mass number': 51,
        'mass': 50.94820847 * u.u,
        'stable': False,
        'half-life': 2772.0 * u.s,
    },

    'Mn-52': {
        'atomic number': 25,
        'mass number': 52,
        'mass': 51.9455639 * u.u,
        'stable': False,
        'half-life': 483062.4 * u.s,
    },

    'Mn-53': {
        'atomic number': 25,
        'mass number': 53,
        'mass': 52.94128889 * u.u,
        'stable': False,
        'half-life': 116760626200000.0 * u.s,
    },

    'Mn-54': {
        'atomic number': 25,
        'mass number': 54,
        'mass': 53.9403576 * u.u,
        'stable': False,
        'half-life': 26959219.200000003 * u.s,
    },

    'Mn-55': {
        'atomic number': 25,
        'mass number': 55,
        'mass': 54.93804391 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Mn-56': {
        'atomic number': 25,
        'mass number': 56,
        'mass': 55.93890369 * u.u,
        'stable': False,
        'half-life': 9284.04 * u.s,
    },

    'Mn-57': {
        'atomic number': 25,
        'mass number': 57,
        'mass': 56.9382861 * u.u,
        'stable': False,
        'half-life': 85.4 * u.s,
    },

    'Mn-58': {
        'atomic number': 25,
        'mass number': 58,
        'mass': 57.9400666 * u.u,
        'stable': False,
        'half-life': 3.0 * u.s,
    },

    'Mn-59': {
        'atomic number': 25,
        'mass number': 59,
        'mass': 58.9403911 * u.u,
        'stable': False,
        'half-life': 4.59 * u.s,
    },

    'Mn-60': {
        'atomic number': 25,
        'mass number': 60,
        'mass': 59.9431366 * u.u,
        'stable': False,
        'half-life': 0.28 * u.s,
    },

    'Mn-61': {
        'atomic number': 25,
        'mass number': 61,
        'mass': 60.9444525 * u.u,
        'stable': False,
        'half-life': 0.709 * u.s,
    },

    'Mn-62': {
        'atomic number': 25,
        'mass number': 62,
        'mass': 61.94795 * u.u,
        'stable': False,
        'half-life': 0.092 * u.s,
    },

    'Mn-63': {
        'atomic number': 25,
        'mass number': 63,
        'mass': 62.9496647 * u.u,
        'stable': False,
        'half-life': 0.275 * u.s,
    },

    'Mn-64': {
        'atomic number': 25,
        'mass number': 64,
        'mass': 63.9538494 * u.u,
        'stable': False,
        'half-life': 0.0888 * u.s,
    },

    'Mn-65': {
        'atomic number': 25,
        'mass number': 65,
        'mass': 64.9560198 * u.u,
        'stable': False,
        'half-life': 0.0919 * u.s,
    },

    'Mn-66': {
        'atomic number': 25,
        'mass number': 66,
        'mass': 65.960547 * u.u,
        'stable': False,
        'half-life': 0.0642 * u.s,
    },

    'Mn-67': {
        'atomic number': 25,
        'mass number': 67,
        'mass': 66.96424 * u.u,
        'stable': False,
        'half-life': 0.0467 * u.s,
    },

    'Mn-68': {
        'atomic number': 25,
        'mass number': 68,
        'mass': 67.96962 * u.u,
        'stable': False,
        'half-life': 0.0337 * u.s,
    },

    'Mn-69': {
        'atomic number': 25,
        'mass number': 69,
        'mass': 68.97366 * u.u,
        'stable': False,
        'half-life': 0.0221 * u.s,
    },

    'Mn-70': {
        'atomic number': 25,
        'mass number': 70,
        'mass': 69.97937 * u.u,
        'stable': False,
        'half-life': 0.0199 * u.s,
    },

    'Mn-71': {
        'atomic number': 25,
        'mass number': 71,
        'mass': 70.98368 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Fe-45': {
        'atomic number': 26,
        'mass number': 45,
        'mass': 45.01442 * u.u,
        'stable': False,
        'half-life': 0.0022 * u.s,
    },

    'Fe-46': {
        'atomic number': 26,
        'mass number': 46,
        'mass': 46.00063 * u.u,
        'stable': False,
        'half-life': 0.013 * u.s,
    },

    'Fe-47': {
        'atomic number': 26,
        'mass number': 47,
        'mass': 46.99185 * u.u,
        'stable': False,
        'half-life': 0.0219 * u.s,
    },

    'Fe-48': {
        'atomic number': 26,
        'mass number': 48,
        'mass': 47.98023 * u.u,
        'stable': False,
        'half-life': 0.0453 * u.s,
    },

    'Fe-49': {
        'atomic number': 26,
        'mass number': 49,
        'mass': 48.973429 * u.u,
        'stable': False,
        'half-life': 0.0647 * u.s,
    },

    'Fe-50': {
        'atomic number': 26,
        'mass number': 50,
        'mass': 49.962975 * u.u,
        'stable': False,
        'half-life': 0.1521 * u.s,
    },

    'Fe-51': {
        'atomic number': 26,
        'mass number': 51,
        'mass': 50.956841 * u.u,
        'stable': False,
        'half-life': 0.3054 * u.s,
    },

    'Fe-52': {
        'atomic number': 26,
        'mass number': 52,
        'mass': 51.9481131 * u.u,
        'stable': False,
        'half-life': 29790.0 * u.s,
    },

    'Fe-53': {
        'atomic number': 26,
        'mass number': 53,
        'mass': 52.9453064 * u.u,
        'stable': False,
        'half-life': 510.6 * u.s,
    },

    'Fe-54': {
        'atomic number': 26,
        'mass number': 54,
        'mass': 53.93960899 * u.u,
        'stable': True,
        'abundance': 0.05845,
    },

    'Fe-55': {
        'atomic number': 26,
        'mass number': 55,
        'mass': 54.93829199 * u.u,
        'stable': False,
        'half-life': 86592204.944 * u.s,
    },

    'Fe-56': {
        'atomic number': 26,
        'mass number': 56,
        'mass': 55.93493633 * u.u,
        'stable': True,
        'abundance': 0.91754,
    },

    'Fe-57': {
        'atomic number': 26,
        'mass number': 57,
        'mass': 56.93539284 * u.u,
        'stable': True,
        'abundance': 0.02119,
    },

    'Fe-58': {
        'atomic number': 26,
        'mass number': 58,
        'mass': 57.93327443 * u.u,
        'stable': True,
        'abundance': 0.00282,
    },

    'Fe-59': {
        'atomic number': 26,
        'mass number': 59,
        'mass': 58.93487434 * u.u,
        'stable': False,
        'half-life': 3845439.36 * u.s,
    },

    'Fe-60': {
        'atomic number': 26,
        'mass number': 60,
        'mass': 59.9340711 * u.u,
        'stable': False,
        'half-life': 82679146120000.0 * u.s,
    },

    'Fe-61': {
        'atomic number': 26,
        'mass number': 61,
        'mass': 60.9367462 * u.u,
        'stable': False,
        'half-life': 358.8 * u.s,
    },

    'Fe-62': {
        'atomic number': 26,
        'mass number': 62,
        'mass': 61.9367918 * u.u,
        'stable': False,
        'half-life': 68.0 * u.s,
    },

    'Fe-63': {
        'atomic number': 26,
        'mass number': 63,
        'mass': 62.9402727 * u.u,
        'stable': False,
        'half-life': 6.1 * u.s,
    },

    'Fe-64': {
        'atomic number': 26,
        'mass number': 64,
        'mass': 63.9409878 * u.u,
        'stable': False,
        'half-life': 2.0 * u.s,
    },

    'Fe-65': {
        'atomic number': 26,
        'mass number': 65,
        'mass': 64.9450115 * u.u,
        'stable': False,
        'half-life': 0.81 * u.s,
    },

    'Fe-66': {
        'atomic number': 26,
        'mass number': 66,
        'mass': 65.94625 * u.u,
        'stable': False,
        'half-life': 0.351 * u.s,
    },

    'Fe-67': {
        'atomic number': 26,
        'mass number': 67,
        'mass': 66.95054 * u.u,
        'stable': False,
        'half-life': 0.394 * u.s,
    },

    'Fe-68': {
        'atomic number': 26,
        'mass number': 68,
        'mass': 67.95295 * u.u,
        'stable': False,
        'half-life': 0.188 * u.s,
    },

    'Fe-69': {
        'atomic number': 26,
        'mass number': 69,
        'mass': 68.95807 * u.u,
        'stable': False,
        'half-life': 0.1082 * u.s,
    },

    'Fe-70': {
        'atomic number': 26,
        'mass number': 70,
        'mass': 69.96102 * u.u,
        'stable': False,
        'half-life': 0.063 * u.s,
    },

    'Fe-71': {
        'atomic number': 26,
        'mass number': 71,
        'mass': 70.96672 * u.u,
        'stable': False,
        'half-life': 0.0337 * u.s,
    },

    'Fe-72': {
        'atomic number': 26,
        'mass number': 72,
        'mass': 71.96983 * u.u,
        'stable': False,
        'half-life': 0.019 * u.s,
    },

    'Fe-73': {
        'atomic number': 26,
        'mass number': 73,
        'mass': 72.97572 * u.u,
        'stable': False,
        'half-life': 0.0129 * u.s,
    },

    'Fe-74': {
        'atomic number': 26,
        'mass number': 74,
        'mass': 73.97935 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'Co-47': {
        'atomic number': 27,
        'mass number': 47,
        'mass': 47.01057 * u.u,
        'stable': False,
    },

    'Co-48': {
        'atomic number': 27,
        'mass number': 48,
        'mass': 48.00093 * u.u,
        'stable': False,
    },

    'Co-49': {
        'atomic number': 27,
        'mass number': 49,
        'mass': 48.98891 * u.u,
        'stable': False,
    },

    'Co-50': {
        'atomic number': 27,
        'mass number': 50,
        'mass': 49.98091 * u.u,
        'stable': False,
        'half-life': 0.0388 * u.s,
    },

    'Co-51': {
        'atomic number': 27,
        'mass number': 51,
        'mass': 50.970647 * u.u,
        'stable': False,
        'half-life': 0.0688 * u.s,
    },

    'Co-52': {
        'atomic number': 27,
        'mass number': 52,
        'mass': 51.96351 * u.u,
        'stable': False,
        'half-life': 0.1111 * u.s,
    },

    'Co-53': {
        'atomic number': 27,
        'mass number': 53,
        'mass': 52.9542041 * u.u,
        'stable': False,
        'half-life': 0.242 * u.s,
    },

    'Co-54': {
        'atomic number': 27,
        'mass number': 54,
        'mass': 53.94845987 * u.u,
        'stable': False,
        'half-life': 0.19328 * u.s,
    },

    'Co-55': {
        'atomic number': 27,
        'mass number': 55,
        'mass': 54.9419972 * u.u,
        'stable': False,
        'half-life': 63108.0 * u.s,
    },

    'Co-56': {
        'atomic number': 27,
        'mass number': 56,
        'mass': 55.9398388 * u.u,
        'stable': False,
        'half-life': 6673190.4 * u.s,
    },

    'Co-57': {
        'atomic number': 27,
        'mass number': 57,
        'mass': 56.93629057 * u.u,
        'stable': False,
        'half-life': 23510304.0 * u.s,
    },

    'Co-58': {
        'atomic number': 27,
        'mass number': 58,
        'mass': 57.9357521 * u.u,
        'stable': False,
        'half-life': 6114528.0 * u.s,
    },

    'Co-59': {
        'atomic number': 27,
        'mass number': 59,
        'mass': 58.93319429 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Co-60': {
        'atomic number': 27,
        'mass number': 60,
        'mass': 59.9338163 * u.u,
        'stable': False,
        'half-life': 166337280.0 * u.s,
    },

    'Co-61': {
        'atomic number': 27,
        'mass number': 61,
        'mass': 60.93247662 * u.u,
        'stable': False,
        'half-life': 5936.4 * u.s,
    },

    'Co-62': {
        'atomic number': 27,
        'mass number': 62,
        'mass': 61.934059 * u.u,
        'stable': False,
        'half-life': 92.4 * u.s,
    },

    'Co-63': {
        'atomic number': 27,
        'mass number': 63,
        'mass': 62.9336 * u.u,
        'stable': False,
        'half-life': 26.9 * u.s,
    },

    'Co-64': {
        'atomic number': 27,
        'mass number': 64,
        'mass': 63.935811 * u.u,
        'stable': False,
        'half-life': 0.3 * u.s,
    },

    'Co-65': {
        'atomic number': 27,
        'mass number': 65,
        'mass': 64.9364621 * u.u,
        'stable': False,
        'half-life': 1.16 * u.s,
    },

    'Co-66': {
        'atomic number': 27,
        'mass number': 66,
        'mass': 65.939443 * u.u,
        'stable': False,
        'half-life': 0.194 * u.s,
    },

    'Co-67': {
        'atomic number': 27,
        'mass number': 67,
        'mass': 66.9406096 * u.u,
        'stable': False,
        'half-life': 0.329 * u.s,
    },

    'Co-68': {
        'atomic number': 27,
        'mass number': 68,
        'mass': 67.94426 * u.u,
        'stable': False,
        'half-life': 0.2 * u.s,
    },

    'Co-69': {
        'atomic number': 27,
        'mass number': 69,
        'mass': 68.94614 * u.u,
        'stable': False,
        'half-life': 0.18 * u.s,
    },

    'Co-70': {
        'atomic number': 27,
        'mass number': 70,
        'mass': 69.94963 * u.u,
        'stable': False,
        'half-life': 0.112 * u.s,
    },

    'Co-71': {
        'atomic number': 27,
        'mass number': 71,
        'mass': 70.95237 * u.u,
        'stable': False,
        'half-life': 0.08 * u.s,
    },

    'Co-72': {
        'atomic number': 27,
        'mass number': 72,
        'mass': 71.95729 * u.u,
        'stable': False,
        'half-life': 0.0525 * u.s,
    },

    'Co-73': {
        'atomic number': 27,
        'mass number': 73,
        'mass': 72.96039 * u.u,
        'stable': False,
        'half-life': 0.0407 * u.s,
    },

    'Co-74': {
        'atomic number': 27,
        'mass number': 74,
        'mass': 73.96515 * u.u,
        'stable': False,
        'half-life': 0.0313 * u.s,
    },

    'Co-75': {
        'atomic number': 27,
        'mass number': 75,
        'mass': 74.96876 * u.u,
        'stable': False,
        'half-life': 0.0265 * u.s,
    },

    'Co-76': {
        'atomic number': 27,
        'mass number': 76,
        'mass': 75.97413 * u.u,
        'stable': False,
        'half-life': 0.023 * u.s,
    },

    'Ni-48': {
        'atomic number': 28,
        'mass number': 48,
        'mass': 48.01769 * u.u,
        'stable': False,
        'half-life': 0.0028 * u.s,
    },

    'Ni-49': {
        'atomic number': 28,
        'mass number': 49,
        'mass': 49.0077 * u.u,
        'stable': False,
        'half-life': 0.0075 * u.s,
    },

    'Ni-50': {
        'atomic number': 28,
        'mass number': 50,
        'mass': 49.99474 * u.u,
        'stable': False,
        'half-life': 0.0185 * u.s,
    },

    'Ni-51': {
        'atomic number': 28,
        'mass number': 51,
        'mass': 50.98611 * u.u,
        'stable': False,
        'half-life': 0.0238 * u.s,
    },

    'Ni-52': {
        'atomic number': 28,
        'mass number': 52,
        'mass': 51.9748 * u.u,
        'stable': False,
        'half-life': 0.0418 * u.s,
    },

    'Ni-53': {
        'atomic number': 28,
        'mass number': 53,
        'mass': 52.96819 * u.u,
        'stable': False,
        'half-life': 0.0552 * u.s,
    },

    'Ni-54': {
        'atomic number': 28,
        'mass number': 54,
        'mass': 53.957892 * u.u,
        'stable': False,
        'half-life': 0.1142 * u.s,
    },

    'Ni-55': {
        'atomic number': 28,
        'mass number': 55,
        'mass': 54.95133063 * u.u,
        'stable': False,
        'half-life': 0.2047 * u.s,
    },

    'Ni-56': {
        'atomic number': 28,
        'mass number': 56,
        'mass': 55.94212855 * u.u,
        'stable': False,
        'half-life': 524880.0 * u.s,
    },

    'Ni-57': {
        'atomic number': 28,
        'mass number': 57,
        'mass': 56.93979218 * u.u,
        'stable': False,
        'half-life': 128160.0 * u.s,
    },

    'Ni-58': {
        'atomic number': 28,
        'mass number': 58,
        'mass': 57.93534241 * u.u,
        'stable': True,
        'abundance': 0.68077,
    },

    'Ni-59': {
        'atomic number': 28,
        'mass number': 59,
        'mass': 58.9343462 * u.u,
        'stable': False,
        'half-life': 2556111006000.0 * u.s,
    },

    'Ni-60': {
        'atomic number': 28,
        'mass number': 60,
        'mass': 59.93078588 * u.u,
        'stable': True,
        'abundance': 0.26223,
    },

    'Ni-61': {
        'atomic number': 28,
        'mass number': 61,
        'mass': 60.93105557 * u.u,
        'stable': True,
        'abundance': 0.011399,
    },

    'Ni-62': {
        'atomic number': 28,
        'mass number': 62,
        'mass': 61.92834537 * u.u,
        'stable': True,
        'abundance': 0.036346,
    },

    'Ni-63': {
        'atomic number': 28,
        'mass number': 63,
        'mass': 62.92966963 * u.u,
        'stable': False,
        'half-life': 3193560911.2 * u.s,
    },

    'Ni-64': {
        'atomic number': 28,
        'mass number': 64,
        'mass': 63.92796682 * u.u,
        'stable': True,
        'abundance': 0.009255,
    },

    'Ni-65': {
        'atomic number': 28,
        'mass number': 65,
        'mass': 64.93008517 * u.u,
        'stable': False,
        'half-life': 9063.0 * u.s,
    },

    'Ni-66': {
        'atomic number': 28,
        'mass number': 66,
        'mass': 65.9291393 * u.u,
        'stable': False,
        'half-life': 196560.0 * u.s,
    },

    'Ni-67': {
        'atomic number': 28,
        'mass number': 67,
        'mass': 66.9315694 * u.u,
        'stable': False,
        'half-life': 21.0 * u.s,
    },

    'Ni-68': {
        'atomic number': 28,
        'mass number': 68,
        'mass': 67.9318688 * u.u,
        'stable': False,
        'half-life': 29.0 * u.s,
    },

    'Ni-69': {
        'atomic number': 28,
        'mass number': 69,
        'mass': 68.9356103 * u.u,
        'stable': False,
        'half-life': 11.4 * u.s,
    },

    'Ni-70': {
        'atomic number': 28,
        'mass number': 70,
        'mass': 69.9364313 * u.u,
        'stable': False,
        'half-life': 6.0 * u.s,
    },

    'Ni-71': {
        'atomic number': 28,
        'mass number': 71,
        'mass': 70.940519 * u.u,
        'stable': False,
        'half-life': 2.56 * u.s,
    },

    'Ni-72': {
        'atomic number': 28,
        'mass number': 72,
        'mass': 71.9417859 * u.u,
        'stable': False,
        'half-life': 1.57 * u.s,
    },

    'Ni-73': {
        'atomic number': 28,
        'mass number': 73,
        'mass': 72.9462067 * u.u,
        'stable': False,
        'half-life': 0.84 * u.s,
    },

    'Ni-74': {
        'atomic number': 28,
        'mass number': 74,
        'mass': 73.94798 * u.u,
        'stable': False,
        'half-life': 0.5077 * u.s,
    },

    'Ni-75': {
        'atomic number': 28,
        'mass number': 75,
        'mass': 74.9525 * u.u,
        'stable': False,
        'half-life': 0.3316 * u.s,
    },

    'Ni-76': {
        'atomic number': 28,
        'mass number': 76,
        'mass': 75.95533 * u.u,
        'stable': False,
        'half-life': 0.2346 * u.s,
    },

    'Ni-77': {
        'atomic number': 28,
        'mass number': 77,
        'mass': 76.96055 * u.u,
        'stable': False,
        'half-life': 0.1589 * u.s,
    },

    'Ni-78': {
        'atomic number': 28,
        'mass number': 78,
        'mass': 77.96336 * u.u,
        'stable': False,
        'half-life': 0.1222 * u.s,
    },

    'Ni-79': {
        'atomic number': 28,
        'mass number': 79,
        'mass': 78.97025 * u.u,
        'stable': False,
        'half-life': 0.044 * u.s,
    },

    'Cu-52': {
        'atomic number': 29,
        'mass number': 52,
        'mass': 51.99671 * u.u,
        'stable': False,
    },

    'Cu-53': {
        'atomic number': 29,
        'mass number': 53,
        'mass': 52.98459 * u.u,
        'stable': False,
    },

    'Cu-54': {
        'atomic number': 29,
        'mass number': 54,
        'mass': 53.97666 * u.u,
        'stable': False,
    },

    'Cu-55': {
        'atomic number': 29,
        'mass number': 55,
        'mass': 54.96604 * u.u,
        'stable': False,
        'half-life': 0.057 * u.s,
    },

    'Cu-56': {
        'atomic number': 29,
        'mass number': 56,
        'mass': 55.95895 * u.u,
        'stable': False,
        'half-life': 0.093 * u.s,
    },

    'Cu-57': {
        'atomic number': 29,
        'mass number': 57,
        'mass': 56.9492125 * u.u,
        'stable': False,
        'half-life': 0.1963 * u.s,
    },

    'Cu-58': {
        'atomic number': 29,
        'mass number': 58,
        'mass': 57.94453305 * u.u,
        'stable': False,
        'half-life': 3.204 * u.s,
    },

    'Cu-59': {
        'atomic number': 29,
        'mass number': 59,
        'mass': 58.93949748 * u.u,
        'stable': False,
        'half-life': 81.5 * u.s,
    },

    'Cu-60': {
        'atomic number': 29,
        'mass number': 60,
        'mass': 59.9373645 * u.u,
        'stable': False,
        'half-life': 1422.0 * u.s,
    },

    'Cu-61': {
        'atomic number': 29,
        'mass number': 61,
        'mass': 60.9334576 * u.u,
        'stable': False,
        'half-life': 12020.4 * u.s,
    },

    'Cu-62': {
        'atomic number': 29,
        'mass number': 62,
        'mass': 61.93259541 * u.u,
        'stable': False,
        'half-life': 580.2 * u.s,
    },

    'Cu-63': {
        'atomic number': 29,
        'mass number': 63,
        'mass': 62.92959772 * u.u,
        'stable': True,
        'abundance': 0.6915,
    },

    'Cu-64': {
        'atomic number': 29,
        'mass number': 64,
        'mass': 63.92976434 * u.u,
        'stable': False,
        'half-life': 45721.44 * u.s,
    },

    'Cu-65': {
        'atomic number': 29,
        'mass number': 65,
        'mass': 64.9277897 * u.u,
        'stable': True,
        'abundance': 0.3085,
    },

    'Cu-66': {
        'atomic number': 29,
        'mass number': 66,
        'mass': 65.92886903 * u.u,
        'stable': False,
        'half-life': 307.2 * u.s,
    },

    'Cu-67': {
        'atomic number': 29,
        'mass number': 67,
        'mass': 66.9277303 * u.u,
        'stable': False,
        'half-life': 222588.0 * u.s,
    },

    'Cu-68': {
        'atomic number': 29,
        'mass number': 68,
        'mass': 67.9296109 * u.u,
        'stable': False,
        'half-life': 30.9 * u.s,
    },

    'Cu-69': {
        'atomic number': 29,
        'mass number': 69,
        'mass': 68.9294293 * u.u,
        'stable': False,
        'half-life': 171.0 * u.s,
    },

    'Cu-70': {
        'atomic number': 29,
        'mass number': 70,
        'mass': 69.9323921 * u.u,
        'stable': False,
        'half-life': 44.5 * u.s,
    },

    'Cu-71': {
        'atomic number': 29,
        'mass number': 71,
        'mass': 70.9326768 * u.u,
        'stable': False,
        'half-life': 19.4 * u.s,
    },

    'Cu-72': {
        'atomic number': 29,
        'mass number': 72,
        'mass': 71.9358203 * u.u,
        'stable': False,
        'half-life': 6.63 * u.s,
    },

    'Cu-73': {
        'atomic number': 29,
        'mass number': 73,
        'mass': 72.9366744 * u.u,
        'stable': False,
        'half-life': 4.2 * u.s,
    },

    'Cu-74': {
        'atomic number': 29,
        'mass number': 74,
        'mass': 73.9398749 * u.u,
        'stable': False,
        'half-life': 1.63 * u.s,
    },

    'Cu-75': {
        'atomic number': 29,
        'mass number': 75,
        'mass': 74.9415226 * u.u,
        'stable': False,
        'half-life': 1.224 * u.s,
    },

    'Cu-76': {
        'atomic number': 29,
        'mass number': 76,
        'mass': 75.945275 * u.u,
        'stable': False,
        'half-life': 0.6377 * u.s,
    },

    'Cu-77': {
        'atomic number': 29,
        'mass number': 77,
        'mass': 76.94792 * u.u,
        'stable': False,
        'half-life': 0.4679 * u.s,
    },

    'Cu-78': {
        'atomic number': 29,
        'mass number': 78,
        'mass': 77.95223 * u.u,
        'stable': False,
        'half-life': 0.3307 * u.s,
    },

    'Cu-79': {
        'atomic number': 29,
        'mass number': 79,
        'mass': 78.95502 * u.u,
        'stable': False,
        'half-life': 0.241 * u.s,
    },

    'Cu-80': {
        'atomic number': 29,
        'mass number': 80,
        'mass': 79.96089 * u.u,
        'stable': False,
        'half-life': 0.1133 * u.s,
    },

    'Cu-81': {
        'atomic number': 29,
        'mass number': 81,
        'mass': 80.96587 * u.u,
        'stable': False,
        'half-life': 0.0732 * u.s,
    },

    'Cu-82': {
        'atomic number': 29,
        'mass number': 82,
        'mass': 81.97244 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'Zn-54': {
        'atomic number': 30,
        'mass number': 54,
        'mass': 53.99204 * u.u,
        'stable': False,
        'half-life': 0.0018 * u.s,
    },

    'Zn-55': {
        'atomic number': 30,
        'mass number': 55,
        'mass': 54.98398 * u.u,
        'stable': False,
        'half-life': 0.0198 * u.s,
    },

    'Zn-56': {
        'atomic number': 30,
        'mass number': 56,
        'mass': 55.97254 * u.u,
        'stable': False,
        'half-life': 0.0329 * u.s,
    },

    'Zn-57': {
        'atomic number': 30,
        'mass number': 57,
        'mass': 56.96506 * u.u,
        'stable': False,
        'half-life': 0.038 * u.s,
    },

    'Zn-58': {
        'atomic number': 30,
        'mass number': 58,
        'mass': 57.954591 * u.u,
        'stable': False,
        'half-life': 0.0867 * u.s,
    },

    'Zn-59': {
        'atomic number': 30,
        'mass number': 59,
        'mass': 58.94931266 * u.u,
        'stable': False,
        'half-life': 0.182 * u.s,
    },

    'Zn-60': {
        'atomic number': 30,
        'mass number': 60,
        'mass': 59.9418421 * u.u,
        'stable': False,
        'half-life': 142.8 * u.s,
    },

    'Zn-61': {
        'atomic number': 30,
        'mass number': 61,
        'mass': 60.939507 * u.u,
        'stable': False,
        'half-life': 89.1 * u.s,
    },

    'Zn-62': {
        'atomic number': 30,
        'mass number': 62,
        'mass': 61.93433397 * u.u,
        'stable': False,
        'half-life': 33094.8 * u.s,
    },

    'Zn-63': {
        'atomic number': 30,
        'mass number': 63,
        'mass': 62.9332115 * u.u,
        'stable': False,
        'half-life': 2308.2 * u.s,
    },

    'Zn-64': {
        'atomic number': 30,
        'mass number': 64,
        'mass': 63.92914201 * u.u,
        'stable': True,
        'abundance': 0.4917,
    },

    'Zn-65': {
        'atomic number': 30,
        'mass number': 65,
        'mass': 64.92924077 * u.u,
        'stable': False,
        'half-life': 21095769.599999998 * u.s,
    },

    'Zn-66': {
        'atomic number': 30,
        'mass number': 66,
        'mass': 65.92603381 * u.u,
        'stable': True,
        'abundance': 0.2773,
    },

    'Zn-67': {
        'atomic number': 30,
        'mass number': 67,
        'mass': 66.92712775 * u.u,
        'stable': True,
        'abundance': 0.0404,
    },

    'Zn-68': {
        'atomic number': 30,
        'mass number': 68,
        'mass': 67.92484455 * u.u,
        'stable': True,
        'abundance': 0.1845,
    },

    'Zn-69': {
        'atomic number': 30,
        'mass number': 69,
        'mass': 68.9265507 * u.u,
        'stable': False,
        'half-life': 3384.0 * u.s,
    },

    'Zn-70': {
        'atomic number': 30,
        'mass number': 70,
        'mass': 69.9253192 * u.u,
        'stable': True,
        'abundance': 0.0061,
    },

    'Zn-71': {
        'atomic number': 30,
        'mass number': 71,
        'mass': 70.9277196 * u.u,
        'stable': False,
        'half-life': 147.0 * u.s,
    },

    'Zn-72': {
        'atomic number': 30,
        'mass number': 72,
        'mass': 71.9268428 * u.u,
        'stable': False,
        'half-life': 167400.0 * u.s,
    },

    'Zn-73': {
        'atomic number': 30,
        'mass number': 73,
        'mass': 72.9295826 * u.u,
        'stable': False,
        'half-life': 23.5 * u.s,
    },

    'Zn-74': {
        'atomic number': 30,
        'mass number': 74,
        'mass': 73.9294073 * u.u,
        'stable': False,
        'half-life': 95.6 * u.s,
    },

    'Zn-75': {
        'atomic number': 30,
        'mass number': 75,
        'mass': 74.9328402 * u.u,
        'stable': False,
        'half-life': 10.2 * u.s,
    },

    'Zn-76': {
        'atomic number': 30,
        'mass number': 76,
        'mass': 75.933115 * u.u,
        'stable': False,
        'half-life': 5.7 * u.s,
    },

    'Zn-77': {
        'atomic number': 30,
        'mass number': 77,
        'mass': 76.9368872 * u.u,
        'stable': False,
        'half-life': 2.08 * u.s,
    },

    'Zn-78': {
        'atomic number': 30,
        'mass number': 78,
        'mass': 77.9382892 * u.u,
        'stable': False,
        'half-life': 1.47 * u.s,
    },

    'Zn-79': {
        'atomic number': 30,
        'mass number': 79,
        'mass': 78.9426381 * u.u,
        'stable': False,
        'half-life': 0.746 * u.s,
    },

    'Zn-80': {
        'atomic number': 30,
        'mass number': 80,
        'mass': 79.9445529 * u.u,
        'stable': False,
        'half-life': 0.5622 * u.s,
    },

    'Zn-81': {
        'atomic number': 30,
        'mass number': 81,
        'mass': 80.9504026 * u.u,
        'stable': False,
        'half-life': 0.3032 * u.s,
    },

    'Zn-82': {
        'atomic number': 30,
        'mass number': 82,
        'mass': 81.95426 * u.u,
        'stable': False,
        'half-life': 0.1779 * u.s,
    },

    'Zn-83': {
        'atomic number': 30,
        'mass number': 83,
        'mass': 82.96056 * u.u,
        'stable': False,
        'half-life': 0.119 * u.s,
    },

    'Zn-84': {
        'atomic number': 30,
        'mass number': 84,
        'mass': 83.96521 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'Zn-85': {
        'atomic number': 30,
        'mass number': 85,
        'mass': 84.97226 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'Ga-56': {
        'atomic number': 31,
        'mass number': 56,
        'mass': 55.99536 * u.u,
        'stable': False,
    },

    'Ga-57': {
        'atomic number': 31,
        'mass number': 57,
        'mass': 56.9832 * u.u,
        'stable': False,
    },

    'Ga-58': {
        'atomic number': 31,
        'mass number': 58,
        'mass': 57.97478 * u.u,
        'stable': False,
    },

    'Ga-59': {
        'atomic number': 31,
        'mass number': 59,
        'mass': 58.96353 * u.u,
        'stable': False,
    },

    'Ga-60': {
        'atomic number': 31,
        'mass number': 60,
        'mass': 59.95729 * u.u,
        'stable': False,
        'half-life': 0.07 * u.s,
    },

    'Ga-61': {
        'atomic number': 31,
        'mass number': 61,
        'mass': 60.949399 * u.u,
        'stable': False,
        'half-life': 0.167 * u.s,
    },

    'Ga-62': {
        'atomic number': 31,
        'mass number': 62,
        'mass': 61.94419025 * u.u,
        'stable': False,
        'half-life': 0.116121 * u.s,
    },

    'Ga-63': {
        'atomic number': 31,
        'mass number': 63,
        'mass': 62.9392942 * u.u,
        'stable': False,
        'half-life': 32.4 * u.s,
    },

    'Ga-64': {
        'atomic number': 31,
        'mass number': 64,
        'mass': 63.9368404 * u.u,
        'stable': False,
        'half-life': 157.62 * u.s,
    },

    'Ga-65': {
        'atomic number': 31,
        'mass number': 65,
        'mass': 64.93273459 * u.u,
        'stable': False,
        'half-life': 912.0 * u.s,
    },

    'Ga-66': {
        'atomic number': 31,
        'mass number': 66,
        'mass': 65.9315894 * u.u,
        'stable': False,
        'half-life': 33494.4 * u.s,
    },

    'Ga-67': {
        'atomic number': 31,
        'mass number': 67,
        'mass': 66.9282025 * u.u,
        'stable': False,
        'half-life': 281797.056 * u.s,
    },

    'Ga-68': {
        'atomic number': 31,
        'mass number': 68,
        'mass': 67.9279805 * u.u,
        'stable': False,
        'half-life': 4070.7 * u.s,
    },

    'Ga-69': {
        'atomic number': 31,
        'mass number': 69,
        'mass': 68.9255735 * u.u,
        'stable': True,
        'abundance': 0.60108,
    },

    'Ga-70': {
        'atomic number': 31,
        'mass number': 70,
        'mass': 69.9260219 * u.u,
        'stable': False,
        'half-life': 1268.4 * u.s,
    },

    'Ga-71': {
        'atomic number': 31,
        'mass number': 71,
        'mass': 70.92470258 * u.u,
        'stable': True,
        'abundance': 0.39892,
    },

    'Ga-72': {
        'atomic number': 31,
        'mass number': 72,
        'mass': 71.92636747 * u.u,
        'stable': False,
        'half-life': 50490.0 * u.s,
    },

    'Ga-73': {
        'atomic number': 31,
        'mass number': 73,
        'mass': 72.9251747 * u.u,
        'stable': False,
        'half-life': 17496.0 * u.s,
    },

    'Ga-74': {
        'atomic number': 31,
        'mass number': 74,
        'mass': 73.9269457 * u.u,
        'stable': False,
        'half-life': 487.2 * u.s,
    },

    'Ga-75': {
        'atomic number': 31,
        'mass number': 75,
        'mass': 74.9265002 * u.u,
        'stable': False,
        'half-life': 126.0 * u.s,
    },

    'Ga-76': {
        'atomic number': 31,
        'mass number': 76,
        'mass': 75.9288276 * u.u,
        'stable': False,
        'half-life': 32.6 * u.s,
    },

    'Ga-77': {
        'atomic number': 31,
        'mass number': 77,
        'mass': 76.9291543 * u.u,
        'stable': False,
        'half-life': 13.2 * u.s,
    },

    'Ga-78': {
        'atomic number': 31,
        'mass number': 78,
        'mass': 77.9316088 * u.u,
        'stable': False,
        'half-life': 5.09 * u.s,
    },

    'Ga-79': {
        'atomic number': 31,
        'mass number': 79,
        'mass': 78.9328523 * u.u,
        'stable': False,
        'half-life': 2.848 * u.s,
    },

    'Ga-80': {
        'atomic number': 31,
        'mass number': 80,
        'mass': 79.9364208 * u.u,
        'stable': False,
        'half-life': 1.9 * u.s,
    },

    'Ga-81': {
        'atomic number': 31,
        'mass number': 81,
        'mass': 80.9381338 * u.u,
        'stable': False,
        'half-life': 1.217 * u.s,
    },

    'Ga-82': {
        'atomic number': 31,
        'mass number': 82,
        'mass': 81.9431765 * u.u,
        'stable': False,
        'half-life': 0.599 * u.s,
    },

    'Ga-83': {
        'atomic number': 31,
        'mass number': 83,
        'mass': 82.9471203 * u.u,
        'stable': False,
        'half-life': 0.3081 * u.s,
    },

    'Ga-84': {
        'atomic number': 31,
        'mass number': 84,
        'mass': 83.95246 * u.u,
        'stable': False,
        'half-life': 0.085 * u.s,
    },

    'Ga-85': {
        'atomic number': 31,
        'mass number': 85,
        'mass': 84.95699 * u.u,
        'stable': False,
        'half-life': 0.0922 * u.s,
    },

    'Ga-86': {
        'atomic number': 31,
        'mass number': 86,
        'mass': 85.96301 * u.u,
        'stable': False,
        'half-life': 0.047 * u.s,
    },

    'Ga-87': {
        'atomic number': 31,
        'mass number': 87,
        'mass': 86.96824 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Ge-58': {
        'atomic number': 32,
        'mass number': 58,
        'mass': 57.99172 * u.u,
        'stable': False,
    },

    'Ge-59': {
        'atomic number': 32,
        'mass number': 59,
        'mass': 58.98249 * u.u,
        'stable': False,
        'half-life': '8# ms',
    },

    'Ge-60': {
        'atomic number': 32,
        'mass number': 60,
        'mass': 59.97036 * u.u,
        'stable': False,
        'half-life': '30# ms',
    },

    'Ge-61': {
        'atomic number': 32,
        'mass number': 61,
        'mass': 60.96379 * u.u,
        'stable': False,
        'half-life': 0.044 * u.s,
    },

    'Ge-62': {
        'atomic number': 32,
        'mass number': 62,
        'mass': 61.95502 * u.u,
        'stable': False,
        'half-life': 0.129 * u.s,
    },

    'Ge-63': {
        'atomic number': 32,
        'mass number': 63,
        'mass': 62.949628 * u.u,
        'stable': False,
        'half-life': 0.142 * u.s,
    },

    'Ge-64': {
        'atomic number': 32,
        'mass number': 64,
        'mass': 63.9416899 * u.u,
        'stable': False,
        'half-life': 63.7 * u.s,
    },

    'Ge-65': {
        'atomic number': 32,
        'mass number': 65,
        'mass': 64.9393681 * u.u,
        'stable': False,
        'half-life': 30.9 * u.s,
    },

    'Ge-66': {
        'atomic number': 32,
        'mass number': 66,
        'mass': 65.9338621 * u.u,
        'stable': False,
        'half-life': 8136.0 * u.s,
    },

    'Ge-67': {
        'atomic number': 32,
        'mass number': 67,
        'mass': 66.9327339 * u.u,
        'stable': False,
        'half-life': 1134.0 * u.s,
    },

    'Ge-68': {
        'atomic number': 32,
        'mass number': 68,
        'mass': 67.9280953 * u.u,
        'stable': False,
        'half-life': 23408352.0 * u.s,
    },

    'Ge-69': {
        'atomic number': 32,
        'mass number': 69,
        'mass': 68.9279645 * u.u,
        'stable': False,
        'half-life': 140580.0 * u.s,
    },

    'Ge-70': {
        'atomic number': 32,
        'mass number': 70,
        'mass': 69.92424875 * u.u,
        'stable': True,
        'abundance': 0.2057,
    },

    'Ge-71': {
        'atomic number': 32,
        'mass number': 71,
        'mass': 70.92495233 * u.u,
        'stable': False,
        'half-life': 987552.0 * u.s,
    },

    'Ge-72': {
        'atomic number': 32,
        'mass number': 72,
        'mass': 71.922075826 * u.u,
        'stable': True,
        'abundance': 0.2745,
    },

    'Ge-73': {
        'atomic number': 32,
        'mass number': 73,
        'mass': 72.923458956 * u.u,
        'stable': True,
        'abundance': 0.0775,
    },

    'Ge-74': {
        'atomic number': 32,
        'mass number': 74,
        'mass': 73.921177761 * u.u,
        'stable': True,
        'abundance': 0.365,
    },

    'Ge-75': {
        'atomic number': 32,
        'mass number': 75,
        'mass': 74.92285837 * u.u,
        'stable': False,
        'half-life': 4966.8 * u.s,
    },

    'Ge-76': {
        'atomic number': 32,
        'mass number': 76,
        'mass': 75.921402726 * u.u,
        'stable': False,
        'abundance': 0.0773,
    },

    'Ge-77': {
        'atomic number': 32,
        'mass number': 77,
        'mass': 76.923549843 * u.u,
        'stable': False,
        'half-life': 40359.6 * u.s,
    },

    'Ge-78': {
        'atomic number': 32,
        'mass number': 78,
        'mass': 77.9228529 * u.u,
        'stable': False,
        'half-life': 5280.0 * u.s,
    },

    'Ge-79': {
        'atomic number': 32,
        'mass number': 79,
        'mass': 78.92536 * u.u,
        'stable': False,
        'half-life': 18.98 * u.s,
    },

    'Ge-80': {
        'atomic number': 32,
        'mass number': 80,
        'mass': 79.9253508 * u.u,
        'stable': False,
        'half-life': 29.5 * u.s,
    },

    'Ge-81': {
        'atomic number': 32,
        'mass number': 81,
        'mass': 80.9288329 * u.u,
        'stable': False,
        'half-life': 8.0 * u.s,
    },

    'Ge-82': {
        'atomic number': 32,
        'mass number': 82,
        'mass': 81.929774 * u.u,
        'stable': False,
        'half-life': 4.56 * u.s,
    },

    'Ge-83': {
        'atomic number': 32,
        'mass number': 83,
        'mass': 82.9345391 * u.u,
        'stable': False,
        'half-life': 1.85 * u.s,
    },

    'Ge-84': {
        'atomic number': 32,
        'mass number': 84,
        'mass': 83.9375751 * u.u,
        'stable': False,
        'half-life': 0.951 * u.s,
    },

    'Ge-85': {
        'atomic number': 32,
        'mass number': 85,
        'mass': 84.9429697 * u.u,
        'stable': False,
        'half-life': 0.494 * u.s,
    },

    'Ge-86': {
        'atomic number': 32,
        'mass number': 86,
        'mass': 85.94658 * u.u,
        'stable': False,
        'half-life': 0.226 * u.s,
    },

    'Ge-87': {
        'atomic number': 32,
        'mass number': 87,
        'mass': 86.95268 * u.u,
        'stable': False,
        'half-life': '150# ms',
    },

    'Ge-88': {
        'atomic number': 32,
        'mass number': 88,
        'mass': 87.95691 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Ge-89': {
        'atomic number': 32,
        'mass number': 89,
        'mass': 88.96379 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'Ge-90': {
        'atomic number': 32,
        'mass number': 90,
        'mass': 89.96863 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'As-60': {
        'atomic number': 33,
        'mass number': 60,
        'mass': 59.99388 * u.u,
        'stable': False,
    },

    'As-61': {
        'atomic number': 33,
        'mass number': 61,
        'mass': 60.98112 * u.u,
        'stable': False,
    },

    'As-62': {
        'atomic number': 33,
        'mass number': 62,
        'mass': 61.97361 * u.u,
        'stable': False,
    },

    'As-63': {
        'atomic number': 33,
        'mass number': 63,
        'mass': 62.9639 * u.u,
        'stable': False,
    },

    'As-64': {
        'atomic number': 33,
        'mass number': 64,
        'mass': 63.95743 * u.u,
        'stable': False,
        'half-life': 0.04 * u.s,
    },

    'As-65': {
        'atomic number': 33,
        'mass number': 65,
        'mass': 64.949611 * u.u,
        'stable': False,
        'half-life': 0.17 * u.s,
    },

    'As-66': {
        'atomic number': 33,
        'mass number': 66,
        'mass': 65.9441488 * u.u,
        'stable': False,
        'half-life': 0.09577 * u.s,
    },

    'As-67': {
        'atomic number': 33,
        'mass number': 67,
        'mass': 66.93925111 * u.u,
        'stable': False,
        'half-life': 42.5 * u.s,
    },

    'As-68': {
        'atomic number': 33,
        'mass number': 68,
        'mass': 67.9367741 * u.u,
        'stable': False,
        'half-life': 151.6 * u.s,
    },

    'As-69': {
        'atomic number': 33,
        'mass number': 69,
        'mass': 68.932246 * u.u,
        'stable': False,
        'half-life': 912.0 * u.s,
    },

    'As-70': {
        'atomic number': 33,
        'mass number': 70,
        'mass': 69.930926 * u.u,
        'stable': False,
        'half-life': 3156.0 * u.s,
    },

    'As-71': {
        'atomic number': 33,
        'mass number': 71,
        'mass': 70.9271138 * u.u,
        'stable': False,
        'half-life': 235080.0 * u.s,
    },

    'As-72': {
        'atomic number': 33,
        'mass number': 72,
        'mass': 71.9267523 * u.u,
        'stable': False,
        'half-life': 93600.0 * u.s,
    },

    'As-73': {
        'atomic number': 33,
        'mass number': 73,
        'mass': 72.9238291 * u.u,
        'stable': False,
        'half-life': 6937920.0 * u.s,
    },

    'As-74': {
        'atomic number': 33,
        'mass number': 74,
        'mass': 73.9239286 * u.u,
        'stable': False,
        'half-life': 1535328.0 * u.s,
    },

    'As-75': {
        'atomic number': 33,
        'mass number': 75,
        'mass': 74.92159457 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'As-76': {
        'atomic number': 33,
        'mass number': 76,
        'mass': 75.92239202 * u.u,
        'stable': False,
        'half-life': 93121.92 * u.s,
    },

    'As-77': {
        'atomic number': 33,
        'mass number': 77,
        'mass': 76.9206476 * u.u,
        'stable': False,
        'half-life': 139644.0 * u.s,
    },

    'As-78': {
        'atomic number': 33,
        'mass number': 78,
        'mass': 77.921828 * u.u,
        'stable': False,
        'half-life': 5442.0 * u.s,
    },

    'As-79': {
        'atomic number': 33,
        'mass number': 79,
        'mass': 78.9209484 * u.u,
        'stable': False,
        'half-life': 540.6 * u.s,
    },

    'As-80': {
        'atomic number': 33,
        'mass number': 80,
        'mass': 79.9224746 * u.u,
        'stable': False,
        'half-life': 15.2 * u.s,
    },

    'As-81': {
        'atomic number': 33,
        'mass number': 81,
        'mass': 80.9221323 * u.u,
        'stable': False,
        'half-life': 33.3 * u.s,
    },

    'As-82': {
        'atomic number': 33,
        'mass number': 82,
        'mass': 81.9247412 * u.u,
        'stable': False,
        'half-life': 19.1 * u.s,
    },

    'As-83': {
        'atomic number': 33,
        'mass number': 83,
        'mass': 82.9252069 * u.u,
        'stable': False,
        'half-life': 13.4 * u.s,
    },

    'As-84': {
        'atomic number': 33,
        'mass number': 84,
        'mass': 83.9293033 * u.u,
        'stable': False,
        'half-life': 4.02 * u.s,
    },

    'As-85': {
        'atomic number': 33,
        'mass number': 85,
        'mass': 84.9321637 * u.u,
        'stable': False,
        'half-life': 2.021 * u.s,
    },

    'As-86': {
        'atomic number': 33,
        'mass number': 86,
        'mass': 85.9367015 * u.u,
        'stable': False,
        'half-life': 0.945 * u.s,
    },

    'As-87': {
        'atomic number': 33,
        'mass number': 87,
        'mass': 86.9402917 * u.u,
        'stable': False,
        'half-life': 0.492 * u.s,
    },

    'As-88': {
        'atomic number': 33,
        'mass number': 88,
        'mass': 87.94555 * u.u,
        'stable': False,
        'half-life': 0.27 * u.s,
    },

    'As-89': {
        'atomic number': 33,
        'mass number': 89,
        'mass': 88.94976 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'As-90': {
        'atomic number': 33,
        'mass number': 90,
        'mass': 89.95563 * u.u,
        'stable': False,
        'half-life': '80# ms',
    },

    'As-91': {
        'atomic number': 33,
        'mass number': 91,
        'mass': 90.96039 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'As-92': {
        'atomic number': 33,
        'mass number': 92,
        'mass': 91.96674 * u.u,
        'stable': False,
        'half-life': '30# ms',
    },

    'Se-64': {
        'atomic number': 34,
        'mass number': 64,
        'mass': 63.97109 * u.u,
        'stable': False,
        'half-life': '30# ms',
    },

    'Se-65': {
        'atomic number': 34,
        'mass number': 65,
        'mass': 64.9644 * u.u,
        'stable': False,
        'half-life': 0.033 * u.s,
    },

    'Se-66': {
        'atomic number': 34,
        'mass number': 66,
        'mass': 65.95559 * u.u,
        'stable': False,
        'half-life': 0.033 * u.s,
    },

    'Se-67': {
        'atomic number': 34,
        'mass number': 67,
        'mass': 66.949994 * u.u,
        'stable': False,
        'half-life': 0.133 * u.s,
    },

    'Se-68': {
        'atomic number': 34,
        'mass number': 68,
        'mass': 67.94182524 * u.u,
        'stable': False,
        'half-life': 35.5 * u.s,
    },

    'Se-69': {
        'atomic number': 34,
        'mass number': 69,
        'mass': 68.9394148 * u.u,
        'stable': False,
        'half-life': 27.4 * u.s,
    },

    'Se-70': {
        'atomic number': 34,
        'mass number': 70,
        'mass': 69.9335155 * u.u,
        'stable': False,
        'half-life': 2466.0 * u.s,
    },

    'Se-71': {
        'atomic number': 34,
        'mass number': 71,
        'mass': 70.9322094 * u.u,
        'stable': False,
        'half-life': 284.4 * u.s,
    },

    'Se-72': {
        'atomic number': 34,
        'mass number': 72,
        'mass': 71.9271405 * u.u,
        'stable': False,
        'half-life': 725760.0 * u.s,
    },

    'Se-73': {
        'atomic number': 34,
        'mass number': 73,
        'mass': 72.9267549 * u.u,
        'stable': False,
        'half-life': 25740.0 * u.s,
    },

    'Se-74': {
        'atomic number': 34,
        'mass number': 74,
        'mass': 73.922475934 * u.u,
        'stable': True,
        'abundance': 0.0089,
    },

    'Se-75': {
        'atomic number': 34,
        'mass number': 75,
        'mass': 74.92252287 * u.u,
        'stable': False,
        'half-life': 10351497.6 * u.s,
    },

    'Se-76': {
        'atomic number': 34,
        'mass number': 76,
        'mass': 75.919213704 * u.u,
        'stable': True,
        'abundance': 0.0937,
    },

    'Se-77': {
        'atomic number': 34,
        'mass number': 77,
        'mass': 76.919914154 * u.u,
        'stable': True,
        'abundance': 0.0763,
    },

    'Se-78': {
        'atomic number': 34,
        'mass number': 78,
        'mass': 77.91730928 * u.u,
        'stable': True,
        'abundance': 0.2377,
    },

    'Se-79': {
        'atomic number': 34,
        'mass number': 79,
        'mass': 78.91849929 * u.u,
        'stable': False,
        'half-life': 10319114802000.0 * u.s,
    },

    'Se-80': {
        'atomic number': 34,
        'mass number': 80,
        'mass': 79.9165218 * u.u,
        'stable': True,
        'abundance': 0.4961,
    },

    'Se-81': {
        'atomic number': 34,
        'mass number': 81,
        'mass': 80.917993 * u.u,
        'stable': False,
        'half-life': 1107.0 * u.s,
    },

    'Se-82': {
        'atomic number': 34,
        'mass number': 82,
        'mass': 81.9166995 * u.u,
        'stable': False,
        'abundance': 0.0873,
    },

    'Se-83': {
        'atomic number': 34,
        'mass number': 83,
        'mass': 82.9191186 * u.u,
        'stable': False,
        'half-life': 1335.0 * u.s,
    },

    'Se-84': {
        'atomic number': 34,
        'mass number': 84,
        'mass': 83.9184668 * u.u,
        'stable': False,
        'half-life': 195.6 * u.s,
    },

    'Se-85': {
        'atomic number': 34,
        'mass number': 85,
        'mass': 84.9222608 * u.u,
        'stable': False,
        'half-life': 32.9 * u.s,
    },

    'Se-86': {
        'atomic number': 34,
        'mass number': 86,
        'mass': 85.9243117 * u.u,
        'stable': False,
        'half-life': 14.3 * u.s,
    },

    'Se-87': {
        'atomic number': 34,
        'mass number': 87,
        'mass': 86.9286886 * u.u,
        'stable': False,
        'half-life': 5.5 * u.s,
    },

    'Se-88': {
        'atomic number': 34,
        'mass number': 88,
        'mass': 87.9314175 * u.u,
        'stable': False,
        'half-life': 1.53 * u.s,
    },

    'Se-89': {
        'atomic number': 34,
        'mass number': 89,
        'mass': 88.9366691 * u.u,
        'stable': False,
        'half-life': 0.43 * u.s,
    },

    'Se-90': {
        'atomic number': 34,
        'mass number': 90,
        'mass': 89.9401 * u.u,
        'stable': False,
        'half-life': 0.21 * u.s,
    },

    'Se-91': {
        'atomic number': 34,
        'mass number': 91,
        'mass': 90.94596 * u.u,
        'stable': False,
        'half-life': 0.27 * u.s,
    },

    'Se-92': {
        'atomic number': 34,
        'mass number': 92,
        'mass': 91.94984 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Se-93': {
        'atomic number': 34,
        'mass number': 93,
        'mass': 92.95629 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'Se-94': {
        'atomic number': 34,
        'mass number': 94,
        'mass': 93.96049 * u.u,
        'stable': False,
        'half-life': '20# ms',
    },

    'Se-95': {
        'atomic number': 34,
        'mass number': 95,
        'mass': 94.9673 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Br-67': {
        'atomic number': 35,
        'mass number': 67,
        'mass': 66.96465 * u.u,
        'stable': False,
    },

    'Br-68': {
        'atomic number': 35,
        'mass number': 68,
        'mass': 67.95873 * u.u,
        'stable': False,
    },

    'Br-69': {
        'atomic number': 35,
        'mass number': 69,
        'mass': 68.950497 * u.u,
        'stable': False,
        'half-life': '<24 ns',
    },

    'Br-70': {
        'atomic number': 35,
        'mass number': 70,
        'mass': 69.944792 * u.u,
        'stable': False,
        'half-life': 0.0791 * u.s,
    },

    'Br-71': {
        'atomic number': 35,
        'mass number': 71,
        'mass': 70.9393422 * u.u,
        'stable': False,
        'half-life': 21.4 * u.s,
    },

    'Br-72': {
        'atomic number': 35,
        'mass number': 72,
        'mass': 71.9365886 * u.u,
        'stable': False,
        'half-life': 78.6 * u.s,
    },

    'Br-73': {
        'atomic number': 35,
        'mass number': 73,
        'mass': 72.9316715 * u.u,
        'stable': False,
        'half-life': 204.0 * u.s,
    },

    'Br-74': {
        'atomic number': 35,
        'mass number': 74,
        'mass': 73.9299102 * u.u,
        'stable': False,
        'half-life': 1524.0 * u.s,
    },

    'Br-75': {
        'atomic number': 35,
        'mass number': 75,
        'mass': 74.9258105 * u.u,
        'stable': False,
        'half-life': 5802.0 * u.s,
    },

    'Br-76': {
        'atomic number': 35,
        'mass number': 76,
        'mass': 75.924542 * u.u,
        'stable': False,
        'half-life': 58320.0 * u.s,
    },

    'Br-77': {
        'atomic number': 35,
        'mass number': 77,
        'mass': 76.9213792 * u.u,
        'stable': False,
        'half-life': 205344.0 * u.s,
    },

    'Br-78': {
        'atomic number': 35,
        'mass number': 78,
        'mass': 77.9211459 * u.u,
        'stable': False,
        'half-life': 387.0 * u.s,
    },

    'Br-79': {
        'atomic number': 35,
        'mass number': 79,
        'mass': 78.9183376 * u.u,
        'stable': True,
        'abundance': 0.5069,
    },

    'Br-80': {
        'atomic number': 35,
        'mass number': 80,
        'mass': 79.9185298 * u.u,
        'stable': False,
        'half-life': 1060.8 * u.s,
    },

    'Br-81': {
        'atomic number': 35,
        'mass number': 81,
        'mass': 80.9162897 * u.u,
        'stable': True,
        'abundance': 0.4931,
    },

    'Br-82': {
        'atomic number': 35,
        'mass number': 82,
        'mass': 81.9168032 * u.u,
        'stable': False,
        'half-life': 127015.2 * u.s,
    },

    'Br-83': {
        'atomic number': 35,
        'mass number': 83,
        'mass': 82.9151756 * u.u,
        'stable': False,
        'half-life': 8546.4 * u.s,
    },

    'Br-84': {
        'atomic number': 35,
        'mass number': 84,
        'mass': 83.916496 * u.u,
        'stable': False,
        'half-life': 1905.6 * u.s,
    },

    'Br-85': {
        'atomic number': 35,
        'mass number': 85,
        'mass': 84.9156458 * u.u,
        'stable': False,
        'half-life': 174.0 * u.s,
    },

    'Br-86': {
        'atomic number': 35,
        'mass number': 86,
        'mass': 85.9188054 * u.u,
        'stable': False,
        'half-life': 55.1 * u.s,
    },

    'Br-87': {
        'atomic number': 35,
        'mass number': 87,
        'mass': 86.920674 * u.u,
        'stable': False,
        'half-life': 55.65 * u.s,
    },

    'Br-88': {
        'atomic number': 35,
        'mass number': 88,
        'mass': 87.9240833 * u.u,
        'stable': False,
        'half-life': 16.34 * u.s,
    },

    'Br-89': {
        'atomic number': 35,
        'mass number': 89,
        'mass': 88.9267046 * u.u,
        'stable': False,
        'half-life': 4.357 * u.s,
    },

    'Br-90': {
        'atomic number': 35,
        'mass number': 90,
        'mass': 89.9312928 * u.u,
        'stable': False,
        'half-life': 1.91 * u.s,
    },

    'Br-91': {
        'atomic number': 35,
        'mass number': 91,
        'mass': 90.9343986 * u.u,
        'stable': False,
        'half-life': 0.543 * u.s,
    },

    'Br-92': {
        'atomic number': 35,
        'mass number': 92,
        'mass': 91.9396316 * u.u,
        'stable': False,
        'half-life': 0.314 * u.s,
    },

    'Br-93': {
        'atomic number': 35,
        'mass number': 93,
        'mass': 92.94313 * u.u,
        'stable': False,
        'half-life': 0.152 * u.s,
    },

    'Br-94': {
        'atomic number': 35,
        'mass number': 94,
        'mass': 93.9489 * u.u,
        'stable': False,
        'half-life': 0.07 * u.s,
    },

    'Br-95': {
        'atomic number': 35,
        'mass number': 95,
        'mass': 94.95301 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'Br-96': {
        'atomic number': 35,
        'mass number': 96,
        'mass': 95.95903 * u.u,
        'stable': False,
        'half-life': '20# ms',
    },

    'Br-97': {
        'atomic number': 35,
        'mass number': 97,
        'mass': 96.96344 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Br-98': {
        'atomic number': 35,
        'mass number': 98,
        'mass': 97.96946 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Kr-69': {
        'atomic number': 36,
        'mass number': 69,
        'mass': 68.96518 * u.u,
        'stable': False,
        'half-life': 0.028 * u.s,
    },

    'Kr-70': {
        'atomic number': 36,
        'mass number': 70,
        'mass': 69.95604 * u.u,
        'stable': False,
        'half-life': 0.052 * u.s,
    },

    'Kr-71': {
        'atomic number': 36,
        'mass number': 71,
        'mass': 70.95027 * u.u,
        'stable': False,
        'half-life': 0.1 * u.s,
    },

    'Kr-72': {
        'atomic number': 36,
        'mass number': 72,
        'mass': 71.9420924 * u.u,
        'stable': False,
        'half-life': 17.16 * u.s,
    },

    'Kr-73': {
        'atomic number': 36,
        'mass number': 73,
        'mass': 72.9392892 * u.u,
        'stable': False,
        'half-life': 27.3 * u.s,
    },

    'Kr-74': {
        'atomic number': 36,
        'mass number': 74,
        'mass': 73.933084 * u.u,
        'stable': False,
        'half-life': 690.0 * u.s,
    },

    'Kr-75': {
        'atomic number': 36,
        'mass number': 75,
        'mass': 74.9309457 * u.u,
        'stable': False,
        'half-life': 276.0 * u.s,
    },

    'Kr-76': {
        'atomic number': 36,
        'mass number': 76,
        'mass': 75.9259103 * u.u,
        'stable': False,
        'half-life': 53280.0 * u.s,
    },

    'Kr-77': {
        'atomic number': 36,
        'mass number': 77,
        'mass': 76.92467 * u.u,
        'stable': False,
        'half-life': 4464.0 * u.s,
    },

    'Kr-78': {
        'atomic number': 36,
        'mass number': 78,
        'mass': 77.92036494 * u.u,
        'stable': True,
        'abundance': 0.00355,
    },

    'Kr-79': {
        'atomic number': 36,
        'mass number': 79,
        'mass': 78.9200829 * u.u,
        'stable': False,
        'half-life': 126144.0 * u.s,
    },

    'Kr-80': {
        'atomic number': 36,
        'mass number': 80,
        'mass': 79.91637808 * u.u,
        'stable': True,
        'abundance': 0.02286,
    },

    'Kr-81': {
        'atomic number': 36,
        'mass number': 81,
        'mass': 80.9165912 * u.u,
        'stable': False,
        'half-life': 7226536054000.0 * u.s,
    },

    'Kr-82': {
        'atomic number': 36,
        'mass number': 82,
        'mass': 81.91348273 * u.u,
        'stable': True,
        'abundance': 0.11593,
    },

    'Kr-83': {
        'atomic number': 36,
        'mass number': 83,
        'mass': 82.91412716 * u.u,
        'stable': True,
        'abundance': 0.115,
    },

    'Kr-84': {
        'atomic number': 36,
        'mass number': 84,
        'mass': 83.9114977282 * u.u,
        'stable': True,
        'abundance': 0.56987,
    },

    'Kr-85': {
        'atomic number': 36,
        'mass number': 85,
        'mass': 84.9125273 * u.u,
        'stable': False,
        'half-life': 340044480.0 * u.s,
    },

    'Kr-86': {
        'atomic number': 36,
        'mass number': 86,
        'mass': 85.9106106269 * u.u,
        'stable': True,
        'abundance': 0.17279,
    },

    'Kr-87': {
        'atomic number': 36,
        'mass number': 87,
        'mass': 86.91335476 * u.u,
        'stable': False,
        'half-life': 4578.0 * u.s,
    },

    'Kr-88': {
        'atomic number': 36,
        'mass number': 88,
        'mass': 87.9144479 * u.u,
        'stable': False,
        'half-life': 10170.0 * u.s,
    },

    'Kr-89': {
        'atomic number': 36,
        'mass number': 89,
        'mass': 88.9178355 * u.u,
        'stable': False,
        'half-life': 189.0 * u.s,
    },

    'Kr-90': {
        'atomic number': 36,
        'mass number': 90,
        'mass': 89.9195279 * u.u,
        'stable': False,
        'half-life': 32.32 * u.s,
    },

    'Kr-91': {
        'atomic number': 36,
        'mass number': 91,
        'mass': 90.9238063 * u.u,
        'stable': False,
        'half-life': 8.57 * u.s,
    },

    'Kr-92': {
        'atomic number': 36,
        'mass number': 92,
        'mass': 91.9261731 * u.u,
        'stable': False,
        'half-life': 1.84 * u.s,
    },

    'Kr-93': {
        'atomic number': 36,
        'mass number': 93,
        'mass': 92.9311472 * u.u,
        'stable': False,
        'half-life': 1.286 * u.s,
    },

    'Kr-94': {
        'atomic number': 36,
        'mass number': 94,
        'mass': 93.93414 * u.u,
        'stable': False,
        'half-life': 0.212 * u.s,
    },

    'Kr-95': {
        'atomic number': 36,
        'mass number': 95,
        'mass': 94.939711 * u.u,
        'stable': False,
        'half-life': 0.114 * u.s,
    },

    'Kr-96': {
        'atomic number': 36,
        'mass number': 96,
        'mass': 95.943017 * u.u,
        'stable': False,
        'half-life': 0.08 * u.s,
    },

    'Kr-97': {
        'atomic number': 36,
        'mass number': 97,
        'mass': 96.94909 * u.u,
        'stable': False,
        'half-life': 0.0622 * u.s,
    },

    'Kr-98': {
        'atomic number': 36,
        'mass number': 98,
        'mass': 97.95243 * u.u,
        'stable': False,
        'half-life': 0.0428 * u.s,
    },

    'Kr-99': {
        'atomic number': 36,
        'mass number': 99,
        'mass': 98.95839 * u.u,
        'stable': False,
        'half-life': 0.04 * u.s,
    },

    'Kr-100': {
        'atomic number': 36,
        'mass number': 100,
        'mass': 99.96237 * u.u,
        'stable': False,
        'half-life': 0.012 * u.s,
    },

    'Kr-101': {
        'atomic number': 36,
        'mass number': 101,
        'mass': 100.96873 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Rb-71': {
        'atomic number': 37,
        'mass number': 71,
        'mass': 70.96532 * u.u,
        'stable': False,
    },

    'Rb-72': {
        'atomic number': 37,
        'mass number': 72,
        'mass': 71.95908 * u.u,
        'stable': False,
    },

    'Rb-73': {
        'atomic number': 37,
        'mass number': 73,
        'mass': 72.95053 * u.u,
        'stable': False,
    },

    'Rb-74': {
        'atomic number': 37,
        'mass number': 74,
        'mass': 73.9442659 * u.u,
        'stable': False,
        'half-life': 0.064776 * u.s,
    },

    'Rb-75': {
        'atomic number': 37,
        'mass number': 75,
        'mass': 74.9385732 * u.u,
        'stable': False,
        'half-life': 19.0 * u.s,
    },

    'Rb-76': {
        'atomic number': 37,
        'mass number': 76,
        'mass': 75.935073 * u.u,
        'stable': False,
        'half-life': 36.5 * u.s,
    },

    'Rb-77': {
        'atomic number': 37,
        'mass number': 77,
        'mass': 76.9304016 * u.u,
        'stable': False,
        'half-life': 226.8 * u.s,
    },

    'Rb-78': {
        'atomic number': 37,
        'mass number': 78,
        'mass': 77.9281419 * u.u,
        'stable': False,
        'half-life': 1059.6 * u.s,
    },

    'Rb-79': {
        'atomic number': 37,
        'mass number': 79,
        'mass': 78.9239899 * u.u,
        'stable': False,
        'half-life': 1374.0 * u.s,
    },

    'Rb-80': {
        'atomic number': 37,
        'mass number': 80,
        'mass': 79.9225164 * u.u,
        'stable': False,
        'half-life': 33.4 * u.s,
    },

    'Rb-81': {
        'atomic number': 37,
        'mass number': 81,
        'mass': 80.9189939 * u.u,
        'stable': False,
        'half-life': 16459.2 * u.s,
    },

    'Rb-82': {
        'atomic number': 37,
        'mass number': 82,
        'mass': 81.918209 * u.u,
        'stable': False,
        'half-life': 76.38 * u.s,
    },

    'Rb-83': {
        'atomic number': 37,
        'mass number': 83,
        'mass': 82.9151142 * u.u,
        'stable': False,
        'half-life': 7447680.0 * u.s,
    },

    'Rb-84': {
        'atomic number': 37,
        'mass number': 84,
        'mass': 83.9143752 * u.u,
        'stable': False,
        'half-life': 2835648.0 * u.s,
    },

    'Rb-85': {
        'atomic number': 37,
        'mass number': 85,
        'mass': 84.9117897379 * u.u,
        'stable': True,
        'abundance': 0.7217,
    },

    'Rb-86': {
        'atomic number': 37,
        'mass number': 86,
        'mass': 85.91116743 * u.u,
        'stable': False,
        'half-life': 1610668.8 * u.s,
    },

    'Rb-87': {
        'atomic number': 37,
        'mass number': 87,
        'mass': 86.909180531 * u.u,
        'stable': False,
        'abundance': 0.2783,
    },

    'Rb-88': {
        'atomic number': 37,
        'mass number': 88,
        'mass': 87.91131559 * u.u,
        'stable': False,
        'half-life': 1066.38 * u.s,
    },

    'Rb-89': {
        'atomic number': 37,
        'mass number': 89,
        'mass': 88.9122783 * u.u,
        'stable': False,
        'half-life': 919.2 * u.s,
    },

    'Rb-90': {
        'atomic number': 37,
        'mass number': 90,
        'mass': 89.9147985 * u.u,
        'stable': False,
        'half-life': 158.0 * u.s,
    },

    'Rb-91': {
        'atomic number': 37,
        'mass number': 91,
        'mass': 90.9165372 * u.u,
        'stable': False,
        'half-life': 58.2 * u.s,
    },

    'Rb-92': {
        'atomic number': 37,
        'mass number': 92,
        'mass': 91.9197284 * u.u,
        'stable': False,
        'half-life': 4.48 * u.s,
    },

    'Rb-93': {
        'atomic number': 37,
        'mass number': 93,
        'mass': 92.9220393 * u.u,
        'stable': False,
        'half-life': 5.84 * u.s,
    },

    'Rb-94': {
        'atomic number': 37,
        'mass number': 94,
        'mass': 93.9263948 * u.u,
        'stable': False,
        'half-life': 2.702 * u.s,
    },

    'Rb-95': {
        'atomic number': 37,
        'mass number': 95,
        'mass': 94.92926 * u.u,
        'stable': False,
        'half-life': 0.3777 * u.s,
    },

    'Rb-96': {
        'atomic number': 37,
        'mass number': 96,
        'mass': 95.9341334 * u.u,
        'stable': False,
        'half-life': 0.201 * u.s,
    },

    'Rb-97': {
        'atomic number': 37,
        'mass number': 97,
        'mass': 96.9371771 * u.u,
        'stable': False,
        'half-life': 0.1691 * u.s,
    },

    'Rb-98': {
        'atomic number': 37,
        'mass number': 98,
        'mass': 97.9416869 * u.u,
        'stable': False,
        'half-life': 0.114 * u.s,
    },

    'Rb-99': {
        'atomic number': 37,
        'mass number': 99,
        'mass': 98.94503 * u.u,
        'stable': False,
        'half-life': 0.0564 * u.s,
    },

    'Rb-100': {
        'atomic number': 37,
        'mass number': 100,
        'mass': 99.95003 * u.u,
        'stable': False,
        'half-life': 0.048 * u.s,
    },

    'Rb-101': {
        'atomic number': 37,
        'mass number': 101,
        'mass': 100.95404 * u.u,
        'stable': False,
        'half-life': 0.0318 * u.s,
    },

    'Rb-102': {
        'atomic number': 37,
        'mass number': 102,
        'mass': 101.95952 * u.u,
        'stable': False,
        'half-life': 0.037 * u.s,
    },

    'Rb-103': {
        'atomic number': 37,
        'mass number': 103,
        'mass': 102.96392 * u.u,
        'stable': False,
        'half-life': 0.026 * u.s,
    },

    'Sr-73': {
        'atomic number': 38,
        'mass number': 73,
        'mass': 72.9657 * u.u,
        'stable': False,
        'half-life': '>25 ms',
    },

    'Sr-74': {
        'atomic number': 38,
        'mass number': 74,
        'mass': 73.95617 * u.u,
        'stable': False,
        'half-life': 0.027 * u.s,
    },

    'Sr-75': {
        'atomic number': 38,
        'mass number': 75,
        'mass': 74.94995 * u.u,
        'stable': False,
        'half-life': 0.088 * u.s,
    },

    'Sr-76': {
        'atomic number': 38,
        'mass number': 76,
        'mass': 75.941763 * u.u,
        'stable': False,
        'half-life': 7.89 * u.s,
    },

    'Sr-77': {
        'atomic number': 38,
        'mass number': 77,
        'mass': 76.9379455 * u.u,
        'stable': False,
        'half-life': 9.0 * u.s,
    },

    'Sr-78': {
        'atomic number': 38,
        'mass number': 78,
        'mass': 77.93218 * u.u,
        'stable': False,
        'half-life': 156.1 * u.s,
    },

    'Sr-79': {
        'atomic number': 38,
        'mass number': 79,
        'mass': 78.9297077 * u.u,
        'stable': False,
        'half-life': 135.0 * u.s,
    },

    'Sr-80': {
        'atomic number': 38,
        'mass number': 80,
        'mass': 79.9245175 * u.u,
        'stable': False,
        'half-life': 6378.0 * u.s,
    },

    'Sr-81': {
        'atomic number': 38,
        'mass number': 81,
        'mass': 80.9232114 * u.u,
        'stable': False,
        'half-life': 1338.0 * u.s,
    },

    'Sr-82': {
        'atomic number': 38,
        'mass number': 82,
        'mass': 81.9183999 * u.u,
        'stable': False,
        'half-life': 2191104.0 * u.s,
    },

    'Sr-83': {
        'atomic number': 38,
        'mass number': 83,
        'mass': 82.9175544 * u.u,
        'stable': False,
        'half-life': 116676.0 * u.s,
    },

    'Sr-84': {
        'atomic number': 38,
        'mass number': 84,
        'mass': 83.9134191 * u.u,
        'stable': True,
        'abundance': 0.0056,
    },

    'Sr-85': {
        'atomic number': 38,
        'mass number': 85,
        'mass': 84.912932 * u.u,
        'stable': False,
        'half-life': 5603299.199999999 * u.s,
    },

    'Sr-86': {
        'atomic number': 38,
        'mass number': 86,
        'mass': 85.9092606 * u.u,
        'stable': True,
        'abundance': 0.0986,
    },

    'Sr-87': {
        'atomic number': 38,
        'mass number': 87,
        'mass': 86.9088775 * u.u,
        'stable': True,
        'abundance': 0.07,
    },

    'Sr-88': {
        'atomic number': 38,
        'mass number': 88,
        'mass': 87.9056125 * u.u,
        'stable': True,
        'abundance': 0.8258,
    },

    'Sr-89': {
        'atomic number': 38,
        'mass number': 89,
        'mass': 88.9074511 * u.u,
        'stable': False,
        'half-life': 4368643.2 * u.s,
    },

    'Sr-90': {
        'atomic number': 38,
        'mass number': 90,
        'mass': 89.90773 * u.u,
        'stable': False,
        'half-life': 908523899.54 * u.s,
    },

    'Sr-91': {
        'atomic number': 38,
        'mass number': 91,
        'mass': 90.9101954 * u.u,
        'stable': False,
        'half-life': 34740.0 * u.s,
    },

    'Sr-92': {
        'atomic number': 38,
        'mass number': 92,
        'mass': 91.9110382 * u.u,
        'stable': False,
        'half-life': 9399.6 * u.s,
    },

    'Sr-93': {
        'atomic number': 38,
        'mass number': 93,
        'mass': 92.9140242 * u.u,
        'stable': False,
        'half-life': 445.8 * u.s,
    },

    'Sr-94': {
        'atomic number': 38,
        'mass number': 94,
        'mass': 93.9153556 * u.u,
        'stable': False,
        'half-life': 75.3 * u.s,
    },

    'Sr-95': {
        'atomic number': 38,
        'mass number': 95,
        'mass': 94.9193529 * u.u,
        'stable': False,
        'half-life': 23.9 * u.s,
    },

    'Sr-96': {
        'atomic number': 38,
        'mass number': 96,
        'mass': 95.9217066 * u.u,
        'stable': False,
        'half-life': 1.07 * u.s,
    },

    'Sr-97': {
        'atomic number': 38,
        'mass number': 97,
        'mass': 96.926374 * u.u,
        'stable': False,
        'half-life': 0.429 * u.s,
    },

    'Sr-98': {
        'atomic number': 38,
        'mass number': 98,
        'mass': 97.9286888 * u.u,
        'stable': False,
        'half-life': 0.653 * u.s,
    },

    'Sr-99': {
        'atomic number': 38,
        'mass number': 99,
        'mass': 98.9328907 * u.u,
        'stable': False,
        'half-life': 0.269 * u.s,
    },

    'Sr-100': {
        'atomic number': 38,
        'mass number': 100,
        'mass': 99.93577 * u.u,
        'stable': False,
        'half-life': 0.202 * u.s,
    },

    'Sr-101': {
        'atomic number': 38,
        'mass number': 101,
        'mass': 100.940352 * u.u,
        'stable': False,
        'half-life': 0.1138 * u.s,
    },

    'Sr-102': {
        'atomic number': 38,
        'mass number': 102,
        'mass': 101.943791 * u.u,
        'stable': False,
        'half-life': 0.069 * u.s,
    },

    'Sr-103': {
        'atomic number': 38,
        'mass number': 103,
        'mass': 102.94909 * u.u,
        'stable': False,
        'half-life': 0.053 * u.s,
    },

    'Sr-104': {
        'atomic number': 38,
        'mass number': 104,
        'mass': 103.95265 * u.u,
        'stable': False,
        'half-life': 0.0506 * u.s,
    },

    'Sr-105': {
        'atomic number': 38,
        'mass number': 105,
        'mass': 104.95855 * u.u,
        'stable': False,
        'half-life': 0.039 * u.s,
    },

    'Sr-106': {
        'atomic number': 38,
        'mass number': 106,
        'mass': 105.96265 * u.u,
        'stable': False,
        'half-life': 0.021 * u.s,
    },

    'Sr-107': {
        'atomic number': 38,
        'mass number': 107,
        'mass': 106.96897 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Y-76': {
        'atomic number': 39,
        'mass number': 76,
        'mass': 75.95856 * u.u,
        'stable': False,
        'half-life': '120# us',
    },

    'Y-77': {
        'atomic number': 39,
        'mass number': 77,
        'mass': 76.949781 * u.u,
        'stable': False,
        'half-life': 0.063 * u.s,
    },

    'Y-78': {
        'atomic number': 39,
        'mass number': 78,
        'mass': 77.94361 * u.u,
        'stable': False,
        'half-life': 0.054 * u.s,
    },

    'Y-79': {
        'atomic number': 39,
        'mass number': 79,
        'mass': 78.93735 * u.u,
        'stable': False,
        'half-life': 14.8 * u.s,
    },

    'Y-80': {
        'atomic number': 39,
        'mass number': 80,
        'mass': 79.9343561 * u.u,
        'stable': False,
        'half-life': 30.1 * u.s,
    },

    'Y-81': {
        'atomic number': 39,
        'mass number': 81,
        'mass': 80.9294556 * u.u,
        'stable': False,
        'half-life': 70.4 * u.s,
    },

    'Y-82': {
        'atomic number': 39,
        'mass number': 82,
        'mass': 81.9269314 * u.u,
        'stable': False,
        'half-life': 8.3 * u.s,
    },

    'Y-83': {
        'atomic number': 39,
        'mass number': 83,
        'mass': 82.922485 * u.u,
        'stable': False,
        'half-life': 424.8 * u.s,
    },

    'Y-84': {
        'atomic number': 39,
        'mass number': 84,
        'mass': 83.9206721 * u.u,
        'stable': False,
        'half-life': 2370.0 * u.s,
    },

    'Y-85': {
        'atomic number': 39,
        'mass number': 85,
        'mass': 84.916433 * u.u,
        'stable': False,
        'half-life': 9648.0 * u.s,
    },

    'Y-86': {
        'atomic number': 39,
        'mass number': 86,
        'mass': 85.914886 * u.u,
        'stable': False,
        'half-life': 53064.0 * u.s,
    },

    'Y-87': {
        'atomic number': 39,
        'mass number': 87,
        'mass': 86.9108761 * u.u,
        'stable': False,
        'half-life': 287280.0 * u.s,
    },

    'Y-88': {
        'atomic number': 39,
        'mass number': 88,
        'mass': 87.9095016 * u.u,
        'stable': False,
        'half-life': 9212486.4 * u.s,
    },

    'Y-89': {
        'atomic number': 39,
        'mass number': 89,
        'mass': 88.9058403 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Y-90': {
        'atomic number': 39,
        'mass number': 90,
        'mass': 89.9071439 * u.u,
        'stable': False,
        'half-life': 230400.0 * u.s,
    },

    'Y-91': {
        'atomic number': 39,
        'mass number': 91,
        'mass': 90.9072974 * u.u,
        'stable': False,
        'half-life': 5055264.0 * u.s,
    },

    'Y-92': {
        'atomic number': 39,
        'mass number': 92,
        'mass': 91.9089451 * u.u,
        'stable': False,
        'half-life': 12744.0 * u.s,
    },

    'Y-93': {
        'atomic number': 39,
        'mass number': 93,
        'mass': 92.909578 * u.u,
        'stable': False,
        'half-life': 36648.0 * u.s,
    },

    'Y-94': {
        'atomic number': 39,
        'mass number': 94,
        'mass': 93.9115906 * u.u,
        'stable': False,
        'half-life': 1122.0 * u.s,
    },

    'Y-95': {
        'atomic number': 39,
        'mass number': 95,
        'mass': 94.9128161 * u.u,
        'stable': False,
        'half-life': 618.0 * u.s,
    },

    'Y-96': {
        'atomic number': 39,
        'mass number': 96,
        'mass': 95.9158968 * u.u,
        'stable': False,
        'half-life': 5.34 * u.s,
    },

    'Y-97': {
        'atomic number': 39,
        'mass number': 97,
        'mass': 96.9182741 * u.u,
        'stable': False,
        'half-life': 3.75 * u.s,
    },

    'Y-98': {
        'atomic number': 39,
        'mass number': 98,
        'mass': 97.9223821 * u.u,
        'stable': False,
        'half-life': 0.548 * u.s,
    },

    'Y-99': {
        'atomic number': 39,
        'mass number': 99,
        'mass': 98.924148 * u.u,
        'stable': False,
        'half-life': 1.484 * u.s,
    },

    'Y-100': {
        'atomic number': 39,
        'mass number': 100,
        'mass': 99.927715 * u.u,
        'stable': False,
        'half-life': 0.735 * u.s,
    },

    'Y-101': {
        'atomic number': 39,
        'mass number': 101,
        'mass': 100.9301477 * u.u,
        'stable': False,
        'half-life': 0.426 * u.s,
    },

    'Y-102': {
        'atomic number': 39,
        'mass number': 102,
        'mass': 101.9343277 * u.u,
        'stable': False,
        'half-life': 0.298 * u.s,
    },

    'Y-103': {
        'atomic number': 39,
        'mass number': 103,
        'mass': 102.937243 * u.u,
        'stable': False,
        'half-life': 0.239 * u.s,
    },

    'Y-104': {
        'atomic number': 39,
        'mass number': 104,
        'mass': 103.94196 * u.u,
        'stable': False,
        'half-life': 0.197 * u.s,
    },

    'Y-105': {
        'atomic number': 39,
        'mass number': 105,
        'mass': 104.94544 * u.u,
        'stable': False,
        'half-life': 0.095 * u.s,
    },

    'Y-106': {
        'atomic number': 39,
        'mass number': 106,
        'mass': 105.95056 * u.u,
        'stable': False,
        'half-life': 0.074 * u.s,
    },

    'Y-107': {
        'atomic number': 39,
        'mass number': 107,
        'mass': 106.95452 * u.u,
        'stable': False,
        'half-life': 0.0335 * u.s,
    },

    'Y-108': {
        'atomic number': 39,
        'mass number': 108,
        'mass': 107.95996 * u.u,
        'stable': False,
        'half-life': 0.03 * u.s,
    },

    'Y-109': {
        'atomic number': 39,
        'mass number': 109,
        'mass': 108.96436 * u.u,
        'stable': False,
        'half-life': 0.025 * u.s,
    },

    'Zr-78': {
        'atomic number': 40,
        'mass number': 78,
        'mass': 77.95566 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'Zr-79': {
        'atomic number': 40,
        'mass number': 79,
        'mass': 78.94948 * u.u,
        'stable': False,
        'half-life': 0.056 * u.s,
    },

    'Zr-80': {
        'atomic number': 40,
        'mass number': 80,
        'mass': 79.9404 * u.u,
        'stable': False,
        'half-life': 4.6 * u.s,
    },

    'Zr-81': {
        'atomic number': 40,
        'mass number': 81,
        'mass': 80.93731 * u.u,
        'stable': False,
        'half-life': 5.5 * u.s,
    },

    'Zr-82': {
        'atomic number': 40,
        'mass number': 82,
        'mass': 81.93135 * u.u,
        'stable': False,
        'half-life': 32.0 * u.s,
    },

    'Zr-83': {
        'atomic number': 40,
        'mass number': 83,
        'mass': 82.9292421 * u.u,
        'stable': False,
        'half-life': 42.0 * u.s,
    },

    'Zr-84': {
        'atomic number': 40,
        'mass number': 84,
        'mass': 83.9233269 * u.u,
        'stable': False,
        'half-life': 1548.0 * u.s,
    },

    'Zr-85': {
        'atomic number': 40,
        'mass number': 85,
        'mass': 84.9214444 * u.u,
        'stable': False,
        'half-life': 471.6 * u.s,
    },

    'Zr-86': {
        'atomic number': 40,
        'mass number': 86,
        'mass': 85.9162972 * u.u,
        'stable': False,
        'half-life': 59400.0 * u.s,
    },

    'Zr-87': {
        'atomic number': 40,
        'mass number': 87,
        'mass': 86.914818 * u.u,
        'stable': False,
        'half-life': 6048.0 * u.s,
    },

    'Zr-88': {
        'atomic number': 40,
        'mass number': 88,
        'mass': 87.9102213 * u.u,
        'stable': False,
        'half-life': 7205760.0 * u.s,
    },

    'Zr-89': {
        'atomic number': 40,
        'mass number': 89,
        'mass': 88.9088814 * u.u,
        'stable': False,
        'half-life': 282276.0 * u.s,
    },

    'Zr-90': {
        'atomic number': 40,
        'mass number': 90,
        'mass': 89.9046977 * u.u,
        'stable': True,
        'abundance': 0.5145,
    },

    'Zr-91': {
        'atomic number': 40,
        'mass number': 91,
        'mass': 90.9056396 * u.u,
        'stable': True,
        'abundance': 0.1122,
    },

    'Zr-92': {
        'atomic number': 40,
        'mass number': 92,
        'mass': 91.9050347 * u.u,
        'stable': True,
        'abundance': 0.1715,
    },

    'Zr-93': {
        'atomic number': 40,
        'mass number': 93,
        'mass': 92.9064699 * u.u,
        'stable': False,
        'half-life': 50806650860000.0 * u.s,
    },

    'Zr-94': {
        'atomic number': 40,
        'mass number': 94,
        'mass': 93.9063108 * u.u,
        'stable': True,
        'abundance': 0.1738,
    },

    'Zr-95': {
        'atomic number': 40,
        'mass number': 95,
        'mass': 94.9080385 * u.u,
        'stable': False,
        'half-life': 5532364.8 * u.s,
    },

    'Zr-96': {
        'atomic number': 40,
        'mass number': 96,
        'mass': 95.9082714 * u.u,
        'stable': False,
        'abundance': 0.028,
    },

    'Zr-97': {
        'atomic number': 40,
        'mass number': 97,
        'mass': 96.9109512 * u.u,
        'stable': False,
        'half-life': 60296.4 * u.s,
    },

    'Zr-98': {
        'atomic number': 40,
        'mass number': 98,
        'mass': 97.9127289 * u.u,
        'stable': False,
        'half-life': 30.7 * u.s,
    },

    'Zr-99': {
        'atomic number': 40,
        'mass number': 99,
        'mass': 98.916667 * u.u,
        'stable': False,
        'half-life': 2.1 * u.s,
    },

    'Zr-100': {
        'atomic number': 40,
        'mass number': 100,
        'mass': 99.9180006 * u.u,
        'stable': False,
        'half-life': 7.1 * u.s,
    },

    'Zr-101': {
        'atomic number': 40,
        'mass number': 101,
        'mass': 100.921448 * u.u,
        'stable': False,
        'half-life': 2.3 * u.s,
    },

    'Zr-102': {
        'atomic number': 40,
        'mass number': 102,
        'mass': 101.9231409 * u.u,
        'stable': False,
        'half-life': 2.9 * u.s,
    },

    'Zr-103': {
        'atomic number': 40,
        'mass number': 103,
        'mass': 102.927191 * u.u,
        'stable': False,
        'half-life': 1.38 * u.s,
    },

    'Zr-104': {
        'atomic number': 40,
        'mass number': 104,
        'mass': 103.929436 * u.u,
        'stable': False,
        'half-life': 0.92 * u.s,
    },

    'Zr-105': {
        'atomic number': 40,
        'mass number': 105,
        'mass': 104.934008 * u.u,
        'stable': False,
        'half-life': 0.67 * u.s,
    },

    'Zr-106': {
        'atomic number': 40,
        'mass number': 106,
        'mass': 105.93676 * u.u,
        'stable': False,
        'half-life': 0.1786 * u.s,
    },

    'Zr-107': {
        'atomic number': 40,
        'mass number': 107,
        'mass': 106.94174 * u.u,
        'stable': False,
        'half-life': 0.1457 * u.s,
    },

    'Zr-108': {
        'atomic number': 40,
        'mass number': 108,
        'mass': 107.94487 * u.u,
        'stable': False,
        'half-life': 0.0785 * u.s,
    },

    'Zr-109': {
        'atomic number': 40,
        'mass number': 109,
        'mass': 108.95041 * u.u,
        'stable': False,
        'half-life': 0.056 * u.s,
    },

    'Zr-110': {
        'atomic number': 40,
        'mass number': 110,
        'mass': 109.95396 * u.u,
        'stable': False,
        'half-life': 0.0375 * u.s,
    },

    'Zr-111': {
        'atomic number': 40,
        'mass number': 111,
        'mass': 110.95968 * u.u,
        'stable': False,
        'half-life': 0.024 * u.s,
    },

    'Zr-112': {
        'atomic number': 40,
        'mass number': 112,
        'mass': 111.9637 * u.u,
        'stable': False,
        'half-life': 0.043 * u.s,
    },

    'Nb-81': {
        'atomic number': 41,
        'mass number': 81,
        'mass': 80.9496 * u.u,
        'stable': False,
        'half-life': '<44 ns',
    },

    'Nb-82': {
        'atomic number': 41,
        'mass number': 82,
        'mass': 81.94396 * u.u,
        'stable': False,
        'half-life': 0.05 * u.s,
    },

    'Nb-83': {
        'atomic number': 41,
        'mass number': 83,
        'mass': 82.93729 * u.u,
        'stable': False,
        'half-life': 3.9 * u.s,
    },

    'Nb-84': {
        'atomic number': 41,
        'mass number': 84,
        'mass': 83.93449 * u.u,
        'stable': False,
        'half-life': 9.8 * u.s,
    },

    'Nb-85': {
        'atomic number': 41,
        'mass number': 85,
        'mass': 84.9288458 * u.u,
        'stable': False,
        'half-life': 20.5 * u.s,
    },

    'Nb-86': {
        'atomic number': 41,
        'mass number': 86,
        'mass': 85.9257828 * u.u,
        'stable': False,
        'half-life': 88.0 * u.s,
    },

    'Nb-87': {
        'atomic number': 41,
        'mass number': 87,
        'mass': 86.9206937 * u.u,
        'stable': False,
        'half-life': 222.0 * u.s,
    },

    'Nb-88': {
        'atomic number': 41,
        'mass number': 88,
        'mass': 87.918222 * u.u,
        'stable': False,
        'half-life': 870.0 * u.s,
    },

    'Nb-89': {
        'atomic number': 41,
        'mass number': 89,
        'mass': 88.913445 * u.u,
        'stable': False,
        'half-life': 7308.0 * u.s,
    },

    'Nb-90': {
        'atomic number': 41,
        'mass number': 90,
        'mass': 89.9112584 * u.u,
        'stable': False,
        'half-life': 52560.0 * u.s,
    },

    'Nb-91': {
        'atomic number': 41,
        'mass number': 91,
        'mass': 90.9069897 * u.u,
        'stable': False,
        'half-life': 21458709680.0 * u.s,
    },

    'Nb-92': {
        'atomic number': 41,
        'mass number': 92,
        'mass': 91.9071881 * u.u,
        'stable': False,
        'half-life': 1095025332200000.0 * u.s,
    },

    'Nb-93': {
        'atomic number': 41,
        'mass number': 93,
        'mass': 92.906373 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Nb-94': {
        'atomic number': 41,
        'mass number': 94,
        'mass': 93.9072788 * u.u,
        'stable': False,
        'half-life': 643761290400.0 * u.s,
    },

    'Nb-95': {
        'atomic number': 41,
        'mass number': 95,
        'mass': 94.9068324 * u.u,
        'stable': False,
        'half-life': 3023222.4 * u.s,
    },

    'Nb-96': {
        'atomic number': 41,
        'mass number': 96,
        'mass': 95.9080973 * u.u,
        'stable': False,
        'half-life': 84060.0 * u.s,
    },

    'Nb-97': {
        'atomic number': 41,
        'mass number': 97,
        'mass': 96.9080959 * u.u,
        'stable': False,
        'half-life': 4326.0 * u.s,
    },

    'Nb-98': {
        'atomic number': 41,
        'mass number': 98,
        'mass': 97.9103265 * u.u,
        'stable': False,
        'half-life': 2.86 * u.s,
    },

    'Nb-99': {
        'atomic number': 41,
        'mass number': 99,
        'mass': 98.911613 * u.u,
        'stable': False,
        'half-life': 15.0 * u.s,
    },

    'Nb-100': {
        'atomic number': 41,
        'mass number': 100,
        'mass': 99.9143276 * u.u,
        'stable': False,
        'half-life': 1.5 * u.s,
    },

    'Nb-101': {
        'atomic number': 41,
        'mass number': 101,
        'mass': 100.9153103 * u.u,
        'stable': False,
        'half-life': 7.1 * u.s,
    },

    'Nb-102': {
        'atomic number': 41,
        'mass number': 102,
        'mass': 101.9180772 * u.u,
        'stable': False,
        'half-life': 4.3 * u.s,
    },

    'Nb-103': {
        'atomic number': 41,
        'mass number': 103,
        'mass': 102.9194572 * u.u,
        'stable': False,
        'half-life': 1.5 * u.s,
    },

    'Nb-104': {
        'atomic number': 41,
        'mass number': 104,
        'mass': 103.9228925 * u.u,
        'stable': False,
        'half-life': 4.9 * u.s,
    },

    'Nb-105': {
        'atomic number': 41,
        'mass number': 105,
        'mass': 104.9249465 * u.u,
        'stable': False,
        'half-life': 2.95 * u.s,
    },

    'Nb-106': {
        'atomic number': 41,
        'mass number': 106,
        'mass': 105.9289317 * u.u,
        'stable': False,
        'half-life': 1.05 * u.s,
    },

    'Nb-107': {
        'atomic number': 41,
        'mass number': 107,
        'mass': 106.9315937 * u.u,
        'stable': False,
        'half-life': 0.289 * u.s,
    },

    'Nb-108': {
        'atomic number': 41,
        'mass number': 108,
        'mass': 107.9360748 * u.u,
        'stable': False,
        'half-life': 0.198 * u.s,
    },

    'Nb-109': {
        'atomic number': 41,
        'mass number': 109,
        'mass': 108.93922 * u.u,
        'stable': False,
        'half-life': 0.1069 * u.s,
    },

    'Nb-110': {
        'atomic number': 41,
        'mass number': 110,
        'mass': 109.94403 * u.u,
        'stable': False,
        'half-life': 0.082 * u.s,
    },

    'Nb-111': {
        'atomic number': 41,
        'mass number': 111,
        'mass': 110.94753 * u.u,
        'stable': False,
        'half-life': 0.054 * u.s,
    },

    'Nb-112': {
        'atomic number': 41,
        'mass number': 112,
        'mass': 111.95247 * u.u,
        'stable': False,
        'half-life': 0.038 * u.s,
    },

    'Nb-113': {
        'atomic number': 41,
        'mass number': 113,
        'mass': 112.95651 * u.u,
        'stable': False,
        'half-life': 0.032 * u.s,
    },

    'Nb-114': {
        'atomic number': 41,
        'mass number': 114,
        'mass': 113.96201 * u.u,
        'stable': False,
        'half-life': 0.017 * u.s,
    },

    'Nb-115': {
        'atomic number': 41,
        'mass number': 115,
        'mass': 114.96634 * u.u,
        'stable': False,
        'half-life': 0.023 * u.s,
    },

    'Mo-83': {
        'atomic number': 42,
        'mass number': 83,
        'mass': 82.94988 * u.u,
        'stable': False,
        'half-life': 0.023 * u.s,
    },

    'Mo-84': {
        'atomic number': 42,
        'mass number': 84,
        'mass': 83.94149 * u.u,
        'stable': False,
        'half-life': 2.3 * u.s,
    },

    'Mo-85': {
        'atomic number': 42,
        'mass number': 85,
        'mass': 84.938261 * u.u,
        'stable': False,
        'half-life': 3.2 * u.s,
    },

    'Mo-86': {
        'atomic number': 42,
        'mass number': 86,
        'mass': 85.9311748 * u.u,
        'stable': False,
        'half-life': 19.1 * u.s,
    },

    'Mo-87': {
        'atomic number': 42,
        'mass number': 87,
        'mass': 86.9281962 * u.u,
        'stable': False,
        'half-life': 14.1 * u.s,
    },

    'Mo-88': {
        'atomic number': 42,
        'mass number': 88,
        'mass': 87.9219678 * u.u,
        'stable': False,
        'half-life': 480.0 * u.s,
    },

    'Mo-89': {
        'atomic number': 42,
        'mass number': 89,
        'mass': 88.9194682 * u.u,
        'stable': False,
        'half-life': 126.6 * u.s,
    },

    'Mo-90': {
        'atomic number': 42,
        'mass number': 90,
        'mass': 89.9139309 * u.u,
        'stable': False,
        'half-life': 20016.0 * u.s,
    },

    'Mo-91': {
        'atomic number': 42,
        'mass number': 91,
        'mass': 90.9117453 * u.u,
        'stable': False,
        'half-life': 929.4 * u.s,
    },

    'Mo-92': {
        'atomic number': 42,
        'mass number': 92,
        'mass': 91.90680796 * u.u,
        'stable': True,
        'abundance': 0.1453,
    },

    'Mo-93': {
        'atomic number': 42,
        'mass number': 93,
        'mass': 92.90680958 * u.u,
        'stable': False,
        'half-life': 126227704000.0 * u.s,
    },

    'Mo-94': {
        'atomic number': 42,
        'mass number': 94,
        'mass': 93.9050849 * u.u,
        'stable': True,
        'abundance': 0.0915,
    },

    'Mo-95': {
        'atomic number': 42,
        'mass number': 95,
        'mass': 94.90583877 * u.u,
        'stable': True,
        'abundance': 0.1584,
    },

    'Mo-96': {
        'atomic number': 42,
        'mass number': 96,
        'mass': 95.90467612 * u.u,
        'stable': True,
        'abundance': 0.1667,
    },

    'Mo-97': {
        'atomic number': 42,
        'mass number': 97,
        'mass': 96.90601812 * u.u,
        'stable': True,
        'abundance': 0.096,
    },

    'Mo-98': {
        'atomic number': 42,
        'mass number': 98,
        'mass': 97.90540482 * u.u,
        'stable': True,
        'abundance': 0.2439,
    },

    'Mo-99': {
        'atomic number': 42,
        'mass number': 99,
        'mass': 98.90770851 * u.u,
        'stable': False,
        'half-life': 237326.04 * u.s,
    },

    'Mo-100': {
        'atomic number': 42,
        'mass number': 100,
        'mass': 99.9074718 * u.u,
        'stable': False,
        'abundance': 0.0982,
    },

    'Mo-101': {
        'atomic number': 42,
        'mass number': 101,
        'mass': 100.9103414 * u.u,
        'stable': False,
        'half-life': 876.6 * u.s,
    },

    'Mo-102': {
        'atomic number': 42,
        'mass number': 102,
        'mass': 101.9102834 * u.u,
        'stable': False,
        'half-life': 678.0 * u.s,
    },

    'Mo-103': {
        'atomic number': 42,
        'mass number': 103,
        'mass': 102.913079 * u.u,
        'stable': False,
        'half-life': 67.5 * u.s,
    },

    'Mo-104': {
        'atomic number': 42,
        'mass number': 104,
        'mass': 103.9137344 * u.u,
        'stable': False,
        'half-life': 60.0 * u.s,
    },

    'Mo-105': {
        'atomic number': 42,
        'mass number': 105,
        'mass': 104.916969 * u.u,
        'stable': False,
        'half-life': 35.6 * u.s,
    },

    'Mo-106': {
        'atomic number': 42,
        'mass number': 106,
        'mass': 105.918259 * u.u,
        'stable': False,
        'half-life': 8.73 * u.s,
    },

    'Mo-107': {
        'atomic number': 42,
        'mass number': 107,
        'mass': 106.922106 * u.u,
        'stable': False,
        'half-life': 3.5 * u.s,
    },

    'Mo-108': {
        'atomic number': 42,
        'mass number': 108,
        'mass': 107.924033 * u.u,
        'stable': False,
        'half-life': 1.105 * u.s,
    },

    'Mo-109': {
        'atomic number': 42,
        'mass number': 109,
        'mass': 108.928424 * u.u,
        'stable': False,
        'half-life': 0.7 * u.s,
    },

    'Mo-110': {
        'atomic number': 42,
        'mass number': 110,
        'mass': 109.930704 * u.u,
        'stable': False,
        'half-life': 0.292 * u.s,
    },

    'Mo-111': {
        'atomic number': 42,
        'mass number': 111,
        'mass': 110.935654 * u.u,
        'stable': False,
        'half-life': 0.1936 * u.s,
    },

    'Mo-112': {
        'atomic number': 42,
        'mass number': 112,
        'mass': 111.93831 * u.u,
        'stable': False,
        'half-life': 0.125 * u.s,
    },

    'Mo-113': {
        'atomic number': 42,
        'mass number': 113,
        'mass': 112.94335 * u.u,
        'stable': False,
        'half-life': 0.08 * u.s,
    },

    'Mo-114': {
        'atomic number': 42,
        'mass number': 114,
        'mass': 113.94653 * u.u,
        'stable': False,
        'half-life': 0.058 * u.s,
    },

    'Mo-115': {
        'atomic number': 42,
        'mass number': 115,
        'mass': 114.95196 * u.u,
        'stable': False,
        'half-life': 0.0455 * u.s,
    },

    'Mo-116': {
        'atomic number': 42,
        'mass number': 116,
        'mass': 115.95545 * u.u,
        'stable': False,
        'half-life': 0.032 * u.s,
    },

    'Mo-117': {
        'atomic number': 42,
        'mass number': 117,
        'mass': 116.96117 * u.u,
        'stable': False,
        'half-life': 0.022 * u.s,
    },

    'Tc-85': {
        'atomic number': 43,
        'mass number': 85,
        'mass': 84.95058 * u.u,
        'stable': False,
    },

    'Tc-86': {
        'atomic number': 43,
        'mass number': 86,
        'mass': 85.94493 * u.u,
        'stable': False,
        'half-life': 0.055 * u.s,
    },

    'Tc-87': {
        'atomic number': 43,
        'mass number': 87,
        'mass': 86.9380672 * u.u,
        'stable': False,
        'half-life': 2.2 * u.s,
    },

    'Tc-88': {
        'atomic number': 43,
        'mass number': 88,
        'mass': 87.93378 * u.u,
        'stable': False,
        'half-life': 6.4 * u.s,
    },

    'Tc-89': {
        'atomic number': 43,
        'mass number': 89,
        'mass': 88.9276487 * u.u,
        'stable': False,
        'half-life': 12.8 * u.s,
    },

    'Tc-90': {
        'atomic number': 43,
        'mass number': 90,
        'mass': 89.9240739 * u.u,
        'stable': False,
        'half-life': 49.2 * u.s,
    },

    'Tc-91': {
        'atomic number': 43,
        'mass number': 91,
        'mass': 90.9184254 * u.u,
        'stable': False,
        'half-life': 188.4 * u.s,
    },

    'Tc-92': {
        'atomic number': 43,
        'mass number': 92,
        'mass': 91.9152698 * u.u,
        'stable': False,
        'half-life': 255.0 * u.s,
    },

    'Tc-93': {
        'atomic number': 43,
        'mass number': 93,
        'mass': 92.910246 * u.u,
        'stable': False,
        'half-life': 9900.0 * u.s,
    },

    'Tc-94': {
        'atomic number': 43,
        'mass number': 94,
        'mass': 93.9096536 * u.u,
        'stable': False,
        'half-life': 17580.0 * u.s,
    },

    'Tc-95': {
        'atomic number': 43,
        'mass number': 95,
        'mass': 94.9076536 * u.u,
        'stable': False,
        'half-life': 72000.0 * u.s,
    },

    'Tc-96': {
        'atomic number': 43,
        'mass number': 96,
        'mass': 95.907868 * u.u,
        'stable': False,
        'half-life': 369792.0 * u.s,
    },

    'Tc-97': {
        'atomic number': 43,
        'mass number': 97,
        'mass': 96.9063667 * u.u,
        'stable': False,
        'half-life': 132854658460000.0 * u.s,
    },

    'Tc-98': {
        'atomic number': 43,
        'mass number': 98,
        'mass': 97.9072124 * u.u,
        'stable': False,
        'half-life': 132539089200000.0 * u.s,
    },

    'Tc-99': {
        'atomic number': 43,
        'mass number': 99,
        'mass': 98.9062508 * u.u,
        'stable': False,
        'half-life': 21636.0 * u.s,
    },

    'Tc-100': {
        'atomic number': 43,
        'mass number': 100,
        'mass': 99.9076539 * u.u,
        'stable': False,
        'half-life': 15.46 * u.s,
    },

    'Tc-101': {
        'atomic number': 43,
        'mass number': 101,
        'mass': 100.907309 * u.u,
        'stable': False,
        'half-life': 853.2 * u.s,
    },

    'Tc-102': {
        'atomic number': 43,
        'mass number': 102,
        'mass': 101.9092097 * u.u,
        'stable': False,
        'half-life': 5.28 * u.s,
    },

    'Tc-103': {
        'atomic number': 43,
        'mass number': 103,
        'mass': 102.909176 * u.u,
        'stable': False,
        'half-life': 54.2 * u.s,
    },

    'Tc-104': {
        'atomic number': 43,
        'mass number': 104,
        'mass': 103.911425 * u.u,
        'stable': False,
        'half-life': 1098.0 * u.s,
    },

    'Tc-105': {
        'atomic number': 43,
        'mass number': 105,
        'mass': 104.911655 * u.u,
        'stable': False,
        'half-life': 456.0 * u.s,
    },

    'Tc-106': {
        'atomic number': 43,
        'mass number': 106,
        'mass': 105.914358 * u.u,
        'stable': False,
        'half-life': 35.6 * u.s,
    },

    'Tc-107': {
        'atomic number': 43,
        'mass number': 107,
        'mass': 106.9154606 * u.u,
        'stable': False,
        'half-life': 21.2 * u.s,
    },

    'Tc-108': {
        'atomic number': 43,
        'mass number': 108,
        'mass': 107.9184957 * u.u,
        'stable': False,
        'half-life': 5.17 * u.s,
    },

    'Tc-109': {
        'atomic number': 43,
        'mass number': 109,
        'mass': 108.920256 * u.u,
        'stable': False,
        'half-life': 1.14 * u.s,
    },

    'Tc-110': {
        'atomic number': 43,
        'mass number': 110,
        'mass': 109.923744 * u.u,
        'stable': False,
        'half-life': 0.9 * u.s,
    },

    'Tc-111': {
        'atomic number': 43,
        'mass number': 111,
        'mass': 110.925901 * u.u,
        'stable': False,
        'half-life': 0.35 * u.s,
    },

    'Tc-112': {
        'atomic number': 43,
        'mass number': 112,
        'mass': 111.9299458 * u.u,
        'stable': False,
        'half-life': 0.323 * u.s,
    },

    'Tc-113': {
        'atomic number': 43,
        'mass number': 113,
        'mass': 112.932569 * u.u,
        'stable': False,
        'half-life': 0.152 * u.s,
    },

    'Tc-114': {
        'atomic number': 43,
        'mass number': 114,
        'mass': 113.93691 * u.u,
        'stable': False,
        'half-life': 0.09 * u.s,
    },

    'Tc-115': {
        'atomic number': 43,
        'mass number': 115,
        'mass': 114.93998 * u.u,
        'stable': False,
        'half-life': 0.078 * u.s,
    },

    'Tc-116': {
        'atomic number': 43,
        'mass number': 116,
        'mass': 115.94476 * u.u,
        'stable': False,
        'half-life': 0.057 * u.s,
    },

    'Tc-117': {
        'atomic number': 43,
        'mass number': 117,
        'mass': 116.94806 * u.u,
        'stable': False,
        'half-life': 0.0445 * u.s,
    },

    'Tc-118': {
        'atomic number': 43,
        'mass number': 118,
        'mass': 117.95299 * u.u,
        'stable': False,
        'half-life': 0.03 * u.s,
    },

    'Tc-119': {
        'atomic number': 43,
        'mass number': 119,
        'mass': 118.95666 * u.u,
        'stable': False,
        'half-life': 0.022 * u.s,
    },

    'Tc-120': {
        'atomic number': 43,
        'mass number': 120,
        'mass': 119.96187 * u.u,
        'stable': False,
        'half-life': 0.021 * u.s,
    },

    'Ru-87': {
        'atomic number': 44,
        'mass number': 87,
        'mass': 86.95069 * u.u,
        'stable': False,
        'half-life': '50# ms',
    },

    'Ru-88': {
        'atomic number': 44,
        'mass number': 88,
        'mass': 87.9416 * u.u,
        'stable': False,
        'half-life': 1.3 * u.s,
    },

    'Ru-89': {
        'atomic number': 44,
        'mass number': 89,
        'mass': 88.93762 * u.u,
        'stable': False,
        'half-life': 1.5 * u.s,
    },

    'Ru-90': {
        'atomic number': 44,
        'mass number': 90,
        'mass': 89.9303444 * u.u,
        'stable': False,
        'half-life': 11.0 * u.s,
    },

    'Ru-91': {
        'atomic number': 44,
        'mass number': 91,
        'mass': 90.9267419 * u.u,
        'stable': False,
        'half-life': 8.0 * u.s,
    },

    'Ru-92': {
        'atomic number': 44,
        'mass number': 92,
        'mass': 91.9202344 * u.u,
        'stable': False,
        'half-life': 219.0 * u.s,
    },

    'Ru-93': {
        'atomic number': 44,
        'mass number': 93,
        'mass': 92.9171044 * u.u,
        'stable': False,
        'half-life': 59.7 * u.s,
    },

    'Ru-94': {
        'atomic number': 44,
        'mass number': 94,
        'mass': 93.9113429 * u.u,
        'stable': False,
        'half-life': 3108.0 * u.s,
    },

    'Ru-95': {
        'atomic number': 44,
        'mass number': 95,
        'mass': 94.910406 * u.u,
        'stable': False,
        'half-life': 5914.8 * u.s,
    },

    'Ru-96': {
        'atomic number': 44,
        'mass number': 96,
        'mass': 95.90759025 * u.u,
        'stable': True,
        'abundance': 0.0554,
    },

    'Ru-97': {
        'atomic number': 44,
        'mass number': 97,
        'mass': 96.9075471 * u.u,
        'stable': False,
        'half-life': 245116.8 * u.s,
    },

    'Ru-98': {
        'atomic number': 44,
        'mass number': 98,
        'mass': 97.9052868 * u.u,
        'stable': True,
        'abundance': 0.0187,
    },

    'Ru-99': {
        'atomic number': 44,
        'mass number': 99,
        'mass': 98.9059341 * u.u,
        'stable': True,
        'abundance': 0.1276,
    },

    'Ru-100': {
        'atomic number': 44,
        'mass number': 100,
        'mass': 99.9042143 * u.u,
        'stable': True,
        'abundance': 0.126,
    },

    'Ru-101': {
        'atomic number': 44,
        'mass number': 101,
        'mass': 100.9055769 * u.u,
        'stable': True,
        'abundance': 0.1706,
    },

    'Ru-102': {
        'atomic number': 44,
        'mass number': 102,
        'mass': 101.9043441 * u.u,
        'stable': True,
        'abundance': 0.3155,
    },

    'Ru-103': {
        'atomic number': 44,
        'mass number': 103,
        'mass': 102.9063186 * u.u,
        'stable': False,
        'half-life': 3396384.0 * u.s,
    },

    'Ru-104': {
        'atomic number': 44,
        'mass number': 104,
        'mass': 103.9054275 * u.u,
        'stable': True,
        'abundance': 0.1862,
    },

    'Ru-105': {
        'atomic number': 44,
        'mass number': 105,
        'mass': 104.9077476 * u.u,
        'stable': False,
        'half-life': 15984.0 * u.s,
    },

    'Ru-106': {
        'atomic number': 44,
        'mass number': 106,
        'mass': 105.9073291 * u.u,
        'stable': False,
        'half-life': 32123520.0 * u.s,
    },

    'Ru-107': {
        'atomic number': 44,
        'mass number': 107,
        'mass': 106.909972 * u.u,
        'stable': False,
        'half-life': 225.0 * u.s,
    },

    'Ru-108': {
        'atomic number': 44,
        'mass number': 108,
        'mass': 107.910188 * u.u,
        'stable': False,
        'half-life': 273.0 * u.s,
    },

    'Ru-109': {
        'atomic number': 44,
        'mass number': 109,
        'mass': 108.913326 * u.u,
        'stable': False,
        'half-life': 34.5 * u.s,
    },

    'Ru-110': {
        'atomic number': 44,
        'mass number': 110,
        'mass': 109.9140407 * u.u,
        'stable': False,
        'half-life': 12.04 * u.s,
    },

    'Ru-111': {
        'atomic number': 44,
        'mass number': 111,
        'mass': 110.91757 * u.u,
        'stable': False,
        'half-life': 2.12 * u.s,
    },

    'Ru-112': {
        'atomic number': 44,
        'mass number': 112,
        'mass': 111.918809 * u.u,
        'stable': False,
        'half-life': 1.75 * u.s,
    },

    'Ru-113': {
        'atomic number': 44,
        'mass number': 113,
        'mass': 112.922844 * u.u,
        'stable': False,
        'half-life': 0.8 * u.s,
    },

    'Ru-114': {
        'atomic number': 44,
        'mass number': 114,
        'mass': 113.9246136 * u.u,
        'stable': False,
        'half-life': 0.54 * u.s,
    },

    'Ru-115': {
        'atomic number': 44,
        'mass number': 115,
        'mass': 114.92882 * u.u,
        'stable': False,
        'half-life': 0.318 * u.s,
    },

    'Ru-116': {
        'atomic number': 44,
        'mass number': 116,
        'mass': 115.9312192 * u.u,
        'stable': False,
        'half-life': 0.204 * u.s,
    },

    'Ru-117': {
        'atomic number': 44,
        'mass number': 117,
        'mass': 116.9361 * u.u,
        'stable': False,
        'half-life': 0.151 * u.s,
    },

    'Ru-118': {
        'atomic number': 44,
        'mass number': 118,
        'mass': 117.93853 * u.u,
        'stable': False,
        'half-life': 0.099 * u.s,
    },

    'Ru-119': {
        'atomic number': 44,
        'mass number': 119,
        'mass': 118.94357 * u.u,
        'stable': False,
        'half-life': 0.0695 * u.s,
    },

    'Ru-120': {
        'atomic number': 44,
        'mass number': 120,
        'mass': 119.94631 * u.u,
        'stable': False,
        'half-life': 0.045 * u.s,
    },

    'Ru-121': {
        'atomic number': 44,
        'mass number': 121,
        'mass': 120.95164 * u.u,
        'stable': False,
        'half-life': 0.029 * u.s,
    },

    'Ru-122': {
        'atomic number': 44,
        'mass number': 122,
        'mass': 121.95447 * u.u,
        'stable': False,
        'half-life': 0.025 * u.s,
    },

    'Ru-123': {
        'atomic number': 44,
        'mass number': 123,
        'mass': 122.95989 * u.u,
        'stable': False,
        'half-life': 0.019 * u.s,
    },

    'Ru-124': {
        'atomic number': 44,
        'mass number': 124,
        'mass': 123.96305 * u.u,
        'stable': False,
        'half-life': 0.015 * u.s,
    },

    'Rh-89': {
        'atomic number': 45,
        'mass number': 89,
        'mass': 88.95058 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Rh-90': {
        'atomic number': 45,
        'mass number': 90,
        'mass': 89.94422 * u.u,
        'stable': False,
        'half-life': 0.015 * u.s,
    },

    'Rh-91': {
        'atomic number': 45,
        'mass number': 91,
        'mass': 90.93688 * u.u,
        'stable': False,
        'half-life': 1.6 * u.s,
    },

    'Rh-92': {
        'atomic number': 45,
        'mass number': 92,
        'mass': 91.9323677 * u.u,
        'stable': False,
        'half-life': 4.66 * u.s,
    },

    'Rh-93': {
        'atomic number': 45,
        'mass number': 93,
        'mass': 92.9259128 * u.u,
        'stable': False,
        'half-life': 13.9 * u.s,
    },

    'Rh-94': {
        'atomic number': 45,
        'mass number': 94,
        'mass': 93.9217305 * u.u,
        'stable': False,
        'half-life': 70.6 * u.s,
    },

    'Rh-95': {
        'atomic number': 45,
        'mass number': 95,
        'mass': 94.9158979 * u.u,
        'stable': False,
        'half-life': 301.2 * u.s,
    },

    'Rh-96': {
        'atomic number': 45,
        'mass number': 96,
        'mass': 95.914453 * u.u,
        'stable': False,
        'half-life': 594.0 * u.s,
    },

    'Rh-97': {
        'atomic number': 45,
        'mass number': 97,
        'mass': 96.911329 * u.u,
        'stable': False,
        'half-life': 1842.0 * u.s,
    },

    'Rh-98': {
        'atomic number': 45,
        'mass number': 98,
        'mass': 97.910708 * u.u,
        'stable': False,
        'half-life': 523.2 * u.s,
    },

    'Rh-99': {
        'atomic number': 45,
        'mass number': 99,
        'mass': 98.9081282 * u.u,
        'stable': False,
        'half-life': 1391040.0 * u.s,
    },

    'Rh-100': {
        'atomic number': 45,
        'mass number': 100,
        'mass': 99.908117 * u.u,
        'stable': False,
        'half-life': 74880.0 * u.s,
    },

    'Rh-101': {
        'atomic number': 45,
        'mass number': 101,
        'mass': 100.9061606 * u.u,
        'stable': False,
        'half-life': 104137855.8 * u.s,
    },

    'Rh-102': {
        'atomic number': 45,
        'mass number': 102,
        'mass': 101.9068374 * u.u,
        'stable': False,
        'half-life': 17884800.0 * u.s,
    },

    'Rh-103': {
        'atomic number': 45,
        'mass number': 103,
        'mass': 102.905498 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Rh-104': {
        'atomic number': 45,
        'mass number': 104,
        'mass': 103.9066492 * u.u,
        'stable': False,
        'half-life': 42.3 * u.s,
    },

    'Rh-105': {
        'atomic number': 45,
        'mass number': 105,
        'mass': 104.9056885 * u.u,
        'stable': False,
        'half-life': 127285.2 * u.s,
    },

    'Rh-106': {
        'atomic number': 45,
        'mass number': 106,
        'mass': 105.9072868 * u.u,
        'stable': False,
        'half-life': 30.07 * u.s,
    },

    'Rh-107': {
        'atomic number': 45,
        'mass number': 107,
        'mass': 106.906748 * u.u,
        'stable': False,
        'half-life': 1302.0 * u.s,
    },

    'Rh-108': {
        'atomic number': 45,
        'mass number': 108,
        'mass': 107.908714 * u.u,
        'stable': False,
        'half-life': 16.8 * u.s,
    },

    'Rh-109': {
        'atomic number': 45,
        'mass number': 109,
        'mass': 108.9087488 * u.u,
        'stable': False,
        'half-life': 80.0 * u.s,
    },

    'Rh-110': {
        'atomic number': 45,
        'mass number': 110,
        'mass': 109.911079 * u.u,
        'stable': False,
        'half-life': 3.35 * u.s,
    },

    'Rh-111': {
        'atomic number': 45,
        'mass number': 111,
        'mass': 110.9116423 * u.u,
        'stable': False,
        'half-life': 11.0 * u.s,
    },

    'Rh-112': {
        'atomic number': 45,
        'mass number': 112,
        'mass': 111.914403 * u.u,
        'stable': False,
        'half-life': 3.4 * u.s,
    },

    'Rh-113': {
        'atomic number': 45,
        'mass number': 113,
        'mass': 112.9154393 * u.u,
        'stable': False,
        'half-life': 2.8 * u.s,
    },

    'Rh-114': {
        'atomic number': 45,
        'mass number': 114,
        'mass': 113.918718 * u.u,
        'stable': False,
        'half-life': 1.85 * u.s,
    },

    'Rh-115': {
        'atomic number': 45,
        'mass number': 115,
        'mass': 114.9203116 * u.u,
        'stable': False,
        'half-life': 0.99 * u.s,
    },

    'Rh-116': {
        'atomic number': 45,
        'mass number': 116,
        'mass': 115.924059 * u.u,
        'stable': False,
        'half-life': 0.685 * u.s,
    },

    'Rh-117': {
        'atomic number': 45,
        'mass number': 117,
        'mass': 116.9260354 * u.u,
        'stable': False,
        'half-life': 0.421 * u.s,
    },

    'Rh-118': {
        'atomic number': 45,
        'mass number': 118,
        'mass': 117.93034 * u.u,
        'stable': False,
        'half-life': 0.284 * u.s,
    },

    'Rh-119': {
        'atomic number': 45,
        'mass number': 119,
        'mass': 118.932557 * u.u,
        'stable': False,
        'half-life': 0.19 * u.s,
    },

    'Rh-120': {
        'atomic number': 45,
        'mass number': 120,
        'mass': 119.93686 * u.u,
        'stable': False,
        'half-life': 0.1296 * u.s,
    },

    'Rh-121': {
        'atomic number': 45,
        'mass number': 121,
        'mass': 120.93942 * u.u,
        'stable': False,
        'half-life': 0.076 * u.s,
    },

    'Rh-122': {
        'atomic number': 45,
        'mass number': 122,
        'mass': 121.94399 * u.u,
        'stable': False,
        'half-life': 0.051 * u.s,
    },

    'Rh-123': {
        'atomic number': 45,
        'mass number': 123,
        'mass': 122.94685 * u.u,
        'stable': False,
        'half-life': 0.042 * u.s,
    },

    'Rh-124': {
        'atomic number': 45,
        'mass number': 124,
        'mass': 123.95151 * u.u,
        'stable': False,
        'half-life': 0.03 * u.s,
    },

    'Rh-125': {
        'atomic number': 45,
        'mass number': 125,
        'mass': 124.95469 * u.u,
        'stable': False,
        'half-life': 0.0265 * u.s,
    },

    'Rh-126': {
        'atomic number': 45,
        'mass number': 126,
        'mass': 125.95946 * u.u,
        'stable': False,
        'half-life': 0.019 * u.s,
    },

    'Pd-91': {
        'atomic number': 46,
        'mass number': 91,
        'mass': 90.95032 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Pd-92': {
        'atomic number': 46,
        'mass number': 92,
        'mass': 91.94088 * u.u,
        'stable': False,
        'half-life': 1.1 * u.s,
    },

    'Pd-93': {
        'atomic number': 46,
        'mass number': 93,
        'mass': 92.93651 * u.u,
        'stable': False,
        'half-life': 1.15 * u.s,
    },

    'Pd-94': {
        'atomic number': 46,
        'mass number': 94,
        'mass': 93.9290376 * u.u,
        'stable': False,
        'half-life': 9.0 * u.s,
    },

    'Pd-95': {
        'atomic number': 46,
        'mass number': 95,
        'mass': 94.9248898 * u.u,
        'stable': False,
        'half-life': 7.5 * u.s,
    },

    'Pd-96': {
        'atomic number': 46,
        'mass number': 96,
        'mass': 95.9182151 * u.u,
        'stable': False,
        'half-life': 122.0 * u.s,
    },

    'Pd-97': {
        'atomic number': 46,
        'mass number': 97,
        'mass': 96.916472 * u.u,
        'stable': False,
        'half-life': 186.0 * u.s,
    },

    'Pd-98': {
        'atomic number': 46,
        'mass number': 98,
        'mass': 97.9126983 * u.u,
        'stable': False,
        'half-life': 1062.0 * u.s,
    },

    'Pd-99': {
        'atomic number': 46,
        'mass number': 99,
        'mass': 98.9117748 * u.u,
        'stable': False,
        'half-life': 1284.0 * u.s,
    },

    'Pd-100': {
        'atomic number': 46,
        'mass number': 100,
        'mass': 99.908505 * u.u,
        'stable': False,
        'half-life': 313632.0 * u.s,
    },

    'Pd-101': {
        'atomic number': 46,
        'mass number': 101,
        'mass': 100.9082864 * u.u,
        'stable': False,
        'half-life': 30492.0 * u.s,
    },

    'Pd-102': {
        'atomic number': 46,
        'mass number': 102,
        'mass': 101.9056022 * u.u,
        'stable': True,
        'abundance': 0.0102,
    },

    'Pd-103': {
        'atomic number': 46,
        'mass number': 103,
        'mass': 102.9060809 * u.u,
        'stable': False,
        'half-life': 1468022.4 * u.s,
    },

    'Pd-104': {
        'atomic number': 46,
        'mass number': 104,
        'mass': 103.9040305 * u.u,
        'stable': True,
        'abundance': 0.1114,
    },

    'Pd-105': {
        'atomic number': 46,
        'mass number': 105,
        'mass': 104.9050796 * u.u,
        'stable': True,
        'abundance': 0.2233,
    },

    'Pd-106': {
        'atomic number': 46,
        'mass number': 106,
        'mass': 105.9034804 * u.u,
        'stable': True,
        'abundance': 0.2733,
    },

    'Pd-107': {
        'atomic number': 46,
        'mass number': 107,
        'mass': 106.9051282 * u.u,
        'stable': False,
        'half-life': 205120019000000.0 * u.s,
    },

    'Pd-108': {
        'atomic number': 46,
        'mass number': 108,
        'mass': 107.9038916 * u.u,
        'stable': True,
        'abundance': 0.2646,
    },

    'Pd-109': {
        'atomic number': 46,
        'mass number': 109,
        'mass': 108.9059504 * u.u,
        'stable': False,
        'half-life': 49324.32 * u.s,
    },

    'Pd-110': {
        'atomic number': 46,
        'mass number': 110,
        'mass': 109.9051722 * u.u,
        'stable': True,
        'abundance': 0.1172,
    },

    'Pd-111': {
        'atomic number': 46,
        'mass number': 111,
        'mass': 110.90768968 * u.u,
        'stable': False,
        'half-life': 1404.0 * u.s,
    },

    'Pd-112': {
        'atomic number': 46,
        'mass number': 112,
        'mass': 111.9073297 * u.u,
        'stable': False,
        'half-life': 75744.0 * u.s,
    },

    'Pd-113': {
        'atomic number': 46,
        'mass number': 113,
        'mass': 112.910261 * u.u,
        'stable': False,
        'half-life': 93.0 * u.s,
    },

    'Pd-114': {
        'atomic number': 46,
        'mass number': 114,
        'mass': 113.9103686 * u.u,
        'stable': False,
        'half-life': 145.2 * u.s,
    },

    'Pd-115': {
        'atomic number': 46,
        'mass number': 115,
        'mass': 114.913659 * u.u,
        'stable': False,
        'half-life': 25.0 * u.s,
    },

    'Pd-116': {
        'atomic number': 46,
        'mass number': 116,
        'mass': 115.914297 * u.u,
        'stable': False,
        'half-life': 11.8 * u.s,
    },

    'Pd-117': {
        'atomic number': 46,
        'mass number': 117,
        'mass': 116.9179547 * u.u,
        'stable': False,
        'half-life': 4.3 * u.s,
    },

    'Pd-118': {
        'atomic number': 46,
        'mass number': 118,
        'mass': 117.9190667 * u.u,
        'stable': False,
        'half-life': 1.9 * u.s,
    },

    'Pd-119': {
        'atomic number': 46,
        'mass number': 119,
        'mass': 118.9233402 * u.u,
        'stable': False,
        'half-life': 0.92 * u.s,
    },

    'Pd-120': {
        'atomic number': 46,
        'mass number': 120,
        'mass': 119.9245511 * u.u,
        'stable': False,
        'half-life': 0.492 * u.s,
    },

    'Pd-121': {
        'atomic number': 46,
        'mass number': 121,
        'mass': 120.9289503 * u.u,
        'stable': False,
        'half-life': 0.29 * u.s,
    },

    'Pd-122': {
        'atomic number': 46,
        'mass number': 122,
        'mass': 121.930632 * u.u,
        'stable': False,
        'half-life': 0.195 * u.s,
    },

    'Pd-123': {
        'atomic number': 46,
        'mass number': 123,
        'mass': 122.93514 * u.u,
        'stable': False,
        'half-life': 0.108 * u.s,
    },

    'Pd-124': {
        'atomic number': 46,
        'mass number': 124,
        'mass': 123.93714 * u.u,
        'stable': False,
        'half-life': 0.088 * u.s,
    },

    'Pd-125': {
        'atomic number': 46,
        'mass number': 125,
        'mass': 124.94179 * u.u,
        'stable': False,
        'half-life': 0.057 * u.s,
    },

    'Pd-126': {
        'atomic number': 46,
        'mass number': 126,
        'mass': 125.94416 * u.u,
        'stable': False,
        'half-life': 0.0486 * u.s,
    },

    'Pd-127': {
        'atomic number': 46,
        'mass number': 127,
        'mass': 126.94907 * u.u,
        'stable': False,
        'half-life': 0.038 * u.s,
    },

    'Pd-128': {
        'atomic number': 46,
        'mass number': 128,
        'mass': 127.95183 * u.u,
        'stable': False,
        'half-life': 0.035 * u.s,
    },

    'Ag-93': {
        'atomic number': 47,
        'mass number': 93,
        'mass': 92.95033 * u.u,
        'stable': False,
        'half-life': '20# ms',
    },

    'Ag-94': {
        'atomic number': 47,
        'mass number': 94,
        'mass': 93.94373 * u.u,
        'stable': False,
        'half-life': 0.037 * u.s,
    },

    'Ag-95': {
        'atomic number': 47,
        'mass number': 95,
        'mass': 94.93602 * u.u,
        'stable': False,
        'half-life': 1.76 * u.s,
    },

    'Ag-96': {
        'atomic number': 47,
        'mass number': 96,
        'mass': 95.930744 * u.u,
        'stable': False,
        'half-life': 4.44 * u.s,
    },

    'Ag-97': {
        'atomic number': 47,
        'mass number': 97,
        'mass': 96.92397 * u.u,
        'stable': False,
        'half-life': 25.5 * u.s,
    },

    'Ag-98': {
        'atomic number': 47,
        'mass number': 98,
        'mass': 97.92156 * u.u,
        'stable': False,
        'half-life': 47.5 * u.s,
    },

    'Ag-99': {
        'atomic number': 47,
        'mass number': 99,
        'mass': 98.9176458 * u.u,
        'stable': False,
        'half-life': 124.2 * u.s,
    },

    'Ag-100': {
        'atomic number': 47,
        'mass number': 100,
        'mass': 99.9161154 * u.u,
        'stable': False,
        'half-life': 120.6 * u.s,
    },

    'Ag-101': {
        'atomic number': 47,
        'mass number': 101,
        'mass': 100.912684 * u.u,
        'stable': False,
        'half-life': 666.0 * u.s,
    },

    'Ag-102': {
        'atomic number': 47,
        'mass number': 102,
        'mass': 101.9117047 * u.u,
        'stable': False,
        'half-life': 774.0 * u.s,
    },

    'Ag-103': {
        'atomic number': 47,
        'mass number': 103,
        'mass': 102.9089631 * u.u,
        'stable': False,
        'half-life': 3942.0 * u.s,
    },

    'Ag-104': {
        'atomic number': 47,
        'mass number': 104,
        'mass': 103.9086239 * u.u,
        'stable': False,
        'half-life': 4152.0 * u.s,
    },

    'Ag-105': {
        'atomic number': 47,
        'mass number': 105,
        'mass': 104.9065256 * u.u,
        'stable': False,
        'half-life': 3567456.0 * u.s,
    },

    'Ag-106': {
        'atomic number': 47,
        'mass number': 106,
        'mass': 105.9066636 * u.u,
        'stable': False,
        'half-life': 1437.6 * u.s,
    },

    'Ag-107': {
        'atomic number': 47,
        'mass number': 107,
        'mass': 106.9050916 * u.u,
        'stable': True,
        'abundance': 0.51839,
    },

    'Ag-108': {
        'atomic number': 47,
        'mass number': 108,
        'mass': 107.9059503 * u.u,
        'stable': False,
        'half-life': 142.92 * u.s,
    },

    'Ag-109': {
        'atomic number': 47,
        'mass number': 109,
        'mass': 108.9047553 * u.u,
        'stable': True,
        'abundance': 0.48161,
    },

    'Ag-110': {
        'atomic number': 47,
        'mass number': 110,
        'mass': 109.9061102 * u.u,
        'stable': False,
        'half-life': 24.56 * u.s,
    },

    'Ag-111': {
        'atomic number': 47,
        'mass number': 111,
        'mass': 110.9052959 * u.u,
        'stable': False,
        'half-life': 642211.2 * u.s,
    },

    'Ag-112': {
        'atomic number': 47,
        'mass number': 112,
        'mass': 111.9070486 * u.u,
        'stable': False,
        'half-life': 11268.0 * u.s,
    },

    'Ag-113': {
        'atomic number': 47,
        'mass number': 113,
        'mass': 112.906573 * u.u,
        'stable': False,
        'half-life': 19332.0 * u.s,
    },

    'Ag-114': {
        'atomic number': 47,
        'mass number': 114,
        'mass': 113.908823 * u.u,
        'stable': False,
        'half-life': 4.6 * u.s,
    },

    'Ag-115': {
        'atomic number': 47,
        'mass number': 115,
        'mass': 114.908767 * u.u,
        'stable': False,
        'half-life': 1200.0 * u.s,
    },

    'Ag-116': {
        'atomic number': 47,
        'mass number': 116,
        'mass': 115.9113868 * u.u,
        'stable': False,
        'half-life': 229.8 * u.s,
    },

    'Ag-117': {
        'atomic number': 47,
        'mass number': 117,
        'mass': 116.911774 * u.u,
        'stable': False,
        'half-life': 73.6 * u.s,
    },

    'Ag-118': {
        'atomic number': 47,
        'mass number': 118,
        'mass': 117.9145955 * u.u,
        'stable': False,
        'half-life': 3.76 * u.s,
    },

    'Ag-119': {
        'atomic number': 47,
        'mass number': 119,
        'mass': 118.91557 * u.u,
        'stable': False,
        'half-life': 6.0 * u.s,
    },

    'Ag-120': {
        'atomic number': 47,
        'mass number': 120,
        'mass': 119.9187848 * u.u,
        'stable': False,
        'half-life': 1.52 * u.s,
    },

    'Ag-121': {
        'atomic number': 47,
        'mass number': 121,
        'mass': 120.920125 * u.u,
        'stable': False,
        'half-life': 0.78 * u.s,
    },

    'Ag-122': {
        'atomic number': 47,
        'mass number': 122,
        'mass': 121.923664 * u.u,
        'stable': False,
        'half-life': 0.529 * u.s,
    },

    'Ag-123': {
        'atomic number': 47,
        'mass number': 123,
        'mass': 122.925337 * u.u,
        'stable': False,
        'half-life': 0.3 * u.s,
    },

    'Ag-124': {
        'atomic number': 47,
        'mass number': 124,
        'mass': 123.92893 * u.u,
        'stable': False,
        'half-life': 0.1779 * u.s,
    },

    'Ag-125': {
        'atomic number': 47,
        'mass number': 125,
        'mass': 124.93105 * u.u,
        'stable': False,
        'half-life': 0.159 * u.s,
    },

    'Ag-126': {
        'atomic number': 47,
        'mass number': 126,
        'mass': 125.93475 * u.u,
        'stable': False,
        'half-life': 0.0993 * u.s,
    },

    'Ag-127': {
        'atomic number': 47,
        'mass number': 127,
        'mass': 126.93711 * u.u,
        'stable': False,
        'half-life': 0.089 * u.s,
    },

    'Ag-128': {
        'atomic number': 47,
        'mass number': 128,
        'mass': 127.94106 * u.u,
        'stable': False,
        'half-life': 0.059 * u.s,
    },

    'Ag-129': {
        'atomic number': 47,
        'mass number': 129,
        'mass': 128.94395 * u.u,
        'stable': False,
        'half-life': 0.0499 * u.s,
    },

    'Ag-130': {
        'atomic number': 47,
        'mass number': 130,
        'mass': 129.9507 * u.u,
        'stable': False,
        'half-life': 0.0406 * u.s,
    },

    'Cd-95': {
        'atomic number': 48,
        'mass number': 95,
        'mass': 94.94994 * u.u,
        'stable': False,
        'half-life': 0.09 * u.s,
    },

    'Cd-96': {
        'atomic number': 48,
        'mass number': 96,
        'mass': 95.94034 * u.u,
        'stable': False,
        'half-life': 0.88 * u.s,
    },

    'Cd-97': {
        'atomic number': 48,
        'mass number': 97,
        'mass': 96.9351 * u.u,
        'stable': False,
        'half-life': 1.1 * u.s,
    },

    'Cd-98': {
        'atomic number': 48,
        'mass number': 98,
        'mass': 97.927389 * u.u,
        'stable': False,
        'half-life': 9.2 * u.s,
    },

    'Cd-99': {
        'atomic number': 48,
        'mass number': 99,
        'mass': 98.9249258 * u.u,
        'stable': False,
        'half-life': 16.0 * u.s,
    },

    'Cd-100': {
        'atomic number': 48,
        'mass number': 100,
        'mass': 99.9203488 * u.u,
        'stable': False,
        'half-life': 49.1 * u.s,
    },

    'Cd-101': {
        'atomic number': 48,
        'mass number': 101,
        'mass': 100.9185862 * u.u,
        'stable': False,
        'half-life': 81.6 * u.s,
    },

    'Cd-102': {
        'atomic number': 48,
        'mass number': 102,
        'mass': 101.914482 * u.u,
        'stable': False,
        'half-life': 330.0 * u.s,
    },

    'Cd-103': {
        'atomic number': 48,
        'mass number': 103,
        'mass': 102.9134165 * u.u,
        'stable': False,
        'half-life': 438.0 * u.s,
    },

    'Cd-104': {
        'atomic number': 48,
        'mass number': 104,
        'mass': 103.9098564 * u.u,
        'stable': False,
        'half-life': 3462.0 * u.s,
    },

    'Cd-105': {
        'atomic number': 48,
        'mass number': 105,
        'mass': 104.9094639 * u.u,
        'stable': False,
        'half-life': 3330.0 * u.s,
    },

    'Cd-106': {
        'atomic number': 48,
        'mass number': 106,
        'mass': 105.9064599 * u.u,
        'stable': True,
        'abundance': 0.0125,
    },

    'Cd-107': {
        'atomic number': 48,
        'mass number': 107,
        'mass': 106.9066121 * u.u,
        'stable': False,
        'half-life': 23400.0 * u.s,
    },

    'Cd-108': {
        'atomic number': 48,
        'mass number': 108,
        'mass': 107.9041834 * u.u,
        'stable': True,
        'abundance': 0.0089,
    },

    'Cd-109': {
        'atomic number': 48,
        'mass number': 109,
        'mass': 108.9049867 * u.u,
        'stable': False,
        'half-life': 40025664.0 * u.s,
    },

    'Cd-110': {
        'atomic number': 48,
        'mass number': 110,
        'mass': 109.90300661 * u.u,
        'stable': True,
        'abundance': 0.1249,
    },

    'Cd-111': {
        'atomic number': 48,
        'mass number': 111,
        'mass': 110.90418287 * u.u,
        'stable': True,
        'abundance': 0.128,
    },

    'Cd-112': {
        'atomic number': 48,
        'mass number': 112,
        'mass': 111.90276287 * u.u,
        'stable': True,
        'abundance': 0.2413,
    },

    'Cd-113': {
        'atomic number': 48,
        'mass number': 113,
        'mass': 112.90440813 * u.u,
        'stable': False,
        'abundance': 0.1222,
    },

    'Cd-114': {
        'atomic number': 48,
        'mass number': 114,
        'mass': 113.90336509 * u.u,
        'stable': True,
        'abundance': 0.2873,
    },

    'Cd-115': {
        'atomic number': 48,
        'mass number': 115,
        'mass': 114.90543751 * u.u,
        'stable': False,
        'half-life': 192456.0 * u.s,
    },

    'Cd-116': {
        'atomic number': 48,
        'mass number': 116,
        'mass': 115.90476315 * u.u,
        'stable': False,
        'abundance': 0.0749,
    },

    'Cd-117': {
        'atomic number': 48,
        'mass number': 117,
        'mass': 116.907226 * u.u,
        'stable': False,
        'half-life': 8964.0 * u.s,
    },

    'Cd-118': {
        'atomic number': 48,
        'mass number': 118,
        'mass': 117.906922 * u.u,
        'stable': False,
        'half-life': 3018.0 * u.s,
    },

    'Cd-119': {
        'atomic number': 48,
        'mass number': 119,
        'mass': 118.909847 * u.u,
        'stable': False,
        'half-life': 161.4 * u.s,
    },

    'Cd-120': {
        'atomic number': 48,
        'mass number': 120,
        'mass': 119.9098681 * u.u,
        'stable': False,
        'half-life': 50.8 * u.s,
    },

    'Cd-121': {
        'atomic number': 48,
        'mass number': 121,
        'mass': 120.9129637 * u.u,
        'stable': False,
        'half-life': 13.5 * u.s,
    },

    'Cd-122': {
        'atomic number': 48,
        'mass number': 122,
        'mass': 121.9134591 * u.u,
        'stable': False,
        'half-life': 5.24 * u.s,
    },

    'Cd-123': {
        'atomic number': 48,
        'mass number': 123,
        'mass': 122.9168925 * u.u,
        'stable': False,
        'half-life': 2.1 * u.s,
    },

    'Cd-124': {
        'atomic number': 48,
        'mass number': 124,
        'mass': 123.9176574 * u.u,
        'stable': False,
        'half-life': 1.25 * u.s,
    },

    'Cd-125': {
        'atomic number': 48,
        'mass number': 125,
        'mass': 124.9212576 * u.u,
        'stable': False,
        'half-life': 0.68 * u.s,
    },

    'Cd-126': {
        'atomic number': 48,
        'mass number': 126,
        'mass': 125.9224291 * u.u,
        'stable': False,
        'half-life': 0.513 * u.s,
    },

    'Cd-127': {
        'atomic number': 48,
        'mass number': 127,
        'mass': 126.926472 * u.u,
        'stable': False,
        'half-life': 0.33 * u.s,
    },

    'Cd-128': {
        'atomic number': 48,
        'mass number': 128,
        'mass': 127.9278129 * u.u,
        'stable': False,
        'half-life': 0.246 * u.s,
    },

    'Cd-129': {
        'atomic number': 48,
        'mass number': 129,
        'mass': 128.93182 * u.u,
        'stable': False,
        'half-life': 0.1515 * u.s,
    },

    'Cd-130': {
        'atomic number': 48,
        'mass number': 130,
        'mass': 129.93394 * u.u,
        'stable': False,
        'half-life': 0.1268 * u.s,
    },

    'Cd-131': {
        'atomic number': 48,
        'mass number': 131,
        'mass': 130.9406 * u.u,
        'stable': False,
        'half-life': 0.098 * u.s,
    },

    'Cd-132': {
        'atomic number': 48,
        'mass number': 132,
        'mass': 131.94604 * u.u,
        'stable': False,
        'half-life': 0.082 * u.s,
    },

    'Cd-133': {
        'atomic number': 48,
        'mass number': 133,
        'mass': 132.95285 * u.u,
        'stable': False,
        'half-life': 0.061 * u.s,
    },

    'In-97': {
        'atomic number': 49,
        'mass number': 97,
        'mass': 96.94934 * u.u,
        'stable': False,
        'half-life': 0.05 * u.s,
    },

    'In-98': {
        'atomic number': 49,
        'mass number': 98,
        'mass': 97.94214 * u.u,
        'stable': False,
        'half-life': 0.037 * u.s,
    },

    'In-99': {
        'atomic number': 49,
        'mass number': 99,
        'mass': 98.93411 * u.u,
        'stable': False,
        'half-life': 3.1 * u.s,
    },

    'In-100': {
        'atomic number': 49,
        'mass number': 100,
        'mass': 99.93096 * u.u,
        'stable': False,
        'half-life': 5.83 * u.s,
    },

    'In-101': {
        'atomic number': 49,
        'mass number': 101,
        'mass': 100.92634 * u.u,
        'stable': False,
        'half-life': 15.1 * u.s,
    },

    'In-102': {
        'atomic number': 49,
        'mass number': 102,
        'mass': 101.9241071 * u.u,
        'stable': False,
        'half-life': 23.3 * u.s,
    },

    'In-103': {
        'atomic number': 49,
        'mass number': 103,
        'mass': 102.9198819 * u.u,
        'stable': False,
        'half-life': 60.0 * u.s,
    },

    'In-104': {
        'atomic number': 49,
        'mass number': 104,
        'mass': 103.9182145 * u.u,
        'stable': False,
        'half-life': 108.0 * u.s,
    },

    'In-105': {
        'atomic number': 49,
        'mass number': 105,
        'mass': 104.914502 * u.u,
        'stable': False,
        'half-life': 304.2 * u.s,
    },

    'In-106': {
        'atomic number': 49,
        'mass number': 106,
        'mass': 105.913464 * u.u,
        'stable': False,
        'half-life': 372.0 * u.s,
    },

    'In-107': {
        'atomic number': 49,
        'mass number': 107,
        'mass': 106.91029 * u.u,
        'stable': False,
        'half-life': 1944.0 * u.s,
    },

    'In-108': {
        'atomic number': 49,
        'mass number': 108,
        'mass': 107.9096935 * u.u,
        'stable': False,
        'half-life': 3480.0 * u.s,
    },

    'In-109': {
        'atomic number': 49,
        'mass number': 109,
        'mass': 108.9071514 * u.u,
        'stable': False,
        'half-life': 15001.2 * u.s,
    },

    'In-110': {
        'atomic number': 49,
        'mass number': 110,
        'mass': 109.90717 * u.u,
        'stable': False,
        'half-life': 17712.0 * u.s,
    },

    'In-111': {
        'atomic number': 49,
        'mass number': 111,
        'mass': 110.9051085 * u.u,
        'stable': False,
        'half-life': 242332.128 * u.s,
    },

    'In-112': {
        'atomic number': 49,
        'mass number': 112,
        'mass': 111.9055377 * u.u,
        'stable': False,
        'half-life': 892.8 * u.s,
    },

    'In-113': {
        'atomic number': 49,
        'mass number': 113,
        'mass': 112.90406184 * u.u,
        'stable': True,
        'abundance': 0.0429,
    },

    'In-114': {
        'atomic number': 49,
        'mass number': 114,
        'mass': 113.90491791 * u.u,
        'stable': False,
        'half-life': 71.9 * u.s,
    },

    'In-115': {
        'atomic number': 49,
        'mass number': 115,
        'mass': 114.903878776 * u.u,
        'stable': False,
        'abundance': 0.9571,
    },

    'In-116': {
        'atomic number': 49,
        'mass number': 116,
        'mass': 115.90525999 * u.u,
        'stable': False,
        'half-life': 14.1 * u.s,
    },

    'In-117': {
        'atomic number': 49,
        'mass number': 117,
        'mass': 116.9045157 * u.u,
        'stable': False,
        'half-life': 2592.0 * u.s,
    },

    'In-118': {
        'atomic number': 49,
        'mass number': 118,
        'mass': 117.9063566 * u.u,
        'stable': False,
        'half-life': 5.0 * u.s,
    },

    'In-119': {
        'atomic number': 49,
        'mass number': 119,
        'mass': 118.9058507 * u.u,
        'stable': False,
        'half-life': 144.0 * u.s,
    },

    'In-120': {
        'atomic number': 49,
        'mass number': 120,
        'mass': 119.907967 * u.u,
        'stable': False,
        'half-life': 3.08 * u.s,
    },

    'In-121': {
        'atomic number': 49,
        'mass number': 121,
        'mass': 120.907851 * u.u,
        'stable': False,
        'half-life': 23.1 * u.s,
    },

    'In-122': {
        'atomic number': 49,
        'mass number': 122,
        'mass': 121.910281 * u.u,
        'stable': False,
        'half-life': 1.5 * u.s,
    },

    'In-123': {
        'atomic number': 49,
        'mass number': 123,
        'mass': 122.910434 * u.u,
        'stable': False,
        'half-life': 6.17 * u.s,
    },

    'In-124': {
        'atomic number': 49,
        'mass number': 124,
        'mass': 123.913182 * u.u,
        'stable': False,
        'half-life': 3.12 * u.s,
    },

    'In-125': {
        'atomic number': 49,
        'mass number': 125,
        'mass': 124.913605 * u.u,
        'stable': False,
        'half-life': 2.36 * u.s,
    },

    'In-126': {
        'atomic number': 49,
        'mass number': 126,
        'mass': 125.916507 * u.u,
        'stable': False,
        'half-life': 1.53 * u.s,
    },

    'In-127': {
        'atomic number': 49,
        'mass number': 127,
        'mass': 126.917446 * u.u,
        'stable': False,
        'half-life': 1.09 * u.s,
    },

    'In-128': {
        'atomic number': 49,
        'mass number': 128,
        'mass': 127.9204 * u.u,
        'stable': False,
        'half-life': 0.816 * u.s,
    },

    'In-129': {
        'atomic number': 49,
        'mass number': 129,
        'mass': 128.9218053 * u.u,
        'stable': False,
        'half-life': 0.57 * u.s,
    },

    'In-130': {
        'atomic number': 49,
        'mass number': 130,
        'mass': 129.924977 * u.u,
        'stable': False,
        'half-life': 0.284 * u.s,
    },

    'In-131': {
        'atomic number': 49,
        'mass number': 131,
        'mass': 130.9269715 * u.u,
        'stable': False,
        'half-life': 0.261 * u.s,
    },

    'In-132': {
        'atomic number': 49,
        'mass number': 132,
        'mass': 131.933001 * u.u,
        'stable': False,
        'half-life': 0.198 * u.s,
    },

    'In-133': {
        'atomic number': 49,
        'mass number': 133,
        'mass': 132.93831 * u.u,
        'stable': False,
        'half-life': 0.165 * u.s,
    },

    'In-134': {
        'atomic number': 49,
        'mass number': 134,
        'mass': 133.94454 * u.u,
        'stable': False,
        'half-life': 0.14 * u.s,
    },

    'In-135': {
        'atomic number': 49,
        'mass number': 135,
        'mass': 134.95005 * u.u,
        'stable': False,
        'half-life': 0.101 * u.s,
    },

    'Sn-99': {
        'atomic number': 50,
        'mass number': 99,
        'mass': 98.94853 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Sn-100': {
        'atomic number': 50,
        'mass number': 100,
        'mass': 99.9385 * u.u,
        'stable': False,
        'half-life': 1.16 * u.s,
    },

    'Sn-101': {
        'atomic number': 50,
        'mass number': 101,
        'mass': 100.93526 * u.u,
        'stable': False,
        'half-life': 1.97 * u.s,
    },

    'Sn-102': {
        'atomic number': 50,
        'mass number': 102,
        'mass': 101.93029 * u.u,
        'stable': False,
        'half-life': 3.8 * u.s,
    },

    'Sn-103': {
        'atomic number': 50,
        'mass number': 103,
        'mass': 102.928105 * u.u,
        'stable': False,
        'half-life': 7.0 * u.s,
    },

    'Sn-104': {
        'atomic number': 50,
        'mass number': 104,
        'mass': 103.9231052 * u.u,
        'stable': False,
        'half-life': 20.8 * u.s,
    },

    'Sn-105': {
        'atomic number': 50,
        'mass number': 105,
        'mass': 104.9212684 * u.u,
        'stable': False,
        'half-life': 34.0 * u.s,
    },

    'Sn-106': {
        'atomic number': 50,
        'mass number': 106,
        'mass': 105.9169574 * u.u,
        'stable': False,
        'half-life': 115.2 * u.s,
    },

    'Sn-107': {
        'atomic number': 50,
        'mass number': 107,
        'mass': 106.9157137 * u.u,
        'stable': False,
        'half-life': 174.0 * u.s,
    },

    'Sn-108': {
        'atomic number': 50,
        'mass number': 108,
        'mass': 107.9118943 * u.u,
        'stable': False,
        'half-life': 618.0 * u.s,
    },

    'Sn-109': {
        'atomic number': 50,
        'mass number': 109,
        'mass': 108.9112921 * u.u,
        'stable': False,
        'half-life': 1080.0 * u.s,
    },

    'Sn-110': {
        'atomic number': 50,
        'mass number': 110,
        'mass': 109.907845 * u.u,
        'stable': False,
        'half-life': 14954.4 * u.s,
    },

    'Sn-111': {
        'atomic number': 50,
        'mass number': 111,
        'mass': 110.9077401 * u.u,
        'stable': False,
        'half-life': 2118.0 * u.s,
    },

    'Sn-112': {
        'atomic number': 50,
        'mass number': 112,
        'mass': 111.90482387 * u.u,
        'stable': True,
        'abundance': 0.0097,
    },

    'Sn-113': {
        'atomic number': 50,
        'mass number': 113,
        'mass': 112.9051757 * u.u,
        'stable': False,
        'half-life': 9942825.6 * u.s,
    },

    'Sn-114': {
        'atomic number': 50,
        'mass number': 114,
        'mass': 113.9027827 * u.u,
        'stable': True,
        'abundance': 0.0066,
    },

    'Sn-115': {
        'atomic number': 50,
        'mass number': 115,
        'mass': 114.903344699 * u.u,
        'stable': True,
        'abundance': 0.0034,
    },

    'Sn-116': {
        'atomic number': 50,
        'mass number': 116,
        'mass': 115.9017428 * u.u,
        'stable': True,
        'abundance': 0.1454,
    },

    'Sn-117': {
        'atomic number': 50,
        'mass number': 117,
        'mass': 116.90295398 * u.u,
        'stable': True,
        'abundance': 0.0768,
    },

    'Sn-118': {
        'atomic number': 50,
        'mass number': 118,
        'mass': 117.90160657 * u.u,
        'stable': True,
        'abundance': 0.2422,
    },

    'Sn-119': {
        'atomic number': 50,
        'mass number': 119,
        'mass': 118.90331117 * u.u,
        'stable': True,
        'abundance': 0.0859,
    },

    'Sn-120': {
        'atomic number': 50,
        'mass number': 120,
        'mass': 119.90220163 * u.u,
        'stable': True,
        'abundance': 0.3258,
    },

    'Sn-121': {
        'atomic number': 50,
        'mass number': 121,
        'mass': 120.9042426 * u.u,
        'stable': False,
        'half-life': 97308.0 * u.s,
    },

    'Sn-122': {
        'atomic number': 50,
        'mass number': 122,
        'mass': 121.9034438 * u.u,
        'stable': True,
        'abundance': 0.0463,
    },

    'Sn-123': {
        'atomic number': 50,
        'mass number': 123,
        'mass': 122.9057252 * u.u,
        'stable': False,
        'half-life': 11162880.0 * u.s,
    },

    'Sn-124': {
        'atomic number': 50,
        'mass number': 124,
        'mass': 123.9052766 * u.u,
        'stable': True,
        'abundance': 0.0579,
    },

    'Sn-125': {
        'atomic number': 50,
        'mass number': 125,
        'mass': 124.9077864 * u.u,
        'stable': False,
        'half-life': 832896.0 * u.s,
    },

    'Sn-126': {
        'atomic number': 50,
        'mass number': 126,
        'mass': 125.907659 * u.u,
        'stable': False,
        'half-life': 7258092980000.0 * u.s,
    },

    'Sn-127': {
        'atomic number': 50,
        'mass number': 127,
        'mass': 126.91039 * u.u,
        'stable': False,
        'half-life': 7560.0 * u.s,
    },

    'Sn-128': {
        'atomic number': 50,
        'mass number': 128,
        'mass': 127.910507 * u.u,
        'stable': False,
        'half-life': 3544.2 * u.s,
    },

    'Sn-129': {
        'atomic number': 50,
        'mass number': 129,
        'mass': 128.913465 * u.u,
        'stable': False,
        'half-life': 133.8 * u.s,
    },

    'Sn-130': {
        'atomic number': 50,
        'mass number': 130,
        'mass': 129.9139738 * u.u,
        'stable': False,
        'half-life': 223.2 * u.s,
    },

    'Sn-131': {
        'atomic number': 50,
        'mass number': 131,
        'mass': 130.917045 * u.u,
        'stable': False,
        'half-life': 56.0 * u.s,
    },

    'Sn-132': {
        'atomic number': 50,
        'mass number': 132,
        'mass': 131.9178267 * u.u,
        'stable': False,
        'half-life': 39.7 * u.s,
    },

    'Sn-133': {
        'atomic number': 50,
        'mass number': 133,
        'mass': 132.9239134 * u.u,
        'stable': False,
        'half-life': 1.46 * u.s,
    },

    'Sn-134': {
        'atomic number': 50,
        'mass number': 134,
        'mass': 133.9286821 * u.u,
        'stable': False,
        'half-life': 0.89 * u.s,
    },

    'Sn-135': {
        'atomic number': 50,
        'mass number': 135,
        'mass': 134.9349086 * u.u,
        'stable': False,
        'half-life': 0.515 * u.s,
    },

    'Sn-136': {
        'atomic number': 50,
        'mass number': 136,
        'mass': 135.93999 * u.u,
        'stable': False,
        'half-life': 0.35 * u.s,
    },

    'Sn-137': {
        'atomic number': 50,
        'mass number': 137,
        'mass': 136.94655 * u.u,
        'stable': False,
        'half-life': 0.273 * u.s,
    },

    'Sn-138': {
        'atomic number': 50,
        'mass number': 138,
        'mass': 137.95184 * u.u,
        'stable': False,
        'half-life': 0.15 * u.s,
    },

    'Sb-103': {
        'atomic number': 51,
        'mass number': 103,
        'mass': 102.93969 * u.u,
        'stable': False,
    },

    'Sb-104': {
        'atomic number': 51,
        'mass number': 104,
        'mass': 103.93648 * u.u,
        'stable': False,
        'half-life': 0.47 * u.s,
    },

    'Sb-105': {
        'atomic number': 51,
        'mass number': 105,
        'mass': 104.931276 * u.u,
        'stable': False,
        'half-life': 1.12 * u.s,
    },

    'Sb-106': {
        'atomic number': 51,
        'mass number': 106,
        'mass': 105.928638 * u.u,
        'stable': False,
        'half-life': 0.6 * u.s,
    },

    'Sb-107': {
        'atomic number': 51,
        'mass number': 107,
        'mass': 106.9241506 * u.u,
        'stable': False,
        'half-life': 4.0 * u.s,
    },

    'Sb-108': {
        'atomic number': 51,
        'mass number': 108,
        'mass': 107.9222267 * u.u,
        'stable': False,
        'half-life': 7.4 * u.s,
    },

    'Sb-109': {
        'atomic number': 51,
        'mass number': 109,
        'mass': 108.9181411 * u.u,
        'stable': False,
        'half-life': 17.0 * u.s,
    },

    'Sb-110': {
        'atomic number': 51,
        'mass number': 110,
        'mass': 109.9168543 * u.u,
        'stable': False,
        'half-life': 23.6 * u.s,
    },

    'Sb-111': {
        'atomic number': 51,
        'mass number': 111,
        'mass': 110.9132182 * u.u,
        'stable': False,
        'half-life': 75.0 * u.s,
    },

    'Sb-112': {
        'atomic number': 51,
        'mass number': 112,
        'mass': 111.9124 * u.u,
        'stable': False,
        'half-life': 53.5 * u.s,
    },

    'Sb-113': {
        'atomic number': 51,
        'mass number': 113,
        'mass': 112.909375 * u.u,
        'stable': False,
        'half-life': 400.2 * u.s,
    },

    'Sb-114': {
        'atomic number': 51,
        'mass number': 114,
        'mass': 113.90929 * u.u,
        'stable': False,
        'half-life': 209.4 * u.s,
    },

    'Sb-115': {
        'atomic number': 51,
        'mass number': 115,
        'mass': 114.906598 * u.u,
        'stable': False,
        'half-life': 1926.0 * u.s,
    },

    'Sb-116': {
        'atomic number': 51,
        'mass number': 116,
        'mass': 115.9067931 * u.u,
        'stable': False,
        'half-life': 948.0 * u.s,
    },

    'Sb-117': {
        'atomic number': 51,
        'mass number': 117,
        'mass': 116.9048415 * u.u,
        'stable': False,
        'half-life': 10080.0 * u.s,
    },

    'Sb-118': {
        'atomic number': 51,
        'mass number': 118,
        'mass': 117.9055321 * u.u,
        'stable': False,
        'half-life': 216.0 * u.s,
    },

    'Sb-119': {
        'atomic number': 51,
        'mass number': 119,
        'mass': 118.9039455 * u.u,
        'stable': False,
        'half-life': 137484.0 * u.s,
    },

    'Sb-120': {
        'atomic number': 51,
        'mass number': 120,
        'mass': 119.9050794 * u.u,
        'stable': False,
        'half-life': 953.4 * u.s,
    },

    'Sb-121': {
        'atomic number': 51,
        'mass number': 121,
        'mass': 120.903812 * u.u,
        'stable': True,
        'abundance': 0.5721,
    },

    'Sb-122': {
        'atomic number': 51,
        'mass number': 122,
        'mass': 121.9051699 * u.u,
        'stable': False,
        'half-life': 235336.32 * u.s,
    },

    'Sb-123': {
        'atomic number': 51,
        'mass number': 123,
        'mass': 122.9042132 * u.u,
        'stable': True,
        'abundance': 0.4279,
    },

    'Sb-124': {
        'atomic number': 51,
        'mass number': 124,
        'mass': 123.905935 * u.u,
        'stable': False,
        'half-life': 5201280.0 * u.s,
    },

    'Sb-125': {
        'atomic number': 51,
        'mass number': 125,
        'mass': 124.905253 * u.u,
        'stable': False,
        'half-life': 87053184.0 * u.s,
    },

    'Sb-126': {
        'atomic number': 51,
        'mass number': 126,
        'mass': 125.907253 * u.u,
        'stable': False,
        'half-life': 1067040.0 * u.s,
    },

    'Sb-127': {
        'atomic number': 51,
        'mass number': 127,
        'mass': 126.9069243 * u.u,
        'stable': False,
        'half-life': 332640.0 * u.s,
    },

    'Sb-128': {
        'atomic number': 51,
        'mass number': 128,
        'mass': 127.909146 * u.u,
        'stable': False,
        'half-life': 32580.0 * u.s,
    },

    'Sb-129': {
        'atomic number': 51,
        'mass number': 129,
        'mass': 128.909147 * u.u,
        'stable': False,
        'half-life': 15717.6 * u.s,
    },

    'Sb-130': {
        'atomic number': 51,
        'mass number': 130,
        'mass': 129.911662 * u.u,
        'stable': False,
        'half-life': 2370.0 * u.s,
    },

    'Sb-131': {
        'atomic number': 51,
        'mass number': 131,
        'mass': 130.9119888 * u.u,
        'stable': False,
        'half-life': 1381.8 * u.s,
    },

    'Sb-132': {
        'atomic number': 51,
        'mass number': 132,
        'mass': 131.9145077 * u.u,
        'stable': False,
        'half-life': 167.4 * u.s,
    },

    'Sb-133': {
        'atomic number': 51,
        'mass number': 133,
        'mass': 132.9152732 * u.u,
        'stable': False,
        'half-life': 140.4 * u.s,
    },

    'Sb-134': {
        'atomic number': 51,
        'mass number': 134,
        'mass': 133.9205357 * u.u,
        'stable': False,
        'half-life': 0.78 * u.s,
    },

    'Sb-135': {
        'atomic number': 51,
        'mass number': 135,
        'mass': 134.9251851 * u.u,
        'stable': False,
        'half-life': 1.679 * u.s,
    },

    'Sb-136': {
        'atomic number': 51,
        'mass number': 136,
        'mass': 135.9307459 * u.u,
        'stable': False,
        'half-life': 0.923 * u.s,
    },

    'Sb-137': {
        'atomic number': 51,
        'mass number': 137,
        'mass': 136.93555 * u.u,
        'stable': False,
        'half-life': 0.484 * u.s,
    },

    'Sb-138': {
        'atomic number': 51,
        'mass number': 138,
        'mass': 137.94145 * u.u,
        'stable': False,
        'half-life': 0.348 * u.s,
    },

    'Sb-139': {
        'atomic number': 51,
        'mass number': 139,
        'mass': 138.94655 * u.u,
        'stable': False,
        'half-life': 0.093 * u.s,
    },

    'Sb-140': {
        'atomic number': 51,
        'mass number': 140,
        'mass': 139.95283 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Te-105': {
        'atomic number': 52,
        'mass number': 105,
        'mass': 104.9433 * u.u,
        'stable': False,
        'half-life': 6.33e-07 * u.s,
    },

    'Te-106': {
        'atomic number': 52,
        'mass number': 106,
        'mass': 105.9375 * u.u,
        'stable': False,
        'half-life': 7.8e-05 * u.s,
    },

    'Te-107': {
        'atomic number': 52,
        'mass number': 107,
        'mass': 106.935012 * u.u,
        'stable': False,
        'half-life': 0.0031 * u.s,
    },

    'Te-108': {
        'atomic number': 52,
        'mass number': 108,
        'mass': 107.9293805 * u.u,
        'stable': False,
        'half-life': 2.1 * u.s,
    },

    'Te-109': {
        'atomic number': 52,
        'mass number': 109,
        'mass': 108.9273045 * u.u,
        'stable': False,
        'half-life': 4.6 * u.s,
    },

    'Te-110': {
        'atomic number': 52,
        'mass number': 110,
        'mass': 109.9224581 * u.u,
        'stable': False,
        'half-life': 18.6 * u.s,
    },

    'Te-111': {
        'atomic number': 52,
        'mass number': 111,
        'mass': 110.9210006 * u.u,
        'stable': False,
        'half-life': 26.2 * u.s,
    },

    'Te-112': {
        'atomic number': 52,
        'mass number': 112,
        'mass': 111.9167279 * u.u,
        'stable': False,
        'half-life': 120.0 * u.s,
    },

    'Te-113': {
        'atomic number': 52,
        'mass number': 113,
        'mass': 112.915891 * u.u,
        'stable': False,
        'half-life': 102.0 * u.s,
    },

    'Te-114': {
        'atomic number': 52,
        'mass number': 114,
        'mass': 113.912089 * u.u,
        'stable': False,
        'half-life': 912.0 * u.s,
    },

    'Te-115': {
        'atomic number': 52,
        'mass number': 115,
        'mass': 114.911902 * u.u,
        'stable': False,
        'half-life': 348.0 * u.s,
    },

    'Te-116': {
        'atomic number': 52,
        'mass number': 116,
        'mass': 115.90846 * u.u,
        'stable': False,
        'half-life': 8964.0 * u.s,
    },

    'Te-117': {
        'atomic number': 52,
        'mass number': 117,
        'mass': 116.908646 * u.u,
        'stable': False,
        'half-life': 3720.0 * u.s,
    },

    'Te-118': {
        'atomic number': 52,
        'mass number': 118,
        'mass': 117.905854 * u.u,
        'stable': False,
        'half-life': 518400.0 * u.s,
    },

    'Te-119': {
        'atomic number': 52,
        'mass number': 119,
        'mass': 118.9064071 * u.u,
        'stable': False,
        'half-life': 57780.0 * u.s,
    },

    'Te-120': {
        'atomic number': 52,
        'mass number': 120,
        'mass': 119.9040593 * u.u,
        'stable': True,
        'abundance': 0.0009,
    },

    'Te-121': {
        'atomic number': 52,
        'mass number': 121,
        'mass': 120.904944 * u.u,
        'stable': False,
        'half-life': 1656288.0 * u.s,
    },

    'Te-122': {
        'atomic number': 52,
        'mass number': 122,
        'mass': 121.9030435 * u.u,
        'stable': True,
        'abundance': 0.0255,
    },

    'Te-123': {
        'atomic number': 52,
        'mass number': 123,
        'mass': 122.9042698 * u.u,
        'stable': True,
        'abundance': 0.0089,
    },

    'Te-124': {
        'atomic number': 52,
        'mass number': 124,
        'mass': 123.9028171 * u.u,
        'stable': True,
        'abundance': 0.0474,
    },

    'Te-125': {
        'atomic number': 52,
        'mass number': 125,
        'mass': 124.9044299 * u.u,
        'stable': True,
        'abundance': 0.0707,
    },

    'Te-126': {
        'atomic number': 52,
        'mass number': 126,
        'mass': 125.9033109 * u.u,
        'stable': True,
        'abundance': 0.1884,
    },

    'Te-127': {
        'atomic number': 52,
        'mass number': 127,
        'mass': 126.9052257 * u.u,
        'stable': False,
        'half-life': 33660.0 * u.s,
    },

    'Te-128': {
        'atomic number': 52,
        'mass number': 128,
        'mass': 127.90446128 * u.u,
        'stable': False,
        'abundance': 0.3174,
    },

    'Te-129': {
        'atomic number': 52,
        'mass number': 129,
        'mass': 128.90659646 * u.u,
        'stable': False,
        'half-life': 4176.0 * u.s,
    },

    'Te-130': {
        'atomic number': 52,
        'mass number': 130,
        'mass': 129.906222748 * u.u,
        'stable': False,
        'abundance': 0.3408,
    },

    'Te-131': {
        'atomic number': 52,
        'mass number': 131,
        'mass': 130.908522213 * u.u,
        'stable': False,
        'half-life': 1500.0 * u.s,
    },

    'Te-132': {
        'atomic number': 52,
        'mass number': 132,
        'mass': 131.9085467 * u.u,
        'stable': False,
        'half-life': 276825.6 * u.s,
    },

    'Te-133': {
        'atomic number': 52,
        'mass number': 133,
        'mass': 132.9109688 * u.u,
        'stable': False,
        'half-life': 750.0 * u.s,
    },

    'Te-134': {
        'atomic number': 52,
        'mass number': 134,
        'mass': 133.911394 * u.u,
        'stable': False,
        'half-life': 2508.0 * u.s,
    },

    'Te-135': {
        'atomic number': 52,
        'mass number': 135,
        'mass': 134.9165557 * u.u,
        'stable': False,
        'half-life': 19.0 * u.s,
    },

    'Te-136': {
        'atomic number': 52,
        'mass number': 136,
        'mass': 135.9201006 * u.u,
        'stable': False,
        'half-life': 17.63 * u.s,
    },

    'Te-137': {
        'atomic number': 52,
        'mass number': 137,
        'mass': 136.9255989 * u.u,
        'stable': False,
        'half-life': 2.49 * u.s,
    },

    'Te-138': {
        'atomic number': 52,
        'mass number': 138,
        'mass': 137.9294722 * u.u,
        'stable': False,
        'half-life': 1.4 * u.s,
    },

    'Te-139': {
        'atomic number': 52,
        'mass number': 139,
        'mass': 138.9353672 * u.u,
        'stable': False,
        'half-life': '500# ms',
    },

    'Te-140': {
        'atomic number': 52,
        'mass number': 140,
        'mass': 139.939499 * u.u,
        'stable': False,
        'half-life': '300# ms',
    },

    'Te-141': {
        'atomic number': 52,
        'mass number': 141,
        'mass': 140.9458 * u.u,
        'stable': False,
        'half-life': '150# ms',
    },

    'Te-142': {
        'atomic number': 52,
        'mass number': 142,
        'mass': 141.95022 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Te-143': {
        'atomic number': 52,
        'mass number': 143,
        'mass': 142.95676 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'I-107': {
        'atomic number': 53,
        'mass number': 107,
        'mass': 106.94678 * u.u,
        'stable': False,
        'half-life': '20# us',
    },

    'I-108': {
        'atomic number': 53,
        'mass number': 108,
        'mass': 107.94348 * u.u,
        'stable': False,
        'half-life': 0.036 * u.s,
    },

    'I-109': {
        'atomic number': 53,
        'mass number': 109,
        'mass': 108.9380853 * u.u,
        'stable': False,
        'half-life': 0.000103 * u.s,
    },

    'I-110': {
        'atomic number': 53,
        'mass number': 110,
        'mass': 109.935089 * u.u,
        'stable': False,
        'half-life': 0.664 * u.s,
    },

    'I-111': {
        'atomic number': 53,
        'mass number': 111,
        'mass': 110.9302692 * u.u,
        'stable': False,
        'half-life': 2.5 * u.s,
    },

    'I-112': {
        'atomic number': 53,
        'mass number': 112,
        'mass': 111.928005 * u.u,
        'stable': False,
        'half-life': 3.34 * u.s,
    },

    'I-113': {
        'atomic number': 53,
        'mass number': 113,
        'mass': 112.9236501 * u.u,
        'stable': False,
        'half-life': 6.6 * u.s,
    },

    'I-114': {
        'atomic number': 53,
        'mass number': 114,
        'mass': 113.92185 * u.u,
        'stable': False,
        'half-life': 2.1 * u.s,
    },

    'I-115': {
        'atomic number': 53,
        'mass number': 115,
        'mass': 114.918048 * u.u,
        'stable': False,
        'half-life': 78.0 * u.s,
    },

    'I-116': {
        'atomic number': 53,
        'mass number': 116,
        'mass': 115.91681 * u.u,
        'stable': False,
        'half-life': 2.91 * u.s,
    },

    'I-117': {
        'atomic number': 53,
        'mass number': 117,
        'mass': 116.913648 * u.u,
        'stable': False,
        'half-life': 133.2 * u.s,
    },

    'I-118': {
        'atomic number': 53,
        'mass number': 118,
        'mass': 117.913074 * u.u,
        'stable': False,
        'half-life': 822.0 * u.s,
    },

    'I-119': {
        'atomic number': 53,
        'mass number': 119,
        'mass': 118.910074 * u.u,
        'stable': False,
        'half-life': 1146.0 * u.s,
    },

    'I-120': {
        'atomic number': 53,
        'mass number': 120,
        'mass': 119.910087 * u.u,
        'stable': False,
        'half-life': 4900.2 * u.s,
    },

    'I-121': {
        'atomic number': 53,
        'mass number': 121,
        'mass': 120.9074051 * u.u,
        'stable': False,
        'half-life': 7632.0 * u.s,
    },

    'I-122': {
        'atomic number': 53,
        'mass number': 122,
        'mass': 121.9075888 * u.u,
        'stable': False,
        'half-life': 217.8 * u.s,
    },

    'I-123': {
        'atomic number': 53,
        'mass number': 123,
        'mass': 122.9055885 * u.u,
        'stable': False,
        'half-life': 47604.6 * u.s,
    },

    'I-124': {
        'atomic number': 53,
        'mass number': 124,
        'mass': 123.906209 * u.u,
        'stable': False,
        'half-life': 360806.4 * u.s,
    },

    'I-125': {
        'atomic number': 53,
        'mass number': 125,
        'mass': 124.9046294 * u.u,
        'stable': False,
        'half-life': 5139936.0 * u.s,
    },

    'I-126': {
        'atomic number': 53,
        'mass number': 126,
        'mass': 125.9056233 * u.u,
        'stable': False,
        'half-life': 1117152.0 * u.s,
    },

    'I-127': {
        'atomic number': 53,
        'mass number': 127,
        'mass': 126.9044719 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'I-128': {
        'atomic number': 53,
        'mass number': 128,
        'mass': 127.9058086 * u.u,
        'stable': False,
        'half-life': 1499.4 * u.s,
    },

    'I-129': {
        'atomic number': 53,
        'mass number': 129,
        'mass': 128.9049837 * u.u,
        'stable': False,
        'half-life': 495443738200000.0 * u.s,
    },

    'I-130': {
        'atomic number': 53,
        'mass number': 130,
        'mass': 129.9066702 * u.u,
        'stable': False,
        'half-life': 44496.0 * u.s,
    },

    'I-131': {
        'atomic number': 53,
        'mass number': 131,
        'mass': 130.9061263 * u.u,
        'stable': False,
        'half-life': 692902.0800000001 * u.s,
    },

    'I-132': {
        'atomic number': 53,
        'mass number': 132,
        'mass': 131.9079935 * u.u,
        'stable': False,
        'half-life': 8262.0 * u.s,
    },

    'I-133': {
        'atomic number': 53,
        'mass number': 133,
        'mass': 132.907797 * u.u,
        'stable': False,
        'half-life': 74988.0 * u.s,
    },

    'I-134': {
        'atomic number': 53,
        'mass number': 134,
        'mass': 133.9097588 * u.u,
        'stable': False,
        'half-life': 3150.0 * u.s,
    },

    'I-135': {
        'atomic number': 53,
        'mass number': 135,
        'mass': 134.9100488 * u.u,
        'stable': False,
        'half-life': 23688.0 * u.s,
    },

    'I-136': {
        'atomic number': 53,
        'mass number': 136,
        'mass': 135.914604 * u.u,
        'stable': False,
        'half-life': 83.4 * u.s,
    },

    'I-137': {
        'atomic number': 53,
        'mass number': 137,
        'mass': 136.9180282 * u.u,
        'stable': False,
        'half-life': 24.13 * u.s,
    },

    'I-138': {
        'atomic number': 53,
        'mass number': 138,
        'mass': 137.9227264 * u.u,
        'stable': False,
        'half-life': 6.23 * u.s,
    },

    'I-139': {
        'atomic number': 53,
        'mass number': 139,
        'mass': 138.926506 * u.u,
        'stable': False,
        'half-life': 2.282 * u.s,
    },

    'I-140': {
        'atomic number': 53,
        'mass number': 140,
        'mass': 139.93173 * u.u,
        'stable': False,
        'half-life': 0.86 * u.s,
    },

    'I-141': {
        'atomic number': 53,
        'mass number': 141,
        'mass': 140.93569 * u.u,
        'stable': False,
        'half-life': 0.43 * u.s,
    },

    'I-142': {
        'atomic number': 53,
        'mass number': 142,
        'mass': 141.9412 * u.u,
        'stable': False,
        'half-life': 0.222 * u.s,
    },

    'I-143': {
        'atomic number': 53,
        'mass number': 143,
        'mass': 142.94565 * u.u,
        'stable': False,
        'half-life': 0.13 * u.s,
    },

    'I-144': {
        'atomic number': 53,
        'mass number': 144,
        'mass': 143.95139 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'I-145': {
        'atomic number': 53,
        'mass number': 145,
        'mass': 144.95605 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Xe-109': {
        'atomic number': 54,
        'mass number': 109,
        'mass': 108.95043 * u.u,
        'stable': False,
        'half-life': 0.013 * u.s,
    },

    'Xe-110': {
        'atomic number': 54,
        'mass number': 110,
        'mass': 109.94426 * u.u,
        'stable': False,
        'half-life': 0.093 * u.s,
    },

    'Xe-111': {
        'atomic number': 54,
        'mass number': 111,
        'mass': 110.941607 * u.u,
        'stable': False,
        'half-life': 0.74 * u.s,
    },

    'Xe-112': {
        'atomic number': 54,
        'mass number': 112,
        'mass': 111.935559 * u.u,
        'stable': False,
        'half-life': 2.7 * u.s,
    },

    'Xe-113': {
        'atomic number': 54,
        'mass number': 113,
        'mass': 112.9332217 * u.u,
        'stable': False,
        'half-life': 2.74 * u.s,
    },

    'Xe-114': {
        'atomic number': 54,
        'mass number': 114,
        'mass': 113.92798 * u.u,
        'stable': False,
        'half-life': 10.0 * u.s,
    },

    'Xe-115': {
        'atomic number': 54,
        'mass number': 115,
        'mass': 114.926294 * u.u,
        'stable': False,
        'half-life': 18.0 * u.s,
    },

    'Xe-116': {
        'atomic number': 54,
        'mass number': 116,
        'mass': 115.921581 * u.u,
        'stable': False,
        'half-life': 59.0 * u.s,
    },

    'Xe-117': {
        'atomic number': 54,
        'mass number': 117,
        'mass': 116.920359 * u.u,
        'stable': False,
        'half-life': 61.0 * u.s,
    },

    'Xe-118': {
        'atomic number': 54,
        'mass number': 118,
        'mass': 117.916179 * u.u,
        'stable': False,
        'half-life': 228.0 * u.s,
    },

    'Xe-119': {
        'atomic number': 54,
        'mass number': 119,
        'mass': 118.915411 * u.u,
        'stable': False,
        'half-life': 348.0 * u.s,
    },

    'Xe-120': {
        'atomic number': 54,
        'mass number': 120,
        'mass': 119.911784 * u.u,
        'stable': False,
        'half-life': 2760.0 * u.s,
    },

    'Xe-121': {
        'atomic number': 54,
        'mass number': 121,
        'mass': 120.911453 * u.u,
        'stable': False,
        'half-life': 2406.0 * u.s,
    },

    'Xe-122': {
        'atomic number': 54,
        'mass number': 122,
        'mass': 121.908368 * u.u,
        'stable': False,
        'half-life': 72360.0 * u.s,
    },

    'Xe-123': {
        'atomic number': 54,
        'mass number': 123,
        'mass': 122.908482 * u.u,
        'stable': False,
        'half-life': 7488.0 * u.s,
    },

    'Xe-124': {
        'atomic number': 54,
        'mass number': 124,
        'mass': 123.905892 * u.u,
        'stable': True,
        'abundance': 0.000952,
    },

    'Xe-125': {
        'atomic number': 54,
        'mass number': 125,
        'mass': 124.9063944 * u.u,
        'stable': False,
        'half-life': 60840.0 * u.s,
    },

    'Xe-126': {
        'atomic number': 54,
        'mass number': 126,
        'mass': 125.9042983 * u.u,
        'stable': True,
        'abundance': 0.00089,
    },

    'Xe-127': {
        'atomic number': 54,
        'mass number': 127,
        'mass': 126.9051829 * u.u,
        'stable': False,
        'half-life': 3140173.44 * u.s,
    },

    'Xe-128': {
        'atomic number': 54,
        'mass number': 128,
        'mass': 127.903531 * u.u,
        'stable': True,
        'abundance': 0.019102,
    },

    'Xe-129': {
        'atomic number': 54,
        'mass number': 129,
        'mass': 128.9047808611 * u.u,
        'stable': True,
        'abundance': 0.264006,
    },

    'Xe-130': {
        'atomic number': 54,
        'mass number': 130,
        'mass': 129.903509349 * u.u,
        'stable': True,
        'abundance': 0.04071,
    },

    'Xe-131': {
        'atomic number': 54,
        'mass number': 131,
        'mass': 130.90508406 * u.u,
        'stable': True,
        'abundance': 0.212324,
    },

    'Xe-132': {
        'atomic number': 54,
        'mass number': 132,
        'mass': 131.9041550856 * u.u,
        'stable': True,
        'abundance': 0.269086,
    },

    'Xe-133': {
        'atomic number': 54,
        'mass number': 133,
        'mass': 132.9059108 * u.u,
        'stable': False,
        'half-life': 453381.408 * u.s,
    },

    'Xe-134': {
        'atomic number': 54,
        'mass number': 134,
        'mass': 133.90539466 * u.u,
        'stable': True,
        'abundance': 0.104357,
    },

    'Xe-135': {
        'atomic number': 54,
        'mass number': 135,
        'mass': 134.9072278 * u.u,
        'stable': False,
        'half-life': 32904.0 * u.s,
    },

    'Xe-136': {
        'atomic number': 54,
        'mass number': 136,
        'mass': 135.907214484 * u.u,
        'stable': False,
        'abundance': 0.088573,
    },

    'Xe-137': {
        'atomic number': 54,
        'mass number': 137,
        'mass': 136.91155778 * u.u,
        'stable': False,
        'half-life': 229.08 * u.s,
    },

    'Xe-138': {
        'atomic number': 54,
        'mass number': 138,
        'mass': 137.9141463 * u.u,
        'stable': False,
        'half-life': 848.4 * u.s,
    },

    'Xe-139': {
        'atomic number': 54,
        'mass number': 139,
        'mass': 138.9187922 * u.u,
        'stable': False,
        'half-life': 39.68 * u.s,
    },

    'Xe-140': {
        'atomic number': 54,
        'mass number': 140,
        'mass': 139.9216458 * u.u,
        'stable': False,
        'half-life': 13.6 * u.s,
    },

    'Xe-141': {
        'atomic number': 54,
        'mass number': 141,
        'mass': 140.9267872 * u.u,
        'stable': False,
        'half-life': 1.73 * u.s,
    },

    'Xe-142': {
        'atomic number': 54,
        'mass number': 142,
        'mass': 141.9299731 * u.u,
        'stable': False,
        'half-life': 1.23 * u.s,
    },

    'Xe-143': {
        'atomic number': 54,
        'mass number': 143,
        'mass': 142.9353696 * u.u,
        'stable': False,
        'half-life': 0.511 * u.s,
    },

    'Xe-144': {
        'atomic number': 54,
        'mass number': 144,
        'mass': 143.9389451 * u.u,
        'stable': False,
        'half-life': 0.388 * u.s,
    },

    'Xe-145': {
        'atomic number': 54,
        'mass number': 145,
        'mass': 144.94472 * u.u,
        'stable': False,
        'half-life': 0.188 * u.s,
    },

    'Xe-146': {
        'atomic number': 54,
        'mass number': 146,
        'mass': 145.948518 * u.u,
        'stable': False,
        'half-life': 0.146 * u.s,
    },

    'Xe-147': {
        'atomic number': 54,
        'mass number': 147,
        'mass': 146.95426 * u.u,
        'stable': False,
        'half-life': 0.13 * u.s,
    },

    'Xe-148': {
        'atomic number': 54,
        'mass number': 148,
        'mass': 147.95813 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Cs-112': {
        'atomic number': 55,
        'mass number': 112,
        'mass': 111.950309 * u.u,
        'stable': False,
        'half-life': 0.00049 * u.s,
    },

    'Cs-113': {
        'atomic number': 55,
        'mass number': 113,
        'mass': 112.9444291 * u.u,
        'stable': False,
        'half-life': 1.77e-05 * u.s,
    },

    'Cs-114': {
        'atomic number': 55,
        'mass number': 114,
        'mass': 113.941296 * u.u,
        'stable': False,
        'half-life': 0.57 * u.s,
    },

    'Cs-115': {
        'atomic number': 55,
        'mass number': 115,
        'mass': 114.93591 * u.u,
        'stable': False,
        'half-life': 1.4 * u.s,
    },

    'Cs-116': {
        'atomic number': 55,
        'mass number': 116,
        'mass': 115.93337 * u.u,
        'stable': False,
        'half-life': 0.7 * u.s,
    },

    'Cs-117': {
        'atomic number': 55,
        'mass number': 117,
        'mass': 116.928617 * u.u,
        'stable': False,
        'half-life': 8.4 * u.s,
    },

    'Cs-118': {
        'atomic number': 55,
        'mass number': 118,
        'mass': 117.92656 * u.u,
        'stable': False,
        'half-life': 14.0 * u.s,
    },

    'Cs-119': {
        'atomic number': 55,
        'mass number': 119,
        'mass': 118.922377 * u.u,
        'stable': False,
        'half-life': 43.0 * u.s,
    },

    'Cs-120': {
        'atomic number': 55,
        'mass number': 120,
        'mass': 119.920677 * u.u,
        'stable': False,
        'half-life': 60.4 * u.s,
    },

    'Cs-121': {
        'atomic number': 55,
        'mass number': 121,
        'mass': 120.917227 * u.u,
        'stable': False,
        'half-life': 155.0 * u.s,
    },

    'Cs-122': {
        'atomic number': 55,
        'mass number': 122,
        'mass': 121.916108 * u.u,
        'stable': False,
        'half-life': 21.18 * u.s,
    },

    'Cs-123': {
        'atomic number': 55,
        'mass number': 123,
        'mass': 122.912996 * u.u,
        'stable': False,
        'half-life': 352.8 * u.s,
    },

    'Cs-124': {
        'atomic number': 55,
        'mass number': 124,
        'mass': 123.9122578 * u.u,
        'stable': False,
        'half-life': 30.9 * u.s,
    },

    'Cs-125': {
        'atomic number': 55,
        'mass number': 125,
        'mass': 124.909728 * u.u,
        'stable': False,
        'half-life': 2802.0 * u.s,
    },

    'Cs-126': {
        'atomic number': 55,
        'mass number': 126,
        'mass': 125.909446 * u.u,
        'stable': False,
        'half-life': 98.4 * u.s,
    },

    'Cs-127': {
        'atomic number': 55,
        'mass number': 127,
        'mass': 126.9074174 * u.u,
        'stable': False,
        'half-life': 22500.0 * u.s,
    },

    'Cs-128': {
        'atomic number': 55,
        'mass number': 128,
        'mass': 127.9077487 * u.u,
        'stable': False,
        'half-life': 218.4 * u.s,
    },

    'Cs-129': {
        'atomic number': 55,
        'mass number': 129,
        'mass': 128.9060657 * u.u,
        'stable': False,
        'half-life': 115416.0 * u.s,
    },

    'Cs-130': {
        'atomic number': 55,
        'mass number': 130,
        'mass': 129.9067093 * u.u,
        'stable': False,
        'half-life': 1752.6 * u.s,
    },

    'Cs-131': {
        'atomic number': 55,
        'mass number': 131,
        'mass': 130.9054649 * u.u,
        'stable': False,
        'half-life': 837129.6 * u.s,
    },

    'Cs-132': {
        'atomic number': 55,
        'mass number': 132,
        'mass': 131.9064339 * u.u,
        'stable': False,
        'half-life': 559872.0 * u.s,
    },

    'Cs-133': {
        'atomic number': 55,
        'mass number': 133,
        'mass': 132.905451961 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Cs-134': {
        'atomic number': 55,
        'mass number': 134,
        'mass': 133.906718503 * u.u,
        'stable': False,
        'half-life': 65135232.0 * u.s,
    },

    'Cs-135': {
        'atomic number': 55,
        'mass number': 135,
        'mass': 134.905977 * u.u,
        'stable': False,
        'half-life': 41970711580000.0 * u.s,
    },

    'Cs-136': {
        'atomic number': 55,
        'mass number': 136,
        'mass': 135.9073114 * u.u,
        'stable': False,
        'half-life': 1137024.0 * u.s,
    },

    'Cs-137': {
        'atomic number': 55,
        'mass number': 137,
        'mass': 136.90708923 * u.u,
        'stable': False,
        'half-life': 951981119.9999999 * u.s,
    },

    'Cs-138': {
        'atomic number': 55,
        'mass number': 138,
        'mass': 137.9110171 * u.u,
        'stable': False,
        'half-life': 2004.6 * u.s,
    },

    'Cs-139': {
        'atomic number': 55,
        'mass number': 139,
        'mass': 138.9133638 * u.u,
        'stable': False,
        'half-life': 556.2 * u.s,
    },

    'Cs-140': {
        'atomic number': 55,
        'mass number': 140,
        'mass': 139.9172831 * u.u,
        'stable': False,
        'half-life': 63.7 * u.s,
    },

    'Cs-141': {
        'atomic number': 55,
        'mass number': 141,
        'mass': 140.9200455 * u.u,
        'stable': False,
        'half-life': 24.84 * u.s,
    },

    'Cs-142': {
        'atomic number': 55,
        'mass number': 142,
        'mass': 141.924296 * u.u,
        'stable': False,
        'half-life': 1.684 * u.s,
    },

    'Cs-143': {
        'atomic number': 55,
        'mass number': 143,
        'mass': 142.927349 * u.u,
        'stable': False,
        'half-life': 1.791 * u.s,
    },

    'Cs-144': {
        'atomic number': 55,
        'mass number': 144,
        'mass': 143.932076 * u.u,
        'stable': False,
        'half-life': 0.994 * u.s,
    },

    'Cs-145': {
        'atomic number': 55,
        'mass number': 145,
        'mass': 144.935527 * u.u,
        'stable': False,
        'half-life': 0.582 * u.s,
    },

    'Cs-146': {
        'atomic number': 55,
        'mass number': 146,
        'mass': 145.940344 * u.u,
        'stable': False,
        'half-life': 0.323 * u.s,
    },

    'Cs-147': {
        'atomic number': 55,
        'mass number': 147,
        'mass': 146.944156 * u.u,
        'stable': False,
        'half-life': 0.23 * u.s,
    },

    'Cs-148': {
        'atomic number': 55,
        'mass number': 148,
        'mass': 147.94923 * u.u,
        'stable': False,
        'half-life': 0.145 * u.s,
    },

    'Cs-149': {
        'atomic number': 55,
        'mass number': 149,
        'mass': 148.95302 * u.u,
        'stable': False,
        'half-life': 0.113 * u.s,
    },

    'Cs-150': {
        'atomic number': 55,
        'mass number': 150,
        'mass': 149.95833 * u.u,
        'stable': False,
        'half-life': 0.0844 * u.s,
    },

    'Cs-151': {
        'atomic number': 55,
        'mass number': 151,
        'mass': 150.96258 * u.u,
        'stable': False,
        'half-life': 0.069 * u.s,
    },

    'Ba-114': {
        'atomic number': 56,
        'mass number': 114,
        'mass': 113.95066 * u.u,
        'stable': False,
        'half-life': 0.46 * u.s,
    },

    'Ba-115': {
        'atomic number': 56,
        'mass number': 115,
        'mass': 114.94737 * u.u,
        'stable': False,
        'half-life': 0.45 * u.s,
    },

    'Ba-116': {
        'atomic number': 56,
        'mass number': 116,
        'mass': 115.94128 * u.u,
        'stable': False,
        'half-life': 1.3 * u.s,
    },

    'Ba-117': {
        'atomic number': 56,
        'mass number': 117,
        'mass': 116.93814 * u.u,
        'stable': False,
        'half-life': 1.75 * u.s,
    },

    'Ba-118': {
        'atomic number': 56,
        'mass number': 118,
        'mass': 117.93306 * u.u,
        'stable': False,
        'half-life': 5.2 * u.s,
    },

    'Ba-119': {
        'atomic number': 56,
        'mass number': 119,
        'mass': 118.93066 * u.u,
        'stable': False,
        'half-life': 5.4 * u.s,
    },

    'Ba-120': {
        'atomic number': 56,
        'mass number': 120,
        'mass': 119.92605 * u.u,
        'stable': False,
        'half-life': 24.0 * u.s,
    },

    'Ba-121': {
        'atomic number': 56,
        'mass number': 121,
        'mass': 120.92405 * u.u,
        'stable': False,
        'half-life': 29.7 * u.s,
    },

    'Ba-122': {
        'atomic number': 56,
        'mass number': 122,
        'mass': 121.919904 * u.u,
        'stable': False,
        'half-life': 117.0 * u.s,
    },

    'Ba-123': {
        'atomic number': 56,
        'mass number': 123,
        'mass': 122.918781 * u.u,
        'stable': False,
        'half-life': 162.0 * u.s,
    },

    'Ba-124': {
        'atomic number': 56,
        'mass number': 124,
        'mass': 123.915094 * u.u,
        'stable': False,
        'half-life': 660.0 * u.s,
    },

    'Ba-125': {
        'atomic number': 56,
        'mass number': 125,
        'mass': 124.914472 * u.u,
        'stable': False,
        'half-life': 198.0 * u.s,
    },

    'Ba-126': {
        'atomic number': 56,
        'mass number': 126,
        'mass': 125.91125 * u.u,
        'stable': False,
        'half-life': 6000.0 * u.s,
    },

    'Ba-127': {
        'atomic number': 56,
        'mass number': 127,
        'mass': 126.911091 * u.u,
        'stable': False,
        'half-life': 762.0 * u.s,
    },

    'Ba-128': {
        'atomic number': 56,
        'mass number': 128,
        'mass': 127.908342 * u.u,
        'stable': False,
        'half-life': 209952.0 * u.s,
    },

    'Ba-129': {
        'atomic number': 56,
        'mass number': 129,
        'mass': 128.908681 * u.u,
        'stable': False,
        'half-life': 8028.0 * u.s,
    },

    'Ba-130': {
        'atomic number': 56,
        'mass number': 130,
        'mass': 129.9063207 * u.u,
        'stable': False,
        'abundance': 0.00106,
    },

    'Ba-131': {
        'atomic number': 56,
        'mass number': 131,
        'mass': 130.906941 * u.u,
        'stable': False,
        'half-life': 995328.0 * u.s,
    },

    'Ba-132': {
        'atomic number': 56,
        'mass number': 132,
        'mass': 131.9050611 * u.u,
        'stable': True,
        'abundance': 0.00101,
    },

    'Ba-133': {
        'atomic number': 56,
        'mass number': 133,
        'mass': 132.9060074 * u.u,
        'stable': False,
        'half-life': 333046080.0 * u.s,
    },

    'Ba-134': {
        'atomic number': 56,
        'mass number': 134,
        'mass': 133.90450818 * u.u,
        'stable': True,
        'abundance': 0.02417,
    },

    'Ba-135': {
        'atomic number': 56,
        'mass number': 135,
        'mass': 134.90568838 * u.u,
        'stable': True,
        'abundance': 0.06592,
    },

    'Ba-136': {
        'atomic number': 56,
        'mass number': 136,
        'mass': 135.90457573 * u.u,
        'stable': True,
        'abundance': 0.07854,
    },

    'Ba-137': {
        'atomic number': 56,
        'mass number': 137,
        'mass': 136.90582714 * u.u,
        'stable': True,
        'abundance': 0.11232,
    },

    'Ba-138': {
        'atomic number': 56,
        'mass number': 138,
        'mass': 137.905247 * u.u,
        'stable': True,
        'abundance': 0.71698,
    },

    'Ba-139': {
        'atomic number': 56,
        'mass number': 139,
        'mass': 138.9088411 * u.u,
        'stable': False,
        'half-life': 4987.8 * u.s,
    },

    'Ba-140': {
        'atomic number': 56,
        'mass number': 140,
        'mass': 139.9106057 * u.u,
        'stable': False,
        'half-life': 1101833.28 * u.s,
    },

    'Ba-141': {
        'atomic number': 56,
        'mass number': 141,
        'mass': 140.9144033 * u.u,
        'stable': False,
        'half-life': 1096.2 * u.s,
    },

    'Ba-142': {
        'atomic number': 56,
        'mass number': 142,
        'mass': 141.9164324 * u.u,
        'stable': False,
        'half-life': 636.0 * u.s,
    },

    'Ba-143': {
        'atomic number': 56,
        'mass number': 143,
        'mass': 142.9206253 * u.u,
        'stable': False,
        'half-life': 14.5 * u.s,
    },

    'Ba-144': {
        'atomic number': 56,
        'mass number': 144,
        'mass': 143.9229549 * u.u,
        'stable': False,
        'half-life': 11.5 * u.s,
    },

    'Ba-145': {
        'atomic number': 56,
        'mass number': 145,
        'mass': 144.9275184 * u.u,
        'stable': False,
        'half-life': 4.31 * u.s,
    },

    'Ba-146': {
        'atomic number': 56,
        'mass number': 146,
        'mass': 145.930284 * u.u,
        'stable': False,
        'half-life': 2.22 * u.s,
    },

    'Ba-147': {
        'atomic number': 56,
        'mass number': 147,
        'mass': 146.935304 * u.u,
        'stable': False,
        'half-life': 0.894 * u.s,
    },

    'Ba-148': {
        'atomic number': 56,
        'mass number': 148,
        'mass': 147.938171 * u.u,
        'stable': False,
        'half-life': 0.62 * u.s,
    },

    'Ba-149': {
        'atomic number': 56,
        'mass number': 149,
        'mass': 148.94308 * u.u,
        'stable': False,
        'half-life': 0.348 * u.s,
    },

    'Ba-150': {
        'atomic number': 56,
        'mass number': 150,
        'mass': 149.94605 * u.u,
        'stable': False,
        'half-life': 0.259 * u.s,
    },

    'Ba-151': {
        'atomic number': 56,
        'mass number': 151,
        'mass': 150.95127 * u.u,
        'stable': False,
        'half-life': 0.167 * u.s,
    },

    'Ba-152': {
        'atomic number': 56,
        'mass number': 152,
        'mass': 151.95481 * u.u,
        'stable': False,
        'half-life': 0.139 * u.s,
    },

    'Ba-153': {
        'atomic number': 56,
        'mass number': 153,
        'mass': 152.96036 * u.u,
        'stable': False,
        'half-life': 0.116 * u.s,
    },

    'La-116': {
        'atomic number': 57,
        'mass number': 116,
        'mass': 115.9563 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'La-117': {
        'atomic number': 57,
        'mass number': 117,
        'mass': 116.94999 * u.u,
        'stable': False,
        'half-life': 0.0217 * u.s,
    },

    'La-118': {
        'atomic number': 57,
        'mass number': 118,
        'mass': 117.94673 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'La-119': {
        'atomic number': 57,
        'mass number': 119,
        'mass': 118.94099 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'La-120': {
        'atomic number': 57,
        'mass number': 120,
        'mass': 119.93807 * u.u,
        'stable': False,
        'half-life': 2.8 * u.s,
    },

    'La-121': {
        'atomic number': 57,
        'mass number': 121,
        'mass': 120.93315 * u.u,
        'stable': False,
        'half-life': 5.3 * u.s,
    },

    'La-122': {
        'atomic number': 57,
        'mass number': 122,
        'mass': 121.93071 * u.u,
        'stable': False,
        'half-life': 8.6 * u.s,
    },

    'La-123': {
        'atomic number': 57,
        'mass number': 123,
        'mass': 122.9263 * u.u,
        'stable': False,
        'half-life': 17.0 * u.s,
    },

    'La-124': {
        'atomic number': 57,
        'mass number': 124,
        'mass': 123.924574 * u.u,
        'stable': False,
        'half-life': 29.21 * u.s,
    },

    'La-125': {
        'atomic number': 57,
        'mass number': 125,
        'mass': 124.920816 * u.u,
        'stable': False,
        'half-life': 64.8 * u.s,
    },

    'La-126': {
        'atomic number': 57,
        'mass number': 126,
        'mass': 125.919513 * u.u,
        'stable': False,
        'half-life': 54.0 * u.s,
    },

    'La-127': {
        'atomic number': 57,
        'mass number': 127,
        'mass': 126.916375 * u.u,
        'stable': False,
        'half-life': 306.0 * u.s,
    },

    'La-128': {
        'atomic number': 57,
        'mass number': 128,
        'mass': 127.915592 * u.u,
        'stable': False,
        'half-life': 310.8 * u.s,
    },

    'La-129': {
        'atomic number': 57,
        'mass number': 129,
        'mass': 128.912694 * u.u,
        'stable': False,
        'half-life': 696.0 * u.s,
    },

    'La-130': {
        'atomic number': 57,
        'mass number': 130,
        'mass': 129.912369 * u.u,
        'stable': False,
        'half-life': 522.0 * u.s,
    },

    'La-131': {
        'atomic number': 57,
        'mass number': 131,
        'mass': 130.91007 * u.u,
        'stable': False,
        'half-life': 3540.0 * u.s,
    },

    'La-132': {
        'atomic number': 57,
        'mass number': 132,
        'mass': 131.910119 * u.u,
        'stable': False,
        'half-life': 17280.0 * u.s,
    },

    'La-133': {
        'atomic number': 57,
        'mass number': 133,
        'mass': 132.908218 * u.u,
        'stable': False,
        'half-life': 14083.2 * u.s,
    },

    'La-134': {
        'atomic number': 57,
        'mass number': 134,
        'mass': 133.908514 * u.u,
        'stable': False,
        'half-life': 387.0 * u.s,
    },

    'La-135': {
        'atomic number': 57,
        'mass number': 135,
        'mass': 134.906984 * u.u,
        'stable': False,
        'half-life': 70200.0 * u.s,
    },

    'La-136': {
        'atomic number': 57,
        'mass number': 136,
        'mass': 135.907635 * u.u,
        'stable': False,
        'half-life': 592.2 * u.s,
    },

    'La-137': {
        'atomic number': 57,
        'mass number': 137,
        'mass': 136.9064504 * u.u,
        'stable': False,
        'half-life': 1893415560000.0 * u.s,
    },

    'La-138': {
        'atomic number': 57,
        'mass number': 138,
        'mass': 137.9071149 * u.u,
        'stable': False,
        'abundance': 0.0008881,
    },

    'La-139': {
        'atomic number': 57,
        'mass number': 139,
        'mass': 138.9063563 * u.u,
        'stable': True,
        'abundance': 0.9991119,
    },

    'La-140': {
        'atomic number': 57,
        'mass number': 140,
        'mass': 139.9094806 * u.u,
        'stable': False,
        'half-life': 145054.8 * u.s,
    },

    'La-141': {
        'atomic number': 57,
        'mass number': 141,
        'mass': 140.910966 * u.u,
        'stable': False,
        'half-life': 14112.0 * u.s,
    },

    'La-142': {
        'atomic number': 57,
        'mass number': 142,
        'mass': 141.9140909 * u.u,
        'stable': False,
        'half-life': 5466.0 * u.s,
    },

    'La-143': {
        'atomic number': 57,
        'mass number': 143,
        'mass': 142.9160795 * u.u,
        'stable': False,
        'half-life': 852.0 * u.s,
    },

    'La-144': {
        'atomic number': 57,
        'mass number': 144,
        'mass': 143.919646 * u.u,
        'stable': False,
        'half-life': 40.8 * u.s,
    },

    'La-145': {
        'atomic number': 57,
        'mass number': 145,
        'mass': 144.921808 * u.u,
        'stable': False,
        'half-life': 24.8 * u.s,
    },

    'La-146': {
        'atomic number': 57,
        'mass number': 146,
        'mass': 145.925875 * u.u,
        'stable': False,
        'half-life': 6.27 * u.s,
    },

    'La-147': {
        'atomic number': 57,
        'mass number': 147,
        'mass': 146.928418 * u.u,
        'stable': False,
        'half-life': 4.06 * u.s,
    },

    'La-148': {
        'atomic number': 57,
        'mass number': 148,
        'mass': 147.932679 * u.u,
        'stable': False,
        'half-life': 1.35 * u.s,
    },

    'La-149': {
        'atomic number': 57,
        'mass number': 149,
        'mass': 148.93535 * u.u,
        'stable': False,
        'half-life': 1.07 * u.s,
    },

    'La-150': {
        'atomic number': 57,
        'mass number': 150,
        'mass': 149.93947 * u.u,
        'stable': False,
        'half-life': 0.504 * u.s,
    },

    'La-151': {
        'atomic number': 57,
        'mass number': 151,
        'mass': 150.94232 * u.u,
        'stable': False,
        'half-life': 0.465 * u.s,
    },

    'La-152': {
        'atomic number': 57,
        'mass number': 152,
        'mass': 151.94682 * u.u,
        'stable': False,
        'half-life': 0.287 * u.s,
    },

    'La-153': {
        'atomic number': 57,
        'mass number': 153,
        'mass': 152.95036 * u.u,
        'stable': False,
        'half-life': 0.245 * u.s,
    },

    'La-154': {
        'atomic number': 57,
        'mass number': 154,
        'mass': 153.95517 * u.u,
        'stable': False,
        'half-life': 0.161 * u.s,
    },

    'La-155': {
        'atomic number': 57,
        'mass number': 155,
        'mass': 154.95901 * u.u,
        'stable': False,
        'half-life': 0.101 * u.s,
    },

    'Ce-119': {
        'atomic number': 58,
        'mass number': 119,
        'mass': 118.95271 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'Ce-120': {
        'atomic number': 58,
        'mass number': 120,
        'mass': 119.94654 * u.u,
        'stable': False,
        'half-life': '250# ms',
    },

    'Ce-121': {
        'atomic number': 58,
        'mass number': 121,
        'mass': 120.94335 * u.u,
        'stable': False,
        'half-life': 1.1 * u.s,
    },

    'Ce-122': {
        'atomic number': 58,
        'mass number': 122,
        'mass': 121.93787 * u.u,
        'stable': False,
        'half-life': '2# s',
    },

    'Ce-123': {
        'atomic number': 58,
        'mass number': 123,
        'mass': 122.93528 * u.u,
        'stable': False,
        'half-life': 3.8 * u.s,
    },

    'Ce-124': {
        'atomic number': 58,
        'mass number': 124,
        'mass': 123.93031 * u.u,
        'stable': False,
        'half-life': 9.1 * u.s,
    },

    'Ce-125': {
        'atomic number': 58,
        'mass number': 125,
        'mass': 124.92844 * u.u,
        'stable': False,
        'half-life': 9.7 * u.s,
    },

    'Ce-126': {
        'atomic number': 58,
        'mass number': 126,
        'mass': 125.923971 * u.u,
        'stable': False,
        'half-life': 51.0 * u.s,
    },

    'Ce-127': {
        'atomic number': 58,
        'mass number': 127,
        'mass': 126.922727 * u.u,
        'stable': False,
        'half-life': 34.0 * u.s,
    },

    'Ce-128': {
        'atomic number': 58,
        'mass number': 128,
        'mass': 127.918911 * u.u,
        'stable': False,
        'half-life': 235.8 * u.s,
    },

    'Ce-129': {
        'atomic number': 58,
        'mass number': 129,
        'mass': 128.918102 * u.u,
        'stable': False,
        'half-life': 210.0 * u.s,
    },

    'Ce-130': {
        'atomic number': 58,
        'mass number': 130,
        'mass': 129.914736 * u.u,
        'stable': False,
        'half-life': 1374.0 * u.s,
    },

    'Ce-131': {
        'atomic number': 58,
        'mass number': 131,
        'mass': 130.914429 * u.u,
        'stable': False,
        'half-life': 618.0 * u.s,
    },

    'Ce-132': {
        'atomic number': 58,
        'mass number': 132,
        'mass': 131.911464 * u.u,
        'stable': False,
        'half-life': 12636.0 * u.s,
    },

    'Ce-133': {
        'atomic number': 58,
        'mass number': 133,
        'mass': 132.91152 * u.u,
        'stable': False,
        'half-life': 5820.0 * u.s,
    },

    'Ce-134': {
        'atomic number': 58,
        'mass number': 134,
        'mass': 133.908928 * u.u,
        'stable': False,
        'half-life': 273024.0 * u.s,
    },

    'Ce-135': {
        'atomic number': 58,
        'mass number': 135,
        'mass': 134.909161 * u.u,
        'stable': False,
        'half-life': 63720.0 * u.s,
    },

    'Ce-136': {
        'atomic number': 58,
        'mass number': 136,
        'mass': 135.90712921 * u.u,
        'stable': True,
        'abundance': 0.00185,
    },

    'Ce-137': {
        'atomic number': 58,
        'mass number': 137,
        'mass': 136.90776236 * u.u,
        'stable': False,
        'half-life': 32400.0 * u.s,
    },

    'Ce-138': {
        'atomic number': 58,
        'mass number': 138,
        'mass': 137.905991 * u.u,
        'stable': True,
        'abundance': 0.00251,
    },

    'Ce-139': {
        'atomic number': 58,
        'mass number': 139,
        'mass': 138.9066551 * u.u,
        'stable': False,
        'half-life': 11900217.600000001 * u.s,
    },

    'Ce-140': {
        'atomic number': 58,
        'mass number': 140,
        'mass': 139.9054431 * u.u,
        'stable': True,
        'abundance': 0.8845,
    },

    'Ce-141': {
        'atomic number': 58,
        'mass number': 141,
        'mass': 140.9082807 * u.u,
        'stable': False,
        'half-life': 2808864.0 * u.s,
    },

    'Ce-142': {
        'atomic number': 58,
        'mass number': 142,
        'mass': 141.9092504 * u.u,
        'stable': True,
        'abundance': 0.11114,
    },

    'Ce-143': {
        'atomic number': 58,
        'mass number': 143,
        'mass': 142.9123921 * u.u,
        'stable': False,
        'half-life': 118940.4 * u.s,
    },

    'Ce-144': {
        'atomic number': 58,
        'mass number': 144,
        'mass': 143.9136529 * u.u,
        'stable': False,
        'half-life': 24583737.599999998 * u.s,
    },

    'Ce-145': {
        'atomic number': 58,
        'mass number': 145,
        'mass': 144.917265 * u.u,
        'stable': False,
        'half-life': 180.6 * u.s,
    },

    'Ce-146': {
        'atomic number': 58,
        'mass number': 146,
        'mass': 145.918802 * u.u,
        'stable': False,
        'half-life': 811.2 * u.s,
    },

    'Ce-147': {
        'atomic number': 58,
        'mass number': 147,
        'mass': 146.9226899 * u.u,
        'stable': False,
        'half-life': 56.4 * u.s,
    },

    'Ce-148': {
        'atomic number': 58,
        'mass number': 148,
        'mass': 147.924424 * u.u,
        'stable': False,
        'half-life': 56.8 * u.s,
    },

    'Ce-149': {
        'atomic number': 58,
        'mass number': 149,
        'mass': 148.928427 * u.u,
        'stable': False,
        'half-life': 4.94 * u.s,
    },

    'Ce-150': {
        'atomic number': 58,
        'mass number': 150,
        'mass': 149.930384 * u.u,
        'stable': False,
        'half-life': 6.05 * u.s,
    },

    'Ce-151': {
        'atomic number': 58,
        'mass number': 151,
        'mass': 150.934272 * u.u,
        'stable': False,
        'half-life': 1.76 * u.s,
    },

    'Ce-152': {
        'atomic number': 58,
        'mass number': 152,
        'mass': 151.9366 * u.u,
        'stable': False,
        'half-life': 1.42 * u.s,
    },

    'Ce-153': {
        'atomic number': 58,
        'mass number': 153,
        'mass': 152.94093 * u.u,
        'stable': False,
        'half-life': 0.865 * u.s,
    },

    'Ce-154': {
        'atomic number': 58,
        'mass number': 154,
        'mass': 153.9438 * u.u,
        'stable': False,
        'half-life': 0.722 * u.s,
    },

    'Ce-155': {
        'atomic number': 58,
        'mass number': 155,
        'mass': 154.94855 * u.u,
        'stable': False,
        'half-life': 0.313 * u.s,
    },

    'Ce-156': {
        'atomic number': 58,
        'mass number': 156,
        'mass': 155.95183 * u.u,
        'stable': False,
        'half-life': 0.233 * u.s,
    },

    'Ce-157': {
        'atomic number': 58,
        'mass number': 157,
        'mass': 156.95705 * u.u,
        'stable': False,
        'half-life': 0.175 * u.s,
    },

    'Pr-121': {
        'atomic number': 59,
        'mass number': 121,
        'mass': 120.95532 * u.u,
        'stable': False,
        'half-life': 0.012 * u.s,
    },

    'Pr-122': {
        'atomic number': 59,
        'mass number': 122,
        'mass': 121.95175 * u.u,
        'stable': False,
        'half-life': '500# ms',
    },

    'Pr-123': {
        'atomic number': 59,
        'mass number': 123,
        'mass': 122.94596 * u.u,
        'stable': False,
        'half-life': '800# ms',
    },

    'Pr-124': {
        'atomic number': 59,
        'mass number': 124,
        'mass': 123.94294 * u.u,
        'stable': False,
        'half-life': 1.2 * u.s,
    },

    'Pr-125': {
        'atomic number': 59,
        'mass number': 125,
        'mass': 124.9377 * u.u,
        'stable': False,
        'half-life': 3.3 * u.s,
    },

    'Pr-126': {
        'atomic number': 59,
        'mass number': 126,
        'mass': 125.93524 * u.u,
        'stable': False,
        'half-life': 3.12 * u.s,
    },

    'Pr-127': {
        'atomic number': 59,
        'mass number': 127,
        'mass': 126.93071 * u.u,
        'stable': False,
        'half-life': 4.2 * u.s,
    },

    'Pr-128': {
        'atomic number': 59,
        'mass number': 128,
        'mass': 127.928791 * u.u,
        'stable': False,
        'half-life': 2.85 * u.s,
    },

    'Pr-129': {
        'atomic number': 59,
        'mass number': 129,
        'mass': 128.925095 * u.u,
        'stable': False,
        'half-life': 30.0 * u.s,
    },

    'Pr-130': {
        'atomic number': 59,
        'mass number': 130,
        'mass': 129.92359 * u.u,
        'stable': False,
        'half-life': 40.0 * u.s,
    },

    'Pr-131': {
        'atomic number': 59,
        'mass number': 131,
        'mass': 130.920235 * u.u,
        'stable': False,
        'half-life': 90.0 * u.s,
    },

    'Pr-132': {
        'atomic number': 59,
        'mass number': 132,
        'mass': 131.919255 * u.u,
        'stable': False,
        'half-life': 89.4 * u.s,
    },

    'Pr-133': {
        'atomic number': 59,
        'mass number': 133,
        'mass': 132.916331 * u.u,
        'stable': False,
        'half-life': 390.0 * u.s,
    },

    'Pr-134': {
        'atomic number': 59,
        'mass number': 134,
        'mass': 133.915697 * u.u,
        'stable': False,
        'half-life': 1020.0 * u.s,
    },

    'Pr-135': {
        'atomic number': 59,
        'mass number': 135,
        'mass': 134.913112 * u.u,
        'stable': False,
        'half-life': 1440.0 * u.s,
    },

    'Pr-136': {
        'atomic number': 59,
        'mass number': 136,
        'mass': 135.912677 * u.u,
        'stable': False,
        'half-life': 786.0 * u.s,
    },

    'Pr-137': {
        'atomic number': 59,
        'mass number': 137,
        'mass': 136.9106792 * u.u,
        'stable': False,
        'half-life': 4608.0 * u.s,
    },

    'Pr-138': {
        'atomic number': 59,
        'mass number': 138,
        'mass': 137.910754 * u.u,
        'stable': False,
        'half-life': 87.0 * u.s,
    },

    'Pr-139': {
        'atomic number': 59,
        'mass number': 139,
        'mass': 138.9089408 * u.u,
        'stable': False,
        'half-life': 15876.0 * u.s,
    },

    'Pr-140': {
        'atomic number': 59,
        'mass number': 140,
        'mass': 139.9090803 * u.u,
        'stable': False,
        'half-life': 203.4 * u.s,
    },

    'Pr-141': {
        'atomic number': 59,
        'mass number': 141,
        'mass': 140.9076576 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Pr-142': {
        'atomic number': 59,
        'mass number': 142,
        'mass': 141.9100496 * u.u,
        'stable': False,
        'half-life': 68832.0 * u.s,
    },

    'Pr-143': {
        'atomic number': 59,
        'mass number': 143,
        'mass': 142.9108228 * u.u,
        'stable': False,
        'half-life': 1172448.0 * u.s,
    },

    'Pr-144': {
        'atomic number': 59,
        'mass number': 144,
        'mass': 143.9133109 * u.u,
        'stable': False,
        'half-life': 1036.8 * u.s,
    },

    'Pr-145': {
        'atomic number': 59,
        'mass number': 145,
        'mass': 144.9145182 * u.u,
        'stable': False,
        'half-life': 21542.4 * u.s,
    },

    'Pr-146': {
        'atomic number': 59,
        'mass number': 146,
        'mass': 145.91768 * u.u,
        'stable': False,
        'half-life': 1449.0 * u.s,
    },

    'Pr-147': {
        'atomic number': 59,
        'mass number': 147,
        'mass': 146.919008 * u.u,
        'stable': False,
        'half-life': 804.0 * u.s,
    },

    'Pr-148': {
        'atomic number': 59,
        'mass number': 148,
        'mass': 147.92213 * u.u,
        'stable': False,
        'half-life': 137.4 * u.s,
    },

    'Pr-149': {
        'atomic number': 59,
        'mass number': 149,
        'mass': 148.923736 * u.u,
        'stable': False,
        'half-life': 135.6 * u.s,
    },

    'Pr-150': {
        'atomic number': 59,
        'mass number': 150,
        'mass': 149.9266765 * u.u,
        'stable': False,
        'half-life': 6.19 * u.s,
    },

    'Pr-151': {
        'atomic number': 59,
        'mass number': 151,
        'mass': 150.928309 * u.u,
        'stable': False,
        'half-life': 18.9 * u.s,
    },

    'Pr-152': {
        'atomic number': 59,
        'mass number': 152,
        'mass': 151.931553 * u.u,
        'stable': False,
        'half-life': 3.57 * u.s,
    },

    'Pr-153': {
        'atomic number': 59,
        'mass number': 153,
        'mass': 152.933904 * u.u,
        'stable': False,
        'half-life': 4.28 * u.s,
    },

    'Pr-154': {
        'atomic number': 59,
        'mass number': 154,
        'mass': 153.93753 * u.u,
        'stable': False,
        'half-life': 2.3 * u.s,
    },

    'Pr-155': {
        'atomic number': 59,
        'mass number': 155,
        'mass': 154.940509 * u.u,
        'stable': False,
        'half-life': 1.47 * u.s,
    },

    'Pr-156': {
        'atomic number': 59,
        'mass number': 156,
        'mass': 155.94464 * u.u,
        'stable': False,
        'half-life': 0.444 * u.s,
    },

    'Pr-157': {
        'atomic number': 59,
        'mass number': 157,
        'mass': 156.94789 * u.u,
        'stable': False,
        'half-life': 0.307 * u.s,
    },

    'Pr-158': {
        'atomic number': 59,
        'mass number': 158,
        'mass': 157.95241 * u.u,
        'stable': False,
        'half-life': 0.181 * u.s,
    },

    'Pr-159': {
        'atomic number': 59,
        'mass number': 159,
        'mass': 158.95589 * u.u,
        'stable': False,
        'half-life': 0.134 * u.s,
    },

    'Nd-124': {
        'atomic number': 60,
        'mass number': 124,
        'mass': 123.9522 * u.u,
        'stable': False,
        'half-life': '500# ms',
    },

    'Nd-125': {
        'atomic number': 60,
        'mass number': 125,
        'mass': 124.9489 * u.u,
        'stable': False,
        'half-life': 0.65 * u.s,
    },

    'Nd-126': {
        'atomic number': 60,
        'mass number': 126,
        'mass': 125.94311 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Nd-127': {
        'atomic number': 60,
        'mass number': 127,
        'mass': 126.94038 * u.u,
        'stable': False,
        'half-life': 1.8 * u.s,
    },

    'Nd-128': {
        'atomic number': 60,
        'mass number': 128,
        'mass': 127.93525 * u.u,
        'stable': False,
        'half-life': '5# s',
    },

    'Nd-129': {
        'atomic number': 60,
        'mass number': 129,
        'mass': 128.9331 * u.u,
        'stable': False,
        'half-life': 6.8 * u.s,
    },

    'Nd-130': {
        'atomic number': 60,
        'mass number': 130,
        'mass': 129.928506 * u.u,
        'stable': False,
        'half-life': 21.0 * u.s,
    },

    'Nd-131': {
        'atomic number': 60,
        'mass number': 131,
        'mass': 130.927248 * u.u,
        'stable': False,
        'half-life': 25.4 * u.s,
    },

    'Nd-132': {
        'atomic number': 60,
        'mass number': 132,
        'mass': 131.923321 * u.u,
        'stable': False,
        'half-life': 93.6 * u.s,
    },

    'Nd-133': {
        'atomic number': 60,
        'mass number': 133,
        'mass': 132.922348 * u.u,
        'stable': False,
        'half-life': 70.0 * u.s,
    },

    'Nd-134': {
        'atomic number': 60,
        'mass number': 134,
        'mass': 133.91879 * u.u,
        'stable': False,
        'half-life': 510.0 * u.s,
    },

    'Nd-135': {
        'atomic number': 60,
        'mass number': 135,
        'mass': 134.918181 * u.u,
        'stable': False,
        'half-life': 744.0 * u.s,
    },

    'Nd-136': {
        'atomic number': 60,
        'mass number': 136,
        'mass': 135.914976 * u.u,
        'stable': False,
        'half-life': 3042.0 * u.s,
    },

    'Nd-137': {
        'atomic number': 60,
        'mass number': 137,
        'mass': 136.914562 * u.u,
        'stable': False,
        'half-life': 2310.0 * u.s,
    },

    'Nd-138': {
        'atomic number': 60,
        'mass number': 138,
        'mass': 137.91195 * u.u,
        'stable': False,
        'half-life': 18144.0 * u.s,
    },

    'Nd-139': {
        'atomic number': 60,
        'mass number': 139,
        'mass': 138.911954 * u.u,
        'stable': False,
        'half-life': 1782.0 * u.s,
    },

    'Nd-140': {
        'atomic number': 60,
        'mass number': 140,
        'mass': 139.90955 * u.u,
        'stable': False,
        'half-life': 291168.0 * u.s,
    },

    'Nd-141': {
        'atomic number': 60,
        'mass number': 141,
        'mass': 140.9096147 * u.u,
        'stable': False,
        'half-life': 8964.0 * u.s,
    },

    'Nd-142': {
        'atomic number': 60,
        'mass number': 142,
        'mass': 141.907729 * u.u,
        'stable': True,
        'abundance': 0.27152,
    },

    'Nd-143': {
        'atomic number': 60,
        'mass number': 143,
        'mass': 142.90982 * u.u,
        'stable': True,
        'abundance': 0.12174,
    },

    'Nd-144': {
        'atomic number': 60,
        'mass number': 144,
        'mass': 143.910093 * u.u,
        'stable': False,
        'abundance': 0.23798,
    },

    'Nd-145': {
        'atomic number': 60,
        'mass number': 145,
        'mass': 144.9125793 * u.u,
        'stable': True,
        'abundance': 0.08293,
    },

    'Nd-146': {
        'atomic number': 60,
        'mass number': 146,
        'mass': 145.9131226 * u.u,
        'stable': True,
        'abundance': 0.17189,
    },

    'Nd-147': {
        'atomic number': 60,
        'mass number': 147,
        'mass': 146.9161061 * u.u,
        'stable': False,
        'half-life': 948672.0 * u.s,
    },

    'Nd-148': {
        'atomic number': 60,
        'mass number': 148,
        'mass': 147.9168993 * u.u,
        'stable': True,
        'abundance': 0.05756,
    },

    'Nd-149': {
        'atomic number': 60,
        'mass number': 149,
        'mass': 148.9201548 * u.u,
        'stable': False,
        'half-life': 6220.8 * u.s,
    },

    'Nd-150': {
        'atomic number': 60,
        'mass number': 150,
        'mass': 149.9209022 * u.u,
        'stable': False,
        'abundance': 0.05638,
    },

    'Nd-151': {
        'atomic number': 60,
        'mass number': 151,
        'mass': 150.9238403 * u.u,
        'stable': False,
        'half-life': 746.4 * u.s,
    },

    'Nd-152': {
        'atomic number': 60,
        'mass number': 152,
        'mass': 151.924692 * u.u,
        'stable': False,
        'half-life': 684.0 * u.s,
    },

    'Nd-153': {
        'atomic number': 60,
        'mass number': 153,
        'mass': 152.927718 * u.u,
        'stable': False,
        'half-life': 31.6 * u.s,
    },

    'Nd-154': {
        'atomic number': 60,
        'mass number': 154,
        'mass': 153.92948 * u.u,
        'stable': False,
        'half-life': 25.9 * u.s,
    },

    'Nd-155': {
        'atomic number': 60,
        'mass number': 155,
        'mass': 154.9331357 * u.u,
        'stable': False,
        'half-life': 8.9 * u.s,
    },

    'Nd-156': {
        'atomic number': 60,
        'mass number': 156,
        'mass': 155.93508 * u.u,
        'stable': False,
        'half-life': 5.06 * u.s,
    },

    'Nd-157': {
        'atomic number': 60,
        'mass number': 157,
        'mass': 156.939386 * u.u,
        'stable': False,
        'half-life': 1.15 * u.s,
    },

    'Nd-158': {
        'atomic number': 60,
        'mass number': 158,
        'mass': 157.94197 * u.u,
        'stable': False,
        'half-life': 0.81 * u.s,
    },

    'Nd-159': {
        'atomic number': 60,
        'mass number': 159,
        'mass': 158.94653 * u.u,
        'stable': False,
        'half-life': 0.5 * u.s,
    },

    'Nd-160': {
        'atomic number': 60,
        'mass number': 160,
        'mass': 159.9494 * u.u,
        'stable': False,
        'half-life': 0.439 * u.s,
    },

    'Nd-161': {
        'atomic number': 60,
        'mass number': 161,
        'mass': 160.95428 * u.u,
        'stable': False,
        'half-life': 0.215 * u.s,
    },

    'Pm-126': {
        'atomic number': 61,
        'mass number': 126,
        'mass': 125.95792 * u.u,
        'stable': False,
        'half-life': '500# ms',
    },

    'Pm-127': {
        'atomic number': 61,
        'mass number': 127,
        'mass': 126.95192 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Pm-128': {
        'atomic number': 61,
        'mass number': 128,
        'mass': 127.9487 * u.u,
        'stable': False,
        'half-life': 1.0 * u.s,
    },

    'Pm-129': {
        'atomic number': 61,
        'mass number': 129,
        'mass': 128.94323 * u.u,
        'stable': False,
        'half-life': 2.4 * u.s,
    },

    'Pm-130': {
        'atomic number': 61,
        'mass number': 130,
        'mass': 129.94053 * u.u,
        'stable': False,
        'half-life': 2.6 * u.s,
    },

    'Pm-131': {
        'atomic number': 61,
        'mass number': 131,
        'mass': 130.93567 * u.u,
        'stable': False,
        'half-life': 6.3 * u.s,
    },

    'Pm-132': {
        'atomic number': 61,
        'mass number': 132,
        'mass': 131.93384 * u.u,
        'stable': False,
        'half-life': 6.2 * u.s,
    },

    'Pm-133': {
        'atomic number': 61,
        'mass number': 133,
        'mass': 132.929782 * u.u,
        'stable': False,
        'half-life': 13.5 * u.s,
    },

    'Pm-134': {
        'atomic number': 61,
        'mass number': 134,
        'mass': 133.928353 * u.u,
        'stable': False,
        'half-life': 22.0 * u.s,
    },

    'Pm-135': {
        'atomic number': 61,
        'mass number': 135,
        'mass': 134.924823 * u.u,
        'stable': False,
        'half-life': 49.0 * u.s,
    },

    'Pm-136': {
        'atomic number': 61,
        'mass number': 136,
        'mass': 135.923585 * u.u,
        'stable': False,
        'half-life': 107.0 * u.s,
    },

    'Pm-137': {
        'atomic number': 61,
        'mass number': 137,
        'mass': 136.92048 * u.u,
        'stable': False,
        'half-life': '2# m',
    },

    'Pm-138': {
        'atomic number': 61,
        'mass number': 138,
        'mass': 137.919548 * u.u,
        'stable': False,
        'half-life': 10.0 * u.s,
    },

    'Pm-139': {
        'atomic number': 61,
        'mass number': 139,
        'mass': 138.9168 * u.u,
        'stable': False,
        'half-life': 249.0 * u.s,
    },

    'Pm-140': {
        'atomic number': 61,
        'mass number': 140,
        'mass': 139.91604 * u.u,
        'stable': False,
        'half-life': 9.2 * u.s,
    },

    'Pm-141': {
        'atomic number': 61,
        'mass number': 141,
        'mass': 140.913555 * u.u,
        'stable': False,
        'half-life': 1254.0 * u.s,
    },

    'Pm-142': {
        'atomic number': 61,
        'mass number': 142,
        'mass': 141.91289 * u.u,
        'stable': False,
        'half-life': 40.5 * u.s,
    },

    'Pm-143': {
        'atomic number': 61,
        'mass number': 143,
        'mass': 142.9109383 * u.u,
        'stable': False,
        'half-life': 22896000.0 * u.s,
    },

    'Pm-144': {
        'atomic number': 61,
        'mass number': 144,
        'mass': 143.9125964 * u.u,
        'stable': False,
        'half-life': 31363200.0 * u.s,
    },

    'Pm-145': {
        'atomic number': 61,
        'mass number': 145,
        'mass': 144.9127559 * u.u,
        'stable': False,
        'half-life': 558557590.2 * u.s,
    },

    'Pm-146': {
        'atomic number': 61,
        'mass number': 146,
        'mass': 145.9147024 * u.u,
        'stable': False,
        'half-life': 174509800.78 * u.s,
    },

    'Pm-147': {
        'atomic number': 61,
        'mass number': 147,
        'mass': 146.915145 * u.u,
        'stable': False,
        'half-life': 82786439.6684 * u.s,
    },

    'Pm-148': {
        'atomic number': 61,
        'mass number': 148,
        'mass': 147.9174819 * u.u,
        'stable': False,
        'half-life': 463795.2 * u.s,
    },

    'Pm-149': {
        'atomic number': 61,
        'mass number': 149,
        'mass': 148.9183423 * u.u,
        'stable': False,
        'half-life': 191088.0 * u.s,
    },

    'Pm-150': {
        'atomic number': 61,
        'mass number': 150,
        'mass': 149.920991 * u.u,
        'stable': False,
        'half-life': 9712.8 * u.s,
    },

    'Pm-151': {
        'atomic number': 61,
        'mass number': 151,
        'mass': 150.9212175 * u.u,
        'stable': False,
        'half-life': 102240.0 * u.s,
    },

    'Pm-152': {
        'atomic number': 61,
        'mass number': 152,
        'mass': 151.923506 * u.u,
        'stable': False,
        'half-life': 247.2 * u.s,
    },

    'Pm-153': {
        'atomic number': 61,
        'mass number': 153,
        'mass': 152.9241567 * u.u,
        'stable': False,
        'half-life': 315.0 * u.s,
    },

    'Pm-154': {
        'atomic number': 61,
        'mass number': 154,
        'mass': 153.926472 * u.u,
        'stable': False,
        'half-life': 160.8 * u.s,
    },

    'Pm-155': {
        'atomic number': 61,
        'mass number': 155,
        'mass': 154.928137 * u.u,
        'stable': False,
        'half-life': 41.5 * u.s,
    },

    'Pm-156': {
        'atomic number': 61,
        'mass number': 156,
        'mass': 155.9311175 * u.u,
        'stable': False,
        'half-life': 27.2 * u.s,
    },

    'Pm-157': {
        'atomic number': 61,
        'mass number': 157,
        'mass': 156.9331214 * u.u,
        'stable': False,
        'half-life': 10.56 * u.s,
    },

    'Pm-158': {
        'atomic number': 61,
        'mass number': 158,
        'mass': 157.936565 * u.u,
        'stable': False,
        'half-life': 4.8 * u.s,
    },

    'Pm-159': {
        'atomic number': 61,
        'mass number': 159,
        'mass': 158.939287 * u.u,
        'stable': False,
        'half-life': 1.49 * u.s,
    },

    'Pm-160': {
        'atomic number': 61,
        'mass number': 160,
        'mass': 159.9431 * u.u,
        'stable': False,
        'half-life': 0.725 * u.s,
    },

    'Pm-161': {
        'atomic number': 61,
        'mass number': 161,
        'mass': 160.94607 * u.u,
        'stable': False,
        'half-life': 1.05 * u.s,
    },

    'Pm-162': {
        'atomic number': 61,
        'mass number': 162,
        'mass': 161.95022 * u.u,
        'stable': False,
        'half-life': 0.63 * u.s,
    },

    'Pm-163': {
        'atomic number': 61,
        'mass number': 163,
        'mass': 162.95357 * u.u,
        'stable': False,
        'half-life': 0.43 * u.s,
    },

    'Sm-128': {
        'atomic number': 62,
        'mass number': 128,
        'mass': 127.95842 * u.u,
        'stable': False,
        'half-life': '500# ms',
    },

    'Sm-129': {
        'atomic number': 62,
        'mass number': 129,
        'mass': 128.95476 * u.u,
        'stable': False,
        'half-life': 0.55 * u.s,
    },

    'Sm-130': {
        'atomic number': 62,
        'mass number': 130,
        'mass': 129.949 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Sm-131': {
        'atomic number': 62,
        'mass number': 131,
        'mass': 130.94618 * u.u,
        'stable': False,
        'half-life': 1.2 * u.s,
    },

    'Sm-132': {
        'atomic number': 62,
        'mass number': 132,
        'mass': 131.94087 * u.u,
        'stable': False,
        'half-life': 4.0 * u.s,
    },

    'Sm-133': {
        'atomic number': 62,
        'mass number': 133,
        'mass': 132.93856 * u.u,
        'stable': False,
        'half-life': 2.89 * u.s,
    },

    'Sm-134': {
        'atomic number': 62,
        'mass number': 134,
        'mass': 133.93411 * u.u,
        'stable': False,
        'half-life': 9.5 * u.s,
    },

    'Sm-135': {
        'atomic number': 62,
        'mass number': 135,
        'mass': 134.93252 * u.u,
        'stable': False,
        'half-life': 10.3 * u.s,
    },

    'Sm-136': {
        'atomic number': 62,
        'mass number': 136,
        'mass': 135.928276 * u.u,
        'stable': False,
        'half-life': 47.0 * u.s,
    },

    'Sm-137': {
        'atomic number': 62,
        'mass number': 137,
        'mass': 136.926971 * u.u,
        'stable': False,
        'half-life': 45.0 * u.s,
    },

    'Sm-138': {
        'atomic number': 62,
        'mass number': 138,
        'mass': 137.923244 * u.u,
        'stable': False,
        'half-life': 186.0 * u.s,
    },

    'Sm-139': {
        'atomic number': 62,
        'mass number': 139,
        'mass': 138.922297 * u.u,
        'stable': False,
        'half-life': 154.2 * u.s,
    },

    'Sm-140': {
        'atomic number': 62,
        'mass number': 140,
        'mass': 139.918995 * u.u,
        'stable': False,
        'half-life': 889.2 * u.s,
    },

    'Sm-141': {
        'atomic number': 62,
        'mass number': 141,
        'mass': 140.9184816 * u.u,
        'stable': False,
        'half-life': 612.0 * u.s,
    },

    'Sm-142': {
        'atomic number': 62,
        'mass number': 142,
        'mass': 141.9152044 * u.u,
        'stable': False,
        'half-life': 4349.4 * u.s,
    },

    'Sm-143': {
        'atomic number': 62,
        'mass number': 143,
        'mass': 142.9146353 * u.u,
        'stable': False,
        'half-life': 525.0 * u.s,
    },

    'Sm-144': {
        'atomic number': 62,
        'mass number': 144,
        'mass': 143.9120065 * u.u,
        'stable': True,
        'abundance': 0.0307,
    },

    'Sm-145': {
        'atomic number': 62,
        'mass number': 145,
        'mass': 144.9134173 * u.u,
        'stable': False,
        'half-life': 29376000.0 * u.s,
    },

    'Sm-146': {
        'atomic number': 62,
        'mass number': 146,
        'mass': 145.913047 * u.u,
        'stable': False,
        'half-life': 2145870968000000.0 * u.s,
    },

    'Sm-147': {
        'atomic number': 62,
        'mass number': 147,
        'mass': 146.9149044 * u.u,
        'stable': False,
        'abundance': 0.1499,
    },

    'Sm-148': {
        'atomic number': 62,
        'mass number': 148,
        'mass': 147.9148292 * u.u,
        'stable': False,
        'abundance': 0.1124,
    },

    'Sm-149': {
        'atomic number': 62,
        'mass number': 149,
        'mass': 148.9171921 * u.u,
        'stable': True,
        'abundance': 0.1382,
    },

    'Sm-150': {
        'atomic number': 62,
        'mass number': 150,
        'mass': 149.9172829 * u.u,
        'stable': True,
        'abundance': 0.0738,
    },

    'Sm-151': {
        'atomic number': 62,
        'mass number': 151,
        'mass': 150.9199398 * u.u,
        'stable': False,
        'half-life': 2840123340.0 * u.s,
    },

    'Sm-152': {
        'atomic number': 62,
        'mass number': 152,
        'mass': 151.9197397 * u.u,
        'stable': True,
        'abundance': 0.2675,
    },

    'Sm-153': {
        'atomic number': 62,
        'mass number': 153,
        'mass': 152.9221047 * u.u,
        'stable': False,
        'half-life': 166627.08 * u.s,
    },

    'Sm-154': {
        'atomic number': 62,
        'mass number': 154,
        'mass': 153.9222169 * u.u,
        'stable': True,
        'abundance': 0.2275,
    },

    'Sm-155': {
        'atomic number': 62,
        'mass number': 155,
        'mass': 154.9246477 * u.u,
        'stable': False,
        'half-life': 1338.0 * u.s,
    },

    'Sm-156': {
        'atomic number': 62,
        'mass number': 156,
        'mass': 155.925536 * u.u,
        'stable': False,
        'half-life': 33840.0 * u.s,
    },

    'Sm-157': {
        'atomic number': 62,
        'mass number': 157,
        'mass': 156.9284187 * u.u,
        'stable': False,
        'half-life': 481.8 * u.s,
    },

    'Sm-158': {
        'atomic number': 62,
        'mass number': 158,
        'mass': 157.929951 * u.u,
        'stable': False,
        'half-life': 318.0 * u.s,
    },

    'Sm-159': {
        'atomic number': 62,
        'mass number': 159,
        'mass': 158.9332172 * u.u,
        'stable': False,
        'half-life': 11.37 * u.s,
    },

    'Sm-160': {
        'atomic number': 62,
        'mass number': 160,
        'mass': 159.9353353 * u.u,
        'stable': False,
        'half-life': 9.6 * u.s,
    },

    'Sm-161': {
        'atomic number': 62,
        'mass number': 161,
        'mass': 160.9391602 * u.u,
        'stable': False,
        'half-life': 4.8 * u.s,
    },

    'Sm-162': {
        'atomic number': 62,
        'mass number': 162,
        'mass': 161.94146 * u.u,
        'stable': False,
        'half-life': 2.7 * u.s,
    },

    'Sm-163': {
        'atomic number': 62,
        'mass number': 163,
        'mass': 162.94555 * u.u,
        'stable': False,
        'half-life': 1.3 * u.s,
    },

    'Sm-164': {
        'atomic number': 62,
        'mass number': 164,
        'mass': 163.94836 * u.u,
        'stable': False,
        'half-life': 1.43 * u.s,
    },

    'Sm-165': {
        'atomic number': 62,
        'mass number': 165,
        'mass': 164.95297 * u.u,
        'stable': False,
        'half-life': 0.98 * u.s,
    },

    'Eu-130': {
        'atomic number': 63,
        'mass number': 130,
        'mass': 129.96369 * u.u,
        'stable': False,
        'half-life': 0.001 * u.s,
    },

    'Eu-131': {
        'atomic number': 63,
        'mass number': 131,
        'mass': 130.95784 * u.u,
        'stable': False,
        'half-life': 0.0178 * u.s,
    },

    'Eu-132': {
        'atomic number': 63,
        'mass number': 132,
        'mass': 131.95467 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Eu-133': {
        'atomic number': 63,
        'mass number': 133,
        'mass': 132.94929 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'Eu-134': {
        'atomic number': 63,
        'mass number': 134,
        'mass': 133.9464 * u.u,
        'stable': False,
        'half-life': 0.5 * u.s,
    },

    'Eu-135': {
        'atomic number': 63,
        'mass number': 135,
        'mass': 134.94187 * u.u,
        'stable': False,
        'half-life': 1.5 * u.s,
    },

    'Eu-136': {
        'atomic number': 63,
        'mass number': 136,
        'mass': 135.93962 * u.u,
        'stable': False,
        'half-life': 3.3 * u.s,
    },

    'Eu-137': {
        'atomic number': 63,
        'mass number': 137,
        'mass': 136.93546 * u.u,
        'stable': False,
        'half-life': 8.4 * u.s,
    },

    'Eu-138': {
        'atomic number': 63,
        'mass number': 138,
        'mass': 137.933709 * u.u,
        'stable': False,
        'half-life': 12.1 * u.s,
    },

    'Eu-139': {
        'atomic number': 63,
        'mass number': 139,
        'mass': 138.929792 * u.u,
        'stable': False,
        'half-life': 17.9 * u.s,
    },

    'Eu-140': {
        'atomic number': 63,
        'mass number': 140,
        'mass': 139.928088 * u.u,
        'stable': False,
        'half-life': 1.51 * u.s,
    },

    'Eu-141': {
        'atomic number': 63,
        'mass number': 141,
        'mass': 140.924932 * u.u,
        'stable': False,
        'half-life': 40.7 * u.s,
    },

    'Eu-142': {
        'atomic number': 63,
        'mass number': 142,
        'mass': 141.923442 * u.u,
        'stable': False,
        'half-life': 2.36 * u.s,
    },

    'Eu-143': {
        'atomic number': 63,
        'mass number': 143,
        'mass': 142.920299 * u.u,
        'stable': False,
        'half-life': 155.4 * u.s,
    },

    'Eu-144': {
        'atomic number': 63,
        'mass number': 144,
        'mass': 143.91882 * u.u,
        'stable': False,
        'half-life': 10.2 * u.s,
    },

    'Eu-145': {
        'atomic number': 63,
        'mass number': 145,
        'mass': 144.9162726 * u.u,
        'stable': False,
        'half-life': 512352.0 * u.s,
    },

    'Eu-146': {
        'atomic number': 63,
        'mass number': 146,
        'mass': 145.917211 * u.u,
        'stable': False,
        'half-life': 398304.0 * u.s,
    },

    'Eu-147': {
        'atomic number': 63,
        'mass number': 147,
        'mass': 146.9167527 * u.u,
        'stable': False,
        'half-life': 2082240.0 * u.s,
    },

    'Eu-148': {
        'atomic number': 63,
        'mass number': 148,
        'mass': 147.918089 * u.u,
        'stable': False,
        'half-life': 4708800.0 * u.s,
    },

    'Eu-149': {
        'atomic number': 63,
        'mass number': 149,
        'mass': 148.9179378 * u.u,
        'stable': False,
        'half-life': 8043840.0 * u.s,
    },

    'Eu-150': {
        'atomic number': 63,
        'mass number': 150,
        'mass': 149.9197077 * u.u,
        'stable': False,
        'half-life': 1164450569.4 * u.s,
    },

    'Eu-151': {
        'atomic number': 63,
        'mass number': 151,
        'mass': 150.9198578 * u.u,
        'stable': False,
        'abundance': 0.4781,
    },

    'Eu-152': {
        'atomic number': 63,
        'mass number': 152,
        'mass': 151.9217522 * u.u,
        'stable': False,
        'half-life': 427438080.0 * u.s,
    },

    'Eu-153': {
        'atomic number': 63,
        'mass number': 153,
        'mass': 152.921238 * u.u,
        'stable': True,
        'abundance': 0.5219,
    },

    'Eu-154': {
        'atomic number': 63,
        'mass number': 154,
        'mass': 153.922987 * u.u,
        'stable': False,
        'half-life': 271745280.0 * u.s,
    },

    'Eu-155': {
        'atomic number': 63,
        'mass number': 155,
        'mass': 154.9229011 * u.u,
        'stable': False,
        'half-life': 150254784.0 * u.s,
    },

    'Eu-156': {
        'atomic number': 63,
        'mass number': 156,
        'mass': 155.9247605 * u.u,
        'stable': False,
        'half-life': 1312416.0 * u.s,
    },

    'Eu-157': {
        'atomic number': 63,
        'mass number': 157,
        'mass': 156.9254334 * u.u,
        'stable': False,
        'half-life': 54648.0 * u.s,
    },

    'Eu-158': {
        'atomic number': 63,
        'mass number': 158,
        'mass': 157.927799 * u.u,
        'stable': False,
        'half-life': 2754.0 * u.s,
    },

    'Eu-159': {
        'atomic number': 63,
        'mass number': 159,
        'mass': 158.9291001 * u.u,
        'stable': False,
        'half-life': 1086.0 * u.s,
    },

    'Eu-160': {
        'atomic number': 63,
        'mass number': 160,
        'mass': 159.931851 * u.u,
        'stable': False,
        'half-life': 42.4 * u.s,
    },

    'Eu-161': {
        'atomic number': 63,
        'mass number': 161,
        'mass': 160.933664 * u.u,
        'stable': False,
        'half-life': 26.2 * u.s,
    },

    'Eu-162': {
        'atomic number': 63,
        'mass number': 162,
        'mass': 161.936989 * u.u,
        'stable': False,
        'half-life': '~11 s',
    },

    'Eu-163': {
        'atomic number': 63,
        'mass number': 163,
        'mass': 162.939196 * u.u,
        'stable': False,
        'half-life': 7.7 * u.s,
    },

    'Eu-164': {
        'atomic number': 63,
        'mass number': 164,
        'mass': 163.94274 * u.u,
        'stable': False,
        'half-life': 4.15 * u.s,
    },

    'Eu-165': {
        'atomic number': 63,
        'mass number': 165,
        'mass': 164.94559 * u.u,
        'stable': False,
        'half-life': 2.53 * u.s,
    },

    'Eu-166': {
        'atomic number': 63,
        'mass number': 166,
        'mass': 165.94962 * u.u,
        'stable': False,
        'half-life': 1.24 * u.s,
    },

    'Eu-167': {
        'atomic number': 63,
        'mass number': 167,
        'mass': 166.95289 * u.u,
        'stable': False,
        'half-life': 1.33 * u.s,
    },

    'Gd-133': {
        'atomic number': 64,
        'mass number': 133,
        'mass': 132.96133 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Gd-134': {
        'atomic number': 64,
        'mass number': 134,
        'mass': 133.95566 * u.u,
        'stable': False,
        'half-life': '400# ms',
    },

    'Gd-135': {
        'atomic number': 64,
        'mass number': 135,
        'mass': 134.95245 * u.u,
        'stable': False,
        'half-life': 1.1 * u.s,
    },

    'Gd-136': {
        'atomic number': 64,
        'mass number': 136,
        'mass': 135.9473 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Gd-137': {
        'atomic number': 64,
        'mass number': 137,
        'mass': 136.94502 * u.u,
        'stable': False,
        'half-life': 2.2 * u.s,
    },

    'Gd-138': {
        'atomic number': 64,
        'mass number': 138,
        'mass': 137.94025 * u.u,
        'stable': False,
        'half-life': 4.7 * u.s,
    },

    'Gd-139': {
        'atomic number': 64,
        'mass number': 139,
        'mass': 138.93813 * u.u,
        'stable': False,
        'half-life': 5.7 * u.s,
    },

    'Gd-140': {
        'atomic number': 64,
        'mass number': 140,
        'mass': 139.933674 * u.u,
        'stable': False,
        'half-life': 15.8 * u.s,
    },

    'Gd-141': {
        'atomic number': 64,
        'mass number': 141,
        'mass': 140.932126 * u.u,
        'stable': False,
        'half-life': 14.0 * u.s,
    },

    'Gd-142': {
        'atomic number': 64,
        'mass number': 142,
        'mass': 141.928116 * u.u,
        'stable': False,
        'half-life': 70.2 * u.s,
    },

    'Gd-143': {
        'atomic number': 64,
        'mass number': 143,
        'mass': 142.92675 * u.u,
        'stable': False,
        'half-life': 39.0 * u.s,
    },

    'Gd-144': {
        'atomic number': 64,
        'mass number': 144,
        'mass': 143.922963 * u.u,
        'stable': False,
        'half-life': 268.2 * u.s,
    },

    'Gd-145': {
        'atomic number': 64,
        'mass number': 145,
        'mass': 144.921713 * u.u,
        'stable': False,
        'half-life': 1380.0 * u.s,
    },

    'Gd-146': {
        'atomic number': 64,
        'mass number': 146,
        'mass': 145.9183188 * u.u,
        'stable': False,
        'half-life': 4170528.0 * u.s,
    },

    'Gd-147': {
        'atomic number': 64,
        'mass number': 147,
        'mass': 146.9191014 * u.u,
        'stable': False,
        'half-life': 137016.0 * u.s,
    },

    'Gd-148': {
        'atomic number': 64,
        'mass number': 148,
        'mass': 147.9181215 * u.u,
        'stable': False,
        'half-life': 2237386053.4 * u.s,
    },

    'Gd-149': {
        'atomic number': 64,
        'mass number': 149,
        'mass': 148.9193481 * u.u,
        'stable': False,
        'half-life': 801792.0 * u.s,
    },

    'Gd-150': {
        'atomic number': 64,
        'mass number': 150,
        'mass': 149.9186644 * u.u,
        'stable': False,
        'half-life': 56486897540000.0 * u.s,
    },

    'Gd-151': {
        'atomic number': 64,
        'mass number': 151,
        'mass': 150.920356 * u.u,
        'stable': False,
        'half-life': 10704960.0 * u.s,
    },

    'Gd-152': {
        'atomic number': 64,
        'mass number': 152,
        'mass': 151.9197995 * u.u,
        'stable': False,
        'abundance': 0.002,
    },

    'Gd-153': {
        'atomic number': 64,
        'mass number': 153,
        'mass': 152.921758 * u.u,
        'stable': False,
        'half-life': 20690380.8 * u.s,
    },

    'Gd-154': {
        'atomic number': 64,
        'mass number': 154,
        'mass': 153.9208741 * u.u,
        'stable': True,
        'abundance': 0.0218,
    },

    'Gd-155': {
        'atomic number': 64,
        'mass number': 155,
        'mass': 154.9226305 * u.u,
        'stable': True,
        'abundance': 0.148,
    },

    'Gd-156': {
        'atomic number': 64,
        'mass number': 156,
        'mass': 155.9221312 * u.u,
        'stable': True,
        'abundance': 0.2047,
    },

    'Gd-157': {
        'atomic number': 64,
        'mass number': 157,
        'mass': 156.9239686 * u.u,
        'stable': True,
        'abundance': 0.1565,
    },

    'Gd-158': {
        'atomic number': 64,
        'mass number': 158,
        'mass': 157.9241123 * u.u,
        'stable': True,
        'abundance': 0.2484,
    },

    'Gd-159': {
        'atomic number': 64,
        'mass number': 159,
        'mass': 158.926397 * u.u,
        'stable': False,
        'half-life': 66524.4 * u.s,
    },

    'Gd-160': {
        'atomic number': 64,
        'mass number': 160,
        'mass': 159.9270624 * u.u,
        'stable': True,
        'abundance': 0.2186,
    },

    'Gd-161': {
        'atomic number': 64,
        'mass number': 161,
        'mass': 160.9296775 * u.u,
        'stable': False,
        'half-life': 218.76 * u.s,
    },

    'Gd-162': {
        'atomic number': 64,
        'mass number': 162,
        'mass': 161.930993 * u.u,
        'stable': False,
        'half-life': 504.0 * u.s,
    },

    'Gd-163': {
        'atomic number': 64,
        'mass number': 163,
        'mass': 162.9341769 * u.u,
        'stable': False,
        'half-life': 68.0 * u.s,
    },

    'Gd-164': {
        'atomic number': 64,
        'mass number': 164,
        'mass': 163.93583 * u.u,
        'stable': False,
        'half-life': 45.0 * u.s,
    },

    'Gd-165': {
        'atomic number': 64,
        'mass number': 165,
        'mass': 164.93936 * u.u,
        'stable': False,
        'half-life': 11.0 * u.s,
    },

    'Gd-166': {
        'atomic number': 64,
        'mass number': 166,
        'mass': 165.94146 * u.u,
        'stable': False,
        'half-life': 5.1 * u.s,
    },

    'Gd-167': {
        'atomic number': 64,
        'mass number': 167,
        'mass': 166.94545 * u.u,
        'stable': False,
        'half-life': 4.2 * u.s,
    },

    'Gd-168': {
        'atomic number': 64,
        'mass number': 168,
        'mass': 167.94808 * u.u,
        'stable': False,
        'half-life': 3.03 * u.s,
    },

    'Gd-169': {
        'atomic number': 64,
        'mass number': 169,
        'mass': 168.9526 * u.u,
        'stable': False,
        'half-life': 0.75 * u.s,
    },

    'Tb-135': {
        'atomic number': 65,
        'mass number': 135,
        'mass': 134.96476 * u.u,
        'stable': False,
        'half-life': 0.00101 * u.s,
    },

    'Tb-136': {
        'atomic number': 65,
        'mass number': 136,
        'mass': 135.96129 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'Tb-137': {
        'atomic number': 65,
        'mass number': 137,
        'mass': 136.95602 * u.u,
        'stable': False,
        'half-life': '600# ms',
    },

    'Tb-138': {
        'atomic number': 65,
        'mass number': 138,
        'mass': 137.95312 * u.u,
        'stable': False,
        'half-life': '800# ms',
    },

    'Tb-139': {
        'atomic number': 65,
        'mass number': 139,
        'mass': 138.94833 * u.u,
        'stable': False,
        'half-life': 1.6 * u.s,
    },

    'Tb-140': {
        'atomic number': 65,
        'mass number': 140,
        'mass': 139.94581 * u.u,
        'stable': False,
        'half-life': 2.32 * u.s,
    },

    'Tb-141': {
        'atomic number': 65,
        'mass number': 141,
        'mass': 140.94145 * u.u,
        'stable': False,
        'half-life': 3.5 * u.s,
    },

    'Tb-142': {
        'atomic number': 65,
        'mass number': 142,
        'mass': 141.93928 * u.u,
        'stable': False,
        'half-life': 0.597 * u.s,
    },

    'Tb-143': {
        'atomic number': 65,
        'mass number': 143,
        'mass': 142.935137 * u.u,
        'stable': False,
        'half-life': 12.0 * u.s,
    },

    'Tb-144': {
        'atomic number': 65,
        'mass number': 144,
        'mass': 143.933045 * u.u,
        'stable': False,
        'half-life': '~1 s',
    },

    'Tb-145': {
        'atomic number': 65,
        'mass number': 145,
        'mass': 144.92882 * u.u,
        'stable': False,
        'half-life': 30.9 * u.s,
    },

    'Tb-146': {
        'atomic number': 65,
        'mass number': 146,
        'mass': 145.927253 * u.u,
        'stable': False,
        'half-life': 8.0 * u.s,
    },

    'Tb-147': {
        'atomic number': 65,
        'mass number': 147,
        'mass': 146.9240548 * u.u,
        'stable': False,
        'half-life': 5904.0 * u.s,
    },

    'Tb-148': {
        'atomic number': 65,
        'mass number': 148,
        'mass': 147.924282 * u.u,
        'stable': False,
        'half-life': 3600.0 * u.s,
    },

    'Tb-149': {
        'atomic number': 65,
        'mass number': 149,
        'mass': 148.9232535 * u.u,
        'stable': False,
        'half-life': 14824.8 * u.s,
    },

    'Tb-150': {
        'atomic number': 65,
        'mass number': 150,
        'mass': 149.9236649 * u.u,
        'stable': False,
        'half-life': 12528.0 * u.s,
    },

    'Tb-151': {
        'atomic number': 65,
        'mass number': 151,
        'mass': 150.9231096 * u.u,
        'stable': False,
        'half-life': 63392.4 * u.s,
    },

    'Tb-152': {
        'atomic number': 65,
        'mass number': 152,
        'mass': 151.924083 * u.u,
        'stable': False,
        'half-life': 63000.0 * u.s,
    },

    'Tb-153': {
        'atomic number': 65,
        'mass number': 153,
        'mass': 152.9234424 * u.u,
        'stable': False,
        'half-life': 202176.0 * u.s,
    },

    'Tb-154': {
        'atomic number': 65,
        'mass number': 154,
        'mass': 153.924685 * u.u,
        'stable': False,
        'half-life': 77400.0 * u.s,
    },

    'Tb-155': {
        'atomic number': 65,
        'mass number': 155,
        'mass': 154.923511 * u.u,
        'stable': False,
        'half-life': 459648.0 * u.s,
    },

    'Tb-156': {
        'atomic number': 65,
        'mass number': 156,
        'mass': 155.9247552 * u.u,
        'stable': False,
        'half-life': 462240.0 * u.s,
    },

    'Tb-157': {
        'atomic number': 65,
        'mass number': 157,
        'mass': 156.924033 * u.u,
        'stable': False,
        'half-life': 2240541746.0 * u.s,
    },

    'Tb-158': {
        'atomic number': 65,
        'mass number': 158,
        'mass': 157.9254209 * u.u,
        'stable': False,
        'half-life': 5680246680.0 * u.s,
    },

    'Tb-159': {
        'atomic number': 65,
        'mass number': 159,
        'mass': 158.9253547 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Tb-160': {
        'atomic number': 65,
        'mass number': 160,
        'mass': 159.9271756 * u.u,
        'stable': False,
        'half-life': 6246720.0 * u.s,
    },

    'Tb-161': {
        'atomic number': 65,
        'mass number': 161,
        'mass': 160.9275778 * u.u,
        'stable': False,
        'half-life': 595296.0 * u.s,
    },

    'Tb-162': {
        'atomic number': 65,
        'mass number': 162,
        'mass': 161.929495 * u.u,
        'stable': False,
        'half-life': 456.0 * u.s,
    },

    'Tb-163': {
        'atomic number': 65,
        'mass number': 163,
        'mass': 162.9306547 * u.u,
        'stable': False,
        'half-life': 1170.0 * u.s,
    },

    'Tb-164': {
        'atomic number': 65,
        'mass number': 164,
        'mass': 163.93336 * u.u,
        'stable': False,
        'half-life': 180.0 * u.s,
    },

    'Tb-165': {
        'atomic number': 65,
        'mass number': 165,
        'mass': 164.93498 * u.u,
        'stable': False,
        'half-life': 126.6 * u.s,
    },

    'Tb-166': {
        'atomic number': 65,
        'mass number': 166,
        'mass': 165.93786 * u.u,
        'stable': False,
        'half-life': 27.1 * u.s,
    },

    'Tb-167': {
        'atomic number': 65,
        'mass number': 167,
        'mass': 166.93996 * u.u,
        'stable': False,
        'half-life': 18.9 * u.s,
    },

    'Tb-168': {
        'atomic number': 65,
        'mass number': 168,
        'mass': 167.9434 * u.u,
        'stable': False,
        'half-life': 9.4 * u.s,
    },

    'Tb-169': {
        'atomic number': 65,
        'mass number': 169,
        'mass': 168.94597 * u.u,
        'stable': False,
        'half-life': 5.13 * u.s,
    },

    'Tb-170': {
        'atomic number': 65,
        'mass number': 170,
        'mass': 169.94984 * u.u,
        'stable': False,
        'half-life': 0.96 * u.s,
    },

    'Tb-171': {
        'atomic number': 65,
        'mass number': 171,
        'mass': 170.95273 * u.u,
        'stable': False,
        'half-life': 1.23 * u.s,
    },

    'Dy-138': {
        'atomic number': 66,
        'mass number': 138,
        'mass': 137.9625 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'Dy-139': {
        'atomic number': 66,
        'mass number': 139,
        'mass': 138.95959 * u.u,
        'stable': False,
        'half-life': 0.6 * u.s,
    },

    'Dy-140': {
        'atomic number': 66,
        'mass number': 140,
        'mass': 139.95402 * u.u,
        'stable': False,
        'half-life': '700# ms',
    },

    'Dy-141': {
        'atomic number': 66,
        'mass number': 141,
        'mass': 140.95128 * u.u,
        'stable': False,
        'half-life': 0.9 * u.s,
    },

    'Dy-142': {
        'atomic number': 66,
        'mass number': 142,
        'mass': 141.94619 * u.u,
        'stable': False,
        'half-life': 2.3 * u.s,
    },

    'Dy-143': {
        'atomic number': 66,
        'mass number': 143,
        'mass': 142.943994 * u.u,
        'stable': False,
        'half-life': 5.6 * u.s,
    },

    'Dy-144': {
        'atomic number': 66,
        'mass number': 144,
        'mass': 143.9392695 * u.u,
        'stable': False,
        'half-life': 9.1 * u.s,
    },

    'Dy-145': {
        'atomic number': 66,
        'mass number': 145,
        'mass': 144.937474 * u.u,
        'stable': False,
        'half-life': 9.5 * u.s,
    },

    'Dy-146': {
        'atomic number': 66,
        'mass number': 146,
        'mass': 145.9328445 * u.u,
        'stable': False,
        'half-life': 33.2 * u.s,
    },

    'Dy-147': {
        'atomic number': 66,
        'mass number': 147,
        'mass': 146.9310827 * u.u,
        'stable': False,
        'half-life': 67.0 * u.s,
    },

    'Dy-148': {
        'atomic number': 66,
        'mass number': 148,
        'mass': 147.927157 * u.u,
        'stable': False,
        'half-life': 198.0 * u.s,
    },

    'Dy-149': {
        'atomic number': 66,
        'mass number': 149,
        'mass': 148.927322 * u.u,
        'stable': False,
        'half-life': 252.0 * u.s,
    },

    'Dy-150': {
        'atomic number': 66,
        'mass number': 150,
        'mass': 149.9255933 * u.u,
        'stable': False,
        'half-life': 430.2 * u.s,
    },

    'Dy-151': {
        'atomic number': 66,
        'mass number': 151,
        'mass': 150.9261916 * u.u,
        'stable': False,
        'half-life': 1074.0 * u.s,
    },

    'Dy-152': {
        'atomic number': 66,
        'mass number': 152,
        'mass': 151.9247253 * u.u,
        'stable': False,
        'half-life': 8568.0 * u.s,
    },

    'Dy-153': {
        'atomic number': 66,
        'mass number': 153,
        'mass': 152.9257724 * u.u,
        'stable': False,
        'half-life': 23040.0 * u.s,
    },

    'Dy-154': {
        'atomic number': 66,
        'mass number': 154,
        'mass': 153.9244293 * u.u,
        'stable': False,
        'half-life': 94670778000000.0 * u.s,
    },

    'Dy-155': {
        'atomic number': 66,
        'mass number': 155,
        'mass': 154.925759 * u.u,
        'stable': False,
        'half-life': 35640.0 * u.s,
    },

    'Dy-156': {
        'atomic number': 66,
        'mass number': 156,
        'mass': 155.9242847 * u.u,
        'stable': True,
        'abundance': 0.00056,
    },

    'Dy-157': {
        'atomic number': 66,
        'mass number': 157,
        'mass': 156.9254707 * u.u,
        'stable': False,
        'half-life': 29304.0 * u.s,
    },

    'Dy-158': {
        'atomic number': 66,
        'mass number': 158,
        'mass': 157.9244159 * u.u,
        'stable': True,
        'abundance': 0.00095,
    },

    'Dy-159': {
        'atomic number': 66,
        'mass number': 159,
        'mass': 158.925747 * u.u,
        'stable': False,
        'half-life': 12476160.0 * u.s,
    },

    'Dy-160': {
        'atomic number': 66,
        'mass number': 160,
        'mass': 159.9252046 * u.u,
        'stable': True,
        'abundance': 0.02329,
    },

    'Dy-161': {
        'atomic number': 66,
        'mass number': 161,
        'mass': 160.9269405 * u.u,
        'stable': True,
        'abundance': 0.18889,
    },

    'Dy-162': {
        'atomic number': 66,
        'mass number': 162,
        'mass': 161.9268056 * u.u,
        'stable': True,
        'abundance': 0.25475,
    },

    'Dy-163': {
        'atomic number': 66,
        'mass number': 163,
        'mass': 162.9287383 * u.u,
        'stable': True,
        'abundance': 0.24896,
    },

    'Dy-164': {
        'atomic number': 66,
        'mass number': 164,
        'mass': 163.9291819 * u.u,
        'stable': True,
        'abundance': 0.2826,
    },

    'Dy-165': {
        'atomic number': 66,
        'mass number': 165,
        'mass': 164.9317105 * u.u,
        'stable': False,
        'half-life': 8402.4 * u.s,
    },

    'Dy-166': {
        'atomic number': 66,
        'mass number': 166,
        'mass': 165.9328139 * u.u,
        'stable': False,
        'half-life': 293760.0 * u.s,
    },

    'Dy-167': {
        'atomic number': 66,
        'mass number': 167,
        'mass': 166.935661 * u.u,
        'stable': False,
        'half-life': 372.0 * u.s,
    },

    'Dy-168': {
        'atomic number': 66,
        'mass number': 168,
        'mass': 167.93713 * u.u,
        'stable': False,
        'half-life': 522.0 * u.s,
    },

    'Dy-169': {
        'atomic number': 66,
        'mass number': 169,
        'mass': 168.94031 * u.u,
        'stable': False,
        'half-life': 39.0 * u.s,
    },

    'Dy-170': {
        'atomic number': 66,
        'mass number': 170,
        'mass': 169.94239 * u.u,
        'stable': False,
        'half-life': 54.9 * u.s,
    },

    'Dy-171': {
        'atomic number': 66,
        'mass number': 171,
        'mass': 170.94612 * u.u,
        'stable': False,
        'half-life': 4.07 * u.s,
    },

    'Dy-172': {
        'atomic number': 66,
        'mass number': 172,
        'mass': 171.94846 * u.u,
        'stable': False,
        'half-life': 3.4 * u.s,
    },

    'Dy-173': {
        'atomic number': 66,
        'mass number': 173,
        'mass': 172.95283 * u.u,
        'stable': False,
        'half-life': 1.43 * u.s,
    },

    'Ho-140': {
        'atomic number': 67,
        'mass number': 140,
        'mass': 139.96859 * u.u,
        'stable': False,
        'half-life': 0.006 * u.s,
    },

    'Ho-141': {
        'atomic number': 67,
        'mass number': 141,
        'mass': 140.96311 * u.u,
        'stable': False,
        'half-life': 0.0041 * u.s,
    },

    'Ho-142': {
        'atomic number': 67,
        'mass number': 142,
        'mass': 141.96001 * u.u,
        'stable': False,
        'half-life': 0.4 * u.s,
    },

    'Ho-143': {
        'atomic number': 67,
        'mass number': 143,
        'mass': 142.95486 * u.u,
        'stable': False,
        'half-life': '300# ms',
    },

    'Ho-144': {
        'atomic number': 67,
        'mass number': 144,
        'mass': 143.9521097 * u.u,
        'stable': False,
        'half-life': 0.7 * u.s,
    },

    'Ho-145': {
        'atomic number': 67,
        'mass number': 145,
        'mass': 144.9472674 * u.u,
        'stable': False,
        'half-life': 2.4 * u.s,
    },

    'Ho-146': {
        'atomic number': 67,
        'mass number': 146,
        'mass': 145.9449935 * u.u,
        'stable': False,
        'half-life': 2.8 * u.s,
    },

    'Ho-147': {
        'atomic number': 67,
        'mass number': 147,
        'mass': 146.9401423 * u.u,
        'stable': False,
        'half-life': 5.8 * u.s,
    },

    'Ho-148': {
        'atomic number': 67,
        'mass number': 148,
        'mass': 147.937744 * u.u,
        'stable': False,
        'half-life': 2.2 * u.s,
    },

    'Ho-149': {
        'atomic number': 67,
        'mass number': 149,
        'mass': 148.933803 * u.u,
        'stable': False,
        'half-life': 21.1 * u.s,
    },

    'Ho-150': {
        'atomic number': 67,
        'mass number': 150,
        'mass': 149.933498 * u.u,
        'stable': False,
        'half-life': 76.8 * u.s,
    },

    'Ho-151': {
        'atomic number': 67,
        'mass number': 151,
        'mass': 150.9316983 * u.u,
        'stable': False,
        'half-life': 35.2 * u.s,
    },

    'Ho-152': {
        'atomic number': 67,
        'mass number': 152,
        'mass': 151.931724 * u.u,
        'stable': False,
        'half-life': 161.8 * u.s,
    },

    'Ho-153': {
        'atomic number': 67,
        'mass number': 153,
        'mass': 152.9302064 * u.u,
        'stable': False,
        'half-life': 120.6 * u.s,
    },

    'Ho-154': {
        'atomic number': 67,
        'mass number': 154,
        'mass': 153.9306068 * u.u,
        'stable': False,
        'half-life': 705.6 * u.s,
    },

    'Ho-155': {
        'atomic number': 67,
        'mass number': 155,
        'mass': 154.929104 * u.u,
        'stable': False,
        'half-life': 2880.0 * u.s,
    },

    'Ho-156': {
        'atomic number': 67,
        'mass number': 156,
        'mass': 155.929706 * u.u,
        'stable': False,
        'half-life': 3360.0 * u.s,
    },

    'Ho-157': {
        'atomic number': 67,
        'mass number': 157,
        'mass': 156.928254 * u.u,
        'stable': False,
        'half-life': 756.0 * u.s,
    },

    'Ho-158': {
        'atomic number': 67,
        'mass number': 158,
        'mass': 157.928946 * u.u,
        'stable': False,
        'half-life': 678.0 * u.s,
    },

    'Ho-159': {
        'atomic number': 67,
        'mass number': 159,
        'mass': 158.9277197 * u.u,
        'stable': False,
        'half-life': 1983.0 * u.s,
    },

    'Ho-160': {
        'atomic number': 67,
        'mass number': 160,
        'mass': 159.928737 * u.u,
        'stable': False,
        'half-life': 1536.0 * u.s,
    },

    'Ho-161': {
        'atomic number': 67,
        'mass number': 161,
        'mass': 160.9278615 * u.u,
        'stable': False,
        'half-life': 8928.0 * u.s,
    },

    'Ho-162': {
        'atomic number': 67,
        'mass number': 162,
        'mass': 161.9291023 * u.u,
        'stable': False,
        'half-life': 900.0 * u.s,
    },

    'Ho-163': {
        'atomic number': 67,
        'mass number': 163,
        'mass': 162.928741 * u.u,
        'stable': False,
        'half-life': 144215151820.0 * u.s,
    },

    'Ho-164': {
        'atomic number': 67,
        'mass number': 164,
        'mass': 163.9302403 * u.u,
        'stable': False,
        'half-life': 1740.0 * u.s,
    },

    'Ho-165': {
        'atomic number': 67,
        'mass number': 165,
        'mass': 164.9303288 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Ho-166': {
        'atomic number': 67,
        'mass number': 166,
        'mass': 165.9322909 * u.u,
        'stable': False,
        'half-life': 96458.40000000001 * u.s,
    },

    'Ho-167': {
        'atomic number': 67,
        'mass number': 167,
        'mass': 166.9331385 * u.u,
        'stable': False,
        'half-life': 11160.0 * u.s,
    },

    'Ho-168': {
        'atomic number': 67,
        'mass number': 168,
        'mass': 167.935522 * u.u,
        'stable': False,
        'half-life': 179.4 * u.s,
    },

    'Ho-169': {
        'atomic number': 67,
        'mass number': 169,
        'mass': 168.936878 * u.u,
        'stable': False,
        'half-life': 283.2 * u.s,
    },

    'Ho-170': {
        'atomic number': 67,
        'mass number': 170,
        'mass': 169.939625 * u.u,
        'stable': False,
        'half-life': 165.6 * u.s,
    },

    'Ho-171': {
        'atomic number': 67,
        'mass number': 171,
        'mass': 170.94147 * u.u,
        'stable': False,
        'half-life': 53.0 * u.s,
    },

    'Ho-172': {
        'atomic number': 67,
        'mass number': 172,
        'mass': 171.94473 * u.u,
        'stable': False,
        'half-life': 25.0 * u.s,
    },

    'Ho-173': {
        'atomic number': 67,
        'mass number': 173,
        'mass': 172.94702 * u.u,
        'stable': False,
        'half-life': 6.9 * u.s,
    },

    'Ho-174': {
        'atomic number': 67,
        'mass number': 174,
        'mass': 173.95095 * u.u,
        'stable': False,
        'half-life': 3.2 * u.s,
    },

    'Ho-175': {
        'atomic number': 67,
        'mass number': 175,
        'mass': 174.95362 * u.u,
        'stable': False,
        'half-life': 1.88 * u.s,
    },

    'Er-142': {
        'atomic number': 68,
        'mass number': 142,
        'mass': 141.9701 * u.u,
        'stable': False,
        'half-life': '10# us',
    },

    'Er-143': {
        'atomic number': 68,
        'mass number': 143,
        'mass': 142.96662 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'Er-144': {
        'atomic number': 68,
        'mass number': 144,
        'mass': 143.9607 * u.u,
        'stable': False,
        'half-life': '400# ms',
    },

    'Er-145': {
        'atomic number': 68,
        'mass number': 145,
        'mass': 144.95805 * u.u,
        'stable': False,
        'half-life': 0.9 * u.s,
    },

    'Er-146': {
        'atomic number': 68,
        'mass number': 146,
        'mass': 145.9524184 * u.u,
        'stable': False,
        'half-life': 1.7 * u.s,
    },

    'Er-147': {
        'atomic number': 68,
        'mass number': 147,
        'mass': 146.949964 * u.u,
        'stable': False,
        'half-life': 3.2 * u.s,
    },

    'Er-148': {
        'atomic number': 68,
        'mass number': 148,
        'mass': 147.944735 * u.u,
        'stable': False,
        'half-life': 4.6 * u.s,
    },

    'Er-149': {
        'atomic number': 68,
        'mass number': 149,
        'mass': 148.942306 * u.u,
        'stable': False,
        'half-life': 4.0 * u.s,
    },

    'Er-150': {
        'atomic number': 68,
        'mass number': 150,
        'mass': 149.937916 * u.u,
        'stable': False,
        'half-life': 18.5 * u.s,
    },

    'Er-151': {
        'atomic number': 68,
        'mass number': 151,
        'mass': 150.937449 * u.u,
        'stable': False,
        'half-life': 23.5 * u.s,
    },

    'Er-152': {
        'atomic number': 68,
        'mass number': 152,
        'mass': 151.935057 * u.u,
        'stable': False,
        'half-life': 10.3 * u.s,
    },

    'Er-153': {
        'atomic number': 68,
        'mass number': 153,
        'mass': 152.93508 * u.u,
        'stable': False,
        'half-life': 37.1 * u.s,
    },

    'Er-154': {
        'atomic number': 68,
        'mass number': 154,
        'mass': 153.9327908 * u.u,
        'stable': False,
        'half-life': 223.8 * u.s,
    },

    'Er-155': {
        'atomic number': 68,
        'mass number': 155,
        'mass': 154.9332159 * u.u,
        'stable': False,
        'half-life': 318.0 * u.s,
    },

    'Er-156': {
        'atomic number': 68,
        'mass number': 156,
        'mass': 155.931067 * u.u,
        'stable': False,
        'half-life': 1170.0 * u.s,
    },

    'Er-157': {
        'atomic number': 68,
        'mass number': 157,
        'mass': 156.931949 * u.u,
        'stable': False,
        'half-life': 1119.0 * u.s,
    },

    'Er-158': {
        'atomic number': 68,
        'mass number': 158,
        'mass': 157.929893 * u.u,
        'stable': False,
        'half-life': 8244.0 * u.s,
    },

    'Er-159': {
        'atomic number': 68,
        'mass number': 159,
        'mass': 158.9306918 * u.u,
        'stable': False,
        'half-life': 2160.0 * u.s,
    },

    'Er-160': {
        'atomic number': 68,
        'mass number': 160,
        'mass': 159.929077 * u.u,
        'stable': False,
        'half-life': 102888.0 * u.s,
    },

    'Er-161': {
        'atomic number': 68,
        'mass number': 161,
        'mass': 160.9300046 * u.u,
        'stable': False,
        'half-life': 11556.0 * u.s,
    },

    'Er-162': {
        'atomic number': 68,
        'mass number': 162,
        'mass': 161.9287884 * u.u,
        'stable': True,
        'abundance': 0.00139,
    },

    'Er-163': {
        'atomic number': 68,
        'mass number': 163,
        'mass': 162.9300408 * u.u,
        'stable': False,
        'half-life': 4500.0 * u.s,
    },

    'Er-164': {
        'atomic number': 68,
        'mass number': 164,
        'mass': 163.9292088 * u.u,
        'stable': True,
        'abundance': 0.01601,
    },

    'Er-165': {
        'atomic number': 68,
        'mass number': 165,
        'mass': 164.9307345 * u.u,
        'stable': False,
        'half-life': 37296.0 * u.s,
    },

    'Er-166': {
        'atomic number': 68,
        'mass number': 166,
        'mass': 165.9302995 * u.u,
        'stable': True,
        'abundance': 0.33503,
    },

    'Er-167': {
        'atomic number': 68,
        'mass number': 167,
        'mass': 166.9320546 * u.u,
        'stable': True,
        'abundance': 0.22869,
    },

    'Er-168': {
        'atomic number': 68,
        'mass number': 168,
        'mass': 167.9323767 * u.u,
        'stable': True,
        'abundance': 0.26978,
    },

    'Er-169': {
        'atomic number': 68,
        'mass number': 169,
        'mass': 168.9345968 * u.u,
        'stable': False,
        'half-life': 811468.8 * u.s,
    },

    'Er-170': {
        'atomic number': 68,
        'mass number': 170,
        'mass': 169.9354702 * u.u,
        'stable': True,
        'abundance': 0.1491,
    },

    'Er-171': {
        'atomic number': 68,
        'mass number': 171,
        'mass': 170.9380357 * u.u,
        'stable': False,
        'half-life': 27057.6 * u.s,
    },

    'Er-172': {
        'atomic number': 68,
        'mass number': 172,
        'mass': 171.9393619 * u.u,
        'stable': False,
        'half-life': 177480.0 * u.s,
    },

    'Er-173': {
        'atomic number': 68,
        'mass number': 173,
        'mass': 172.9424 * u.u,
        'stable': False,
        'half-life': 86.04 * u.s,
    },

    'Er-174': {
        'atomic number': 68,
        'mass number': 174,
        'mass': 173.94423 * u.u,
        'stable': False,
        'half-life': 192.0 * u.s,
    },

    'Er-175': {
        'atomic number': 68,
        'mass number': 175,
        'mass': 174.94777 * u.u,
        'stable': False,
        'half-life': 72.0 * u.s,
    },

    'Er-176': {
        'atomic number': 68,
        'mass number': 176,
        'mass': 175.94994 * u.u,
        'stable': False,
        'half-life': '20# s',
    },

    'Er-177': {
        'atomic number': 68,
        'mass number': 177,
        'mass': 176.95399 * u.u,
        'stable': False,
        'half-life': '3# s',
    },

    'Tm-144': {
        'atomic number': 69,
        'mass number': 144,
        'mass': 143.97628 * u.u,
        'stable': False,
        'half-life': 2.3e-06 * u.s,
    },

    'Tm-145': {
        'atomic number': 69,
        'mass number': 145,
        'mass': 144.97039 * u.u,
        'stable': False,
        'half-life': 3.17e-06 * u.s,
    },

    'Tm-146': {
        'atomic number': 69,
        'mass number': 146,
        'mass': 145.96684 * u.u,
        'stable': False,
        'half-life': 0.155 * u.s,
    },

    'Tm-147': {
        'atomic number': 69,
        'mass number': 147,
        'mass': 146.9613799 * u.u,
        'stable': False,
        'half-life': 0.58 * u.s,
    },

    'Tm-148': {
        'atomic number': 69,
        'mass number': 148,
        'mass': 147.958384 * u.u,
        'stable': False,
        'half-life': 0.7 * u.s,
    },

    'Tm-149': {
        'atomic number': 69,
        'mass number': 149,
        'mass': 148.95289 * u.u,
        'stable': False,
        'half-life': 0.9 * u.s,
    },

    'Tm-150': {
        'atomic number': 69,
        'mass number': 150,
        'mass': 149.95009 * u.u,
        'stable': False,
        'half-life': '3# s',
    },

    'Tm-151': {
        'atomic number': 69,
        'mass number': 151,
        'mass': 150.945488 * u.u,
        'stable': False,
        'half-life': 4.17 * u.s,
    },

    'Tm-152': {
        'atomic number': 69,
        'mass number': 152,
        'mass': 151.944422 * u.u,
        'stable': False,
        'half-life': 8.0 * u.s,
    },

    'Tm-153': {
        'atomic number': 69,
        'mass number': 153,
        'mass': 152.94204 * u.u,
        'stable': False,
        'half-life': 1.48 * u.s,
    },

    'Tm-154': {
        'atomic number': 69,
        'mass number': 154,
        'mass': 153.94157 * u.u,
        'stable': False,
        'half-life': 8.1 * u.s,
    },

    'Tm-155': {
        'atomic number': 69,
        'mass number': 155,
        'mass': 154.93921 * u.u,
        'stable': False,
        'half-life': 21.6 * u.s,
    },

    'Tm-156': {
        'atomic number': 69,
        'mass number': 156,
        'mass': 155.938992 * u.u,
        'stable': False,
        'half-life': 83.8 * u.s,
    },

    'Tm-157': {
        'atomic number': 69,
        'mass number': 157,
        'mass': 156.936944 * u.u,
        'stable': False,
        'half-life': 217.8 * u.s,
    },

    'Tm-158': {
        'atomic number': 69,
        'mass number': 158,
        'mass': 157.93698 * u.u,
        'stable': False,
        'half-life': 238.8 * u.s,
    },

    'Tm-159': {
        'atomic number': 69,
        'mass number': 159,
        'mass': 158.934975 * u.u,
        'stable': False,
        'half-life': 547.8 * u.s,
    },

    'Tm-160': {
        'atomic number': 69,
        'mass number': 160,
        'mass': 159.935263 * u.u,
        'stable': False,
        'half-life': 564.0 * u.s,
    },

    'Tm-161': {
        'atomic number': 69,
        'mass number': 161,
        'mass': 160.933549 * u.u,
        'stable': False,
        'half-life': 1812.0 * u.s,
    },

    'Tm-162': {
        'atomic number': 69,
        'mass number': 162,
        'mass': 161.934002 * u.u,
        'stable': False,
        'half-life': 1302.0 * u.s,
    },

    'Tm-163': {
        'atomic number': 69,
        'mass number': 163,
        'mass': 162.9326592 * u.u,
        'stable': False,
        'half-life': 6516.0 * u.s,
    },

    'Tm-164': {
        'atomic number': 69,
        'mass number': 164,
        'mass': 163.933544 * u.u,
        'stable': False,
        'half-life': 120.0 * u.s,
    },

    'Tm-165': {
        'atomic number': 69,
        'mass number': 165,
        'mass': 164.9324431 * u.u,
        'stable': False,
        'half-life': 108216.0 * u.s,
    },

    'Tm-166': {
        'atomic number': 69,
        'mass number': 166,
        'mass': 165.933561 * u.u,
        'stable': False,
        'half-life': 27720.0 * u.s,
    },

    'Tm-167': {
        'atomic number': 69,
        'mass number': 167,
        'mass': 166.9328562 * u.u,
        'stable': False,
        'half-life': 799200.0 * u.s,
    },

    'Tm-168': {
        'atomic number': 69,
        'mass number': 168,
        'mass': 167.9341774 * u.u,
        'stable': False,
        'half-life': 8043840.0 * u.s,
    },

    'Tm-169': {
        'atomic number': 69,
        'mass number': 169,
        'mass': 168.9342179 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Tm-170': {
        'atomic number': 69,
        'mass number': 170,
        'mass': 169.935806 * u.u,
        'stable': False,
        'half-life': 11111040.0 * u.s,
    },

    'Tm-171': {
        'atomic number': 69,
        'mass number': 171,
        'mass': 170.9364339 * u.u,
        'stable': False,
        'half-life': 60589297.92 * u.s,
    },

    'Tm-172': {
        'atomic number': 69,
        'mass number': 172,
        'mass': 171.9384055 * u.u,
        'stable': False,
        'half-life': 228960.0 * u.s,
    },

    'Tm-173': {
        'atomic number': 69,
        'mass number': 173,
        'mass': 172.9396084 * u.u,
        'stable': False,
        'half-life': 29664.0 * u.s,
    },

    'Tm-174': {
        'atomic number': 69,
        'mass number': 174,
        'mass': 173.942173 * u.u,
        'stable': False,
        'half-life': 324.0 * u.s,
    },

    'Tm-175': {
        'atomic number': 69,
        'mass number': 175,
        'mass': 174.943841 * u.u,
        'stable': False,
        'half-life': 912.0 * u.s,
    },

    'Tm-176': {
        'atomic number': 69,
        'mass number': 176,
        'mass': 175.947 * u.u,
        'stable': False,
        'half-life': 111.0 * u.s,
    },

    'Tm-177': {
        'atomic number': 69,
        'mass number': 177,
        'mass': 176.94904 * u.u,
        'stable': False,
        'half-life': 90.0 * u.s,
    },

    'Tm-178': {
        'atomic number': 69,
        'mass number': 178,
        'mass': 177.95264 * u.u,
        'stable': False,
        'half-life': '30# s',
    },

    'Tm-179': {
        'atomic number': 69,
        'mass number': 179,
        'mass': 178.95534 * u.u,
        'stable': False,
        'half-life': '20# s',
    },

    'Yb-148': {
        'atomic number': 70,
        'mass number': 148,
        'mass': 147.96758 * u.u,
        'stable': False,
        'half-life': '250# ms',
    },

    'Yb-149': {
        'atomic number': 70,
        'mass number': 149,
        'mass': 148.96436 * u.u,
        'stable': False,
        'half-life': 0.7 * u.s,
    },

    'Yb-150': {
        'atomic number': 70,
        'mass number': 150,
        'mass': 149.95852 * u.u,
        'stable': False,
        'half-life': '700# ms',
    },

    'Yb-151': {
        'atomic number': 70,
        'mass number': 151,
        'mass': 150.9554 * u.u,
        'stable': False,
        'half-life': 1.6 * u.s,
    },

    'Yb-152': {
        'atomic number': 70,
        'mass number': 152,
        'mass': 151.95027 * u.u,
        'stable': False,
        'half-life': 3.03 * u.s,
    },

    'Yb-153': {
        'atomic number': 70,
        'mass number': 153,
        'mass': 152.94932 * u.u,
        'stable': False,
        'half-life': 4.2 * u.s,
    },

    'Yb-154': {
        'atomic number': 70,
        'mass number': 154,
        'mass': 153.946396 * u.u,
        'stable': False,
        'half-life': 0.409 * u.s,
    },

    'Yb-155': {
        'atomic number': 70,
        'mass number': 155,
        'mass': 154.945783 * u.u,
        'stable': False,
        'half-life': 1.793 * u.s,
    },

    'Yb-156': {
        'atomic number': 70,
        'mass number': 156,
        'mass': 155.942825 * u.u,
        'stable': False,
        'half-life': 26.1 * u.s,
    },

    'Yb-157': {
        'atomic number': 70,
        'mass number': 157,
        'mass': 156.942645 * u.u,
        'stable': False,
        'half-life': 38.6 * u.s,
    },

    'Yb-158': {
        'atomic number': 70,
        'mass number': 158,
        'mass': 157.9398705 * u.u,
        'stable': False,
        'half-life': 89.4 * u.s,
    },

    'Yb-159': {
        'atomic number': 70,
        'mass number': 159,
        'mass': 158.940055 * u.u,
        'stable': False,
        'half-life': 100.2 * u.s,
    },

    'Yb-160': {
        'atomic number': 70,
        'mass number': 160,
        'mass': 159.937557 * u.u,
        'stable': False,
        'half-life': 288.0 * u.s,
    },

    'Yb-161': {
        'atomic number': 70,
        'mass number': 161,
        'mass': 160.937907 * u.u,
        'stable': False,
        'half-life': 252.0 * u.s,
    },

    'Yb-162': {
        'atomic number': 70,
        'mass number': 162,
        'mass': 161.935774 * u.u,
        'stable': False,
        'half-life': 1132.2 * u.s,
    },

    'Yb-163': {
        'atomic number': 70,
        'mass number': 163,
        'mass': 162.93634 * u.u,
        'stable': False,
        'half-life': 663.0 * u.s,
    },

    'Yb-164': {
        'atomic number': 70,
        'mass number': 164,
        'mass': 163.934495 * u.u,
        'stable': False,
        'half-life': 4548.0 * u.s,
    },

    'Yb-165': {
        'atomic number': 70,
        'mass number': 165,
        'mass': 164.93527 * u.u,
        'stable': False,
        'half-life': 594.0 * u.s,
    },

    'Yb-166': {
        'atomic number': 70,
        'mass number': 166,
        'mass': 165.9338747 * u.u,
        'stable': False,
        'half-life': 204120.0 * u.s,
    },

    'Yb-167': {
        'atomic number': 70,
        'mass number': 167,
        'mass': 166.934953 * u.u,
        'stable': False,
        'half-life': 1050.0 * u.s,
    },

    'Yb-168': {
        'atomic number': 70,
        'mass number': 168,
        'mass': 167.9338896 * u.u,
        'stable': True,
        'abundance': 0.00123,
    },

    'Yb-169': {
        'atomic number': 70,
        'mass number': 169,
        'mass': 168.9351825 * u.u,
        'stable': False,
        'half-life': 2766070.0799999996 * u.s,
    },

    'Yb-170': {
        'atomic number': 70,
        'mass number': 170,
        'mass': 169.9347664 * u.u,
        'stable': True,
        'abundance': 0.02982,
    },

    'Yb-171': {
        'atomic number': 70,
        'mass number': 171,
        'mass': 170.9363302 * u.u,
        'stable': True,
        'abundance': 0.1409,
    },

    'Yb-172': {
        'atomic number': 70,
        'mass number': 172,
        'mass': 171.9363859 * u.u,
        'stable': True,
        'abundance': 0.2168,
    },

    'Yb-173': {
        'atomic number': 70,
        'mass number': 173,
        'mass': 172.9382151 * u.u,
        'stable': True,
        'abundance': 0.16103,
    },

    'Yb-174': {
        'atomic number': 70,
        'mass number': 174,
        'mass': 173.9388664 * u.u,
        'stable': True,
        'abundance': 0.32026,
    },

    'Yb-175': {
        'atomic number': 70,
        'mass number': 175,
        'mass': 174.9412808 * u.u,
        'stable': False,
        'half-life': 361584.0 * u.s,
    },

    'Yb-176': {
        'atomic number': 70,
        'mass number': 176,
        'mass': 175.9425764 * u.u,
        'stable': True,
        'abundance': 0.12996,
    },

    'Yb-177': {
        'atomic number': 70,
        'mass number': 177,
        'mass': 176.9452656 * u.u,
        'stable': False,
        'half-life': 6879.6 * u.s,
    },

    'Yb-178': {
        'atomic number': 70,
        'mass number': 178,
        'mass': 177.946651 * u.u,
        'stable': False,
        'half-life': 4440.0 * u.s,
    },

    'Yb-179': {
        'atomic number': 70,
        'mass number': 179,
        'mass': 178.95004 * u.u,
        'stable': False,
        'half-life': 480.0 * u.s,
    },

    'Yb-180': {
        'atomic number': 70,
        'mass number': 180,
        'mass': 179.95212 * u.u,
        'stable': False,
        'half-life': 144.0 * u.s,
    },

    'Yb-181': {
        'atomic number': 70,
        'mass number': 181,
        'mass': 180.95589 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Lu-150': {
        'atomic number': 71,
        'mass number': 150,
        'mass': 149.97355 * u.u,
        'stable': False,
        'half-life': 0.045 * u.s,
    },

    'Lu-151': {
        'atomic number': 71,
        'mass number': 151,
        'mass': 150.96768 * u.u,
        'stable': False,
        'half-life': 0.0784 * u.s,
    },

    'Lu-152': {
        'atomic number': 71,
        'mass number': 152,
        'mass': 151.96412 * u.u,
        'stable': False,
        'half-life': 0.65 * u.s,
    },

    'Lu-153': {
        'atomic number': 71,
        'mass number': 153,
        'mass': 152.95875 * u.u,
        'stable': False,
        'half-life': 0.9 * u.s,
    },

    'Lu-154': {
        'atomic number': 71,
        'mass number': 154,
        'mass': 153.95736 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Lu-155': {
        'atomic number': 71,
        'mass number': 155,
        'mass': 154.954321 * u.u,
        'stable': False,
        'half-life': 0.0686 * u.s,
    },

    'Lu-156': {
        'atomic number': 71,
        'mass number': 156,
        'mass': 155.953033 * u.u,
        'stable': False,
        'half-life': 0.494 * u.s,
    },

    'Lu-157': {
        'atomic number': 71,
        'mass number': 157,
        'mass': 156.950127 * u.u,
        'stable': False,
        'half-life': 6.8 * u.s,
    },

    'Lu-158': {
        'atomic number': 71,
        'mass number': 158,
        'mass': 157.949316 * u.u,
        'stable': False,
        'half-life': 10.6 * u.s,
    },

    'Lu-159': {
        'atomic number': 71,
        'mass number': 159,
        'mass': 158.946636 * u.u,
        'stable': False,
        'half-life': 12.1 * u.s,
    },

    'Lu-160': {
        'atomic number': 71,
        'mass number': 160,
        'mass': 159.946033 * u.u,
        'stable': False,
        'half-life': 36.1 * u.s,
    },

    'Lu-161': {
        'atomic number': 71,
        'mass number': 161,
        'mass': 160.943572 * u.u,
        'stable': False,
        'half-life': 77.0 * u.s,
    },

    'Lu-162': {
        'atomic number': 71,
        'mass number': 162,
        'mass': 161.943283 * u.u,
        'stable': False,
        'half-life': 82.2 * u.s,
    },

    'Lu-163': {
        'atomic number': 71,
        'mass number': 163,
        'mass': 162.941179 * u.u,
        'stable': False,
        'half-life': 238.2 * u.s,
    },

    'Lu-164': {
        'atomic number': 71,
        'mass number': 164,
        'mass': 163.941339 * u.u,
        'stable': False,
        'half-life': 188.4 * u.s,
    },

    'Lu-165': {
        'atomic number': 71,
        'mass number': 165,
        'mass': 164.939407 * u.u,
        'stable': False,
        'half-life': 644.4 * u.s,
    },

    'Lu-166': {
        'atomic number': 71,
        'mass number': 166,
        'mass': 165.939859 * u.u,
        'stable': False,
        'half-life': 159.0 * u.s,
    },

    'Lu-167': {
        'atomic number': 71,
        'mass number': 167,
        'mass': 166.93827 * u.u,
        'stable': False,
        'half-life': 3090.0 * u.s,
    },

    'Lu-168': {
        'atomic number': 71,
        'mass number': 168,
        'mass': 167.938736 * u.u,
        'stable': False,
        'half-life': 330.0 * u.s,
    },

    'Lu-169': {
        'atomic number': 71,
        'mass number': 169,
        'mass': 168.9376441 * u.u,
        'stable': False,
        'half-life': 122616.0 * u.s,
    },

    'Lu-170': {
        'atomic number': 71,
        'mass number': 170,
        'mass': 169.938478 * u.u,
        'stable': False,
        'half-life': 173836.8 * u.s,
    },

    'Lu-171': {
        'atomic number': 71,
        'mass number': 171,
        'mass': 170.937917 * u.u,
        'stable': False,
        'half-life': 711936.0 * u.s,
    },

    'Lu-172': {
        'atomic number': 71,
        'mass number': 172,
        'mass': 171.9390891 * u.u,
        'stable': False,
        'half-life': 578880.0 * u.s,
    },

    'Lu-173': {
        'atomic number': 71,
        'mass number': 173,
        'mass': 172.938934 * u.u,
        'stable': False,
        'half-life': 43232988.62 * u.s,
    },

    'Lu-174': {
        'atomic number': 71,
        'mass number': 174,
        'mass': 173.9403409 * u.u,
        'stable': False,
        'half-life': 104453425.06 * u.s,
    },

    'Lu-175': {
        'atomic number': 71,
        'mass number': 175,
        'mass': 174.9407752 * u.u,
        'stable': True,
        'abundance': 0.97401,
    },

    'Lu-176': {
        'atomic number': 71,
        'mass number': 176,
        'mass': 175.9426897 * u.u,
        'stable': False,
        'abundance': 0.02599,
    },

    'Lu-177': {
        'atomic number': 71,
        'mass number': 177,
        'mass': 176.9437615 * u.u,
        'stable': False,
        'half-life': 573696.0 * u.s,
    },

    'Lu-178': {
        'atomic number': 71,
        'mass number': 178,
        'mass': 177.945958 * u.u,
        'stable': False,
        'half-life': 1704.0 * u.s,
    },

    'Lu-179': {
        'atomic number': 71,
        'mass number': 179,
        'mass': 178.9473309 * u.u,
        'stable': False,
        'half-life': 16524.0 * u.s,
    },

    'Lu-180': {
        'atomic number': 71,
        'mass number': 180,
        'mass': 179.949888 * u.u,
        'stable': False,
        'half-life': 342.0 * u.s,
    },

    'Lu-181': {
        'atomic number': 71,
        'mass number': 181,
        'mass': 180.95191 * u.u,
        'stable': False,
        'half-life': 210.0 * u.s,
    },

    'Lu-182': {
        'atomic number': 71,
        'mass number': 182,
        'mass': 181.95504 * u.u,
        'stable': False,
        'half-life': 120.0 * u.s,
    },

    'Lu-183': {
        'atomic number': 71,
        'mass number': 183,
        'mass': 182.957363 * u.u,
        'stable': False,
        'half-life': 58.0 * u.s,
    },

    'Lu-184': {
        'atomic number': 71,
        'mass number': 184,
        'mass': 183.96091 * u.u,
        'stable': False,
        'half-life': 20.0 * u.s,
    },

    'Lu-185': {
        'atomic number': 71,
        'mass number': 185,
        'mass': 184.96362 * u.u,
        'stable': False,
        'half-life': '6# s',
    },

    'Hf-153': {
        'atomic number': 72,
        'mass number': 153,
        'mass': 152.97069 * u.u,
        'stable': False,
        'half-life': '400# ms',
    },

    'Hf-154': {
        'atomic number': 72,
        'mass number': 154,
        'mass': 153.96486 * u.u,
        'stable': False,
        'half-life': 2.0 * u.s,
    },

    'Hf-155': {
        'atomic number': 72,
        'mass number': 155,
        'mass': 154.96311 * u.u,
        'stable': False,
        'half-life': 0.84 * u.s,
    },

    'Hf-156': {
        'atomic number': 72,
        'mass number': 156,
        'mass': 155.95935 * u.u,
        'stable': False,
        'half-life': 0.023 * u.s,
    },

    'Hf-157': {
        'atomic number': 72,
        'mass number': 157,
        'mass': 156.95824 * u.u,
        'stable': False,
        'half-life': 0.115 * u.s,
    },

    'Hf-158': {
        'atomic number': 72,
        'mass number': 158,
        'mass': 157.954801 * u.u,
        'stable': False,
        'half-life': 0.99 * u.s,
    },

    'Hf-159': {
        'atomic number': 72,
        'mass number': 159,
        'mass': 158.953996 * u.u,
        'stable': False,
        'half-life': 5.2 * u.s,
    },

    'Hf-160': {
        'atomic number': 72,
        'mass number': 160,
        'mass': 159.950691 * u.u,
        'stable': False,
        'half-life': 13.6 * u.s,
    },

    'Hf-161': {
        'atomic number': 72,
        'mass number': 161,
        'mass': 160.950278 * u.u,
        'stable': False,
        'half-life': 18.4 * u.s,
    },

    'Hf-162': {
        'atomic number': 72,
        'mass number': 162,
        'mass': 161.9472148 * u.u,
        'stable': False,
        'half-life': 39.4 * u.s,
    },

    'Hf-163': {
        'atomic number': 72,
        'mass number': 163,
        'mass': 162.947113 * u.u,
        'stable': False,
        'half-life': 40.0 * u.s,
    },

    'Hf-164': {
        'atomic number': 72,
        'mass number': 164,
        'mass': 163.944371 * u.u,
        'stable': False,
        'half-life': 111.0 * u.s,
    },

    'Hf-165': {
        'atomic number': 72,
        'mass number': 165,
        'mass': 164.944567 * u.u,
        'stable': False,
        'half-life': 76.0 * u.s,
    },

    'Hf-166': {
        'atomic number': 72,
        'mass number': 166,
        'mass': 165.94218 * u.u,
        'stable': False,
        'half-life': 406.2 * u.s,
    },

    'Hf-167': {
        'atomic number': 72,
        'mass number': 167,
        'mass': 166.9426 * u.u,
        'stable': False,
        'half-life': 123.0 * u.s,
    },

    'Hf-168': {
        'atomic number': 72,
        'mass number': 168,
        'mass': 167.940568 * u.u,
        'stable': False,
        'half-life': 1557.0 * u.s,
    },

    'Hf-169': {
        'atomic number': 72,
        'mass number': 169,
        'mass': 168.941259 * u.u,
        'stable': False,
        'half-life': 194.4 * u.s,
    },

    'Hf-170': {
        'atomic number': 72,
        'mass number': 170,
        'mass': 169.939609 * u.u,
        'stable': False,
        'half-life': 57636.0 * u.s,
    },

    'Hf-171': {
        'atomic number': 72,
        'mass number': 171,
        'mass': 170.940492 * u.u,
        'stable': False,
        'half-life': 43560.0 * u.s,
    },

    'Hf-172': {
        'atomic number': 72,
        'mass number': 172,
        'mass': 171.93945 * u.u,
        'stable': False,
        'half-life': 59011451.62 * u.s,
    },

    'Hf-173': {
        'atomic number': 72,
        'mass number': 173,
        'mass': 172.940513 * u.u,
        'stable': False,
        'half-life': 84960.0 * u.s,
    },

    'Hf-174': {
        'atomic number': 72,
        'mass number': 174,
        'mass': 173.9400461 * u.u,
        'stable': False,
        'abundance': 0.0016,
    },

    'Hf-175': {
        'atomic number': 72,
        'mass number': 175,
        'mass': 174.9415092 * u.u,
        'stable': False,
        'half-life': 6104160.0 * u.s,
    },

    'Hf-176': {
        'atomic number': 72,
        'mass number': 176,
        'mass': 175.9414076 * u.u,
        'stable': True,
        'abundance': 0.0526,
    },

    'Hf-177': {
        'atomic number': 72,
        'mass number': 177,
        'mass': 176.9432277 * u.u,
        'stable': True,
        'abundance': 0.186,
    },

    'Hf-178': {
        'atomic number': 72,
        'mass number': 178,
        'mass': 177.9437058 * u.u,
        'stable': True,
        'abundance': 0.2728,
    },

    'Hf-179': {
        'atomic number': 72,
        'mass number': 179,
        'mass': 178.9458232 * u.u,
        'stable': True,
        'abundance': 0.1362,
    },

    'Hf-180': {
        'atomic number': 72,
        'mass number': 180,
        'mass': 179.946557 * u.u,
        'stable': True,
        'abundance': 0.3508,
    },

    'Hf-181': {
        'atomic number': 72,
        'mass number': 181,
        'mass': 180.9491083 * u.u,
        'stable': False,
        'half-life': 3662496.0 * u.s,
    },

    'Hf-182': {
        'atomic number': 72,
        'mass number': 182,
        'mass': 181.9505612 * u.u,
        'stable': False,
        'half-life': 280856641400000.0 * u.s,
    },

    'Hf-183': {
        'atomic number': 72,
        'mass number': 183,
        'mass': 182.95353 * u.u,
        'stable': False,
        'half-life': 3664.8 * u.s,
    },

    'Hf-184': {
        'atomic number': 72,
        'mass number': 184,
        'mass': 183.955446 * u.u,
        'stable': False,
        'half-life': 14832.0 * u.s,
    },

    'Hf-185': {
        'atomic number': 72,
        'mass number': 185,
        'mass': 184.958862 * u.u,
        'stable': False,
        'half-life': 210.0 * u.s,
    },

    'Hf-186': {
        'atomic number': 72,
        'mass number': 186,
        'mass': 185.960897 * u.u,
        'stable': False,
        'half-life': 156.0 * u.s,
    },

    'Hf-187': {
        'atomic number': 72,
        'mass number': 187,
        'mass': 186.96477 * u.u,
        'stable': False,
        'half-life': '30# s',
    },

    'Hf-188': {
        'atomic number': 72,
        'mass number': 188,
        'mass': 187.96685 * u.u,
        'stable': False,
        'half-life': '20# s',
    },

    'Hf-189': {
        'atomic number': 72,
        'mass number': 189,
        'mass': 188.97084 * u.u,
        'stable': False,
        'half-life': '2# s',
    },

    'Ta-155': {
        'atomic number': 73,
        'mass number': 155,
        'mass': 154.97424 * u.u,
        'stable': False,
        'half-life': 0.0032 * u.s,
    },

    'Ta-156': {
        'atomic number': 73,
        'mass number': 156,
        'mass': 155.97203 * u.u,
        'stable': False,
        'half-life': 0.106 * u.s,
    },

    'Ta-157': {
        'atomic number': 73,
        'mass number': 157,
        'mass': 156.96818 * u.u,
        'stable': False,
        'half-life': 0.0101 * u.s,
    },

    'Ta-158': {
        'atomic number': 73,
        'mass number': 158,
        'mass': 157.96654 * u.u,
        'stable': False,
        'half-life': 0.049 * u.s,
    },

    'Ta-159': {
        'atomic number': 73,
        'mass number': 159,
        'mass': 158.963023 * u.u,
        'stable': False,
        'half-life': 1.04 * u.s,
    },

    'Ta-160': {
        'atomic number': 73,
        'mass number': 160,
        'mass': 159.961488 * u.u,
        'stable': False,
        'half-life': 1.7 * u.s,
    },

    'Ta-161': {
        'atomic number': 73,
        'mass number': 161,
        'mass': 160.958452 * u.u,
        'stable': False,
        'half-life': '3# s',
    },

    'Ta-162': {
        'atomic number': 73,
        'mass number': 162,
        'mass': 161.957294 * u.u,
        'stable': False,
        'half-life': 3.57 * u.s,
    },

    'Ta-163': {
        'atomic number': 73,
        'mass number': 163,
        'mass': 162.954337 * u.u,
        'stable': False,
        'half-life': 10.6 * u.s,
    },

    'Ta-164': {
        'atomic number': 73,
        'mass number': 164,
        'mass': 163.953534 * u.u,
        'stable': False,
        'half-life': 14.2 * u.s,
    },

    'Ta-165': {
        'atomic number': 73,
        'mass number': 165,
        'mass': 164.950781 * u.u,
        'stable': False,
        'half-life': 31.0 * u.s,
    },

    'Ta-166': {
        'atomic number': 73,
        'mass number': 166,
        'mass': 165.950512 * u.u,
        'stable': False,
        'half-life': 34.4 * u.s,
    },

    'Ta-167': {
        'atomic number': 73,
        'mass number': 167,
        'mass': 166.948093 * u.u,
        'stable': False,
        'half-life': 79.8 * u.s,
    },

    'Ta-168': {
        'atomic number': 73,
        'mass number': 168,
        'mass': 167.948047 * u.u,
        'stable': False,
        'half-life': 120.0 * u.s,
    },

    'Ta-169': {
        'atomic number': 73,
        'mass number': 169,
        'mass': 168.946011 * u.u,
        'stable': False,
        'half-life': 294.0 * u.s,
    },

    'Ta-170': {
        'atomic number': 73,
        'mass number': 170,
        'mass': 169.946175 * u.u,
        'stable': False,
        'half-life': 405.6 * u.s,
    },

    'Ta-171': {
        'atomic number': 73,
        'mass number': 171,
        'mass': 170.944476 * u.u,
        'stable': False,
        'half-life': 1398.0 * u.s,
    },

    'Ta-172': {
        'atomic number': 73,
        'mass number': 172,
        'mass': 171.944895 * u.u,
        'stable': False,
        'half-life': 2208.0 * u.s,
    },

    'Ta-173': {
        'atomic number': 73,
        'mass number': 173,
        'mass': 172.94375 * u.u,
        'stable': False,
        'half-life': 11304.0 * u.s,
    },

    'Ta-174': {
        'atomic number': 73,
        'mass number': 174,
        'mass': 173.944454 * u.u,
        'stable': False,
        'half-life': 4104.0 * u.s,
    },

    'Ta-175': {
        'atomic number': 73,
        'mass number': 175,
        'mass': 174.943737 * u.u,
        'stable': False,
        'half-life': 37800.0 * u.s,
    },

    'Ta-176': {
        'atomic number': 73,
        'mass number': 176,
        'mass': 175.944857 * u.u,
        'stable': False,
        'half-life': 29124.0 * u.s,
    },

    'Ta-177': {
        'atomic number': 73,
        'mass number': 177,
        'mass': 176.9444795 * u.u,
        'stable': False,
        'half-life': 203616.0 * u.s,
    },

    'Ta-178': {
        'atomic number': 73,
        'mass number': 178,
        'mass': 177.945678 * u.u,
        'stable': False,
        'half-life': 8496.0 * u.s,
    },

    'Ta-179': {
        'atomic number': 73,
        'mass number': 179,
        'mass': 178.9459366 * u.u,
        'stable': False,
        'half-life': 57433605.32 * u.s,
    },

    'Ta-180': {
        'atomic number': 73,
        'mass number': 180,
        'mass': 179.9474648 * u.u,
        'stable': True,
        'abundance': 0.0001201,
    },

    'Ta-181': {
        'atomic number': 73,
        'mass number': 181,
        'mass': 180.9479958 * u.u,
        'stable': True,
        'abundance': 0.9998799,
    },

    'Ta-182': {
        'atomic number': 73,
        'mass number': 182,
        'mass': 181.9501519 * u.u,
        'stable': False,
        'half-life': 9913536.0 * u.s,
    },

    'Ta-183': {
        'atomic number': 73,
        'mass number': 183,
        'mass': 182.9513726 * u.u,
        'stable': False,
        'half-life': 440640.0 * u.s,
    },

    'Ta-184': {
        'atomic number': 73,
        'mass number': 184,
        'mass': 183.954008 * u.u,
        'stable': False,
        'half-life': 31320.0 * u.s,
    },

    'Ta-185': {
        'atomic number': 73,
        'mass number': 185,
        'mass': 184.955559 * u.u,
        'stable': False,
        'half-life': 2964.0 * u.s,
    },

    'Ta-186': {
        'atomic number': 73,
        'mass number': 186,
        'mass': 185.958551 * u.u,
        'stable': False,
        'half-life': 630.0 * u.s,
    },

    'Ta-187': {
        'atomic number': 73,
        'mass number': 187,
        'mass': 186.960386 * u.u,
        'stable': False,
        'half-life': 138.0 * u.s,
    },

    'Ta-188': {
        'atomic number': 73,
        'mass number': 188,
        'mass': 187.963916 * u.u,
        'stable': False,
        'half-life': 19.6 * u.s,
    },

    'Ta-189': {
        'atomic number': 73,
        'mass number': 189,
        'mass': 188.96583 * u.u,
        'stable': False,
        'half-life': '3# s',
    },

    'Ta-190': {
        'atomic number': 73,
        'mass number': 190,
        'mass': 189.96939 * u.u,
        'stable': False,
        'half-life': 5.3 * u.s,
    },

    'Ta-191': {
        'atomic number': 73,
        'mass number': 191,
        'mass': 190.97156 * u.u,
        'stable': False,
        'half-life': '3# s',
    },

    'Ta-192': {
        'atomic number': 73,
        'mass number': 192,
        'mass': 191.97514 * u.u,
        'stable': False,
        'half-life': 2.2 * u.s,
    },

    'W-157': {
        'atomic number': 74,
        'mass number': 157,
        'mass': 156.97884 * u.u,
        'stable': False,
        'half-life': 0.275 * u.s,
    },

    'W-158': {
        'atomic number': 74,
        'mass number': 158,
        'mass': 157.97456 * u.u,
        'stable': False,
        'half-life': 0.00125 * u.s,
    },

    'W-159': {
        'atomic number': 74,
        'mass number': 159,
        'mass': 158.97264 * u.u,
        'stable': False,
        'half-life': 0.0082 * u.s,
    },

    'W-160': {
        'atomic number': 74,
        'mass number': 160,
        'mass': 159.96846 * u.u,
        'stable': False,
        'half-life': 0.09 * u.s,
    },

    'W-161': {
        'atomic number': 74,
        'mass number': 161,
        'mass': 160.9672 * u.u,
        'stable': False,
        'half-life': 0.409 * u.s,
    },

    'W-162': {
        'atomic number': 74,
        'mass number': 162,
        'mass': 161.963499 * u.u,
        'stable': False,
        'half-life': 1.19 * u.s,
    },

    'W-163': {
        'atomic number': 74,
        'mass number': 163,
        'mass': 162.962524 * u.u,
        'stable': False,
        'half-life': 2.63 * u.s,
    },

    'W-164': {
        'atomic number': 74,
        'mass number': 164,
        'mass': 163.958961 * u.u,
        'stable': False,
        'half-life': 6.3 * u.s,
    },

    'W-165': {
        'atomic number': 74,
        'mass number': 165,
        'mass': 164.958281 * u.u,
        'stable': False,
        'half-life': 5.1 * u.s,
    },

    'W-166': {
        'atomic number': 74,
        'mass number': 166,
        'mass': 165.955031 * u.u,
        'stable': False,
        'half-life': 19.2 * u.s,
    },

    'W-167': {
        'atomic number': 74,
        'mass number': 167,
        'mass': 166.954805 * u.u,
        'stable': False,
        'half-life': 19.9 * u.s,
    },

    'W-168': {
        'atomic number': 74,
        'mass number': 168,
        'mass': 167.951806 * u.u,
        'stable': False,
        'half-life': 50.9 * u.s,
    },

    'W-169': {
        'atomic number': 74,
        'mass number': 169,
        'mass': 168.951779 * u.u,
        'stable': False,
        'half-life': 74.0 * u.s,
    },

    'W-170': {
        'atomic number': 74,
        'mass number': 170,
        'mass': 169.949232 * u.u,
        'stable': False,
        'half-life': 145.2 * u.s,
    },

    'W-171': {
        'atomic number': 74,
        'mass number': 171,
        'mass': 170.949451 * u.u,
        'stable': False,
        'half-life': 142.8 * u.s,
    },

    'W-172': {
        'atomic number': 74,
        'mass number': 172,
        'mass': 171.947292 * u.u,
        'stable': False,
        'half-life': 396.0 * u.s,
    },

    'W-173': {
        'atomic number': 74,
        'mass number': 173,
        'mass': 172.947689 * u.u,
        'stable': False,
        'half-life': 456.0 * u.s,
    },

    'W-174': {
        'atomic number': 74,
        'mass number': 174,
        'mass': 173.946079 * u.u,
        'stable': False,
        'half-life': 1992.0 * u.s,
    },

    'W-175': {
        'atomic number': 74,
        'mass number': 175,
        'mass': 174.946717 * u.u,
        'stable': False,
        'half-life': 2112.0 * u.s,
    },

    'W-176': {
        'atomic number': 74,
        'mass number': 176,
        'mass': 175.945634 * u.u,
        'stable': False,
        'half-life': 9000.0 * u.s,
    },

    'W-177': {
        'atomic number': 74,
        'mass number': 177,
        'mass': 176.946643 * u.u,
        'stable': False,
        'half-life': 7920.0 * u.s,
    },

    'W-178': {
        'atomic number': 74,
        'mass number': 178,
        'mass': 177.945883 * u.u,
        'stable': False,
        'half-life': 1866240.0 * u.s,
    },

    'W-179': {
        'atomic number': 74,
        'mass number': 179,
        'mass': 178.947077 * u.u,
        'stable': False,
        'half-life': 2223.0 * u.s,
    },

    'W-180': {
        'atomic number': 74,
        'mass number': 180,
        'mass': 179.9467108 * u.u,
        'stable': False,
        'abundance': 0.0012,
    },

    'W-181': {
        'atomic number': 74,
        'mass number': 181,
        'mass': 180.9481978 * u.u,
        'stable': False,
        'half-life': 10462608.0 * u.s,
    },

    'W-182': {
        'atomic number': 74,
        'mass number': 182,
        'mass': 181.94820394 * u.u,
        'stable': True,
        'abundance': 0.265,
    },

    'W-183': {
        'atomic number': 74,
        'mass number': 183,
        'mass': 182.95022275 * u.u,
        'stable': True,
        'abundance': 0.1431,
    },

    'W-184': {
        'atomic number': 74,
        'mass number': 184,
        'mass': 183.95093092 * u.u,
        'stable': True,
        'abundance': 0.3064,
    },

    'W-185': {
        'atomic number': 74,
        'mass number': 185,
        'mass': 184.95341897 * u.u,
        'stable': False,
        'half-life': 6488640.0 * u.s,
    },

    'W-186': {
        'atomic number': 74,
        'mass number': 186,
        'mass': 185.9543628 * u.u,
        'stable': True,
        'abundance': 0.2843,
    },

    'W-187': {
        'atomic number': 74,
        'mass number': 187,
        'mass': 186.9571588 * u.u,
        'stable': False,
        'half-life': 86400.0 * u.s,
    },

    'W-188': {
        'atomic number': 74,
        'mass number': 188,
        'mass': 187.9584862 * u.u,
        'stable': False,
        'half-life': 6029251.2 * u.s,
    },

    'W-189': {
        'atomic number': 74,
        'mass number': 189,
        'mass': 188.961763 * u.u,
        'stable': False,
        'half-life': 642.0 * u.s,
    },

    'W-190': {
        'atomic number': 74,
        'mass number': 190,
        'mass': 189.963091 * u.u,
        'stable': False,
        'half-life': 1800.0 * u.s,
    },

    'W-191': {
        'atomic number': 74,
        'mass number': 191,
        'mass': 190.966531 * u.u,
        'stable': False,
        'half-life': '45# s',
    },

    'W-192': {
        'atomic number': 74,
        'mass number': 192,
        'mass': 191.96817 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'W-193': {
        'atomic number': 74,
        'mass number': 193,
        'mass': 192.97178 * u.u,
        'stable': False,
        'half-life': '3# s',
    },

    'W-194': {
        'atomic number': 74,
        'mass number': 194,
        'mass': 193.97367 * u.u,
        'stable': False,
        'half-life': '5# s',
    },

    'Re-159': {
        'atomic number': 75,
        'mass number': 159,
        'mass': 158.98418 * u.u,
        'stable': False,
        'half-life': '40# us',
    },

    'Re-160': {
        'atomic number': 75,
        'mass number': 160,
        'mass': 159.98182 * u.u,
        'stable': False,
        'half-life': 0.000611 * u.s,
    },

    'Re-161': {
        'atomic number': 75,
        'mass number': 161,
        'mass': 160.97757 * u.u,
        'stable': False,
        'half-life': 0.00044 * u.s,
    },

    'Re-162': {
        'atomic number': 75,
        'mass number': 162,
        'mass': 161.97584 * u.u,
        'stable': False,
        'half-life': 0.107 * u.s,
    },

    'Re-163': {
        'atomic number': 75,
        'mass number': 163,
        'mass': 162.97208 * u.u,
        'stable': False,
        'half-life': 0.39 * u.s,
    },

    'Re-164': {
        'atomic number': 75,
        'mass number': 164,
        'mass': 163.970453 * u.u,
        'stable': False,
        'half-life': 0.719 * u.s,
    },

    'Re-165': {
        'atomic number': 75,
        'mass number': 165,
        'mass': 164.967103 * u.u,
        'stable': False,
        'half-life': 2.62 * u.s,
    },

    'Re-166': {
        'atomic number': 75,
        'mass number': 166,
        'mass': 165.965761 * u.u,
        'stable': False,
        'half-life': 2.25 * u.s,
    },

    'Re-167': {
        'atomic number': 75,
        'mass number': 167,
        'mass': 166.962595 * u.u,
        'stable': False,
        'half-life': 3.4 * u.s,
    },

    'Re-168': {
        'atomic number': 75,
        'mass number': 168,
        'mass': 167.961573 * u.u,
        'stable': False,
        'half-life': 4.4 * u.s,
    },

    'Re-169': {
        'atomic number': 75,
        'mass number': 169,
        'mass': 168.958766 * u.u,
        'stable': False,
        'half-life': 8.1 * u.s,
    },

    'Re-170': {
        'atomic number': 75,
        'mass number': 170,
        'mass': 169.95822 * u.u,
        'stable': False,
        'half-life': 9.2 * u.s,
    },

    'Re-171': {
        'atomic number': 75,
        'mass number': 171,
        'mass': 170.955716 * u.u,
        'stable': False,
        'half-life': 15.2 * u.s,
    },

    'Re-172': {
        'atomic number': 75,
        'mass number': 172,
        'mass': 171.95542 * u.u,
        'stable': False,
        'half-life': 15.0 * u.s,
    },

    'Re-173': {
        'atomic number': 75,
        'mass number': 173,
        'mass': 172.953243 * u.u,
        'stable': False,
        'half-life': 120.0 * u.s,
    },

    'Re-174': {
        'atomic number': 75,
        'mass number': 174,
        'mass': 173.953115 * u.u,
        'stable': False,
        'half-life': 144.0 * u.s,
    },

    'Re-175': {
        'atomic number': 75,
        'mass number': 175,
        'mass': 174.951381 * u.u,
        'stable': False,
        'half-life': 353.4 * u.s,
    },

    'Re-176': {
        'atomic number': 75,
        'mass number': 176,
        'mass': 175.951623 * u.u,
        'stable': False,
        'half-life': 318.0 * u.s,
    },

    'Re-177': {
        'atomic number': 75,
        'mass number': 177,
        'mass': 176.950328 * u.u,
        'stable': False,
        'half-life': 840.0 * u.s,
    },

    'Re-178': {
        'atomic number': 75,
        'mass number': 178,
        'mass': 177.950989 * u.u,
        'stable': False,
        'half-life': 792.0 * u.s,
    },

    'Re-179': {
        'atomic number': 75,
        'mass number': 179,
        'mass': 178.949989 * u.u,
        'stable': False,
        'half-life': 1170.0 * u.s,
    },

    'Re-180': {
        'atomic number': 75,
        'mass number': 180,
        'mass': 179.950792 * u.u,
        'stable': False,
        'half-life': 147.6 * u.s,
    },

    'Re-181': {
        'atomic number': 75,
        'mass number': 181,
        'mass': 180.950058 * u.u,
        'stable': False,
        'half-life': 71640.0 * u.s,
    },

    'Re-182': {
        'atomic number': 75,
        'mass number': 182,
        'mass': 181.95121 * u.u,
        'stable': False,
        'half-life': 231120.0 * u.s,
    },

    'Re-183': {
        'atomic number': 75,
        'mass number': 183,
        'mass': 182.9508196 * u.u,
        'stable': False,
        'half-life': 6048000.0 * u.s,
    },

    'Re-184': {
        'atomic number': 75,
        'mass number': 184,
        'mass': 183.9525228 * u.u,
        'stable': False,
        'half-life': 3058560.0 * u.s,
    },

    'Re-185': {
        'atomic number': 75,
        'mass number': 185,
        'mass': 184.9529545 * u.u,
        'stable': True,
        'abundance': 0.374,
    },

    'Re-186': {
        'atomic number': 75,
        'mass number': 186,
        'mass': 185.9549856 * u.u,
        'stable': False,
        'half-life': 321292.8 * u.s,
    },

    'Re-187': {
        'atomic number': 75,
        'mass number': 187,
        'mass': 186.9557501 * u.u,
        'stable': False,
        'abundance': 0.626,
    },

    'Re-188': {
        'atomic number': 75,
        'mass number': 188,
        'mass': 187.9581115 * u.u,
        'stable': False,
        'half-life': 61203.600000000006 * u.s,
    },

    'Re-189': {
        'atomic number': 75,
        'mass number': 189,
        'mass': 188.959226 * u.u,
        'stable': False,
        'half-life': 87480.0 * u.s,
    },

    'Re-190': {
        'atomic number': 75,
        'mass number': 190,
        'mass': 189.961744 * u.u,
        'stable': False,
        'half-life': 186.0 * u.s,
    },

    'Re-191': {
        'atomic number': 75,
        'mass number': 191,
        'mass': 190.963122 * u.u,
        'stable': False,
        'half-life': 588.0 * u.s,
    },

    'Re-192': {
        'atomic number': 75,
        'mass number': 192,
        'mass': 191.966088 * u.u,
        'stable': False,
        'half-life': 16.0 * u.s,
    },

    'Re-193': {
        'atomic number': 75,
        'mass number': 193,
        'mass': 192.967541 * u.u,
        'stable': False,
        'half-life': '20# s',
    },

    'Re-194': {
        'atomic number': 75,
        'mass number': 194,
        'mass': 193.97076 * u.u,
        'stable': False,
        'half-life': 5.0 * u.s,
    },

    'Re-195': {
        'atomic number': 75,
        'mass number': 195,
        'mass': 194.97254 * u.u,
        'stable': False,
        'half-life': 6.0 * u.s,
    },

    'Re-196': {
        'atomic number': 75,
        'mass number': 196,
        'mass': 195.9758 * u.u,
        'stable': False,
        'half-life': 2.4 * u.s,
    },

    'Re-197': {
        'atomic number': 75,
        'mass number': 197,
        'mass': 196.97799 * u.u,
        'stable': False,
        'half-life': '300# ms',
    },

    'Re-198': {
        'atomic number': 75,
        'mass number': 198,
        'mass': 197.9816 * u.u,
        'stable': False,
        'half-life': '300# ms',
    },

    'Os-161': {
        'atomic number': 76,
        'mass number': 161,
        'mass': 160.98903 * u.u,
        'stable': False,
        'half-life': 0.00064 * u.s,
    },

    'Os-162': {
        'atomic number': 76,
        'mass number': 162,
        'mass': 161.98443 * u.u,
        'stable': False,
        'half-life': 0.0021 * u.s,
    },

    'Os-163': {
        'atomic number': 76,
        'mass number': 163,
        'mass': 162.98241 * u.u,
        'stable': False,
        'half-life': 0.0055 * u.s,
    },

    'Os-164': {
        'atomic number': 76,
        'mass number': 164,
        'mass': 163.97802 * u.u,
        'stable': False,
        'half-life': 0.021 * u.s,
    },

    'Os-165': {
        'atomic number': 76,
        'mass number': 165,
        'mass': 164.9766 * u.u,
        'stable': False,
        'half-life': 0.071 * u.s,
    },

    'Os-166': {
        'atomic number': 76,
        'mass number': 166,
        'mass': 165.972692 * u.u,
        'stable': False,
        'half-life': 0.213 * u.s,
    },

    'Os-167': {
        'atomic number': 76,
        'mass number': 167,
        'mass': 166.971549 * u.u,
        'stable': False,
        'half-life': 0.839 * u.s,
    },

    'Os-168': {
        'atomic number': 76,
        'mass number': 168,
        'mass': 167.967808 * u.u,
        'stable': False,
        'half-life': 2.1 * u.s,
    },

    'Os-169': {
        'atomic number': 76,
        'mass number': 169,
        'mass': 168.967018 * u.u,
        'stable': False,
        'half-life': 3.46 * u.s,
    },

    'Os-170': {
        'atomic number': 76,
        'mass number': 170,
        'mass': 169.963578 * u.u,
        'stable': False,
        'half-life': 7.37 * u.s,
    },

    'Os-171': {
        'atomic number': 76,
        'mass number': 171,
        'mass': 170.963174 * u.u,
        'stable': False,
        'half-life': 8.3 * u.s,
    },

    'Os-172': {
        'atomic number': 76,
        'mass number': 172,
        'mass': 171.960017 * u.u,
        'stable': False,
        'half-life': 19.2 * u.s,
    },

    'Os-173': {
        'atomic number': 76,
        'mass number': 173,
        'mass': 172.959808 * u.u,
        'stable': False,
        'half-life': 22.4 * u.s,
    },

    'Os-174': {
        'atomic number': 76,
        'mass number': 174,
        'mass': 173.957064 * u.u,
        'stable': False,
        'half-life': 44.0 * u.s,
    },

    'Os-175': {
        'atomic number': 76,
        'mass number': 175,
        'mass': 174.956945 * u.u,
        'stable': False,
        'half-life': 84.0 * u.s,
    },

    'Os-176': {
        'atomic number': 76,
        'mass number': 176,
        'mass': 175.954806 * u.u,
        'stable': False,
        'half-life': 216.0 * u.s,
    },

    'Os-177': {
        'atomic number': 76,
        'mass number': 177,
        'mass': 176.954966 * u.u,
        'stable': False,
        'half-life': 180.0 * u.s,
    },

    'Os-178': {
        'atomic number': 76,
        'mass number': 178,
        'mass': 177.953254 * u.u,
        'stable': False,
        'half-life': 300.0 * u.s,
    },

    'Os-179': {
        'atomic number': 76,
        'mass number': 179,
        'mass': 178.953817 * u.u,
        'stable': False,
        'half-life': 390.0 * u.s,
    },

    'Os-180': {
        'atomic number': 76,
        'mass number': 180,
        'mass': 179.952375 * u.u,
        'stable': False,
        'half-life': 1290.0 * u.s,
    },

    'Os-181': {
        'atomic number': 76,
        'mass number': 181,
        'mass': 180.953247 * u.u,
        'stable': False,
        'half-life': 6300.0 * u.s,
    },

    'Os-182': {
        'atomic number': 76,
        'mass number': 182,
        'mass': 181.95211 * u.u,
        'stable': False,
        'half-life': 78624.0 * u.s,
    },

    'Os-183': {
        'atomic number': 76,
        'mass number': 183,
        'mass': 182.953125 * u.u,
        'stable': False,
        'half-life': 46800.0 * u.s,
    },

    'Os-184': {
        'atomic number': 76,
        'mass number': 184,
        'mass': 183.9524885 * u.u,
        'stable': True,
        'abundance': 0.0002,
    },

    'Os-185': {
        'atomic number': 76,
        'mass number': 185,
        'mass': 184.9540417 * u.u,
        'stable': False,
        'half-life': 8030880.0 * u.s,
    },

    'Os-186': {
        'atomic number': 76,
        'mass number': 186,
        'mass': 185.953835 * u.u,
        'stable': False,
        'abundance': 0.0159,
    },

    'Os-187': {
        'atomic number': 76,
        'mass number': 187,
        'mass': 186.9557474 * u.u,
        'stable': True,
        'abundance': 0.0196,
    },

    'Os-188': {
        'atomic number': 76,
        'mass number': 188,
        'mass': 187.9558352 * u.u,
        'stable': True,
        'abundance': 0.1324,
    },

    'Os-189': {
        'atomic number': 76,
        'mass number': 189,
        'mass': 188.9581442 * u.u,
        'stable': True,
        'abundance': 0.1615,
    },

    'Os-190': {
        'atomic number': 76,
        'mass number': 190,
        'mass': 189.9584437 * u.u,
        'stable': True,
        'abundance': 0.2626,
    },

    'Os-191': {
        'atomic number': 76,
        'mass number': 191,
        'mass': 190.9609264 * u.u,
        'stable': False,
        'half-life': 1295136.0 * u.s,
    },

    'Os-192': {
        'atomic number': 76,
        'mass number': 192,
        'mass': 191.961477 * u.u,
        'stable': True,
        'abundance': 0.4078,
    },

    'Os-193': {
        'atomic number': 76,
        'mass number': 193,
        'mass': 192.9641479 * u.u,
        'stable': False,
        'half-life': 107388.0 * u.s,
    },

    'Os-194': {
        'atomic number': 76,
        'mass number': 194,
        'mass': 193.9651772 * u.u,
        'stable': False,
        'half-life': 189341556.0 * u.s,
    },

    'Os-195': {
        'atomic number': 76,
        'mass number': 195,
        'mass': 194.968318 * u.u,
        'stable': False,
        'half-life': 390.0 * u.s,
    },

    'Os-196': {
        'atomic number': 76,
        'mass number': 196,
        'mass': 195.969641 * u.u,
        'stable': False,
        'half-life': 2094.0 * u.s,
    },

    'Os-197': {
        'atomic number': 76,
        'mass number': 197,
        'mass': 196.97283 * u.u,
        'stable': False,
        'half-life': 168.0 * u.s,
    },

    'Os-198': {
        'atomic number': 76,
        'mass number': 198,
        'mass': 197.97441 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Os-199': {
        'atomic number': 76,
        'mass number': 199,
        'mass': 198.97801 * u.u,
        'stable': False,
        'half-life': 6.0 * u.s,
    },

    'Os-200': {
        'atomic number': 76,
        'mass number': 200,
        'mass': 199.97984 * u.u,
        'stable': False,
        'half-life': 7.0 * u.s,
    },

    'Os-201': {
        'atomic number': 76,
        'mass number': 201,
        'mass': 200.98364 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Os-202': {
        'atomic number': 76,
        'mass number': 202,
        'mass': 201.98595 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'Ir-164': {
        'atomic number': 77,
        'mass number': 164,
        'mass': 163.99191 * u.u,
        'stable': False,
        'half-life': '1# ms',
    },

    'Ir-165': {
        'atomic number': 77,
        'mass number': 165,
        'mass': 164.9875 * u.u,
        'stable': False,
        'half-life': '50# ns',
    },

    'Ir-166': {
        'atomic number': 77,
        'mass number': 166,
        'mass': 165.98566 * u.u,
        'stable': False,
        'half-life': 0.0105 * u.s,
    },

    'Ir-167': {
        'atomic number': 77,
        'mass number': 167,
        'mass': 166.981666 * u.u,
        'stable': False,
        'half-life': 0.0293 * u.s,
    },

    'Ir-168': {
        'atomic number': 77,
        'mass number': 168,
        'mass': 167.979907 * u.u,
        'stable': False,
        'half-life': 0.23 * u.s,
    },

    'Ir-169': {
        'atomic number': 77,
        'mass number': 169,
        'mass': 168.976298 * u.u,
        'stable': False,
        'half-life': 0.353 * u.s,
    },

    'Ir-170': {
        'atomic number': 77,
        'mass number': 170,
        'mass': 169.974922 * u.u,
        'stable': False,
        'half-life': 0.91 * u.s,
    },

    'Ir-171': {
        'atomic number': 77,
        'mass number': 171,
        'mass': 170.97164 * u.u,
        'stable': False,
        'half-life': 3.1 * u.s,
    },

    'Ir-172': {
        'atomic number': 77,
        'mass number': 172,
        'mass': 171.970607 * u.u,
        'stable': False,
        'half-life': 4.4 * u.s,
    },

    'Ir-173': {
        'atomic number': 77,
        'mass number': 173,
        'mass': 172.967506 * u.u,
        'stable': False,
        'half-life': 9.0 * u.s,
    },

    'Ir-174': {
        'atomic number': 77,
        'mass number': 174,
        'mass': 173.966861 * u.u,
        'stable': False,
        'half-life': 7.9 * u.s,
    },

    'Ir-175': {
        'atomic number': 77,
        'mass number': 175,
        'mass': 174.96415 * u.u,
        'stable': False,
        'half-life': 9.0 * u.s,
    },

    'Ir-176': {
        'atomic number': 77,
        'mass number': 176,
        'mass': 175.96365 * u.u,
        'stable': False,
        'half-life': 8.7 * u.s,
    },

    'Ir-177': {
        'atomic number': 77,
        'mass number': 177,
        'mass': 176.961301 * u.u,
        'stable': False,
        'half-life': 30.0 * u.s,
    },

    'Ir-178': {
        'atomic number': 77,
        'mass number': 178,
        'mass': 177.961082 * u.u,
        'stable': False,
        'half-life': 12.0 * u.s,
    },

    'Ir-179': {
        'atomic number': 77,
        'mass number': 179,
        'mass': 178.95912 * u.u,
        'stable': False,
        'half-life': 79.0 * u.s,
    },

    'Ir-180': {
        'atomic number': 77,
        'mass number': 180,
        'mass': 179.959229 * u.u,
        'stable': False,
        'half-life': 90.0 * u.s,
    },

    'Ir-181': {
        'atomic number': 77,
        'mass number': 181,
        'mass': 180.957625 * u.u,
        'stable': False,
        'half-life': 294.0 * u.s,
    },

    'Ir-182': {
        'atomic number': 77,
        'mass number': 182,
        'mass': 181.958076 * u.u,
        'stable': False,
        'half-life': 900.0 * u.s,
    },

    'Ir-183': {
        'atomic number': 77,
        'mass number': 183,
        'mass': 182.95684 * u.u,
        'stable': False,
        'half-life': 3480.0 * u.s,
    },

    'Ir-184': {
        'atomic number': 77,
        'mass number': 184,
        'mass': 183.957476 * u.u,
        'stable': False,
        'half-life': 11124.0 * u.s,
    },

    'Ir-185': {
        'atomic number': 77,
        'mass number': 185,
        'mass': 184.956698 * u.u,
        'stable': False,
        'half-life': 51840.0 * u.s,
    },

    'Ir-186': {
        'atomic number': 77,
        'mass number': 186,
        'mass': 185.957944 * u.u,
        'stable': False,
        'half-life': 59904.0 * u.s,
    },

    'Ir-187': {
        'atomic number': 77,
        'mass number': 187,
        'mass': 186.957542 * u.u,
        'stable': False,
        'half-life': 37800.0 * u.s,
    },

    'Ir-188': {
        'atomic number': 77,
        'mass number': 188,
        'mass': 187.958828 * u.u,
        'stable': False,
        'half-life': 149400.0 * u.s,
    },

    'Ir-189': {
        'atomic number': 77,
        'mass number': 189,
        'mass': 188.958715 * u.u,
        'stable': False,
        'half-life': 1140480.0 * u.s,
    },

    'Ir-190': {
        'atomic number': 77,
        'mass number': 190,
        'mass': 189.9605412 * u.u,
        'stable': False,
        'half-life': 1017792.0 * u.s,
    },

    'Ir-191': {
        'atomic number': 77,
        'mass number': 191,
        'mass': 190.9605893 * u.u,
        'stable': True,
        'abundance': 0.373,
    },

    'Ir-192': {
        'atomic number': 77,
        'mass number': 192,
        'mass': 191.9626002 * u.u,
        'stable': False,
        'half-life': 6377184.0 * u.s,
    },

    'Ir-193': {
        'atomic number': 77,
        'mass number': 193,
        'mass': 192.9629216 * u.u,
        'stable': True,
        'abundance': 0.627,
    },

    'Ir-194': {
        'atomic number': 77,
        'mass number': 194,
        'mass': 193.9650735 * u.u,
        'stable': False,
        'half-life': 69408.0 * u.s,
    },

    'Ir-195': {
        'atomic number': 77,
        'mass number': 195,
        'mass': 194.9659747 * u.u,
        'stable': False,
        'half-life': 8244.0 * u.s,
    },

    'Ir-196': {
        'atomic number': 77,
        'mass number': 196,
        'mass': 195.968397 * u.u,
        'stable': False,
        'half-life': 52.0 * u.s,
    },

    'Ir-197': {
        'atomic number': 77,
        'mass number': 197,
        'mass': 196.969655 * u.u,
        'stable': False,
        'half-life': 348.0 * u.s,
    },

    'Ir-198': {
        'atomic number': 77,
        'mass number': 198,
        'mass': 197.97228 * u.u,
        'stable': False,
        'half-life': 8.0 * u.s,
    },

    'Ir-199': {
        'atomic number': 77,
        'mass number': 199,
        'mass': 198.973805 * u.u,
        'stable': False,
        'half-life': 7.0 * u.s,
    },

    'Ir-200': {
        'atomic number': 77,
        'mass number': 200,
        'mass': 199.9768 * u.u,
        'stable': False,
        'half-life': 43.0 * u.s,
    },

    'Ir-201': {
        'atomic number': 77,
        'mass number': 201,
        'mass': 200.97864 * u.u,
        'stable': False,
        'half-life': 21.0 * u.s,
    },

    'Ir-202': {
        'atomic number': 77,
        'mass number': 202,
        'mass': 201.98199 * u.u,
        'stable': False,
        'half-life': 11.0 * u.s,
    },

    'Ir-203': {
        'atomic number': 77,
        'mass number': 203,
        'mass': 202.98423 * u.u,
        'stable': False,
        'half-life': '6# s',
    },

    'Ir-204': {
        'atomic number': 77,
        'mass number': 204,
        'mass': 203.9896 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Pt-166': {
        'atomic number': 78,
        'mass number': 166,
        'mass': 165.99486 * u.u,
        'stable': False,
        'half-life': 0.0003 * u.s,
    },

    'Pt-167': {
        'atomic number': 78,
        'mass number': 167,
        'mass': 166.99269 * u.u,
        'stable': False,
        'half-life': 0.0008 * u.s,
    },

    'Pt-168': {
        'atomic number': 78,
        'mass number': 168,
        'mass': 167.98813 * u.u,
        'stable': False,
        'half-life': 0.00202 * u.s,
    },

    'Pt-169': {
        'atomic number': 78,
        'mass number': 169,
        'mass': 168.98657 * u.u,
        'stable': False,
        'half-life': 0.00699 * u.s,
    },

    'Pt-170': {
        'atomic number': 78,
        'mass number': 170,
        'mass': 169.982496 * u.u,
        'stable': False,
        'half-life': 0.01393 * u.s,
    },

    'Pt-171': {
        'atomic number': 78,
        'mass number': 171,
        'mass': 170.981245 * u.u,
        'stable': False,
        'half-life': 0.0455 * u.s,
    },

    'Pt-172': {
        'atomic number': 78,
        'mass number': 172,
        'mass': 171.977351 * u.u,
        'stable': False,
        'half-life': 0.0976 * u.s,
    },

    'Pt-173': {
        'atomic number': 78,
        'mass number': 173,
        'mass': 172.976443 * u.u,
        'stable': False,
        'half-life': 0.382 * u.s,
    },

    'Pt-174': {
        'atomic number': 78,
        'mass number': 174,
        'mass': 173.97282 * u.u,
        'stable': False,
        'half-life': 0.889 * u.s,
    },

    'Pt-175': {
        'atomic number': 78,
        'mass number': 175,
        'mass': 174.97241 * u.u,
        'stable': False,
        'half-life': 2.43 * u.s,
    },

    'Pt-176': {
        'atomic number': 78,
        'mass number': 176,
        'mass': 175.968938 * u.u,
        'stable': False,
        'half-life': 6.33 * u.s,
    },

    'Pt-177': {
        'atomic number': 78,
        'mass number': 177,
        'mass': 176.96847 * u.u,
        'stable': False,
        'half-life': 10.6 * u.s,
    },

    'Pt-178': {
        'atomic number': 78,
        'mass number': 178,
        'mass': 177.96565 * u.u,
        'stable': False,
        'half-life': 20.7 * u.s,
    },

    'Pt-179': {
        'atomic number': 78,
        'mass number': 179,
        'mass': 178.965359 * u.u,
        'stable': False,
        'half-life': 21.2 * u.s,
    },

    'Pt-180': {
        'atomic number': 78,
        'mass number': 180,
        'mass': 179.963032 * u.u,
        'stable': False,
        'half-life': 56.0 * u.s,
    },

    'Pt-181': {
        'atomic number': 78,
        'mass number': 181,
        'mass': 180.963098 * u.u,
        'stable': False,
        'half-life': 52.0 * u.s,
    },

    'Pt-182': {
        'atomic number': 78,
        'mass number': 182,
        'mass': 181.961172 * u.u,
        'stable': False,
        'half-life': 160.2 * u.s,
    },

    'Pt-183': {
        'atomic number': 78,
        'mass number': 183,
        'mass': 182.961597 * u.u,
        'stable': False,
        'half-life': 390.0 * u.s,
    },

    'Pt-184': {
        'atomic number': 78,
        'mass number': 184,
        'mass': 183.959915 * u.u,
        'stable': False,
        'half-life': 1038.0 * u.s,
    },

    'Pt-185': {
        'atomic number': 78,
        'mass number': 185,
        'mass': 184.960614 * u.u,
        'stable': False,
        'half-life': 4254.0 * u.s,
    },

    'Pt-186': {
        'atomic number': 78,
        'mass number': 186,
        'mass': 185.959351 * u.u,
        'stable': False,
        'half-life': 7488.0 * u.s,
    },

    'Pt-187': {
        'atomic number': 78,
        'mass number': 187,
        'mass': 186.960617 * u.u,
        'stable': False,
        'half-life': 8460.0 * u.s,
    },

    'Pt-188': {
        'atomic number': 78,
        'mass number': 188,
        'mass': 187.9593889 * u.u,
        'stable': False,
        'half-life': 881280.0 * u.s,
    },

    'Pt-189': {
        'atomic number': 78,
        'mass number': 189,
        'mass': 188.960831 * u.u,
        'stable': False,
        'half-life': 39132.0 * u.s,
    },

    'Pt-190': {
        'atomic number': 78,
        'mass number': 190,
        'mass': 189.9599297 * u.u,
        'stable': False,
        'abundance': 0.00012,
    },

    'Pt-191': {
        'atomic number': 78,
        'mass number': 191,
        'mass': 190.9616729 * u.u,
        'stable': False,
        'half-life': 244512.0 * u.s,
    },

    'Pt-192': {
        'atomic number': 78,
        'mass number': 192,
        'mass': 191.9610387 * u.u,
        'stable': True,
        'abundance': 0.00782,
    },

    'Pt-193': {
        'atomic number': 78,
        'mass number': 193,
        'mass': 192.9629824 * u.u,
        'stable': False,
        'half-life': 1577846300.0 * u.s,
    },

    'Pt-194': {
        'atomic number': 78,
        'mass number': 194,
        'mass': 193.9626809 * u.u,
        'stable': True,
        'abundance': 0.3286,
    },

    'Pt-195': {
        'atomic number': 78,
        'mass number': 195,
        'mass': 194.9647917 * u.u,
        'stable': True,
        'abundance': 0.3378,
    },

    'Pt-196': {
        'atomic number': 78,
        'mass number': 196,
        'mass': 195.96495209 * u.u,
        'stable': True,
        'abundance': 0.2521,
    },

    'Pt-197': {
        'atomic number': 78,
        'mass number': 197,
        'mass': 196.96734069 * u.u,
        'stable': False,
        'half-life': 71609.4 * u.s,
    },

    'Pt-198': {
        'atomic number': 78,
        'mass number': 198,
        'mass': 197.9678949 * u.u,
        'stable': True,
        'abundance': 0.07356,
    },

    'Pt-199': {
        'atomic number': 78,
        'mass number': 199,
        'mass': 198.9705952 * u.u,
        'stable': False,
        'half-life': 1848.0 * u.s,
    },

    'Pt-200': {
        'atomic number': 78,
        'mass number': 200,
        'mass': 199.971443 * u.u,
        'stable': False,
        'half-life': 45360.0 * u.s,
    },

    'Pt-201': {
        'atomic number': 78,
        'mass number': 201,
        'mass': 200.974513 * u.u,
        'stable': False,
        'half-life': 150.0 * u.s,
    },

    'Pt-202': {
        'atomic number': 78,
        'mass number': 202,
        'mass': 201.975639 * u.u,
        'stable': False,
        'half-life': 158400.0 * u.s,
    },

    'Pt-203': {
        'atomic number': 78,
        'mass number': 203,
        'mass': 202.97893 * u.u,
        'stable': False,
        'half-life': 22.0 * u.s,
    },

    'Pt-204': {
        'atomic number': 78,
        'mass number': 204,
        'mass': 203.98076 * u.u,
        'stable': False,
        'half-life': 10.3 * u.s,
    },

    'Pt-205': {
        'atomic number': 78,
        'mass number': 205,
        'mass': 204.98608 * u.u,
        'stable': False,
        'half-life': '5# s',
    },

    'Pt-206': {
        'atomic number': 78,
        'mass number': 206,
        'mass': 205.98966 * u.u,
        'stable': False,
        'half-life': '5# s',
    },

    'Au-169': {
        'atomic number': 79,
        'mass number': 169,
        'mass': 168.99808 * u.u,
        'stable': False,
        'half-life': '150# us',
    },

    'Au-170': {
        'atomic number': 79,
        'mass number': 170,
        'mass': 169.99597 * u.u,
        'stable': False,
        'half-life': 0.00029 * u.s,
    },

    'Au-171': {
        'atomic number': 79,
        'mass number': 171,
        'mass': 170.991876 * u.u,
        'stable': False,
        'half-life': 2.23e-05 * u.s,
    },

    'Au-172': {
        'atomic number': 79,
        'mass number': 172,
        'mass': 171.989942 * u.u,
        'stable': False,
        'half-life': 0.028 * u.s,
    },

    'Au-173': {
        'atomic number': 79,
        'mass number': 173,
        'mass': 172.986241 * u.u,
        'stable': False,
        'half-life': 0.0255 * u.s,
    },

    'Au-174': {
        'atomic number': 79,
        'mass number': 174,
        'mass': 173.984717 * u.u,
        'stable': False,
        'half-life': 0.139 * u.s,
    },

    'Au-175': {
        'atomic number': 79,
        'mass number': 175,
        'mass': 174.981304 * u.u,
        'stable': False,
        'half-life': 0.202 * u.s,
    },

    'Au-176': {
        'atomic number': 79,
        'mass number': 176,
        'mass': 175.98025 * u.u,
        'stable': False,
        'half-life': 1.05 * u.s,
    },

    'Au-177': {
        'atomic number': 79,
        'mass number': 177,
        'mass': 176.97687 * u.u,
        'stable': False,
        'half-life': 1.46 * u.s,
    },

    'Au-178': {
        'atomic number': 79,
        'mass number': 178,
        'mass': 177.976032 * u.u,
        'stable': False,
        'half-life': 2.6 * u.s,
    },

    'Au-179': {
        'atomic number': 79,
        'mass number': 179,
        'mass': 178.973174 * u.u,
        'stable': False,
        'half-life': 7.1 * u.s,
    },

    'Au-180': {
        'atomic number': 79,
        'mass number': 180,
        'mass': 179.972523 * u.u,
        'stable': False,
        'half-life': 8.4 * u.s,
    },

    'Au-181': {
        'atomic number': 79,
        'mass number': 181,
        'mass': 180.970079 * u.u,
        'stable': False,
        'half-life': 13.7 * u.s,
    },

    'Au-182': {
        'atomic number': 79,
        'mass number': 182,
        'mass': 181.969618 * u.u,
        'stable': False,
        'half-life': 15.5 * u.s,
    },

    'Au-183': {
        'atomic number': 79,
        'mass number': 183,
        'mass': 182.967591 * u.u,
        'stable': False,
        'half-life': 42.8 * u.s,
    },

    'Au-184': {
        'atomic number': 79,
        'mass number': 184,
        'mass': 183.967452 * u.u,
        'stable': False,
        'half-life': 20.6 * u.s,
    },

    'Au-185': {
        'atomic number': 79,
        'mass number': 185,
        'mass': 184.96579 * u.u,
        'stable': False,
        'half-life': 255.0 * u.s,
    },

    'Au-186': {
        'atomic number': 79,
        'mass number': 186,
        'mass': 185.965953 * u.u,
        'stable': False,
        'half-life': 642.0 * u.s,
    },

    'Au-187': {
        'atomic number': 79,
        'mass number': 187,
        'mass': 186.964543 * u.u,
        'stable': False,
        'half-life': 498.0 * u.s,
    },

    'Au-188': {
        'atomic number': 79,
        'mass number': 188,
        'mass': 187.965349 * u.u,
        'stable': False,
        'half-life': 530.4 * u.s,
    },

    'Au-189': {
        'atomic number': 79,
        'mass number': 189,
        'mass': 188.963948 * u.u,
        'stable': False,
        'half-life': 1722.0 * u.s,
    },

    'Au-190': {
        'atomic number': 79,
        'mass number': 190,
        'mass': 189.964698 * u.u,
        'stable': False,
        'half-life': 2568.0 * u.s,
    },

    'Au-191': {
        'atomic number': 79,
        'mass number': 191,
        'mass': 190.963702 * u.u,
        'stable': False,
        'half-life': 11448.0 * u.s,
    },

    'Au-192': {
        'atomic number': 79,
        'mass number': 192,
        'mass': 191.964814 * u.u,
        'stable': False,
        'half-life': 17784.0 * u.s,
    },

    'Au-193': {
        'atomic number': 79,
        'mass number': 193,
        'mass': 192.9641373 * u.u,
        'stable': False,
        'half-life': 63540.0 * u.s,
    },

    'Au-194': {
        'atomic number': 79,
        'mass number': 194,
        'mass': 193.9654178 * u.u,
        'stable': False,
        'half-life': 136872.0 * u.s,
    },

    'Au-195': {
        'atomic number': 79,
        'mass number': 195,
        'mass': 194.9650352 * u.u,
        'stable': False,
        'half-life': 16078867.200000001 * u.s,
    },

    'Au-196': {
        'atomic number': 79,
        'mass number': 196,
        'mass': 195.9665699 * u.u,
        'stable': False,
        'half-life': 532820.16 * u.s,
    },

    'Au-197': {
        'atomic number': 79,
        'mass number': 197,
        'mass': 196.96656879 * u.u,
        'stable': True,
        'abundance': 1,
    },

    'Au-198': {
        'atomic number': 79,
        'mass number': 198,
        'mass': 197.96824242 * u.u,
        'stable': False,
        'half-life': 232862.688 * u.s,
    },

    'Au-199': {
        'atomic number': 79,
        'mass number': 199,
        'mass': 198.96876528 * u.u,
        'stable': False,
        'half-life': 271209.6 * u.s,
    },

    'Au-200': {
        'atomic number': 79,
        'mass number': 200,
        'mass': 199.970756 * u.u,
        'stable': False,
        'half-life': 2904.0 * u.s,
    },

    'Au-201': {
        'atomic number': 79,
        'mass number': 201,
        'mass': 200.9716575 * u.u,
        'stable': False,
        'half-life': 1560.0 * u.s,
    },

    'Au-202': {
        'atomic number': 79,
        'mass number': 202,
        'mass': 201.973856 * u.u,
        'stable': False,
        'half-life': 28.4 * u.s,
    },

    'Au-203': {
        'atomic number': 79,
        'mass number': 203,
        'mass': 202.9751544 * u.u,
        'stable': False,
        'half-life': 60.0 * u.s,
    },

    'Au-204': {
        'atomic number': 79,
        'mass number': 204,
        'mass': 203.97783 * u.u,
        'stable': False,
        'half-life': 38.3 * u.s,
    },

    'Au-205': {
        'atomic number': 79,
        'mass number': 205,
        'mass': 204.97985 * u.u,
        'stable': False,
        'half-life': 32.5 * u.s,
    },

    'Au-206': {
        'atomic number': 79,
        'mass number': 206,
        'mass': 205.98474 * u.u,
        'stable': False,
        'half-life': 47.0 * u.s,
    },

    'Au-207': {
        'atomic number': 79,
        'mass number': 207,
        'mass': 206.9884 * u.u,
        'stable': False,
        'half-life': '10# s',
    },

    'Au-208': {
        'atomic number': 79,
        'mass number': 208,
        'mass': 207.99345 * u.u,
        'stable': False,
        'half-life': '10# s',
    },

    'Au-209': {
        'atomic number': 79,
        'mass number': 209,
        'mass': 208.99735 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Au-210': {
        'atomic number': 79,
        'mass number': 210,
        'mass': 210.0025 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Hg-171': {
        'atomic number': 80,
        'mass number': 171,
        'mass': 171.00353 * u.u,
        'stable': False,
        'half-life': 7e-05 * u.s,
    },

    'Hg-172': {
        'atomic number': 80,
        'mass number': 172,
        'mass': 171.99881 * u.u,
        'stable': False,
        'half-life': 0.000231 * u.s,
    },

    'Hg-173': {
        'atomic number': 80,
        'mass number': 173,
        'mass': 172.99709 * u.u,
        'stable': False,
        'half-life': 0.0008 * u.s,
    },

    'Hg-174': {
        'atomic number': 80,
        'mass number': 174,
        'mass': 173.992865 * u.u,
        'stable': False,
        'half-life': 0.002 * u.s,
    },

    'Hg-175': {
        'atomic number': 80,
        'mass number': 175,
        'mass': 174.991441 * u.u,
        'stable': False,
        'half-life': 0.0106 * u.s,
    },

    'Hg-176': {
        'atomic number': 80,
        'mass number': 176,
        'mass': 175.987361 * u.u,
        'stable': False,
        'half-life': 0.0203 * u.s,
    },

    'Hg-177': {
        'atomic number': 80,
        'mass number': 177,
        'mass': 176.986277 * u.u,
        'stable': False,
        'half-life': 0.1273 * u.s,
    },

    'Hg-178': {
        'atomic number': 80,
        'mass number': 178,
        'mass': 177.982484 * u.u,
        'stable': False,
        'half-life': 0.2665 * u.s,
    },

    'Hg-179': {
        'atomic number': 80,
        'mass number': 179,
        'mass': 178.981831 * u.u,
        'stable': False,
        'half-life': 1.05 * u.s,
    },

    'Hg-180': {
        'atomic number': 80,
        'mass number': 180,
        'mass': 179.97826 * u.u,
        'stable': False,
        'half-life': 2.59 * u.s,
    },

    'Hg-181': {
        'atomic number': 80,
        'mass number': 181,
        'mass': 180.977819 * u.u,
        'stable': False,
        'half-life': 3.6 * u.s,
    },

    'Hg-182': {
        'atomic number': 80,
        'mass number': 182,
        'mass': 181.974689 * u.u,
        'stable': False,
        'half-life': 10.83 * u.s,
    },

    'Hg-183': {
        'atomic number': 80,
        'mass number': 183,
        'mass': 182.9744448 * u.u,
        'stable': False,
        'half-life': 9.4 * u.s,
    },

    'Hg-184': {
        'atomic number': 80,
        'mass number': 184,
        'mass': 183.971714 * u.u,
        'stable': False,
        'half-life': 30.87 * u.s,
    },

    'Hg-185': {
        'atomic number': 80,
        'mass number': 185,
        'mass': 184.971899 * u.u,
        'stable': False,
        'half-life': 49.1 * u.s,
    },

    'Hg-186': {
        'atomic number': 80,
        'mass number': 186,
        'mass': 185.969362 * u.u,
        'stable': False,
        'half-life': 82.8 * u.s,
    },

    'Hg-187': {
        'atomic number': 80,
        'mass number': 187,
        'mass': 186.969814 * u.u,
        'stable': False,
        'half-life': 114.0 * u.s,
    },

    'Hg-188': {
        'atomic number': 80,
        'mass number': 188,
        'mass': 187.967567 * u.u,
        'stable': False,
        'half-life': 195.0 * u.s,
    },

    'Hg-189': {
        'atomic number': 80,
        'mass number': 189,
        'mass': 188.968195 * u.u,
        'stable': False,
        'half-life': 456.0 * u.s,
    },

    'Hg-190': {
        'atomic number': 80,
        'mass number': 190,
        'mass': 189.966323 * u.u,
        'stable': False,
        'half-life': 1200.0 * u.s,
    },

    'Hg-191': {
        'atomic number': 80,
        'mass number': 191,
        'mass': 190.967157 * u.u,
        'stable': False,
        'half-life': 2940.0 * u.s,
    },

    'Hg-192': {
        'atomic number': 80,
        'mass number': 192,
        'mass': 191.965635 * u.u,
        'stable': False,
        'half-life': 17460.0 * u.s,
    },

    'Hg-193': {
        'atomic number': 80,
        'mass number': 193,
        'mass': 192.966653 * u.u,
        'stable': False,
        'half-life': 13680.0 * u.s,
    },

    'Hg-194': {
        'atomic number': 80,
        'mass number': 194,
        'mass': 193.9654491 * u.u,
        'stable': False,
        'half-life': 14105945922.0 * u.s,
    },

    'Hg-195': {
        'atomic number': 80,
        'mass number': 195,
        'mass': 194.966721 * u.u,
        'stable': False,
        'half-life': 38484.0 * u.s,
    },

    'Hg-196': {
        'atomic number': 80,
        'mass number': 196,
        'mass': 195.9658326 * u.u,
        'stable': True,
        'abundance': 0.0015,
    },

    'Hg-197': {
        'atomic number': 80,
        'mass number': 197,
        'mass': 196.9672128 * u.u,
        'stable': False,
        'half-life': 233784.0 * u.s,
    },

    'Hg-198': {
        'atomic number': 80,
        'mass number': 198,
        'mass': 197.9667686 * u.u,
        'stable': True,
        'abundance': 0.0997,
    },

    'Hg-199': {
        'atomic number': 80,
        'mass number': 199,
        'mass': 198.96828064 * u.u,
        'stable': True,
        'abundance': 0.1687,
    },

    'Hg-200': {
        'atomic number': 80,
        'mass number': 200,
        'mass': 199.96832659 * u.u,
        'stable': True,
        'abundance': 0.231,
    },

    'Hg-201': {
        'atomic number': 80,
        'mass number': 201,
        'mass': 200.97030284 * u.u,
        'stable': True,
        'abundance': 0.1318,
    },

    'Hg-202': {
        'atomic number': 80,
        'mass number': 202,
        'mass': 201.9706434 * u.u,
        'stable': True,
        'abundance': 0.2986,
    },

    'Hg-203': {
        'atomic number': 80,
        'mass number': 203,
        'mass': 202.9728728 * u.u,
        'stable': False,
        'half-life': 4027881.6 * u.s,
    },

    'Hg-204': {
        'atomic number': 80,
        'mass number': 204,
        'mass': 203.97349398 * u.u,
        'stable': True,
        'abundance': 0.0687,
    },

    'Hg-205': {
        'atomic number': 80,
        'mass number': 205,
        'mass': 204.9760734 * u.u,
        'stable': False,
        'half-life': 308.4 * u.s,
    },

    'Hg-206': {
        'atomic number': 80,
        'mass number': 206,
        'mass': 205.977514 * u.u,
        'stable': False,
        'half-life': 499.2 * u.s,
    },

    'Hg-207': {
        'atomic number': 80,
        'mass number': 207,
        'mass': 206.9823 * u.u,
        'stable': False,
        'half-life': 174.0 * u.s,
    },

    'Hg-208': {
        'atomic number': 80,
        'mass number': 208,
        'mass': 207.985759 * u.u,
        'stable': False,
        'half-life': 2520.0 * u.s,
    },

    'Hg-209': {
        'atomic number': 80,
        'mass number': 209,
        'mass': 208.99072 * u.u,
        'stable': False,
        'half-life': 38.0 * u.s,
    },

    'Hg-210': {
        'atomic number': 80,
        'mass number': 210,
        'mass': 209.99424 * u.u,
        'stable': False,
        'half-life': 64.0 * u.s,
    },

    'Hg-211': {
        'atomic number': 80,
        'mass number': 211,
        'mass': 210.99933 * u.u,
        'stable': False,
        'half-life': 26.0 * u.s,
    },

    'Hg-212': {
        'atomic number': 80,
        'mass number': 212,
        'mass': 212.00296 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Hg-213': {
        'atomic number': 80,
        'mass number': 213,
        'mass': 213.00823 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Hg-214': {
        'atomic number': 80,
        'mass number': 214,
        'mass': 214.012 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Hg-215': {
        'atomic number': 80,
        'mass number': 215,
        'mass': 215.0174 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Hg-216': {
        'atomic number': 80,
        'mass number': 216,
        'mass': 216.02132 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Tl-176': {
        'atomic number': 81,
        'mass number': 176,
        'mass': 176.000624 * u.u,
        'stable': False,
        'half-life': 0.0062 * u.s,
    },

    'Tl-177': {
        'atomic number': 81,
        'mass number': 177,
        'mass': 176.996431 * u.u,
        'stable': False,
        'half-life': 0.018 * u.s,
    },

    'Tl-178': {
        'atomic number': 81,
        'mass number': 178,
        'mass': 177.99485 * u.u,
        'stable': False,
        'half-life': 0.255 * u.s,
    },

    'Tl-179': {
        'atomic number': 81,
        'mass number': 179,
        'mass': 178.991111 * u.u,
        'stable': False,
        'half-life': 0.265 * u.s,
    },

    'Tl-180': {
        'atomic number': 81,
        'mass number': 180,
        'mass': 179.990057 * u.u,
        'stable': False,
        'half-life': 1.09 * u.s,
    },

    'Tl-181': {
        'atomic number': 81,
        'mass number': 181,
        'mass': 180.98626 * u.u,
        'stable': False,
        'half-life': 3.2 * u.s,
    },

    'Tl-182': {
        'atomic number': 81,
        'mass number': 182,
        'mass': 181.985713 * u.u,
        'stable': False,
        'half-life': 1.9 * u.s,
    },

    'Tl-183': {
        'atomic number': 81,
        'mass number': 183,
        'mass': 182.982193 * u.u,
        'stable': False,
        'half-life': 6.9 * u.s,
    },

    'Tl-184': {
        'atomic number': 81,
        'mass number': 184,
        'mass': 183.981886 * u.u,
        'stable': False,
        'half-life': 9.5 * u.s,
    },

    'Tl-185': {
        'atomic number': 81,
        'mass number': 185,
        'mass': 184.978789 * u.u,
        'stable': False,
        'half-life': 19.5 * u.s,
    },

    'Tl-186': {
        'atomic number': 81,
        'mass number': 186,
        'mass': 185.978651 * u.u,
        'stable': False,
        'half-life': '40# s',
    },

    'Tl-187': {
        'atomic number': 81,
        'mass number': 187,
        'mass': 186.9759063 * u.u,
        'stable': False,
        'half-life': '~51 s',
    },

    'Tl-188': {
        'atomic number': 81,
        'mass number': 188,
        'mass': 187.976021 * u.u,
        'stable': False,
        'half-life': 71.0 * u.s,
    },

    'Tl-189': {
        'atomic number': 81,
        'mass number': 189,
        'mass': 188.973588 * u.u,
        'stable': False,
        'half-life': 138.0 * u.s,
    },

    'Tl-190': {
        'atomic number': 81,
        'mass number': 190,
        'mass': 189.973828 * u.u,
        'stable': False,
        'half-life': 156.0 * u.s,
    },

    'Tl-191': {
        'atomic number': 81,
        'mass number': 191,
        'mass': 190.9717842 * u.u,
        'stable': False,
        'half-life': '20# m',
    },

    'Tl-192': {
        'atomic number': 81,
        'mass number': 192,
        'mass': 191.972225 * u.u,
        'stable': False,
        'half-life': 576.0 * u.s,
    },

    'Tl-193': {
        'atomic number': 81,
        'mass number': 193,
        'mass': 192.970502 * u.u,
        'stable': False,
        'half-life': 1296.0 * u.s,
    },

    'Tl-194': {
        'atomic number': 81,
        'mass number': 194,
        'mass': 193.971081 * u.u,
        'stable': False,
        'half-life': 1980.0 * u.s,
    },

    'Tl-195': {
        'atomic number': 81,
        'mass number': 195,
        'mass': 194.969774 * u.u,
        'stable': False,
        'half-life': 4176.0 * u.s,
    },

    'Tl-196': {
        'atomic number': 81,
        'mass number': 196,
        'mass': 195.970481 * u.u,
        'stable': False,
        'half-life': 6624.0 * u.s,
    },

    'Tl-197': {
        'atomic number': 81,
        'mass number': 197,
        'mass': 196.969576 * u.u,
        'stable': False,
        'half-life': 10224.0 * u.s,
    },

    'Tl-198': {
        'atomic number': 81,
        'mass number': 198,
        'mass': 197.970483 * u.u,
        'stable': False,
        'half-life': 19080.0 * u.s,
    },

    'Tl-199': {
        'atomic number': 81,
        'mass number': 199,
        'mass': 198.969877 * u.u,
        'stable': False,
        'half-life': 26712.0 * u.s,
    },

    'Tl-200': {
        'atomic number': 81,
        'mass number': 200,
        'mass': 199.9709633 * u.u,
        'stable': False,
        'half-life': 93960.0 * u.s,
    },

    'Tl-201': {
        'atomic number': 81,
        'mass number': 201,
        'mass': 200.970822 * u.u,
        'stable': False,
        'half-life': 263139.83999999997 * u.s,
    },

    'Tl-202': {
        'atomic number': 81,
        'mass number': 202,
        'mass': 201.972102 * u.u,
        'stable': False,
        'half-life': 1077062.4 * u.s,
    },

    'Tl-203': {
        'atomic number': 81,
        'mass number': 203,
        'mass': 202.9723446 * u.u,
        'stable': True,
        'abundance': 0.2952,
    },

    'Tl-204': {
        'atomic number': 81,
        'mass number': 204,
        'mass': 203.9738639 * u.u,
        'stable': False,
        'half-life': 119379851.058 * u.s,
    },

    'Tl-205': {
        'atomic number': 81,
        'mass number': 205,
        'mass': 204.9744278 * u.u,
        'stable': True,
        'abundance': 0.7048,
    },

    'Tl-206': {
        'atomic number': 81,
        'mass number': 206,
        'mass': 205.9761106 * u.u,
        'stable': False,
        'half-life': 252.12 * u.s,
    },

    'Tl-207': {
        'atomic number': 81,
        'mass number': 207,
        'mass': 206.9774197 * u.u,
        'stable': False,
        'half-life': 286.2 * u.s,
    },

    'Tl-208': {
        'atomic number': 81,
        'mass number': 208,
        'mass': 207.982019 * u.u,
        'stable': False,
        'half-life': 183.18 * u.s,
    },

    'Tl-209': {
        'atomic number': 81,
        'mass number': 209,
        'mass': 208.9853594 * u.u,
        'stable': False,
        'half-life': 129.72 * u.s,
    },

    'Tl-210': {
        'atomic number': 81,
        'mass number': 210,
        'mass': 209.990074 * u.u,
        'stable': False,
        'half-life': 78.0 * u.s,
    },

    'Tl-211': {
        'atomic number': 81,
        'mass number': 211,
        'mass': 210.993475 * u.u,
        'stable': False,
        'half-life': 80.0 * u.s,
    },

    'Tl-212': {
        'atomic number': 81,
        'mass number': 212,
        'mass': 211.99834 * u.u,
        'stable': False,
        'half-life': 31.0 * u.s,
    },

    'Tl-213': {
        'atomic number': 81,
        'mass number': 213,
        'mass': 213.001915 * u.u,
        'stable': False,
        'half-life': 24.0 * u.s,
    },

    'Tl-214': {
        'atomic number': 81,
        'mass number': 214,
        'mass': 214.00694 * u.u,
        'stable': False,
        'half-life': 11.0 * u.s,
    },

    'Tl-215': {
        'atomic number': 81,
        'mass number': 215,
        'mass': 215.01064 * u.u,
        'stable': False,
        'half-life': 10.0 * u.s,
    },

    'Tl-216': {
        'atomic number': 81,
        'mass number': 216,
        'mass': 216.0158 * u.u,
        'stable': False,
        'half-life': 6.0 * u.s,
    },

    'Tl-217': {
        'atomic number': 81,
        'mass number': 217,
        'mass': 217.01966 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Tl-218': {
        'atomic number': 81,
        'mass number': 218,
        'mass': 218.02479 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'Pb-178': {
        'atomic number': 82,
        'mass number': 178,
        'mass': 178.003831 * u.u,
        'stable': False,
        'half-life': 0.00023 * u.s,
    },

    'Pb-179': {
        'atomic number': 82,
        'mass number': 179,
        'mass': 179.002201 * u.u,
        'stable': False,
        'half-life': 0.0039 * u.s,
    },

    'Pb-180': {
        'atomic number': 82,
        'mass number': 180,
        'mass': 179.997928 * u.u,
        'stable': False,
        'half-life': 0.0041 * u.s,
    },

    'Pb-181': {
        'atomic number': 82,
        'mass number': 181,
        'mass': 180.996653 * u.u,
        'stable': False,
        'half-life': 0.039 * u.s,
    },

    'Pb-182': {
        'atomic number': 82,
        'mass number': 182,
        'mass': 181.992672 * u.u,
        'stable': False,
        'half-life': 0.055 * u.s,
    },

    'Pb-183': {
        'atomic number': 82,
        'mass number': 183,
        'mass': 182.991872 * u.u,
        'stable': False,
        'half-life': 0.535 * u.s,
    },

    'Pb-184': {
        'atomic number': 82,
        'mass number': 184,
        'mass': 183.988136 * u.u,
        'stable': False,
        'half-life': 0.49 * u.s,
    },

    'Pb-185': {
        'atomic number': 82,
        'mass number': 185,
        'mass': 184.98761 * u.u,
        'stable': False,
        'half-life': 6.3 * u.s,
    },

    'Pb-186': {
        'atomic number': 82,
        'mass number': 186,
        'mass': 185.984238 * u.u,
        'stable': False,
        'half-life': 4.82 * u.s,
    },

    'Pb-187': {
        'atomic number': 82,
        'mass number': 187,
        'mass': 186.9839109 * u.u,
        'stable': False,
        'half-life': 15.2 * u.s,
    },

    'Pb-188': {
        'atomic number': 82,
        'mass number': 188,
        'mass': 187.980875 * u.u,
        'stable': False,
        'half-life': 25.1 * u.s,
    },

    'Pb-189': {
        'atomic number': 82,
        'mass number': 189,
        'mass': 188.980807 * u.u,
        'stable': False,
        'half-life': 39.0 * u.s,
    },

    'Pb-190': {
        'atomic number': 82,
        'mass number': 190,
        'mass': 189.978082 * u.u,
        'stable': False,
        'half-life': 71.0 * u.s,
    },

    'Pb-191': {
        'atomic number': 82,
        'mass number': 191,
        'mass': 190.978276 * u.u,
        'stable': False,
        'half-life': 79.8 * u.s,
    },

    'Pb-192': {
        'atomic number': 82,
        'mass number': 192,
        'mass': 191.975775 * u.u,
        'stable': False,
        'half-life': 210.0 * u.s,
    },

    'Pb-193': {
        'atomic number': 82,
        'mass number': 193,
        'mass': 192.976173 * u.u,
        'stable': False,
        'half-life': '5# m',
    },

    'Pb-194': {
        'atomic number': 82,
        'mass number': 194,
        'mass': 193.974012 * u.u,
        'stable': False,
        'half-life': 642.0 * u.s,
    },

    'Pb-195': {
        'atomic number': 82,
        'mass number': 195,
        'mass': 194.974543 * u.u,
        'stable': False,
        'half-life': '~15 m',
    },

    'Pb-196': {
        'atomic number': 82,
        'mass number': 196,
        'mass': 195.972774 * u.u,
        'stable': False,
        'half-life': 2220.0 * u.s,
    },

    'Pb-197': {
        'atomic number': 82,
        'mass number': 197,
        'mass': 196.9734312 * u.u,
        'stable': False,
        'half-life': 486.0 * u.s,
    },

    'Pb-198': {
        'atomic number': 82,
        'mass number': 198,
        'mass': 197.972034 * u.u,
        'stable': False,
        'half-life': 8640.0 * u.s,
    },

    'Pb-199': {
        'atomic number': 82,
        'mass number': 199,
        'mass': 198.972913 * u.u,
        'stable': False,
        'half-life': 5400.0 * u.s,
    },

    'Pb-200': {
        'atomic number': 82,
        'mass number': 200,
        'mass': 199.971819 * u.u,
        'stable': False,
        'half-life': 77400.0 * u.s,
    },

    'Pb-201': {
        'atomic number': 82,
        'mass number': 201,
        'mass': 200.972883 * u.u,
        'stable': False,
        'half-life': 33588.0 * u.s,
    },

    'Pb-202': {
        'atomic number': 82,
        'mass number': 202,
        'mass': 201.972152 * u.u,
        'stable': False,
        'half-life': 1656738615000.0 * u.s,
    },

    'Pb-203': {
        'atomic number': 82,
        'mass number': 203,
        'mass': 202.9733911 * u.u,
        'stable': False,
        'half-life': 186922.80000000002 * u.s,
    },

    'Pb-204': {
        'atomic number': 82,
        'mass number': 204,
        'mass': 203.973044 * u.u,
        'stable': True,
        'abundance': 0.014,
    },

    'Pb-205': {
        'atomic number': 82,
        'mass number': 205,
        'mass': 204.9744822 * u.u,
        'stable': False,
        'half-life': 545934819800000.0 * u.s,
    },

    'Pb-206': {
        'atomic number': 82,
        'mass number': 206,
        'mass': 205.9744657 * u.u,
        'stable': True,
        'abundance': 0.241,
    },

    'Pb-207': {
        'atomic number': 82,
        'mass number': 207,
        'mass': 206.9758973 * u.u,
        'stable': True,
        'abundance': 0.221,
    },

    'Pb-208': {
        'atomic number': 82,
        'mass number': 208,
        'mass': 207.9766525 * u.u,
        'stable': True,
        'abundance': 0.524,
    },

    'Pb-209': {
        'atomic number': 82,
        'mass number': 209,
        'mass': 208.9810905 * u.u,
        'stable': False,
        'half-life': 11642.4 * u.s,
    },

    'Pb-210': {
        'atomic number': 82,
        'mass number': 210,
        'mass': 209.9841889 * u.u,
        'stable': False,
        'half-life': 700563757.2 * u.s,
    },

    'Pb-211': {
        'atomic number': 82,
        'mass number': 211,
        'mass': 210.9887371 * u.u,
        'stable': False,
        'half-life': 2169.84 * u.s,
    },

    'Pb-212': {
        'atomic number': 82,
        'mass number': 212,
        'mass': 211.9918977 * u.u,
        'stable': False,
        'half-life': 38304.0 * u.s,
    },

    'Pb-213': {
        'atomic number': 82,
        'mass number': 213,
        'mass': 212.9965629 * u.u,
        'stable': False,
        'half-life': 612.0 * u.s,
    },

    'Pb-214': {
        'atomic number': 82,
        'mass number': 214,
        'mass': 213.9998059 * u.u,
        'stable': False,
        'half-life': 1623.6 * u.s,
    },

    'Pb-215': {
        'atomic number': 82,
        'mass number': 215,
        'mass': 215.00474 * u.u,
        'stable': False,
        'half-life': 140.4 * u.s,
    },

    'Pb-216': {
        'atomic number': 82,
        'mass number': 216,
        'mass': 216.00803 * u.u,
        'stable': False,
        'half-life': 99.0 * u.s,
    },

    'Pb-217': {
        'atomic number': 82,
        'mass number': 217,
        'mass': 217.01314 * u.u,
        'stable': False,
        'half-life': 20.0 * u.s,
    },

    'Pb-218': {
        'atomic number': 82,
        'mass number': 218,
        'mass': 218.01659 * u.u,
        'stable': False,
        'half-life': 15.0 * u.s,
    },

    'Pb-219': {
        'atomic number': 82,
        'mass number': 219,
        'mass': 219.02177 * u.u,
        'stable': False,
        'half-life': '10# s',
    },

    'Pb-220': {
        'atomic number': 82,
        'mass number': 220,
        'mass': 220.02541 * u.u,
        'stable': False,
        'half-life': '30# s',
    },

    'Bi-184': {
        'atomic number': 83,
        'mass number': 184,
        'mass': 184.001275 * u.u,
        'stable': False,
        'half-life': 0.0066 * u.s,
    },

    'Bi-185': {
        'atomic number': 83,
        'mass number': 185,
        'mass': 184.9976 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'Bi-186': {
        'atomic number': 83,
        'mass number': 186,
        'mass': 185.996644 * u.u,
        'stable': False,
        'half-life': 0.0148 * u.s,
    },

    'Bi-187': {
        'atomic number': 83,
        'mass number': 187,
        'mass': 186.993147 * u.u,
        'stable': False,
        'half-life': 0.037 * u.s,
    },

    'Bi-188': {
        'atomic number': 83,
        'mass number': 188,
        'mass': 187.992287 * u.u,
        'stable': False,
        'half-life': 0.06 * u.s,
    },

    'Bi-189': {
        'atomic number': 83,
        'mass number': 189,
        'mass': 188.989195 * u.u,
        'stable': False,
        'half-life': 0.658 * u.s,
    },

    'Bi-190': {
        'atomic number': 83,
        'mass number': 190,
        'mass': 189.988622 * u.u,
        'stable': False,
        'half-life': 6.3 * u.s,
    },

    'Bi-191': {
        'atomic number': 83,
        'mass number': 191,
        'mass': 190.9857866 * u.u,
        'stable': False,
        'half-life': 11.7 * u.s,
    },

    'Bi-192': {
        'atomic number': 83,
        'mass number': 192,
        'mass': 191.985469 * u.u,
        'stable': False,
        'half-life': 34.6 * u.s,
    },

    'Bi-193': {
        'atomic number': 83,
        'mass number': 193,
        'mass': 192.98296 * u.u,
        'stable': False,
        'half-life': 63.6 * u.s,
    },

    'Bi-194': {
        'atomic number': 83,
        'mass number': 194,
        'mass': 193.982785 * u.u,
        'stable': False,
        'half-life': 95.0 * u.s,
    },

    'Bi-195': {
        'atomic number': 83,
        'mass number': 195,
        'mass': 194.9806488 * u.u,
        'stable': False,
        'half-life': 183.0 * u.s,
    },

    'Bi-196': {
        'atomic number': 83,
        'mass number': 196,
        'mass': 195.980667 * u.u,
        'stable': False,
        'half-life': 306.0 * u.s,
    },

    'Bi-197': {
        'atomic number': 83,
        'mass number': 197,
        'mass': 196.9788651 * u.u,
        'stable': False,
        'half-life': 559.8 * u.s,
    },

    'Bi-198': {
        'atomic number': 83,
        'mass number': 198,
        'mass': 197.979206 * u.u,
        'stable': False,
        'half-life': 618.0 * u.s,
    },

    'Bi-199': {
        'atomic number': 83,
        'mass number': 199,
        'mass': 198.977673 * u.u,
        'stable': False,
        'half-life': 1620.0 * u.s,
    },

    'Bi-200': {
        'atomic number': 83,
        'mass number': 200,
        'mass': 199.978131 * u.u,
        'stable': False,
        'half-life': 2184.0 * u.s,
    },

    'Bi-201': {
        'atomic number': 83,
        'mass number': 201,
        'mass': 200.97701 * u.u,
        'stable': False,
        'half-life': 6180.0 * u.s,
    },

    'Bi-202': {
        'atomic number': 83,
        'mass number': 202,
        'mass': 201.977734 * u.u,
        'stable': False,
        'half-life': 6192.0 * u.s,
    },

    'Bi-203': {
        'atomic number': 83,
        'mass number': 203,
        'mass': 202.976893 * u.u,
        'stable': False,
        'half-life': 42336.0 * u.s,
    },

    'Bi-204': {
        'atomic number': 83,
        'mass number': 204,
        'mass': 203.9778361 * u.u,
        'stable': False,
        'half-life': 40392.0 * u.s,
    },

    'Bi-205': {
        'atomic number': 83,
        'mass number': 205,
        'mass': 204.9773867 * u.u,
        'stable': False,
        'half-life': 1322784.0 * u.s,
    },

    'Bi-206': {
        'atomic number': 83,
        'mass number': 206,
        'mass': 205.9784993 * u.u,
        'stable': False,
        'half-life': 539395.2 * u.s,
    },

    'Bi-207': {
        'atomic number': 83,
        'mass number': 207,
        'mass': 206.978471 * u.u,
        'stable': False,
        'half-life': 995587200.0 * u.s,
    },

    'Bi-208': {
        'atomic number': 83,
        'mass number': 208,
        'mass': 207.9797425 * u.u,
        'stable': False,
        'half-life': 11612948768000.0 * u.s,
    },

    'Bi-209': {
        'atomic number': 83,
        'mass number': 209,
        'mass': 208.9803991 * u.u,
        'stable': False,
        'abundance': 1,
    },

    'Bi-210': {
        'atomic number': 83,
        'mass number': 210,
        'mass': 209.9841207 * u.u,
        'stable': False,
        'half-life': 433036.8 * u.s,
    },

    'Bi-211': {
        'atomic number': 83,
        'mass number': 211,
        'mass': 210.9872697 * u.u,
        'stable': False,
        'half-life': 128.4 * u.s,
    },

    'Bi-212': {
        'atomic number': 83,
        'mass number': 212,
        'mass': 211.991286 * u.u,
        'stable': False,
        'half-life': 3633.0 * u.s,
    },

    'Bi-213': {
        'atomic number': 83,
        'mass number': 213,
        'mass': 212.9943851 * u.u,
        'stable': False,
        'half-life': 2736.6 * u.s,
    },

    'Bi-214': {
        'atomic number': 83,
        'mass number': 214,
        'mass': 213.998712 * u.u,
        'stable': False,
        'half-life': 1194.0 * u.s,
    },

    'Bi-215': {
        'atomic number': 83,
        'mass number': 215,
        'mass': 215.00177 * u.u,
        'stable': False,
        'half-life': 456.0 * u.s,
    },

    'Bi-216': {
        'atomic number': 83,
        'mass number': 216,
        'mass': 216.006306 * u.u,
        'stable': False,
        'half-life': 135.0 * u.s,
    },

    'Bi-217': {
        'atomic number': 83,
        'mass number': 217,
        'mass': 217.009372 * u.u,
        'stable': False,
        'half-life': 98.5 * u.s,
    },

    'Bi-218': {
        'atomic number': 83,
        'mass number': 218,
        'mass': 218.014188 * u.u,
        'stable': False,
        'half-life': 33.0 * u.s,
    },

    'Bi-219': {
        'atomic number': 83,
        'mass number': 219,
        'mass': 219.01748 * u.u,
        'stable': False,
        'half-life': 8.7 * u.s,
    },

    'Bi-220': {
        'atomic number': 83,
        'mass number': 220,
        'mass': 220.02235 * u.u,
        'stable': False,
        'half-life': 9.5 * u.s,
    },

    'Bi-221': {
        'atomic number': 83,
        'mass number': 221,
        'mass': 221.02587 * u.u,
        'stable': False,
        'half-life': '5# s',
    },

    'Bi-222': {
        'atomic number': 83,
        'mass number': 222,
        'mass': 222.03078 * u.u,
        'stable': False,
        'half-life': '2# s',
    },

    'Bi-223': {
        'atomic number': 83,
        'mass number': 223,
        'mass': 223.0345 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Bi-224': {
        'atomic number': 83,
        'mass number': 224,
        'mass': 224.03947 * u.u,
        'stable': False,
        'half-life': '300# ms',
    },

    'Po-186': {
        'atomic number': 84,
        'mass number': 186,
        'mass': 186.004393 * u.u,
        'stable': False,
        'half-life': 3.4e-05 * u.s,
    },

    'Po-187': {
        'atomic number': 84,
        'mass number': 187,
        'mass': 187.003041 * u.u,
        'stable': False,
        'half-life': 0.0014 * u.s,
    },

    'Po-188': {
        'atomic number': 84,
        'mass number': 188,
        'mass': 187.999416 * u.u,
        'stable': False,
        'half-life': 0.000275 * u.s,
    },

    'Po-189': {
        'atomic number': 84,
        'mass number': 189,
        'mass': 188.998473 * u.u,
        'stable': False,
        'half-life': 0.0038 * u.s,
    },

    'Po-190': {
        'atomic number': 84,
        'mass number': 190,
        'mass': 189.995101 * u.u,
        'stable': False,
        'half-life': 0.00246 * u.s,
    },

    'Po-191': {
        'atomic number': 84,
        'mass number': 191,
        'mass': 190.9945585 * u.u,
        'stable': False,
        'half-life': 0.022 * u.s,
    },

    'Po-192': {
        'atomic number': 84,
        'mass number': 192,
        'mass': 191.991336 * u.u,
        'stable': False,
        'half-life': 0.0322 * u.s,
    },

    'Po-193': {
        'atomic number': 84,
        'mass number': 193,
        'mass': 192.991026 * u.u,
        'stable': False,
        'half-life': 0.388 * u.s,
    },

    'Po-194': {
        'atomic number': 84,
        'mass number': 194,
        'mass': 193.988186 * u.u,
        'stable': False,
        'half-life': 0.392 * u.s,
    },

    'Po-195': {
        'atomic number': 84,
        'mass number': 195,
        'mass': 194.988126 * u.u,
        'stable': False,
        'half-life': 4.64 * u.s,
    },

    'Po-196': {
        'atomic number': 84,
        'mass number': 196,
        'mass': 195.985526 * u.u,
        'stable': False,
        'half-life': 5.56 * u.s,
    },

    'Po-197': {
        'atomic number': 84,
        'mass number': 197,
        'mass': 196.98566 * u.u,
        'stable': False,
        'half-life': 53.6 * u.s,
    },

    'Po-198': {
        'atomic number': 84,
        'mass number': 198,
        'mass': 197.983389 * u.u,
        'stable': False,
        'half-life': 105.6 * u.s,
    },

    'Po-199': {
        'atomic number': 84,
        'mass number': 199,
        'mass': 198.983667 * u.u,
        'stable': False,
        'half-life': 328.2 * u.s,
    },

    'Po-200': {
        'atomic number': 84,
        'mass number': 200,
        'mass': 199.981799 * u.u,
        'stable': False,
        'half-life': 690.6 * u.s,
    },

    'Po-201': {
        'atomic number': 84,
        'mass number': 201,
        'mass': 200.9822598 * u.u,
        'stable': False,
        'half-life': 936.0 * u.s,
    },

    'Po-202': {
        'atomic number': 84,
        'mass number': 202,
        'mass': 201.980758 * u.u,
        'stable': False,
        'half-life': 2676.0 * u.s,
    },

    'Po-203': {
        'atomic number': 84,
        'mass number': 203,
        'mass': 202.9814161 * u.u,
        'stable': False,
        'half-life': 2202.0 * u.s,
    },

    'Po-204': {
        'atomic number': 84,
        'mass number': 204,
        'mass': 203.98031 * u.u,
        'stable': False,
        'half-life': 12668.4 * u.s,
    },

    'Po-205': {
        'atomic number': 84,
        'mass number': 205,
        'mass': 204.981203 * u.u,
        'stable': False,
        'half-life': 6264.0 * u.s,
    },

    'Po-206': {
        'atomic number': 84,
        'mass number': 206,
        'mass': 205.980474 * u.u,
        'stable': False,
        'half-life': 760320.0 * u.s,
    },

    'Po-207': {
        'atomic number': 84,
        'mass number': 207,
        'mass': 206.9815938 * u.u,
        'stable': False,
        'half-life': 20880.0 * u.s,
    },

    'Po-208': {
        'atomic number': 84,
        'mass number': 208,
        'mass': 207.9812461 * u.u,
        'stable': False,
        'half-life': 91451971.548 * u.s,
    },

    'Po-209': {
        'atomic number': 84,
        'mass number': 209,
        'mass': 208.9824308 * u.u,
        'stable': False,
        'half-life': 3913058824.0 * u.s,
    },

    'Po-210': {
        'atomic number': 84,
        'mass number': 210,
        'mass': 209.9828741 * u.u,
        'stable': False,
        'half-life': 11955686.4 * u.s,
    },

    'Po-211': {
        'atomic number': 84,
        'mass number': 211,
        'mass': 210.9866536 * u.u,
        'stable': False,
        'half-life': 0.516 * u.s,
    },

    'Po-212': {
        'atomic number': 84,
        'mass number': 212,
        'mass': 211.9888684 * u.u,
        'stable': False,
        'half-life': 2.947e-07 * u.s,
    },

    'Po-213': {
        'atomic number': 84,
        'mass number': 213,
        'mass': 212.9928576 * u.u,
        'stable': False,
        'half-life': 3.708e-06 * u.s,
    },

    'Po-214': {
        'atomic number': 84,
        'mass number': 214,
        'mass': 213.9952017 * u.u,
        'stable': False,
        'half-life': 0.00016372 * u.s,
    },

    'Po-215': {
        'atomic number': 84,
        'mass number': 215,
        'mass': 214.9994201 * u.u,
        'stable': False,
        'half-life': 0.001781 * u.s,
    },

    'Po-216': {
        'atomic number': 84,
        'mass number': 216,
        'mass': 216.0019152 * u.u,
        'stable': False,
        'half-life': 0.145 * u.s,
    },

    'Po-217': {
        'atomic number': 84,
        'mass number': 217,
        'mass': 217.0063182 * u.u,
        'stable': False,
        'half-life': 1.514 * u.s,
    },

    'Po-218': {
        'atomic number': 84,
        'mass number': 218,
        'mass': 218.0089735 * u.u,
        'stable': False,
        'half-life': 185.88 * u.s,
    },

    'Po-219': {
        'atomic number': 84,
        'mass number': 219,
        'mass': 219.013614 * u.u,
        'stable': False,
        'half-life': 618.0 * u.s,
    },

    'Po-220': {
        'atomic number': 84,
        'mass number': 220,
        'mass': 220.016386 * u.u,
        'stable': False,
        'half-life': '40# s',
    },

    'Po-221': {
        'atomic number': 84,
        'mass number': 221,
        'mass': 221.021228 * u.u,
        'stable': False,
        'half-life': 132.0 * u.s,
    },

    'Po-222': {
        'atomic number': 84,
        'mass number': 222,
        'mass': 222.02414 * u.u,
        'stable': False,
        'half-life': 546.0 * u.s,
    },

    'Po-223': {
        'atomic number': 84,
        'mass number': 223,
        'mass': 223.02907 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Po-224': {
        'atomic number': 84,
        'mass number': 224,
        'mass': 224.03211 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Po-225': {
        'atomic number': 84,
        'mass number': 225,
        'mass': 225.03707 * u.u,
        'stable': False,
        'half-life': '20# s',
    },

    'Po-226': {
        'atomic number': 84,
        'mass number': 226,
        'mass': 226.04031 * u.u,
        'stable': False,
        'half-life': '20# s',
    },

    'Po-227': {
        'atomic number': 84,
        'mass number': 227,
        'mass': 227.04539 * u.u,
        'stable': False,
        'half-life': '5# s',
    },

    'At-191': {
        'atomic number': 85,
        'mass number': 191,
        'mass': 191.004148 * u.u,
        'stable': False,
        'half-life': 0.0021 * u.s,
    },

    'At-192': {
        'atomic number': 85,
        'mass number': 192,
        'mass': 192.003152 * u.u,
        'stable': False,
        'half-life': 0.0115 * u.s,
    },

    'At-193': {
        'atomic number': 85,
        'mass number': 193,
        'mass': 192.999927 * u.u,
        'stable': False,
        'half-life': 0.029 * u.s,
    },

    'At-194': {
        'atomic number': 85,
        'mass number': 194,
        'mass': 193.999236 * u.u,
        'stable': False,
        'half-life': 0.286 * u.s,
    },

    'At-195': {
        'atomic number': 85,
        'mass number': 195,
        'mass': 194.9962685 * u.u,
        'stable': False,
        'half-life': 0.29 * u.s,
    },

    'At-196': {
        'atomic number': 85,
        'mass number': 196,
        'mass': 195.9958 * u.u,
        'stable': False,
        'half-life': 0.388 * u.s,
    },

    'At-197': {
        'atomic number': 85,
        'mass number': 197,
        'mass': 196.993189 * u.u,
        'stable': False,
        'half-life': 0.3882 * u.s,
    },

    'At-198': {
        'atomic number': 85,
        'mass number': 198,
        'mass': 197.992784 * u.u,
        'stable': False,
        'half-life': 3.0 * u.s,
    },

    'At-199': {
        'atomic number': 85,
        'mass number': 199,
        'mass': 198.9905277 * u.u,
        'stable': False,
        'half-life': 7.02 * u.s,
    },

    'At-200': {
        'atomic number': 85,
        'mass number': 200,
        'mass': 199.990351 * u.u,
        'stable': False,
        'half-life': 43.2 * u.s,
    },

    'At-201': {
        'atomic number': 85,
        'mass number': 201,
        'mass': 200.9884171 * u.u,
        'stable': False,
        'half-life': 85.2 * u.s,
    },

    'At-202': {
        'atomic number': 85,
        'mass number': 202,
        'mass': 201.98863 * u.u,
        'stable': False,
        'half-life': 184.0 * u.s,
    },

    'At-203': {
        'atomic number': 85,
        'mass number': 203,
        'mass': 202.986943 * u.u,
        'stable': False,
        'half-life': 444.0 * u.s,
    },

    'At-204': {
        'atomic number': 85,
        'mass number': 204,
        'mass': 203.987251 * u.u,
        'stable': False,
        'half-life': 547.2 * u.s,
    },

    'At-205': {
        'atomic number': 85,
        'mass number': 205,
        'mass': 204.986076 * u.u,
        'stable': False,
        'half-life': 2028.0 * u.s,
    },

    'At-206': {
        'atomic number': 85,
        'mass number': 206,
        'mass': 205.986657 * u.u,
        'stable': False,
        'half-life': 1836.0 * u.s,
    },

    'At-207': {
        'atomic number': 85,
        'mass number': 207,
        'mass': 206.9858 * u.u,
        'stable': False,
        'half-life': 6516.0 * u.s,
    },

    'At-208': {
        'atomic number': 85,
        'mass number': 208,
        'mass': 207.9866133 * u.u,
        'stable': False,
        'half-life': 5868.0 * u.s,
    },

    'At-209': {
        'atomic number': 85,
        'mass number': 209,
        'mass': 208.9861702 * u.u,
        'stable': False,
        'half-life': 19512.0 * u.s,
    },

    'At-210': {
        'atomic number': 85,
        'mass number': 210,
        'mass': 209.9871479 * u.u,
        'stable': False,
        'half-life': 29160.0 * u.s,
    },

    'At-211': {
        'atomic number': 85,
        'mass number': 211,
        'mass': 210.9874966 * u.u,
        'stable': False,
        'half-life': 25970.4 * u.s,
    },

    'At-212': {
        'atomic number': 85,
        'mass number': 212,
        'mass': 211.9907377 * u.u,
        'stable': False,
        'half-life': 0.314 * u.s,
    },

    'At-213': {
        'atomic number': 85,
        'mass number': 213,
        'mass': 212.992937 * u.u,
        'stable': False,
        'half-life': 1.25e-07 * u.s,
    },

    'At-214': {
        'atomic number': 85,
        'mass number': 214,
        'mass': 213.9963721 * u.u,
        'stable': False,
        'half-life': 5.58e-07 * u.s,
    },

    'At-215': {
        'atomic number': 85,
        'mass number': 215,
        'mass': 214.9986528 * u.u,
        'stable': False,
        'half-life': 0.0001 * u.s,
    },

    'At-216': {
        'atomic number': 85,
        'mass number': 216,
        'mass': 216.0024236 * u.u,
        'stable': False,
        'half-life': 0.0003 * u.s,
    },

    'At-217': {
        'atomic number': 85,
        'mass number': 217,
        'mass': 217.0047192 * u.u,
        'stable': False,
        'half-life': 0.03262 * u.s,
    },

    'At-218': {
        'atomic number': 85,
        'mass number': 218,
        'mass': 218.008695 * u.u,
        'stable': False,
        'half-life': 1.5 * u.s,
    },

    'At-219': {
        'atomic number': 85,
        'mass number': 219,
        'mass': 219.0111618 * u.u,
        'stable': False,
        'half-life': 56.0 * u.s,
    },

    'At-220': {
        'atomic number': 85,
        'mass number': 220,
        'mass': 220.015433 * u.u,
        'stable': False,
        'half-life': 222.6 * u.s,
    },

    'At-221': {
        'atomic number': 85,
        'mass number': 221,
        'mass': 221.018017 * u.u,
        'stable': False,
        'half-life': 138.0 * u.s,
    },

    'At-222': {
        'atomic number': 85,
        'mass number': 222,
        'mass': 222.022494 * u.u,
        'stable': False,
        'half-life': 54.0 * u.s,
    },

    'At-223': {
        'atomic number': 85,
        'mass number': 223,
        'mass': 223.025151 * u.u,
        'stable': False,
        'half-life': 50.0 * u.s,
    },

    'At-224': {
        'atomic number': 85,
        'mass number': 224,
        'mass': 224.029749 * u.u,
        'stable': False,
        'half-life': 150.0 * u.s,
    },

    'At-225': {
        'atomic number': 85,
        'mass number': 225,
        'mass': 225.03263 * u.u,
        'stable': False,
        'half-life': '2# m',
    },

    'At-226': {
        'atomic number': 85,
        'mass number': 226,
        'mass': 226.03716 * u.u,
        'stable': False,
        'half-life': '20# s',
    },

    'At-227': {
        'atomic number': 85,
        'mass number': 227,
        'mass': 227.04024 * u.u,
        'stable': False,
        'half-life': '20# s',
    },

    'At-228': {
        'atomic number': 85,
        'mass number': 228,
        'mass': 228.04475 * u.u,
        'stable': False,
        'half-life': '5# s',
    },

    'At-229': {
        'atomic number': 85,
        'mass number': 229,
        'mass': 229.04812 * u.u,
        'stable': False,
        'half-life': '5# s',
    },

    'Rn-193': {
        'atomic number': 86,
        'mass number': 193,
        'mass': 193.009708 * u.u,
        'stable': False,
        'half-life': 0.00115 * u.s,
    },

    'Rn-194': {
        'atomic number': 86,
        'mass number': 194,
        'mass': 194.006144 * u.u,
        'stable': False,
        'half-life': 0.00078 * u.s,
    },

    'Rn-195': {
        'atomic number': 86,
        'mass number': 195,
        'mass': 195.005422 * u.u,
        'stable': False,
        'half-life': 0.007 * u.s,
    },

    'Rn-196': {
        'atomic number': 86,
        'mass number': 196,
        'mass': 196.002116 * u.u,
        'stable': False,
        'half-life': 0.0047 * u.s,
    },

    'Rn-197': {
        'atomic number': 86,
        'mass number': 197,
        'mass': 197.001585 * u.u,
        'stable': False,
        'half-life': 0.054 * u.s,
    },

    'Rn-198': {
        'atomic number': 86,
        'mass number': 198,
        'mass': 197.998679 * u.u,
        'stable': False,
        'half-life': 0.065 * u.s,
    },

    'Rn-199': {
        'atomic number': 86,
        'mass number': 199,
        'mass': 198.99839 * u.u,
        'stable': False,
        'half-life': 0.59 * u.s,
    },

    'Rn-200': {
        'atomic number': 86,
        'mass number': 200,
        'mass': 199.99569 * u.u,
        'stable': False,
        'half-life': 1.09 * u.s,
    },

    'Rn-201': {
        'atomic number': 86,
        'mass number': 201,
        'mass': 200.995628 * u.u,
        'stable': False,
        'half-life': 7.0 * u.s,
    },

    'Rn-202': {
        'atomic number': 86,
        'mass number': 202,
        'mass': 201.993264 * u.u,
        'stable': False,
        'half-life': 9.7 * u.s,
    },

    'Rn-203': {
        'atomic number': 86,
        'mass number': 203,
        'mass': 202.993388 * u.u,
        'stable': False,
        'half-life': 44.0 * u.s,
    },

    'Rn-204': {
        'atomic number': 86,
        'mass number': 204,
        'mass': 203.99143 * u.u,
        'stable': False,
        'half-life': 74.52 * u.s,
    },

    'Rn-205': {
        'atomic number': 86,
        'mass number': 205,
        'mass': 204.991719 * u.u,
        'stable': False,
        'half-life': 169.8 * u.s,
    },

    'Rn-206': {
        'atomic number': 86,
        'mass number': 206,
        'mass': 205.990214 * u.u,
        'stable': False,
        'half-life': 340.2 * u.s,
    },

    'Rn-207': {
        'atomic number': 86,
        'mass number': 207,
        'mass': 206.9907303 * u.u,
        'stable': False,
        'half-life': 555.0 * u.s,
    },

    'Rn-208': {
        'atomic number': 86,
        'mass number': 208,
        'mass': 207.989635 * u.u,
        'stable': False,
        'half-life': 1461.0 * u.s,
    },

    'Rn-209': {
        'atomic number': 86,
        'mass number': 209,
        'mass': 208.990415 * u.u,
        'stable': False,
        'half-life': 1728.0 * u.s,
    },

    'Rn-210': {
        'atomic number': 86,
        'mass number': 210,
        'mass': 209.9896891 * u.u,
        'stable': False,
        'half-life': 8640.0 * u.s,
    },

    'Rn-211': {
        'atomic number': 86,
        'mass number': 211,
        'mass': 210.9906011 * u.u,
        'stable': False,
        'half-life': 52560.0 * u.s,
    },

    'Rn-212': {
        'atomic number': 86,
        'mass number': 212,
        'mass': 211.9907039 * u.u,
        'stable': False,
        'half-life': 1434.0 * u.s,
    },

    'Rn-213': {
        'atomic number': 86,
        'mass number': 213,
        'mass': 212.9938831 * u.u,
        'stable': False,
        'half-life': 0.0195 * u.s,
    },

    'Rn-214': {
        'atomic number': 86,
        'mass number': 214,
        'mass': 213.995363 * u.u,
        'stable': False,
        'half-life': 2.7e-07 * u.s,
    },

    'Rn-215': {
        'atomic number': 86,
        'mass number': 215,
        'mass': 214.9987459 * u.u,
        'stable': False,
        'half-life': 2.3e-06 * u.s,
    },

    'Rn-216': {
        'atomic number': 86,
        'mass number': 216,
        'mass': 216.0002719 * u.u,
        'stable': False,
        'half-life': 4.5e-05 * u.s,
    },

    'Rn-217': {
        'atomic number': 86,
        'mass number': 217,
        'mass': 217.003928 * u.u,
        'stable': False,
        'half-life': 0.00054 * u.s,
    },

    'Rn-218': {
        'atomic number': 86,
        'mass number': 218,
        'mass': 218.0056016 * u.u,
        'stable': False,
        'half-life': 0.03375 * u.s,
    },

    'Rn-219': {
        'atomic number': 86,
        'mass number': 219,
        'mass': 219.0094804 * u.u,
        'stable': False,
        'half-life': 3.96 * u.s,
    },

    'Rn-220': {
        'atomic number': 86,
        'mass number': 220,
        'mass': 220.0113941 * u.u,
        'stable': False,
        'half-life': 55.6 * u.s,
    },

    'Rn-221': {
        'atomic number': 86,
        'mass number': 221,
        'mass': 221.0155371 * u.u,
        'stable': False,
        'half-life': 1542.0 * u.s,
    },

    'Rn-222': {
        'atomic number': 86,
        'mass number': 222,
        'mass': 222.0175782 * u.u,
        'stable': False,
        'half-life': 330177.6 * u.s,
    },

    'Rn-223': {
        'atomic number': 86,
        'mass number': 223,
        'mass': 223.0218893 * u.u,
        'stable': False,
        'half-life': 1458.0 * u.s,
    },

    'Rn-224': {
        'atomic number': 86,
        'mass number': 224,
        'mass': 224.024096 * u.u,
        'stable': False,
        'half-life': 6420.0 * u.s,
    },

    'Rn-225': {
        'atomic number': 86,
        'mass number': 225,
        'mass': 225.028486 * u.u,
        'stable': False,
        'half-life': 279.6 * u.s,
    },

    'Rn-226': {
        'atomic number': 86,
        'mass number': 226,
        'mass': 226.030861 * u.u,
        'stable': False,
        'half-life': 444.0 * u.s,
    },

    'Rn-227': {
        'atomic number': 86,
        'mass number': 227,
        'mass': 227.035304 * u.u,
        'stable': False,
        'half-life': 20.2 * u.s,
    },

    'Rn-228': {
        'atomic number': 86,
        'mass number': 228,
        'mass': 228.037835 * u.u,
        'stable': False,
        'half-life': 65.0 * u.s,
    },

    'Rn-229': {
        'atomic number': 86,
        'mass number': 229,
        'mass': 229.042257 * u.u,
        'stable': False,
        'half-life': 11.9 * u.s,
    },

    'Rn-230': {
        'atomic number': 86,
        'mass number': 230,
        'mass': 230.04514 * u.u,
        'stable': False,
        'half-life': '10# s',
    },

    'Rn-231': {
        'atomic number': 86,
        'mass number': 231,
        'mass': 231.04987 * u.u,
        'stable': False,
        'half-life': '300# ms',
    },

    'Fr-199': {
        'atomic number': 87,
        'mass number': 199,
        'mass': 199.007259 * u.u,
        'stable': False,
        'half-life': 0.0066 * u.s,
    },

    'Fr-200': {
        'atomic number': 87,
        'mass number': 200,
        'mass': 200.006586 * u.u,
        'stable': False,
        'half-life': 0.0475 * u.s,
    },

    'Fr-201': {
        'atomic number': 87,
        'mass number': 201,
        'mass': 201.003867 * u.u,
        'stable': False,
        'half-life': 0.0628 * u.s,
    },

    'Fr-202': {
        'atomic number': 87,
        'mass number': 202,
        'mass': 202.00332 * u.u,
        'stable': False,
        'half-life': 0.372 * u.s,
    },

    'Fr-203': {
        'atomic number': 87,
        'mass number': 203,
        'mass': 203.0009407 * u.u,
        'stable': False,
        'half-life': 0.55 * u.s,
    },

    'Fr-204': {
        'atomic number': 87,
        'mass number': 204,
        'mass': 204.000652 * u.u,
        'stable': False,
        'half-life': 1.75 * u.s,
    },

    'Fr-205': {
        'atomic number': 87,
        'mass number': 205,
        'mass': 204.9985939 * u.u,
        'stable': False,
        'half-life': 3.82 * u.s,
    },

    'Fr-206': {
        'atomic number': 87,
        'mass number': 206,
        'mass': 205.998666 * u.u,
        'stable': False,
        'half-life': '~16 s',
    },

    'Fr-207': {
        'atomic number': 87,
        'mass number': 207,
        'mass': 206.996946 * u.u,
        'stable': False,
        'half-life': 14.8 * u.s,
    },

    'Fr-208': {
        'atomic number': 87,
        'mass number': 208,
        'mass': 207.997138 * u.u,
        'stable': False,
        'half-life': 59.1 * u.s,
    },

    'Fr-209': {
        'atomic number': 87,
        'mass number': 209,
        'mass': 208.995955 * u.u,
        'stable': False,
        'half-life': 50.5 * u.s,
    },

    'Fr-210': {
        'atomic number': 87,
        'mass number': 210,
        'mass': 209.996422 * u.u,
        'stable': False,
        'half-life': 190.8 * u.s,
    },

    'Fr-211': {
        'atomic number': 87,
        'mass number': 211,
        'mass': 210.995556 * u.u,
        'stable': False,
        'half-life': 186.0 * u.s,
    },

    'Fr-212': {
        'atomic number': 87,
        'mass number': 212,
        'mass': 211.9962257 * u.u,
        'stable': False,
        'half-life': 1200.0 * u.s,
    },

    'Fr-213': {
        'atomic number': 87,
        'mass number': 213,
        'mass': 212.996186 * u.u,
        'stable': False,
        'half-life': 34.14 * u.s,
    },

    'Fr-214': {
        'atomic number': 87,
        'mass number': 214,
        'mass': 213.9989713 * u.u,
        'stable': False,
        'half-life': 0.00518 * u.s,
    },

    'Fr-215': {
        'atomic number': 87,
        'mass number': 215,
        'mass': 215.0003418 * u.u,
        'stable': False,
        'half-life': 8.6e-08 * u.s,
    },

    'Fr-216': {
        'atomic number': 87,
        'mass number': 216,
        'mass': 216.0031899 * u.u,
        'stable': False,
        'half-life': 7e-07 * u.s,
    },

    'Fr-217': {
        'atomic number': 87,
        'mass number': 217,
        'mass': 217.0046323 * u.u,
        'stable': False,
        'half-life': 1.68e-05 * u.s,
    },

    'Fr-218': {
        'atomic number': 87,
        'mass number': 218,
        'mass': 218.0075787 * u.u,
        'stable': False,
        'half-life': 0.001 * u.s,
    },

    'Fr-219': {
        'atomic number': 87,
        'mass number': 219,
        'mass': 219.0092524 * u.u,
        'stable': False,
        'half-life': 0.02 * u.s,
    },

    'Fr-220': {
        'atomic number': 87,
        'mass number': 220,
        'mass': 220.0123277 * u.u,
        'stable': False,
        'half-life': 27.4 * u.s,
    },

    'Fr-221': {
        'atomic number': 87,
        'mass number': 221,
        'mass': 221.0142552 * u.u,
        'stable': False,
        'half-life': 288.06 * u.s,
    },

    'Fr-222': {
        'atomic number': 87,
        'mass number': 222,
        'mass': 222.017552 * u.u,
        'stable': False,
        'half-life': 852.0 * u.s,
    },

    'Fr-223': {
        'atomic number': 87,
        'mass number': 223,
        'mass': 223.019736 * u.u,
        'stable': False,
        'half-life': 1320.0 * u.s,
    },

    'Fr-224': {
        'atomic number': 87,
        'mass number': 224,
        'mass': 224.023398 * u.u,
        'stable': False,
        'half-life': 199.8 * u.s,
    },

    'Fr-225': {
        'atomic number': 87,
        'mass number': 225,
        'mass': 225.025573 * u.u,
        'stable': False,
        'half-life': 237.0 * u.s,
    },

    'Fr-226': {
        'atomic number': 87,
        'mass number': 226,
        'mass': 226.029566 * u.u,
        'stable': False,
        'half-life': 49.0 * u.s,
    },

    'Fr-227': {
        'atomic number': 87,
        'mass number': 227,
        'mass': 227.031869 * u.u,
        'stable': False,
        'half-life': 148.2 * u.s,
    },

    'Fr-228': {
        'atomic number': 87,
        'mass number': 228,
        'mass': 228.035823 * u.u,
        'stable': False,
        'half-life': 38.0 * u.s,
    },

    'Fr-229': {
        'atomic number': 87,
        'mass number': 229,
        'mass': 229.038298 * u.u,
        'stable': False,
        'half-life': 50.2 * u.s,
    },

    'Fr-230': {
        'atomic number': 87,
        'mass number': 230,
        'mass': 230.042416 * u.u,
        'stable': False,
        'half-life': 19.1 * u.s,
    },

    'Fr-231': {
        'atomic number': 87,
        'mass number': 231,
        'mass': 231.045158 * u.u,
        'stable': False,
        'half-life': 17.6 * u.s,
    },

    'Fr-232': {
        'atomic number': 87,
        'mass number': 232,
        'mass': 232.04937 * u.u,
        'stable': False,
        'half-life': 5.5 * u.s,
    },

    'Fr-233': {
        'atomic number': 87,
        'mass number': 233,
        'mass': 233.05264 * u.u,
        'stable': False,
        'half-life': 0.9 * u.s,
    },

    'Ra-201': {
        'atomic number': 88,
        'mass number': 201,
        'mass': 201.01271 * u.u,
        'stable': False,
        'half-life': 0.02 * u.s,
    },

    'Ra-202': {
        'atomic number': 88,
        'mass number': 202,
        'mass': 202.00976 * u.u,
        'stable': False,
        'half-life': 0.0041 * u.s,
    },

    'Ra-203': {
        'atomic number': 88,
        'mass number': 203,
        'mass': 203.009304 * u.u,
        'stable': False,
        'half-life': 0.036 * u.s,
    },

    'Ra-204': {
        'atomic number': 88,
        'mass number': 204,
        'mass': 204.006492 * u.u,
        'stable': False,
        'half-life': 0.06 * u.s,
    },

    'Ra-205': {
        'atomic number': 88,
        'mass number': 205,
        'mass': 205.006268 * u.u,
        'stable': False,
        'half-life': 0.22 * u.s,
    },

    'Ra-206': {
        'atomic number': 88,
        'mass number': 206,
        'mass': 206.003828 * u.u,
        'stable': False,
        'half-life': 0.24 * u.s,
    },

    'Ra-207': {
        'atomic number': 88,
        'mass number': 207,
        'mass': 207.003799 * u.u,
        'stable': False,
        'half-life': 1.38 * u.s,
    },

    'Ra-208': {
        'atomic number': 88,
        'mass number': 208,
        'mass': 208.001841 * u.u,
        'stable': False,
        'half-life': 1.11 * u.s,
    },

    'Ra-209': {
        'atomic number': 88,
        'mass number': 209,
        'mass': 209.00199 * u.u,
        'stable': False,
        'half-life': 4.71 * u.s,
    },

    'Ra-210': {
        'atomic number': 88,
        'mass number': 210,
        'mass': 210.000494 * u.u,
        'stable': False,
        'half-life': 4.0 * u.s,
    },

    'Ra-211': {
        'atomic number': 88,
        'mass number': 211,
        'mass': 211.0008932 * u.u,
        'stable': False,
        'half-life': 13.2 * u.s,
    },

    'Ra-212': {
        'atomic number': 88,
        'mass number': 212,
        'mass': 211.999787 * u.u,
        'stable': False,
        'half-life': 13.0 * u.s,
    },

    'Ra-213': {
        'atomic number': 88,
        'mass number': 213,
        'mass': 213.000384 * u.u,
        'stable': False,
        'half-life': 163.8 * u.s,
    },

    'Ra-214': {
        'atomic number': 88,
        'mass number': 214,
        'mass': 214.0000997 * u.u,
        'stable': False,
        'half-life': 2.437 * u.s,
    },

    'Ra-215': {
        'atomic number': 88,
        'mass number': 215,
        'mass': 215.0027204 * u.u,
        'stable': False,
        'half-life': 0.00167 * u.s,
    },

    'Ra-216': {
        'atomic number': 88,
        'mass number': 216,
        'mass': 216.0035334 * u.u,
        'stable': False,
        'half-life': 1.82e-07 * u.s,
    },

    'Ra-217': {
        'atomic number': 88,
        'mass number': 217,
        'mass': 217.0063207 * u.u,
        'stable': False,
        'half-life': 1.63e-06 * u.s,
    },

    'Ra-218': {
        'atomic number': 88,
        'mass number': 218,
        'mass': 218.007141 * u.u,
        'stable': False,
        'half-life': 2.52e-05 * u.s,
    },

    'Ra-219': {
        'atomic number': 88,
        'mass number': 219,
        'mass': 219.0100855 * u.u,
        'stable': False,
        'half-life': 0.01 * u.s,
    },

    'Ra-220': {
        'atomic number': 88,
        'mass number': 220,
        'mass': 220.0110259 * u.u,
        'stable': False,
        'half-life': 0.0179 * u.s,
    },

    'Ra-221': {
        'atomic number': 88,
        'mass number': 221,
        'mass': 221.0139177 * u.u,
        'stable': False,
        'half-life': 28.0 * u.s,
    },

    'Ra-222': {
        'atomic number': 88,
        'mass number': 222,
        'mass': 222.0153748 * u.u,
        'stable': False,
        'half-life': 33.6 * u.s,
    },

    'Ra-223': {
        'atomic number': 88,
        'mass number': 223,
        'mass': 223.0185023 * u.u,
        'stable': False,
        'half-life': 988217.28 * u.s,
    },

    'Ra-224': {
        'atomic number': 88,
        'mass number': 224,
        'mass': 224.020212 * u.u,
        'stable': False,
        'half-life': 313796.16 * u.s,
    },

    'Ra-225': {
        'atomic number': 88,
        'mass number': 225,
        'mass': 225.0236119 * u.u,
        'stable': False,
        'half-life': 1287360.0 * u.s,
    },

    'Ra-226': {
        'atomic number': 88,
        'mass number': 226,
        'mass': 226.0254103 * u.u,
        'stable': False,
        'half-life': 50491081600.0 * u.s,
    },

    'Ra-227': {
        'atomic number': 88,
        'mass number': 227,
        'mass': 227.0291783 * u.u,
        'stable': False,
        'half-life': 2532.0 * u.s,
    },

    'Ra-228': {
        'atomic number': 88,
        'mass number': 228,
        'mass': 228.0310707 * u.u,
        'stable': False,
        'half-life': 181452324.5 * u.s,
    },

    'Ra-229': {
        'atomic number': 88,
        'mass number': 229,
        'mass': 229.034942 * u.u,
        'stable': False,
        'half-life': 240.0 * u.s,
    },

    'Ra-230': {
        'atomic number': 88,
        'mass number': 230,
        'mass': 230.037055 * u.u,
        'stable': False,
        'half-life': 5580.0 * u.s,
    },

    'Ra-231': {
        'atomic number': 88,
        'mass number': 231,
        'mass': 231.041027 * u.u,
        'stable': False,
        'half-life': 104.0 * u.s,
    },

    'Ra-232': {
        'atomic number': 88,
        'mass number': 232,
        'mass': 232.0434753 * u.u,
        'stable': False,
        'half-life': 240.0 * u.s,
    },

    'Ra-233': {
        'atomic number': 88,
        'mass number': 233,
        'mass': 233.047582 * u.u,
        'stable': False,
        'half-life': 30.0 * u.s,
    },

    'Ra-234': {
        'atomic number': 88,
        'mass number': 234,
        'mass': 234.050342 * u.u,
        'stable': False,
        'half-life': 30.0 * u.s,
    },

    'Ra-235': {
        'atomic number': 88,
        'mass number': 235,
        'mass': 235.05497 * u.u,
        'stable': False,
        'half-life': '3# s',
    },

    'Ac-206': {
        'atomic number': 89,
        'mass number': 206,
        'mass': 206.014452 * u.u,
        'stable': False,
        'half-life': 0.025 * u.s,
    },

    'Ac-207': {
        'atomic number': 89,
        'mass number': 207,
        'mass': 207.011966 * u.u,
        'stable': False,
        'half-life': 0.031 * u.s,
    },

    'Ac-208': {
        'atomic number': 89,
        'mass number': 208,
        'mass': 208.01155 * u.u,
        'stable': False,
        'half-life': 0.097 * u.s,
    },

    'Ac-209': {
        'atomic number': 89,
        'mass number': 209,
        'mass': 209.009495 * u.u,
        'stable': False,
        'half-life': 0.094 * u.s,
    },

    'Ac-210': {
        'atomic number': 89,
        'mass number': 210,
        'mass': 210.009436 * u.u,
        'stable': False,
        'half-life': 0.35 * u.s,
    },

    'Ac-211': {
        'atomic number': 89,
        'mass number': 211,
        'mass': 211.007732 * u.u,
        'stable': False,
        'half-life': 0.213 * u.s,
    },

    'Ac-212': {
        'atomic number': 89,
        'mass number': 212,
        'mass': 212.007813 * u.u,
        'stable': False,
        'half-life': 0.895 * u.s,
    },

    'Ac-213': {
        'atomic number': 89,
        'mass number': 213,
        'mass': 213.006609 * u.u,
        'stable': False,
        'half-life': 0.738 * u.s,
    },

    'Ac-214': {
        'atomic number': 89,
        'mass number': 214,
        'mass': 214.006918 * u.u,
        'stable': False,
        'half-life': 8.2 * u.s,
    },

    'Ac-215': {
        'atomic number': 89,
        'mass number': 215,
        'mass': 215.006475 * u.u,
        'stable': False,
        'half-life': 0.17 * u.s,
    },

    'Ac-216': {
        'atomic number': 89,
        'mass number': 216,
        'mass': 216.008743 * u.u,
        'stable': False,
        'half-life': 0.00044 * u.s,
    },

    'Ac-217': {
        'atomic number': 89,
        'mass number': 217,
        'mass': 217.009344 * u.u,
        'stable': False,
        'half-life': 6.9e-08 * u.s,
    },

    'Ac-218': {
        'atomic number': 89,
        'mass number': 218,
        'mass': 218.011642 * u.u,
        'stable': False,
        'half-life': 1e-06 * u.s,
    },

    'Ac-219': {
        'atomic number': 89,
        'mass number': 219,
        'mass': 219.012421 * u.u,
        'stable': False,
        'half-life': 1.18e-05 * u.s,
    },

    'Ac-220': {
        'atomic number': 89,
        'mass number': 220,
        'mass': 220.0147549 * u.u,
        'stable': False,
        'half-life': 0.02636 * u.s,
    },

    'Ac-221': {
        'atomic number': 89,
        'mass number': 221,
        'mass': 221.015592 * u.u,
        'stable': False,
        'half-life': 0.052 * u.s,
    },

    'Ac-222': {
        'atomic number': 89,
        'mass number': 222,
        'mass': 222.0178442 * u.u,
        'stable': False,
        'half-life': 5.0 * u.s,
    },

    'Ac-223': {
        'atomic number': 89,
        'mass number': 223,
        'mass': 223.0191377 * u.u,
        'stable': False,
        'half-life': 126.0 * u.s,
    },

    'Ac-224': {
        'atomic number': 89,
        'mass number': 224,
        'mass': 224.0217232 * u.u,
        'stable': False,
        'half-life': 10008.0 * u.s,
    },

    'Ac-225': {
        'atomic number': 89,
        'mass number': 225,
        'mass': 225.02323 * u.u,
        'stable': False,
        'half-life': 857088.0 * u.s,
    },

    'Ac-226': {
        'atomic number': 89,
        'mass number': 226,
        'mass': 226.0260984 * u.u,
        'stable': False,
        'half-life': 105732.0 * u.s,
    },

    'Ac-227': {
        'atomic number': 89,
        'mass number': 227,
        'mass': 227.0277523 * u.u,
        'stable': False,
        'half-life': 687057392.872 * u.s,
    },

    'Ac-228': {
        'atomic number': 89,
        'mass number': 228,
        'mass': 228.0310215 * u.u,
        'stable': False,
        'half-life': 22140.0 * u.s,
    },

    'Ac-229': {
        'atomic number': 89,
        'mass number': 229,
        'mass': 229.032956 * u.u,
        'stable': False,
        'half-life': 3762.0 * u.s,
    },

    'Ac-230': {
        'atomic number': 89,
        'mass number': 230,
        'mass': 230.036327 * u.u,
        'stable': False,
        'half-life': 122.0 * u.s,
    },

    'Ac-231': {
        'atomic number': 89,
        'mass number': 231,
        'mass': 231.038393 * u.u,
        'stable': False,
        'half-life': 450.0 * u.s,
    },

    'Ac-232': {
        'atomic number': 89,
        'mass number': 232,
        'mass': 232.042034 * u.u,
        'stable': False,
        'half-life': 118.8 * u.s,
    },

    'Ac-233': {
        'atomic number': 89,
        'mass number': 233,
        'mass': 233.044346 * u.u,
        'stable': False,
        'half-life': 145.0 * u.s,
    },

    'Ac-234': {
        'atomic number': 89,
        'mass number': 234,
        'mass': 234.048139 * u.u,
        'stable': False,
        'half-life': 45.0 * u.s,
    },

    'Ac-235': {
        'atomic number': 89,
        'mass number': 235,
        'mass': 235.05084 * u.u,
        'stable': False,
        'half-life': 62.0 * u.s,
    },

    'Ac-236': {
        'atomic number': 89,
        'mass number': 236,
        'mass': 236.054988 * u.u,
        'stable': False,
        'half-life': 270.0 * u.s,
    },

    'Ac-237': {
        'atomic number': 89,
        'mass number': 237,
        'mass': 237.05827 * u.u,
        'stable': False,
        'half-life': '4# m',
    },

    'Th-208': {
        'atomic number': 90,
        'mass number': 208,
        'mass': 208.0179 * u.u,
        'stable': False,
        'half-life': 0.0024 * u.s,
    },

    'Th-209': {
        'atomic number': 90,
        'mass number': 209,
        'mass': 209.017753 * u.u,
        'stable': False,
        'half-life': '60# ms',
    },

    'Th-210': {
        'atomic number': 90,
        'mass number': 210,
        'mass': 210.015094 * u.u,
        'stable': False,
        'half-life': 0.016 * u.s,
    },

    'Th-211': {
        'atomic number': 90,
        'mass number': 211,
        'mass': 211.014929 * u.u,
        'stable': False,
        'half-life': 0.048 * u.s,
    },

    'Th-212': {
        'atomic number': 90,
        'mass number': 212,
        'mass': 212.012988 * u.u,
        'stable': False,
        'half-life': 0.0317 * u.s,
    },

    'Th-213': {
        'atomic number': 90,
        'mass number': 213,
        'mass': 213.013009 * u.u,
        'stable': False,
        'half-life': 0.144 * u.s,
    },

    'Th-214': {
        'atomic number': 90,
        'mass number': 214,
        'mass': 214.0115 * u.u,
        'stable': False,
        'half-life': 0.087 * u.s,
    },

    'Th-215': {
        'atomic number': 90,
        'mass number': 215,
        'mass': 215.0117248 * u.u,
        'stable': False,
        'half-life': 1.2 * u.s,
    },

    'Th-216': {
        'atomic number': 90,
        'mass number': 216,
        'mass': 216.011056 * u.u,
        'stable': False,
        'half-life': 0.026 * u.s,
    },

    'Th-217': {
        'atomic number': 90,
        'mass number': 217,
        'mass': 217.013117 * u.u,
        'stable': False,
        'half-life': 0.000247 * u.s,
    },

    'Th-218': {
        'atomic number': 90,
        'mass number': 218,
        'mass': 218.013276 * u.u,
        'stable': False,
        'half-life': 1.17e-07 * u.s,
    },

    'Th-219': {
        'atomic number': 90,
        'mass number': 219,
        'mass': 219.015537 * u.u,
        'stable': False,
        'half-life': 1.021e-06 * u.s,
    },

    'Th-220': {
        'atomic number': 90,
        'mass number': 220,
        'mass': 220.015748 * u.u,
        'stable': False,
        'half-life': 9.7e-06 * u.s,
    },

    'Th-221': {
        'atomic number': 90,
        'mass number': 221,
        'mass': 221.018184 * u.u,
        'stable': False,
        'half-life': 0.00178 * u.s,
    },

    'Th-222': {
        'atomic number': 90,
        'mass number': 222,
        'mass': 222.018469 * u.u,
        'stable': False,
        'half-life': 0.00224 * u.s,
    },

    'Th-223': {
        'atomic number': 90,
        'mass number': 223,
        'mass': 223.0208119 * u.u,
        'stable': False,
        'half-life': 0.6 * u.s,
    },

    'Th-224': {
        'atomic number': 90,
        'mass number': 224,
        'mass': 224.021464 * u.u,
        'stable': False,
        'half-life': 1.04 * u.s,
    },

    'Th-225': {
        'atomic number': 90,
        'mass number': 225,
        'mass': 225.0239514 * u.u,
        'stable': False,
        'half-life': 525.0 * u.s,
    },

    'Th-226': {
        'atomic number': 90,
        'mass number': 226,
        'mass': 226.0249034 * u.u,
        'stable': False,
        'half-life': 1842.0 * u.s,
    },

    'Th-227': {
        'atomic number': 90,
        'mass number': 227,
        'mass': 227.0277042 * u.u,
        'stable': False,
        'half-life': 1615420.8 * u.s,
    },

    'Th-228': {
        'atomic number': 90,
        'mass number': 228,
        'mass': 228.0287413 * u.u,
        'stable': False,
        'half-life': 60359040.0 * u.s,
    },

    'Th-229': {
        'atomic number': 90,
        'mass number': 229,
        'mass': 229.0317627 * u.u,
        'stable': False,
        'half-life': 249930853920.0 * u.s,
    },

    'Th-230': {
        'atomic number': 90,
        'mass number': 230,
        'mass': 230.0331341 * u.u,
        'stable': False,
        'half-life': 2379392220400.0 * u.s,
    },

    'Th-231': {
        'atomic number': 90,
        'mass number': 231,
        'mass': 231.0363046 * u.u,
        'stable': False,
        'half-life': 91872.0 * u.s,
    },

    'Th-232': {
        'atomic number': 90,
        'mass number': 232,
        'mass': 232.0380558 * u.u,
        'stable': False,
        'abundance': 1,
    },

    'Th-233': {
        'atomic number': 90,
        'mass number': 233,
        'mass': 233.0415823 * u.u,
        'stable': False,
        'half-life': 1309.8 * u.s,
    },

    'Th-234': {
        'atomic number': 90,
        'mass number': 234,
        'mass': 234.0436014 * u.u,
        'stable': False,
        'half-life': 2082240.0 * u.s,
    },

    'Th-235': {
        'atomic number': 90,
        'mass number': 235,
        'mass': 235.047255 * u.u,
        'stable': False,
        'half-life': 432.0 * u.s,
    },

    'Th-236': {
        'atomic number': 90,
        'mass number': 236,
        'mass': 236.049657 * u.u,
        'stable': False,
        'half-life': 2238.0 * u.s,
    },

    'Th-237': {
        'atomic number': 90,
        'mass number': 237,
        'mass': 237.053629 * u.u,
        'stable': False,
        'half-life': 288.0 * u.s,
    },

    'Th-238': {
        'atomic number': 90,
        'mass number': 238,
        'mass': 238.0565 * u.u,
        'stable': False,
        'half-life': 564.0 * u.s,
    },

    'Th-239': {
        'atomic number': 90,
        'mass number': 239,
        'mass': 239.06077 * u.u,
        'stable': False,
        'half-life': '2# m',
    },

    'Pa-212': {
        'atomic number': 91,
        'mass number': 212,
        'mass': 212.023203 * u.u,
        'stable': False,
        'half-life': 0.0075 * u.s,
    },

    'Pa-213': {
        'atomic number': 91,
        'mass number': 213,
        'mass': 213.021109 * u.u,
        'stable': False,
        'half-life': 0.007 * u.s,
    },

    'Pa-214': {
        'atomic number': 91,
        'mass number': 214,
        'mass': 214.020918 * u.u,
        'stable': False,
        'half-life': 0.017 * u.s,
    },

    'Pa-215': {
        'atomic number': 91,
        'mass number': 215,
        'mass': 215.019183 * u.u,
        'stable': False,
        'half-life': 0.014 * u.s,
    },

    'Pa-216': {
        'atomic number': 91,
        'mass number': 216,
        'mass': 216.019109 * u.u,
        'stable': False,
        'half-life': 0.105 * u.s,
    },

    'Pa-217': {
        'atomic number': 91,
        'mass number': 217,
        'mass': 217.018325 * u.u,
        'stable': False,
        'half-life': 0.00348 * u.s,
    },

    'Pa-218': {
        'atomic number': 91,
        'mass number': 218,
        'mass': 218.020059 * u.u,
        'stable': False,
        'half-life': 0.000113 * u.s,
    },

    'Pa-219': {
        'atomic number': 91,
        'mass number': 219,
        'mass': 219.019904 * u.u,
        'stable': False,
        'half-life': 5.3e-08 * u.s,
    },

    'Pa-220': {
        'atomic number': 91,
        'mass number': 220,
        'mass': 220.021705 * u.u,
        'stable': False,
        'half-life': 7.8e-07 * u.s,
    },

    'Pa-221': {
        'atomic number': 91,
        'mass number': 221,
        'mass': 221.021875 * u.u,
        'stable': False,
        'half-life': 5.9e-06 * u.s,
    },

    'Pa-222': {
        'atomic number': 91,
        'mass number': 222,
        'mass': 222.023784 * u.u,
        'stable': False,
        'half-life': 0.0032 * u.s,
    },

    'Pa-223': {
        'atomic number': 91,
        'mass number': 223,
        'mass': 223.023963 * u.u,
        'stable': False,
        'half-life': 0.0051 * u.s,
    },

    'Pa-224': {
        'atomic number': 91,
        'mass number': 224,
        'mass': 224.0256176 * u.u,
        'stable': False,
        'half-life': 0.846 * u.s,
    },

    'Pa-225': {
        'atomic number': 91,
        'mass number': 225,
        'mass': 225.026131 * u.u,
        'stable': False,
        'half-life': 1.7 * u.s,
    },

    'Pa-226': {
        'atomic number': 91,
        'mass number': 226,
        'mass': 226.027948 * u.u,
        'stable': False,
        'half-life': 108.0 * u.s,
    },

    'Pa-227': {
        'atomic number': 91,
        'mass number': 227,
        'mass': 227.0288054 * u.u,
        'stable': False,
        'half-life': 2298.0 * u.s,
    },

    'Pa-228': {
        'atomic number': 91,
        'mass number': 228,
        'mass': 228.0310517 * u.u,
        'stable': False,
        'half-life': 79200.0 * u.s,
    },

    'Pa-229': {
        'atomic number': 91,
        'mass number': 229,
        'mass': 229.0320972 * u.u,
        'stable': False,
        'half-life': 129600.0 * u.s,
    },

    'Pa-230': {
        'atomic number': 91,
        'mass number': 230,
        'mass': 230.034541 * u.u,
        'stable': False,
        'half-life': 1503360.0 * u.s,
    },

    'Pa-231': {
        'atomic number': 91,
        'mass number': 231,
        'mass': 231.0358842 * u.u,
        'stable': False,
        'abundance': 1,
    },

    'Pa-232': {
        'atomic number': 91,
        'mass number': 232,
        'mass': 232.0385917 * u.u,
        'stable': False,
        'half-life': 114048.0 * u.s,
    },

    'Pa-233': {
        'atomic number': 91,
        'mass number': 233,
        'mass': 233.0402472 * u.u,
        'stable': False,
        'half-life': 2330640.0 * u.s,
    },

    'Pa-234': {
        'atomic number': 91,
        'mass number': 234,
        'mass': 234.0433072 * u.u,
        'stable': False,
        'half-life': 24120.0 * u.s,
    },

    'Pa-235': {
        'atomic number': 91,
        'mass number': 235,
        'mass': 235.045399 * u.u,
        'stable': False,
        'half-life': 1464.0 * u.s,
    },

    'Pa-236': {
        'atomic number': 91,
        'mass number': 236,
        'mass': 236.048668 * u.u,
        'stable': False,
        'half-life': 546.0 * u.s,
    },

    'Pa-237': {
        'atomic number': 91,
        'mass number': 237,
        'mass': 237.051023 * u.u,
        'stable': False,
        'half-life': 522.0 * u.s,
    },

    'Pa-238': {
        'atomic number': 91,
        'mass number': 238,
        'mass': 238.054637 * u.u,
        'stable': False,
        'half-life': 136.8 * u.s,
    },

    'Pa-239': {
        'atomic number': 91,
        'mass number': 239,
        'mass': 239.05726 * u.u,
        'stable': False,
        'half-life': 6480.0 * u.s,
    },

    'Pa-240': {
        'atomic number': 91,
        'mass number': 240,
        'mass': 240.06098 * u.u,
        'stable': False,
        'half-life': '2# m',
    },

    'Pa-241': {
        'atomic number': 91,
        'mass number': 241,
        'mass': 241.06408 * u.u,
        'stable': False,
        'half-life': '2# m',
    },

    'U-217': {
        'atomic number': 92,
        'mass number': 217,
        'mass': 217.02466 * u.u,
        'stable': False,
        'half-life': 0.0008 * u.s,
    },

    'U-218': {
        'atomic number': 92,
        'mass number': 218,
        'mass': 218.023523 * u.u,
        'stable': False,
        'half-life': 0.00055 * u.s,
    },

    'U-219': {
        'atomic number': 92,
        'mass number': 219,
        'mass': 219.024999 * u.u,
        'stable': False,
        'half-life': 5.5e-05 * u.s,
    },

    'U-220': {
        'atomic number': 92,
        'mass number': 220,
        'mass': 220.02462 * u.u,
        'stable': False,
        'half-life': '60# ns',
    },

    'U-221': {
        'atomic number': 92,
        'mass number': 221,
        'mass': 221.02628 * u.u,
        'stable': False,
        'half-life': 6.6e-07 * u.s,
    },

    'U-222': {
        'atomic number': 92,
        'mass number': 222,
        'mass': 222.026 * u.u,
        'stable': False,
        'half-life': 4.7e-06 * u.s,
    },

    'U-223': {
        'atomic number': 92,
        'mass number': 223,
        'mass': 223.027739 * u.u,
        'stable': False,
        'half-life': 2.1e-05 * u.s,
    },

    'U-224': {
        'atomic number': 92,
        'mass number': 224,
        'mass': 224.027605 * u.u,
        'stable': False,
        'half-life': 0.000396 * u.s,
    },

    'U-225': {
        'atomic number': 92,
        'mass number': 225,
        'mass': 225.029391 * u.u,
        'stable': False,
        'half-life': 0.061 * u.s,
    },

    'U-226': {
        'atomic number': 92,
        'mass number': 226,
        'mass': 226.029339 * u.u,
        'stable': False,
        'half-life': 0.269 * u.s,
    },

    'U-227': {
        'atomic number': 92,
        'mass number': 227,
        'mass': 227.031157 * u.u,
        'stable': False,
        'half-life': 66.0 * u.s,
    },

    'U-228': {
        'atomic number': 92,
        'mass number': 228,
        'mass': 228.031371 * u.u,
        'stable': False,
        'half-life': 546.0 * u.s,
    },

    'U-229': {
        'atomic number': 92,
        'mass number': 229,
        'mass': 229.0335063 * u.u,
        'stable': False,
        'half-life': 3468.0 * u.s,
    },

    'U-230': {
        'atomic number': 92,
        'mass number': 230,
        'mass': 230.0339401 * u.u,
        'stable': False,
        'half-life': 1747872.0 * u.s,
    },

    'U-231': {
        'atomic number': 92,
        'mass number': 231,
        'mass': 231.0362939 * u.u,
        'stable': False,
        'half-life': 362880.0 * u.s,
    },

    'U-232': {
        'atomic number': 92,
        'mass number': 232,
        'mass': 232.0371563 * u.u,
        'stable': False,
        'half-life': 2174272201.4 * u.s,
    },

    'U-233': {
        'atomic number': 92,
        'mass number': 233,
        'mass': 233.0396355 * u.u,
        'stable': False,
        'half-life': 5023862619200.0 * u.s,
    },

    'U-234': {
        'atomic number': 92,
        'mass number': 234,
        'mass': 234.0409523 * u.u,
        'stable': False,
        'abundance': 5.4e-05,
    },

    'U-235': {
        'atomic number': 92,
        'mass number': 235,
        'mass': 235.0439301 * u.u,
        'stable': False,
        'abundance': 0.007204,
    },

    'U-236': {
        'atomic number': 92,
        'mass number': 236,
        'mass': 236.0455682 * u.u,
        'stable': False,
        'half-life': 739063206920000.0 * u.s,
    },

    'U-237': {
        'atomic number': 92,
        'mass number': 237,
        'mass': 237.0487304 * u.u,
        'stable': False,
        'half-life': 583372.8 * u.s,
    },

    'U-238': {
        'atomic number': 92,
        'mass number': 238,
        'mass': 238.0507884 * u.u,
        'stable': False,
        'abundance': 0.992742,
    },

    'U-239': {
        'atomic number': 92,
        'mass number': 239,
        'mass': 239.0542935 * u.u,
        'stable': False,
        'half-life': 1407.0 * u.s,
    },

    'U-240': {
        'atomic number': 92,
        'mass number': 240,
        'mass': 240.0565934 * u.u,
        'stable': False,
        'half-life': 50760.0 * u.s,
    },

    'U-241': {
        'atomic number': 92,
        'mass number': 241,
        'mass': 241.06033 * u.u,
        'stable': False,
        'half-life': '5# m',
    },

    'U-242': {
        'atomic number': 92,
        'mass number': 242,
        'mass': 242.06293 * u.u,
        'stable': False,
        'half-life': 1008.0 * u.s,
    },

    'U-243': {
        'atomic number': 92,
        'mass number': 243,
        'mass': 243.06699 * u.u,
        'stable': False,
        'half-life': '10# m',
    },

    'Np-219': {
        'atomic number': 93,
        'mass number': 219,
        'mass': 219.03143 * u.u,
        'stable': False,
        'half-life': '<5 us',
    },

    'Np-220': {
        'atomic number': 93,
        'mass number': 220,
        'mass': 220.03254 * u.u,
        'stable': False,
        'half-life': '30# ns',
    },

    'Np-221': {
        'atomic number': 93,
        'mass number': 221,
        'mass': 221.03204 * u.u,
        'stable': False,
        'half-life': '30# ns',
    },

    'Np-222': {
        'atomic number': 93,
        'mass number': 222,
        'mass': 222.0333 * u.u,
        'stable': False,
        'half-life': '700# ns',
    },

    'Np-223': {
        'atomic number': 93,
        'mass number': 223,
        'mass': 223.03285 * u.u,
        'stable': False,
        'half-life': '1# us',
    },

    'Np-224': {
        'atomic number': 93,
        'mass number': 224,
        'mass': 224.03422 * u.u,
        'stable': False,
        'half-life': '100# us',
    },

    'Np-225': {
        'atomic number': 93,
        'mass number': 225,
        'mass': 225.033911 * u.u,
        'stable': False,
        'half-life': 0.006 * u.s,
    },

    'Np-226': {
        'atomic number': 93,
        'mass number': 226,
        'mass': 226.035188 * u.u,
        'stable': False,
        'half-life': 0.035 * u.s,
    },

    'Np-227': {
        'atomic number': 93,
        'mass number': 227,
        'mass': 227.034957 * u.u,
        'stable': False,
        'half-life': 0.51 * u.s,
    },

    'Np-228': {
        'atomic number': 93,
        'mass number': 228,
        'mass': 228.036067 * u.u,
        'stable': False,
        'half-life': 61.4 * u.s,
    },

    'Np-229': {
        'atomic number': 93,
        'mass number': 229,
        'mass': 229.036264 * u.u,
        'stable': False,
        'half-life': 240.0 * u.s,
    },

    'Np-230': {
        'atomic number': 93,
        'mass number': 230,
        'mass': 230.037828 * u.u,
        'stable': False,
        'half-life': 276.0 * u.s,
    },

    'Np-231': {
        'atomic number': 93,
        'mass number': 231,
        'mass': 231.038245 * u.u,
        'stable': False,
        'half-life': 2928.0 * u.s,
    },

    'Np-232': {
        'atomic number': 93,
        'mass number': 232,
        'mass': 232.04011 * u.u,
        'stable': False,
        'half-life': 882.0 * u.s,
    },

    'Np-233': {
        'atomic number': 93,
        'mass number': 233,
        'mass': 233.040741 * u.u,
        'stable': False,
        'half-life': 2172.0 * u.s,
    },

    'Np-234': {
        'atomic number': 93,
        'mass number': 234,
        'mass': 234.0428953 * u.u,
        'stable': False,
        'half-life': 380160.0 * u.s,
    },

    'Np-235': {
        'atomic number': 93,
        'mass number': 235,
        'mass': 235.0440635 * u.u,
        'stable': False,
        'half-life': 34223040.0 * u.s,
    },

    'Np-236': {
        'atomic number': 93,
        'mass number': 236,
        'mass': 236.04657 * u.u,
        'stable': False,
        'half-life': 4828209678000.0 * u.s,
    },

    'Np-237': {
        'atomic number': 93,
        'mass number': 237,
        'mass': 237.0481736 * u.u,
        'stable': False,
        'half-life': 67658049344000.0 * u.s,
    },

    'Np-238': {
        'atomic number': 93,
        'mass number': 238,
        'mass': 238.0509466 * u.u,
        'stable': False,
        'half-life': 181353.6 * u.s,
    },

    'Np-239': {
        'atomic number': 93,
        'mass number': 239,
        'mass': 239.0529392 * u.u,
        'stable': False,
        'half-life': 203558.4 * u.s,
    },

    'Np-240': {
        'atomic number': 93,
        'mass number': 240,
        'mass': 240.056165 * u.u,
        'stable': False,
        'half-life': 3714.0 * u.s,
    },

    'Np-241': {
        'atomic number': 93,
        'mass number': 241,
        'mass': 241.058253 * u.u,
        'stable': False,
        'half-life': 834.0 * u.s,
    },

    'Np-242': {
        'atomic number': 93,
        'mass number': 242,
        'mass': 242.06164 * u.u,
        'stable': False,
        'half-life': 132.0 * u.s,
    },

    'Np-243': {
        'atomic number': 93,
        'mass number': 243,
        'mass': 243.06428 * u.u,
        'stable': False,
        'half-life': 111.0 * u.s,
    },

    'Np-244': {
        'atomic number': 93,
        'mass number': 244,
        'mass': 244.06785 * u.u,
        'stable': False,
        'half-life': 137.4 * u.s,
    },

    'Np-245': {
        'atomic number': 93,
        'mass number': 245,
        'mass': 245.0708 * u.u,
        'stable': False,
        'half-life': '2# m',
    },

    'Pu-228': {
        'atomic number': 94,
        'mass number': 228,
        'mass': 228.038732 * u.u,
        'stable': False,
        'half-life': 2.1 * u.s,
    },

    'Pu-229': {
        'atomic number': 94,
        'mass number': 229,
        'mass': 229.040144 * u.u,
        'stable': False,
        'half-life': 91.0 * u.s,
    },

    'Pu-230': {
        'atomic number': 94,
        'mass number': 230,
        'mass': 230.03965 * u.u,
        'stable': False,
        'half-life': 102.0 * u.s,
    },

    'Pu-231': {
        'atomic number': 94,
        'mass number': 231,
        'mass': 231.041102 * u.u,
        'stable': False,
        'half-life': 516.0 * u.s,
    },

    'Pu-232': {
        'atomic number': 94,
        'mass number': 232,
        'mass': 232.041185 * u.u,
        'stable': False,
        'half-life': 2022.0 * u.s,
    },

    'Pu-233': {
        'atomic number': 94,
        'mass number': 233,
        'mass': 233.042998 * u.u,
        'stable': False,
        'half-life': 1254.0 * u.s,
    },

    'Pu-234': {
        'atomic number': 94,
        'mass number': 234,
        'mass': 234.0433174 * u.u,
        'stable': False,
        'half-life': 31680.0 * u.s,
    },

    'Pu-235': {
        'atomic number': 94,
        'mass number': 235,
        'mass': 235.045286 * u.u,
        'stable': False,
        'half-life': 1518.0 * u.s,
    },

    'Pu-236': {
        'atomic number': 94,
        'mass number': 236,
        'mass': 236.0460581 * u.u,
        'stable': False,
        'half-life': 90189694.508 * u.s,
    },

    'Pu-237': {
        'atomic number': 94,
        'mass number': 237,
        'mass': 237.0484098 * u.u,
        'stable': False,
        'half-life': 3943296.0 * u.s,
    },

    'Pu-238': {
        'atomic number': 94,
        'mass number': 238,
        'mass': 238.0495601 * u.u,
        'stable': False,
        'half-life': 2767542410.2 * u.s,
    },

    'Pu-239': {
        'atomic number': 94,
        'mass number': 239,
        'mass': 239.0521636 * u.u,
        'stable': False,
        'half-life': 760837485860.0 * u.s,
    },

    'Pu-240': {
        'atomic number': 94,
        'mass number': 240,
        'mass': 240.0538138 * u.u,
        'stable': False,
        'half-life': 207044991486.0 * u.s,
    },

    'Pu-241': {
        'atomic number': 94,
        'mass number': 241,
        'mass': 241.0568517 * u.u,
        'stable': False,
        'half-life': 452179192.654 * u.s,
    },

    'Pu-242': {
        'atomic number': 94,
        'mass number': 242,
        'mass': 242.0587428 * u.u,
        'stable': False,
        'half-life': 11833847250000.0 * u.s,
    },

    'Pu-243': {
        'atomic number': 94,
        'mass number': 243,
        'mass': 243.0620036 * u.u,
        'stable': False,
        'half-life': 17841.6 * u.s,
    },

    'Pu-244': {
        'atomic number': 94,
        'mass number': 244,
        'mass': 244.0642053 * u.u,
        'stable': False,
        'half-life': 2524554080000000.0 * u.s,
    },

    'Pu-245': {
        'atomic number': 94,
        'mass number': 245,
        'mass': 245.067826 * u.u,
        'stable': False,
        'half-life': 37800.0 * u.s,
    },

    'Pu-246': {
        'atomic number': 94,
        'mass number': 246,
        'mass': 246.070205 * u.u,
        'stable': False,
        'half-life': 936576.0 * u.s,
    },

    'Pu-247': {
        'atomic number': 94,
        'mass number': 247,
        'mass': 247.07419 * u.u,
        'stable': False,
        'half-life': 196128.0 * u.s,
    },

    'Am-230': {
        'atomic number': 95,
        'mass number': 230,
        'mass': 230.04609 * u.u,
        'stable': False,
        'half-life': 40.0 * u.s,
    },

    'Am-231': {
        'atomic number': 95,
        'mass number': 231,
        'mass': 231.04556 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Am-232': {
        'atomic number': 95,
        'mass number': 232,
        'mass': 232.04645 * u.u,
        'stable': False,
        'half-life': 78.6 * u.s,
    },

    'Am-233': {
        'atomic number': 95,
        'mass number': 233,
        'mass': 233.04644 * u.u,
        'stable': False,
        'half-life': 192.0 * u.s,
    },

    'Am-234': {
        'atomic number': 95,
        'mass number': 234,
        'mass': 234.04773 * u.u,
        'stable': False,
        'half-life': 139.2 * u.s,
    },

    'Am-235': {
        'atomic number': 95,
        'mass number': 235,
        'mass': 235.047908 * u.u,
        'stable': False,
        'half-life': 618.0 * u.s,
    },

    'Am-236': {
        'atomic number': 95,
        'mass number': 236,
        'mass': 236.04943 * u.u,
        'stable': False,
        'half-life': 216.0 * u.s,
    },

    'Am-237': {
        'atomic number': 95,
        'mass number': 237,
        'mass': 237.049996 * u.u,
        'stable': False,
        'half-life': 4416.0 * u.s,
    },

    'Am-238': {
        'atomic number': 95,
        'mass number': 238,
        'mass': 238.051985 * u.u,
        'stable': False,
        'half-life': 5880.0 * u.s,
    },

    'Am-239': {
        'atomic number': 95,
        'mass number': 239,
        'mass': 239.0530247 * u.u,
        'stable': False,
        'half-life': 42840.0 * u.s,
    },

    'Am-240': {
        'atomic number': 95,
        'mass number': 240,
        'mass': 240.0553 * u.u,
        'stable': False,
        'half-life': 182880.0 * u.s,
    },

    'Am-241': {
        'atomic number': 95,
        'mass number': 241,
        'mass': 241.0568293 * u.u,
        'stable': False,
        'half-life': 13651526187.6 * u.s,
    },

    'Am-242': {
        'atomic number': 95,
        'mass number': 242,
        'mass': 242.0595494 * u.u,
        'stable': False,
        'half-life': 57672.0 * u.s,
    },

    'Am-243': {
        'atomic number': 95,
        'mass number': 243,
        'mass': 243.0613813 * u.u,
        'stable': False,
        'half-life': 232385203064.0 * u.s,
    },

    'Am-244': {
        'atomic number': 95,
        'mass number': 244,
        'mass': 244.0642851 * u.u,
        'stable': False,
        'half-life': 36360.0 * u.s,
    },

    'Am-245': {
        'atomic number': 95,
        'mass number': 245,
        'mass': 245.0664548 * u.u,
        'stable': False,
        'half-life': 7380.0 * u.s,
    },

    'Am-246': {
        'atomic number': 95,
        'mass number': 246,
        'mass': 246.069775 * u.u,
        'stable': False,
        'half-life': 2340.0 * u.s,
    },

    'Am-247': {
        'atomic number': 95,
        'mass number': 247,
        'mass': 247.07209 * u.u,
        'stable': False,
        'half-life': 1380.0 * u.s,
    },

    'Am-248': {
        'atomic number': 95,
        'mass number': 248,
        'mass': 248.07575 * u.u,
        'stable': False,
        'half-life': '3# m',
    },

    'Am-249': {
        'atomic number': 95,
        'mass number': 249,
        'mass': 249.07848 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Cm-232': {
        'atomic number': 96,
        'mass number': 232,
        'mass': 232.04982 * u.u,
        'stable': False,
        'half-life': '10# s',
    },

    'Cm-233': {
        'atomic number': 96,
        'mass number': 233,
        'mass': 233.05077 * u.u,
        'stable': False,
        'half-life': 27.0 * u.s,
    },

    'Cm-234': {
        'atomic number': 96,
        'mass number': 234,
        'mass': 234.05016 * u.u,
        'stable': False,
        'half-life': 52.0 * u.s,
    },

    'Cm-235': {
        'atomic number': 96,
        'mass number': 235,
        'mass': 235.05154 * u.u,
        'stable': False,
        'half-life': '5# m',
    },

    'Cm-236': {
        'atomic number': 96,
        'mass number': 236,
        'mass': 236.051374 * u.u,
        'stable': False,
        'half-life': 408.0 * u.s,
    },

    'Cm-237': {
        'atomic number': 96,
        'mass number': 237,
        'mass': 237.052869 * u.u,
        'stable': False,
        'half-life': '20# m',
    },

    'Cm-238': {
        'atomic number': 96,
        'mass number': 238,
        'mass': 238.053081 * u.u,
        'stable': False,
        'half-life': 7920.0 * u.s,
    },

    'Cm-239': {
        'atomic number': 96,
        'mass number': 239,
        'mass': 239.05491 * u.u,
        'stable': False,
        'half-life': 9000.0 * u.s,
    },

    'Cm-240': {
        'atomic number': 96,
        'mass number': 240,
        'mass': 240.0555297 * u.u,
        'stable': False,
        'half-life': 2332800.0 * u.s,
    },

    'Cm-241': {
        'atomic number': 96,
        'mass number': 241,
        'mass': 241.0576532 * u.u,
        'stable': False,
        'half-life': 2833920.0 * u.s,
    },

    'Cm-242': {
        'atomic number': 96,
        'mass number': 242,
        'mass': 242.058836 * u.u,
        'stable': False,
        'half-life': 14065920.0 * u.s,
    },

    'Cm-243': {
        'atomic number': 96,
        'mass number': 243,
        'mass': 243.0613893 * u.u,
        'stable': False,
        'half-life': 918306546.6 * u.s,
    },

    'Cm-244': {
        'atomic number': 96,
        'mass number': 244,
        'mass': 244.0627528 * u.u,
        'stable': False,
        'half-life': 571180360.6 * u.s,
    },

    'Cm-245': {
        'atomic number': 96,
        'mass number': 245,
        'mass': 245.0654915 * u.u,
        'stable': False,
        'half-life': 260344639500.0 * u.s,
    },

    'Cm-246': {
        'atomic number': 96,
        'mass number': 246,
        'mass': 246.0672238 * u.u,
        'stable': False,
        'half-life': 148506893756.0 * u.s,
    },

    'Cm-247': {
        'atomic number': 96,
        'mass number': 247,
        'mass': 247.0703541 * u.u,
        'stable': False,
        'half-life': 492288045600000.0 * u.s,
    },

    'Cm-248': {
        'atomic number': 96,
        'mass number': 248,
        'mass': 248.0723499 * u.u,
        'stable': False,
        'half-life': 10981810248000.0 * u.s,
    },

    'Cm-249': {
        'atomic number': 96,
        'mass number': 249,
        'mass': 249.0759548 * u.u,
        'stable': False,
        'half-life': 3849.0 * u.s,
    },

    'Cm-250': {
        'atomic number': 96,
        'mass number': 250,
        'mass': 250.078358 * u.u,
        'stable': False,
        'half-life': '8300# y',
    },

    'Cm-251': {
        'atomic number': 96,
        'mass number': 251,
        'mass': 251.082286 * u.u,
        'stable': False,
        'half-life': 1008.0 * u.s,
    },

    'Cm-252': {
        'atomic number': 96,
        'mass number': 252,
        'mass': 252.08487 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Bk-234': {
        'atomic number': 97,
        'mass number': 234,
        'mass': 234.05727 * u.u,
        'stable': False,
        'half-life': 20.0 * u.s,
    },

    'Bk-235': {
        'atomic number': 97,
        'mass number': 235,
        'mass': 235.05658 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Bk-236': {
        'atomic number': 97,
        'mass number': 236,
        'mass': 236.05748 * u.u,
        'stable': False,
        'half-life': '2# m',
    },

    'Bk-237': {
        'atomic number': 97,
        'mass number': 237,
        'mass': 237.0571 * u.u,
        'stable': False,
        'half-life': '2# m',
    },

    'Bk-238': {
        'atomic number': 97,
        'mass number': 238,
        'mass': 238.0582 * u.u,
        'stable': False,
        'half-life': 144.0 * u.s,
    },

    'Bk-239': {
        'atomic number': 97,
        'mass number': 239,
        'mass': 239.05824 * u.u,
        'stable': False,
        'half-life': '4# m',
    },

    'Bk-240': {
        'atomic number': 97,
        'mass number': 240,
        'mass': 240.05976 * u.u,
        'stable': False,
        'half-life': 288.0 * u.s,
    },

    'Bk-241': {
        'atomic number': 97,
        'mass number': 241,
        'mass': 241.06016 * u.u,
        'stable': False,
        'half-life': 276.0 * u.s,
    },

    'Bk-242': {
        'atomic number': 97,
        'mass number': 242,
        'mass': 242.06198 * u.u,
        'stable': False,
        'half-life': 420.0 * u.s,
    },

    'Bk-243': {
        'atomic number': 97,
        'mass number': 243,
        'mass': 243.0630078 * u.u,
        'stable': False,
        'half-life': 16560.0 * u.s,
    },

    'Bk-244': {
        'atomic number': 97,
        'mass number': 244,
        'mass': 244.065181 * u.u,
        'stable': False,
        'half-life': 18072.0 * u.s,
    },

    'Bk-245': {
        'atomic number': 97,
        'mass number': 245,
        'mass': 245.0663618 * u.u,
        'stable': False,
        'half-life': 427680.0 * u.s,
    },

    'Bk-246': {
        'atomic number': 97,
        'mass number': 246,
        'mass': 246.068673 * u.u,
        'stable': False,
        'half-life': 155520.0 * u.s,
    },

    'Bk-247': {
        'atomic number': 97,
        'mass number': 247,
        'mass': 247.0703073 * u.u,
        'stable': False,
        'half-life': 43548557880.0 * u.s,
    },

    'Bk-248': {
        'atomic number': 97,
        'mass number': 248,
        'mass': 248.073088 * u.u,
        'stable': False,
        'half-life': '>9 y',
    },

    'Bk-249': {
        'atomic number': 97,
        'mass number': 249,
        'mass': 249.0749877 * u.u,
        'stable': False,
        'half-life': 28270080.0 * u.s,
    },

    'Bk-250': {
        'atomic number': 97,
        'mass number': 250,
        'mass': 250.0783167 * u.u,
        'stable': False,
        'half-life': 11563.2 * u.s,
    },

    'Bk-251': {
        'atomic number': 97,
        'mass number': 251,
        'mass': 251.080762 * u.u,
        'stable': False,
        'half-life': 3336.0 * u.s,
    },

    'Bk-252': {
        'atomic number': 97,
        'mass number': 252,
        'mass': 252.08431 * u.u,
        'stable': False,
        'half-life': 108.0 * u.s,
    },

    'Bk-253': {
        'atomic number': 97,
        'mass number': 253,
        'mass': 253.08688 * u.u,
        'stable': False,
        'half-life': '10# m',
    },

    'Bk-254': {
        'atomic number': 97,
        'mass number': 254,
        'mass': 254.0906 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Cf-237': {
        'atomic number': 98,
        'mass number': 237,
        'mass': 237.062198 * u.u,
        'stable': False,
        'half-life': 0.8 * u.s,
    },

    'Cf-238': {
        'atomic number': 98,
        'mass number': 238,
        'mass': 238.06149 * u.u,
        'stable': False,
        'half-life': 0.0211 * u.s,
    },

    'Cf-239': {
        'atomic number': 98,
        'mass number': 239,
        'mass': 239.06253 * u.u,
        'stable': False,
        'half-life': 60.0 * u.s,
    },

    'Cf-240': {
        'atomic number': 98,
        'mass number': 240,
        'mass': 240.062256 * u.u,
        'stable': False,
        'half-life': 40.3 * u.s,
    },

    'Cf-241': {
        'atomic number': 98,
        'mass number': 241,
        'mass': 241.06369 * u.u,
        'stable': False,
        'half-life': 141.0 * u.s,
    },

    'Cf-242': {
        'atomic number': 98,
        'mass number': 242,
        'mass': 242.063754 * u.u,
        'stable': False,
        'half-life': 209.4 * u.s,
    },

    'Cf-243': {
        'atomic number': 98,
        'mass number': 243,
        'mass': 243.06548 * u.u,
        'stable': False,
        'half-life': 642.0 * u.s,
    },

    'Cf-244': {
        'atomic number': 98,
        'mass number': 244,
        'mass': 244.0660008 * u.u,
        'stable': False,
        'half-life': 1164.0 * u.s,
    },

    'Cf-245': {
        'atomic number': 98,
        'mass number': 245,
        'mass': 245.0680487 * u.u,
        'stable': False,
        'half-life': 2700.0 * u.s,
    },

    'Cf-246': {
        'atomic number': 98,
        'mass number': 246,
        'mass': 246.0688055 * u.u,
        'stable': False,
        'half-life': 128520.0 * u.s,
    },

    'Cf-247': {
        'atomic number': 98,
        'mass number': 247,
        'mass': 247.070965 * u.u,
        'stable': False,
        'half-life': 11196.0 * u.s,
    },

    'Cf-248': {
        'atomic number': 98,
        'mass number': 248,
        'mass': 248.0721851 * u.u,
        'stable': False,
        'half-life': 28814400.0 * u.s,
    },

    'Cf-249': {
        'atomic number': 98,
        'mass number': 249,
        'mass': 249.0748539 * u.u,
        'stable': False,
        'half-life': 11076481026.0 * u.s,
    },

    'Cf-250': {
        'atomic number': 98,
        'mass number': 250,
        'mass': 250.0764062 * u.u,
        'stable': False,
        'half-life': 412764592.08 * u.s,
    },

    'Cf-251': {
        'atomic number': 98,
        'mass number': 251,
        'mass': 251.0795886 * u.u,
        'stable': False,
        'half-life': 28401233400.0 * u.s,
    },

    'Cf-252': {
        'atomic number': 98,
        'mass number': 252,
        'mass': 252.0816272 * u.u,
        'stable': False,
        'half-life': 83468069.27 * u.s,
    },

    'Cf-253': {
        'atomic number': 98,
        'mass number': 253,
        'mass': 253.0851345 * u.u,
        'stable': False,
        'half-life': 1538784.0 * u.s,
    },

    'Cf-254': {
        'atomic number': 98,
        'mass number': 254,
        'mass': 254.087324 * u.u,
        'stable': False,
        'half-life': 5227200.0 * u.s,
    },

    'Cf-255': {
        'atomic number': 98,
        'mass number': 255,
        'mass': 255.09105 * u.u,
        'stable': False,
        'half-life': 5100.0 * u.s,
    },

    'Cf-256': {
        'atomic number': 98,
        'mass number': 256,
        'mass': 256.09344 * u.u,
        'stable': False,
        'half-life': 738.0 * u.s,
    },

    'Es-239': {
        'atomic number': 99,
        'mass number': 239,
        'mass': 239.06823 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Es-240': {
        'atomic number': 99,
        'mass number': 240,
        'mass': 240.06892 * u.u,
        'stable': False,
        'half-life': '1# s',
    },

    'Es-241': {
        'atomic number': 99,
        'mass number': 241,
        'mass': 241.06856 * u.u,
        'stable': False,
        'half-life': 10.0 * u.s,
    },

    'Es-242': {
        'atomic number': 99,
        'mass number': 242,
        'mass': 242.06957 * u.u,
        'stable': False,
        'half-life': 17.8 * u.s,
    },

    'Es-243': {
        'atomic number': 99,
        'mass number': 243,
        'mass': 243.06951 * u.u,
        'stable': False,
        'half-life': 21.6 * u.s,
    },

    'Es-244': {
        'atomic number': 99,
        'mass number': 244,
        'mass': 244.07088 * u.u,
        'stable': False,
        'half-life': 37.0 * u.s,
    },

    'Es-245': {
        'atomic number': 99,
        'mass number': 245,
        'mass': 245.07125 * u.u,
        'stable': False,
        'half-life': 66.0 * u.s,
    },

    'Es-246': {
        'atomic number': 99,
        'mass number': 246,
        'mass': 246.0729 * u.u,
        'stable': False,
        'half-life': 450.0 * u.s,
    },

    'Es-247': {
        'atomic number': 99,
        'mass number': 247,
        'mass': 247.073622 * u.u,
        'stable': False,
        'half-life': 273.0 * u.s,
    },

    'Es-248': {
        'atomic number': 99,
        'mass number': 248,
        'mass': 248.075471 * u.u,
        'stable': False,
        'half-life': 1440.0 * u.s,
    },

    'Es-249': {
        'atomic number': 99,
        'mass number': 249,
        'mass': 249.076411 * u.u,
        'stable': False,
        'half-life': 6132.0 * u.s,
    },

    'Es-250': {
        'atomic number': 99,
        'mass number': 250,
        'mass': 250.07861 * u.u,
        'stable': False,
        'half-life': 30960.0 * u.s,
    },

    'Es-251': {
        'atomic number': 99,
        'mass number': 251,
        'mass': 251.0799936 * u.u,
        'stable': False,
        'half-life': 118800.0 * u.s,
    },

    'Es-252': {
        'atomic number': 99,
        'mass number': 252,
        'mass': 252.08298 * u.u,
        'stable': False,
        'half-life': 40754880.0 * u.s,
    },

    'Es-253': {
        'atomic number': 99,
        'mass number': 253,
        'mass': 253.0848257 * u.u,
        'stable': False,
        'half-life': 1768608.0 * u.s,
    },

    'Es-254': {
        'atomic number': 99,
        'mass number': 254,
        'mass': 254.0880222 * u.u,
        'stable': False,
        'half-life': 23820480.0 * u.s,
    },

    'Es-255': {
        'atomic number': 99,
        'mass number': 255,
        'mass': 255.090275 * u.u,
        'stable': False,
        'half-life': 3438720.0 * u.s,
    },

    'Es-256': {
        'atomic number': 99,
        'mass number': 256,
        'mass': 256.0936 * u.u,
        'stable': False,
        'half-life': 1524.0 * u.s,
    },

    'Es-257': {
        'atomic number': 99,
        'mass number': 257,
        'mass': 257.09598 * u.u,
        'stable': False,
        'half-life': 665280.0 * u.s,
    },

    'Es-258': {
        'atomic number': 99,
        'mass number': 258,
        'mass': 258.09952 * u.u,
        'stable': False,
        'half-life': '3# m',
    },

    'Fm-241': {
        'atomic number': 100,
        'mass number': 241,
        'mass': 241.07421 * u.u,
        'stable': False,
        'half-life': 0.00073 * u.s,
    },

    'Fm-242': {
        'atomic number': 100,
        'mass number': 242,
        'mass': 242.07343 * u.u,
        'stable': False,
        'half-life': 0.0008 * u.s,
    },

    'Fm-243': {
        'atomic number': 100,
        'mass number': 243,
        'mass': 243.07446 * u.u,
        'stable': False,
        'half-life': 0.231 * u.s,
    },

    'Fm-244': {
        'atomic number': 100,
        'mass number': 244,
        'mass': 244.07404 * u.u,
        'stable': False,
        'half-life': 0.00312 * u.s,
    },

    'Fm-245': {
        'atomic number': 100,
        'mass number': 245,
        'mass': 245.07535 * u.u,
        'stable': False,
        'half-life': 4.2 * u.s,
    },

    'Fm-246': {
        'atomic number': 100,
        'mass number': 246,
        'mass': 246.07535 * u.u,
        'stable': False,
        'half-life': 1.54 * u.s,
    },

    'Fm-247': {
        'atomic number': 100,
        'mass number': 247,
        'mass': 247.07694 * u.u,
        'stable': False,
        'half-life': 31.0 * u.s,
    },

    'Fm-248': {
        'atomic number': 100,
        'mass number': 248,
        'mass': 248.0771865 * u.u,
        'stable': False,
        'half-life': 34.5 * u.s,
    },

    'Fm-249': {
        'atomic number': 100,
        'mass number': 249,
        'mass': 249.0789275 * u.u,
        'stable': False,
        'half-life': 96.0 * u.s,
    },

    'Fm-250': {
        'atomic number': 100,
        'mass number': 250,
        'mass': 250.079521 * u.u,
        'stable': False,
        'half-life': 1824.0 * u.s,
    },

    'Fm-251': {
        'atomic number': 100,
        'mass number': 251,
        'mass': 251.08154 * u.u,
        'stable': False,
        'half-life': 19080.0 * u.s,
    },

    'Fm-252': {
        'atomic number': 100,
        'mass number': 252,
        'mass': 252.0824671 * u.u,
        'stable': False,
        'half-life': 91404.0 * u.s,
    },

    'Fm-253': {
        'atomic number': 100,
        'mass number': 253,
        'mass': 253.0851846 * u.u,
        'stable': False,
        'half-life': 259200.0 * u.s,
    },

    'Fm-254': {
        'atomic number': 100,
        'mass number': 254,
        'mass': 254.0868544 * u.u,
        'stable': False,
        'half-life': 11664.0 * u.s,
    },

    'Fm-255': {
        'atomic number': 100,
        'mass number': 255,
        'mass': 255.089964 * u.u,
        'stable': False,
        'half-life': 72252.0 * u.s,
    },

    'Fm-256': {
        'atomic number': 100,
        'mass number': 256,
        'mass': 256.0917745 * u.u,
        'stable': False,
        'half-life': 9456.0 * u.s,
    },

    'Fm-257': {
        'atomic number': 100,
        'mass number': 257,
        'mass': 257.0951061 * u.u,
        'stable': False,
        'half-life': 8683200.0 * u.s,
    },

    'Fm-258': {
        'atomic number': 100,
        'mass number': 258,
        'mass': 258.09708 * u.u,
        'stable': False,
        'half-life': 0.00037 * u.s,
    },

    'Fm-259': {
        'atomic number': 100,
        'mass number': 259,
        'mass': 259.1006 * u.u,
        'stable': False,
        'half-life': 1.5 * u.s,
    },

    'Fm-260': {
        'atomic number': 100,
        'mass number': 260,
        'mass': 260.10281 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Md-245': {
        'atomic number': 101,
        'mass number': 245,
        'mass': 245.08081 * u.u,
        'stable': False,
        'half-life': 0.4 * u.s,
    },

    'Md-246': {
        'atomic number': 101,
        'mass number': 246,
        'mass': 246.08171 * u.u,
        'stable': False,
        'half-life': 0.92 * u.s,
    },

    'Md-247': {
        'atomic number': 101,
        'mass number': 247,
        'mass': 247.08152 * u.u,
        'stable': False,
        'half-life': 1.2 * u.s,
    },

    'Md-248': {
        'atomic number': 101,
        'mass number': 248,
        'mass': 248.08282 * u.u,
        'stable': False,
        'half-life': 7.0 * u.s,
    },

    'Md-249': {
        'atomic number': 101,
        'mass number': 249,
        'mass': 249.08291 * u.u,
        'stable': False,
        'half-life': 23.4 * u.s,
    },

    'Md-250': {
        'atomic number': 101,
        'mass number': 250,
        'mass': 250.08441 * u.u,
        'stable': False,
        'half-life': 52.0 * u.s,
    },

    'Md-251': {
        'atomic number': 101,
        'mass number': 251,
        'mass': 251.084774 * u.u,
        'stable': False,
        'half-life': 252.6 * u.s,
    },

    'Md-252': {
        'atomic number': 101,
        'mass number': 252,
        'mass': 252.08643 * u.u,
        'stable': False,
        'half-life': 138.0 * u.s,
    },

    'Md-253': {
        'atomic number': 101,
        'mass number': 253,
        'mass': 253.087144 * u.u,
        'stable': False,
        'half-life': 720.0 * u.s,
    },

    'Md-254': {
        'atomic number': 101,
        'mass number': 254,
        'mass': 254.08959 * u.u,
        'stable': False,
        'half-life': 600.0 * u.s,
    },

    'Md-255': {
        'atomic number': 101,
        'mass number': 255,
        'mass': 255.0910841 * u.u,
        'stable': False,
        'half-life': 1620.0 * u.s,
    },

    'Md-256': {
        'atomic number': 101,
        'mass number': 256,
        'mass': 256.09389 * u.u,
        'stable': False,
        'half-life': '30# m',
    },

    'Md-257': {
        'atomic number': 101,
        'mass number': 257,
        'mass': 257.0955424 * u.u,
        'stable': False,
        'half-life': 19872.0 * u.s,
    },

    'Md-258': {
        'atomic number': 101,
        'mass number': 258,
        'mass': 258.0984315 * u.u,
        'stable': False,
        'half-life': 4449600.0 * u.s,
    },

    'Md-259': {
        'atomic number': 101,
        'mass number': 259,
        'mass': 259.10051 * u.u,
        'stable': False,
        'half-life': 5760.0 * u.s,
    },

    'Md-260': {
        'atomic number': 101,
        'mass number': 260,
        'mass': 260.10365 * u.u,
        'stable': False,
        'half-life': 2401920.0 * u.s,
    },

    'Md-261': {
        'atomic number': 101,
        'mass number': 261,
        'mass': 261.10583 * u.u,
        'stable': False,
        'half-life': '40# m',
    },

    'Md-262': {
        'atomic number': 101,
        'mass number': 262,
        'mass': 262.1091 * u.u,
        'stable': False,
        'half-life': '3# m',
    },

    'No-248': {
        'atomic number': 102,
        'mass number': 248,
        'mass': 248.08655 * u.u,
        'stable': False,
    },

    'No-249': {
        'atomic number': 102,
        'mass number': 249,
        'mass': 249.0878 * u.u,
        'stable': False,
        'half-life': 5.7e-05 * u.s,
    },

    'No-250': {
        'atomic number': 102,
        'mass number': 250,
        'mass': 250.08756 * u.u,
        'stable': False,
        'half-life': 5e-06 * u.s,
    },

    'No-251': {
        'atomic number': 102,
        'mass number': 251,
        'mass': 251.08894 * u.u,
        'stable': False,
        'half-life': 0.8 * u.s,
    },

    'No-252': {
        'atomic number': 102,
        'mass number': 252,
        'mass': 252.088967 * u.u,
        'stable': False,
        'half-life': 2.45 * u.s,
    },

    'No-253': {
        'atomic number': 102,
        'mass number': 253,
        'mass': 253.0905641 * u.u,
        'stable': False,
        'half-life': 93.6 * u.s,
    },

    'No-254': {
        'atomic number': 102,
        'mass number': 254,
        'mass': 254.090956 * u.u,
        'stable': False,
        'half-life': 51.2 * u.s,
    },

    'No-255': {
        'atomic number': 102,
        'mass number': 255,
        'mass': 255.093191 * u.u,
        'stable': False,
        'half-life': 211.2 * u.s,
    },

    'No-256': {
        'atomic number': 102,
        'mass number': 256,
        'mass': 256.0942829 * u.u,
        'stable': False,
        'half-life': 2.91 * u.s,
    },

    'No-257': {
        'atomic number': 102,
        'mass number': 257,
        'mass': 257.0968878 * u.u,
        'stable': False,
        'half-life': 24.5 * u.s,
    },

    'No-258': {
        'atomic number': 102,
        'mass number': 258,
        'mass': 258.09821 * u.u,
        'stable': False,
        'half-life': 0.0012 * u.s,
    },

    'No-259': {
        'atomic number': 102,
        'mass number': 259,
        'mass': 259.10103 * u.u,
        'stable': False,
        'half-life': 3480.0 * u.s,
    },

    'No-260': {
        'atomic number': 102,
        'mass number': 260,
        'mass': 260.10264 * u.u,
        'stable': False,
        'half-life': 0.106 * u.s,
    },

    'No-261': {
        'atomic number': 102,
        'mass number': 261,
        'mass': 261.1057 * u.u,
        'stable': False,
        'half-life': '3# h',
    },

    'No-262': {
        'atomic number': 102,
        'mass number': 262,
        'mass': 262.10746 * u.u,
        'stable': False,
        'half-life': '~5 ms',
    },

    'No-263': {
        'atomic number': 102,
        'mass number': 263,
        'mass': 263.11071 * u.u,
        'stable': False,
        'half-life': '20# m',
    },

    'No-264': {
        'atomic number': 102,
        'mass number': 264,
        'mass': 264.11273 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Lr-251': {
        'atomic number': 103,
        'mass number': 251,
        'mass': 251.09418 * u.u,
        'stable': False,
        'half-life': '150# us',
    },

    'Lr-252': {
        'atomic number': 103,
        'mass number': 252,
        'mass': 252.09526 * u.u,
        'stable': False,
        'half-life': 0.369 * u.s,
    },

    'Lr-253': {
        'atomic number': 103,
        'mass number': 253,
        'mass': 253.09509 * u.u,
        'stable': False,
        'half-life': 0.632 * u.s,
    },

    'Lr-254': {
        'atomic number': 103,
        'mass number': 254,
        'mass': 254.09648 * u.u,
        'stable': False,
        'half-life': 17.1 * u.s,
    },

    'Lr-255': {
        'atomic number': 103,
        'mass number': 255,
        'mass': 255.096562 * u.u,
        'stable': False,
        'half-life': 31.1 * u.s,
    },

    'Lr-256': {
        'atomic number': 103,
        'mass number': 256,
        'mass': 256.098494 * u.u,
        'stable': False,
        'half-life': 27.0 * u.s,
    },

    'Lr-257': {
        'atomic number': 103,
        'mass number': 257,
        'mass': 257.099418 * u.u,
        'stable': False,
        'half-life': 6.0 * u.s,
    },

    'Lr-258': {
        'atomic number': 103,
        'mass number': 258,
        'mass': 258.10176 * u.u,
        'stable': False,
        'half-life': 3.6 * u.s,
    },

    'Lr-259': {
        'atomic number': 103,
        'mass number': 259,
        'mass': 259.102902 * u.u,
        'stable': False,
        'half-life': 6.2 * u.s,
    },

    'Lr-260': {
        'atomic number': 103,
        'mass number': 260,
        'mass': 260.1055 * u.u,
        'stable': False,
        'half-life': 180.0 * u.s,
    },

    'Lr-261': {
        'atomic number': 103,
        'mass number': 261,
        'mass': 261.10688 * u.u,
        'stable': False,
        'half-life': 2340.0 * u.s,
    },

    'Lr-262': {
        'atomic number': 103,
        'mass number': 262,
        'mass': 262.10961 * u.u,
        'stable': False,
        'half-life': '~4 h',
    },

    'Lr-263': {
        'atomic number': 103,
        'mass number': 263,
        'mass': 263.11136 * u.u,
        'stable': False,
        'half-life': '5# h',
    },

    'Lr-264': {
        'atomic number': 103,
        'mass number': 264,
        'mass': 264.1142 * u.u,
        'stable': False,
        'half-life': '10# h',
    },

    'Lr-265': {
        'atomic number': 103,
        'mass number': 265,
        'mass': 265.11619 * u.u,
        'stable': False,
        'half-life': '10# h',
    },

    'Lr-266': {
        'atomic number': 103,
        'mass number': 266,
        'mass': 266.11983 * u.u,
        'stable': False,
        'half-life': 75600.0 * u.s,
    },

    'Rf-253': {
        'atomic number': 104,
        'mass number': 253,
        'mass': 253.10044 * u.u,
        'stable': False,
        'half-life': 0.013 * u.s,
    },

    'Rf-254': {
        'atomic number': 104,
        'mass number': 254,
        'mass': 254.10005 * u.u,
        'stable': False,
        'half-life': 2.32e-05 * u.s,
    },

    'Rf-255': {
        'atomic number': 104,
        'mass number': 255,
        'mass': 255.10127 * u.u,
        'stable': False,
        'half-life': 1.66 * u.s,
    },

    'Rf-256': {
        'atomic number': 104,
        'mass number': 256,
        'mass': 256.101152 * u.u,
        'stable': False,
        'half-life': 0.00667 * u.s,
    },

    'Rf-257': {
        'atomic number': 104,
        'mass number': 257,
        'mass': 257.102918 * u.u,
        'stable': False,
        'half-life': 4.82 * u.s,
    },

    'Rf-258': {
        'atomic number': 104,
        'mass number': 258,
        'mass': 258.103428 * u.u,
        'stable': False,
        'half-life': 0.0138 * u.s,
    },

    'Rf-259': {
        'atomic number': 104,
        'mass number': 259,
        'mass': 259.105596 * u.u,
        'stable': False,
        'half-life': 2.63 * u.s,
    },

    'Rf-260': {
        'atomic number': 104,
        'mass number': 260,
        'mass': 260.10644 * u.u,
        'stable': False,
        'half-life': 0.021 * u.s,
    },

    'Rf-261': {
        'atomic number': 104,
        'mass number': 261,
        'mass': 261.108773 * u.u,
        'stable': False,
        'half-life': 2.2 * u.s,
    },

    'Rf-262': {
        'atomic number': 104,
        'mass number': 262,
        'mass': 262.10992 * u.u,
        'stable': False,
        'half-life': 0.25 * u.s,
    },

    'Rf-263': {
        'atomic number': 104,
        'mass number': 263,
        'mass': 263.11249 * u.u,
        'stable': False,
        'half-life': 660.0 * u.s,
    },

    'Rf-264': {
        'atomic number': 104,
        'mass number': 264,
        'mass': 264.11388 * u.u,
        'stable': False,
        'half-life': '1# h',
    },

    'Rf-265': {
        'atomic number': 104,
        'mass number': 265,
        'mass': 265.11668 * u.u,
        'stable': False,
        'half-life': 96.0 * u.s,
    },

    'Rf-266': {
        'atomic number': 104,
        'mass number': 266,
        'mass': 266.11817 * u.u,
        'stable': False,
        'half-life': '4# h',
    },

    'Rf-267': {
        'atomic number': 104,
        'mass number': 267,
        'mass': 267.12179 * u.u,
        'stable': False,
        'half-life': 9000.0 * u.s,
    },

    'Rf-268': {
        'atomic number': 104,
        'mass number': 268,
        'mass': 268.12397 * u.u,
        'stable': False,
        'half-life': '1# h',
    },

    'Db-255': {
        'atomic number': 105,
        'mass number': 255,
        'mass': 255.10707 * u.u,
        'stable': False,
        'half-life': 1.7 * u.s,
    },

    'Db-256': {
        'atomic number': 105,
        'mass number': 256,
        'mass': 256.10789 * u.u,
        'stable': False,
        'half-life': 1.7 * u.s,
    },

    'Db-257': {
        'atomic number': 105,
        'mass number': 257,
        'mass': 257.10758 * u.u,
        'stable': False,
        'half-life': 2.3 * u.s,
    },

    'Db-258': {
        'atomic number': 105,
        'mass number': 258,
        'mass': 258.10928 * u.u,
        'stable': False,
        'half-life': 4.5 * u.s,
    },

    'Db-259': {
        'atomic number': 105,
        'mass number': 259,
        'mass': 259.109492 * u.u,
        'stable': False,
        'half-life': 0.51 * u.s,
    },

    'Db-260': {
        'atomic number': 105,
        'mass number': 260,
        'mass': 260.1113 * u.u,
        'stable': False,
        'half-life': 1.52 * u.s,
    },

    'Db-261': {
        'atomic number': 105,
        'mass number': 261,
        'mass': 261.11192 * u.u,
        'stable': False,
        'half-life': 4.7 * u.s,
    },

    'Db-262': {
        'atomic number': 105,
        'mass number': 262,
        'mass': 262.11407 * u.u,
        'stable': False,
        'half-life': 34.0 * u.s,
    },

    'Db-263': {
        'atomic number': 105,
        'mass number': 263,
        'mass': 263.11499 * u.u,
        'stable': False,
        'half-life': 29.0 * u.s,
    },

    'Db-264': {
        'atomic number': 105,
        'mass number': 264,
        'mass': 264.11741 * u.u,
        'stable': False,
        'half-life': '3# m',
    },

    'Db-265': {
        'atomic number': 105,
        'mass number': 265,
        'mass': 265.11861 * u.u,
        'stable': False,
        'half-life': '15# m',
    },

    'Db-266': {
        'atomic number': 105,
        'mass number': 266,
        'mass': 266.12103 * u.u,
        'stable': False,
        'half-life': 4800.0 * u.s,
    },

    'Db-267': {
        'atomic number': 105,
        'mass number': 267,
        'mass': 267.12247 * u.u,
        'stable': False,
        'half-life': 6000.0 * u.s,
    },

    'Db-268': {
        'atomic number': 105,
        'mass number': 268,
        'mass': 268.12567 * u.u,
        'stable': False,
        'half-life': 104400.0 * u.s,
    },

    'Db-269': {
        'atomic number': 105,
        'mass number': 269,
        'mass': 269.12791 * u.u,
        'stable': False,
        'half-life': '3# h',
    },

    'Db-270': {
        'atomic number': 105,
        'mass number': 270,
        'mass': 270.13136 * u.u,
        'stable': False,
        'half-life': 7200.0 * u.s,
    },

    'Sg-258': {
        'atomic number': 106,
        'mass number': 258,
        'mass': 258.11298 * u.u,
        'stable': False,
        'half-life': 0.0027 * u.s,
    },

    'Sg-259': {
        'atomic number': 106,
        'mass number': 259,
        'mass': 259.1144 * u.u,
        'stable': False,
        'half-life': 0.402 * u.s,
    },

    'Sg-260': {
        'atomic number': 106,
        'mass number': 260,
        'mass': 260.114384 * u.u,
        'stable': False,
        'half-life': 0.00495 * u.s,
    },

    'Sg-261': {
        'atomic number': 106,
        'mass number': 261,
        'mass': 261.115949 * u.u,
        'stable': False,
        'half-life': 0.183 * u.s,
    },

    'Sg-262': {
        'atomic number': 106,
        'mass number': 262,
        'mass': 262.116337 * u.u,
        'stable': False,
        'half-life': 0.0109 * u.s,
    },

    'Sg-263': {
        'atomic number': 106,
        'mass number': 263,
        'mass': 263.11829 * u.u,
        'stable': False,
        'half-life': 0.94 * u.s,
    },

    'Sg-264': {
        'atomic number': 106,
        'mass number': 264,
        'mass': 264.11893 * u.u,
        'stable': False,
        'half-life': 0.047 * u.s,
    },

    'Sg-265': {
        'atomic number': 106,
        'mass number': 265,
        'mass': 265.12109 * u.u,
        'stable': False,
        'half-life': 9.2 * u.s,
    },

    'Sg-266': {
        'atomic number': 106,
        'mass number': 266,
        'mass': 266.12198 * u.u,
        'stable': False,
        'half-life': 0.39 * u.s,
    },

    'Sg-267': {
        'atomic number': 106,
        'mass number': 267,
        'mass': 267.12436 * u.u,
        'stable': False,
        'half-life': 108.0 * u.s,
    },

    'Sg-268': {
        'atomic number': 106,
        'mass number': 268,
        'mass': 268.12539 * u.u,
        'stable': False,
        'half-life': '2# m',
    },

    'Sg-269': {
        'atomic number': 106,
        'mass number': 269,
        'mass': 269.12863 * u.u,
        'stable': False,
        'half-life': 300.0 * u.s,
    },

    'Sg-270': {
        'atomic number': 106,
        'mass number': 270,
        'mass': 270.13043 * u.u,
        'stable': False,
        'half-life': '3# m',
    },

    'Sg-271': {
        'atomic number': 106,
        'mass number': 271,
        'mass': 271.13393 * u.u,
        'stable': False,
        'half-life': 186.0 * u.s,
    },

    'Sg-272': {
        'atomic number': 106,
        'mass number': 272,
        'mass': 272.13589 * u.u,
        'stable': False,
        'half-life': '4# m',
    },

    'Sg-273': {
        'atomic number': 106,
        'mass number': 273,
        'mass': 273.13958 * u.u,
        'stable': False,
        'half-life': '5# m',
    },

    'Bh-260': {
        'atomic number': 107,
        'mass number': 260,
        'mass': 260.12166 * u.u,
        'stable': False,
        'half-life': 0.041 * u.s,
    },

    'Bh-261': {
        'atomic number': 107,
        'mass number': 261,
        'mass': 261.12145 * u.u,
        'stable': False,
        'half-life': 0.0128 * u.s,
    },

    'Bh-262': {
        'atomic number': 107,
        'mass number': 262,
        'mass': 262.12297 * u.u,
        'stable': False,
        'half-life': 0.084 * u.s,
    },

    'Bh-263': {
        'atomic number': 107,
        'mass number': 263,
        'mass': 263.12292 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'Bh-264': {
        'atomic number': 107,
        'mass number': 264,
        'mass': 264.12459 * u.u,
        'stable': False,
        'half-life': 1.07 * u.s,
    },

    'Bh-265': {
        'atomic number': 107,
        'mass number': 265,
        'mass': 265.12491 * u.u,
        'stable': False,
        'half-life': 1.19 * u.s,
    },

    'Bh-266': {
        'atomic number': 107,
        'mass number': 266,
        'mass': 266.12679 * u.u,
        'stable': False,
        'half-life': 2.5 * u.s,
    },

    'Bh-267': {
        'atomic number': 107,
        'mass number': 267,
        'mass': 267.1275 * u.u,
        'stable': False,
        'half-life': 22.0 * u.s,
    },

    'Bh-268': {
        'atomic number': 107,
        'mass number': 268,
        'mass': 268.12969 * u.u,
        'stable': False,
        'half-life': '25# s',
    },

    'Bh-269': {
        'atomic number': 107,
        'mass number': 269,
        'mass': 269.13042 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Bh-270': {
        'atomic number': 107,
        'mass number': 270,
        'mass': 270.13336 * u.u,
        'stable': False,
        'half-life': 228.0 * u.s,
    },

    'Bh-271': {
        'atomic number': 107,
        'mass number': 271,
        'mass': 271.13526 * u.u,
        'stable': False,
        'half-life': 600.0 * u.s,
    },

    'Bh-272': {
        'atomic number': 107,
        'mass number': 272,
        'mass': 272.13826 * u.u,
        'stable': False,
        'half-life': 11.3 * u.s,
    },

    'Bh-273': {
        'atomic number': 107,
        'mass number': 273,
        'mass': 273.14024 * u.u,
        'stable': False,
        'half-life': '1# m',
    },

    'Bh-274': {
        'atomic number': 107,
        'mass number': 274,
        'mass': 274.14355 * u.u,
        'stable': False,
        'half-life': 60.0 * u.s,
    },

    'Bh-275': {
        'atomic number': 107,
        'mass number': 275,
        'mass': 275.14567 * u.u,
        'stable': False,
        'half-life': '5# m',
    },

    'Hs-263': {
        'atomic number': 108,
        'mass number': 263,
        'mass': 263.12852 * u.u,
        'stable': False,
        'half-life': 0.00076 * u.s,
    },

    'Hs-264': {
        'atomic number': 108,
        'mass number': 264,
        'mass': 264.128357 * u.u,
        'stable': False,
        'half-life': 0.00054 * u.s,
    },

    'Hs-265': {
        'atomic number': 108,
        'mass number': 265,
        'mass': 265.129793 * u.u,
        'stable': False,
        'half-life': 0.00196 * u.s,
    },

    'Hs-266': {
        'atomic number': 108,
        'mass number': 266,
        'mass': 266.130046 * u.u,
        'stable': False,
        'half-life': 0.00302 * u.s,
    },

    'Hs-267': {
        'atomic number': 108,
        'mass number': 267,
        'mass': 267.13167 * u.u,
        'stable': False,
        'half-life': 0.055 * u.s,
    },

    'Hs-268': {
        'atomic number': 108,
        'mass number': 268,
        'mass': 268.13186 * u.u,
        'stable': False,
        'half-life': 1.42 * u.s,
    },

    'Hs-269': {
        'atomic number': 108,
        'mass number': 269,
        'mass': 269.13375 * u.u,
        'stable': False,
        'half-life': 16.0 * u.s,
    },

    'Hs-270': {
        'atomic number': 108,
        'mass number': 270,
        'mass': 270.13429 * u.u,
        'stable': False,
        'half-life': 9.0 * u.s,
    },

    'Hs-271': {
        'atomic number': 108,
        'mass number': 271,
        'mass': 271.13717 * u.u,
        'stable': False,
        'half-life': '10# s',
    },

    'Hs-272': {
        'atomic number': 108,
        'mass number': 272,
        'mass': 272.1385 * u.u,
        'stable': False,
        'half-life': '10# s',
    },

    'Hs-273': {
        'atomic number': 108,
        'mass number': 273,
        'mass': 273.14168 * u.u,
        'stable': False,
        'half-life': 1.06 * u.s,
    },

    'Hs-274': {
        'atomic number': 108,
        'mass number': 274,
        'mass': 274.1433 * u.u,
        'stable': False,
        'half-life': '500# ms',
    },

    'Hs-275': {
        'atomic number': 108,
        'mass number': 275,
        'mass': 275.14667 * u.u,
        'stable': False,
        'half-life': 0.29 * u.s,
    },

    'Hs-276': {
        'atomic number': 108,
        'mass number': 276,
        'mass': 276.14846 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Hs-277': {
        'atomic number': 108,
        'mass number': 277,
        'mass': 277.1519 * u.u,
        'stable': False,
        'half-life': 0.011 * u.s,
    },

    'Mt-265': {
        'atomic number': 109,
        'mass number': 265,
        'mass': 265.136 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'Mt-266': {
        'atomic number': 109,
        'mass number': 266,
        'mass': 266.13737 * u.u,
        'stable': False,
        'half-life': 0.0012 * u.s,
    },

    'Mt-267': {
        'atomic number': 109,
        'mass number': 267,
        'mass': 267.13719 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Mt-268': {
        'atomic number': 109,
        'mass number': 268,
        'mass': 268.13865 * u.u,
        'stable': False,
        'half-life': 0.027 * u.s,
    },

    'Mt-269': {
        'atomic number': 109,
        'mass number': 269,
        'mass': 269.13882 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Mt-270': {
        'atomic number': 109,
        'mass number': 270,
        'mass': 270.14033 * u.u,
        'stable': False,
        'half-life': 0.0063 * u.s,
    },

    'Mt-271': {
        'atomic number': 109,
        'mass number': 271,
        'mass': 271.14074 * u.u,
        'stable': False,
        'half-life': '400# ms',
    },

    'Mt-272': {
        'atomic number': 109,
        'mass number': 272,
        'mass': 272.14341 * u.u,
        'stable': False,
        'half-life': '400# ms',
    },

    'Mt-273': {
        'atomic number': 109,
        'mass number': 273,
        'mass': 273.1444 * u.u,
        'stable': False,
        'half-life': '800# ms',
    },

    'Mt-274': {
        'atomic number': 109,
        'mass number': 274,
        'mass': 274.14724 * u.u,
        'stable': False,
        'half-life': 0.85 * u.s,
    },

    'Mt-275': {
        'atomic number': 109,
        'mass number': 275,
        'mass': 275.14882 * u.u,
        'stable': False,
        'half-life': 0.117 * u.s,
    },

    'Mt-276': {
        'atomic number': 109,
        'mass number': 276,
        'mass': 276.15159 * u.u,
        'stable': False,
        'half-life': 0.63 * u.s,
    },

    'Mt-277': {
        'atomic number': 109,
        'mass number': 277,
        'mass': 277.15327 * u.u,
        'stable': False,
        'half-life': 9.0 * u.s,
    },

    'Mt-278': {
        'atomic number': 109,
        'mass number': 278,
        'mass': 278.15631 * u.u,
        'stable': False,
        'half-life': 7.0 * u.s,
    },

    'Mt-279': {
        'atomic number': 109,
        'mass number': 279,
        'mass': 279.15808 * u.u,
        'stable': False,
        'half-life': '30# s',
    },

    'Ds-267': {
        'atomic number': 110,
        'mass number': 267,
        'mass': 267.14377 * u.u,
        'stable': False,
        'half-life': 1e-05 * u.s,
    },

    'Ds-268': {
        'atomic number': 110,
        'mass number': 268,
        'mass': 268.14348 * u.u,
        'stable': False,
        'half-life': '100# us',
    },

    'Ds-269': {
        'atomic number': 110,
        'mass number': 269,
        'mass': 269.144752 * u.u,
        'stable': False,
        'half-life': 0.00023 * u.s,
    },

    'Ds-270': {
        'atomic number': 110,
        'mass number': 270,
        'mass': 270.144584 * u.u,
        'stable': False,
        'half-life': 0.000205 * u.s,
    },

    'Ds-271': {
        'atomic number': 110,
        'mass number': 271,
        'mass': 271.14595 * u.u,
        'stable': False,
        'half-life': 0.09 * u.s,
    },

    'Ds-272': {
        'atomic number': 110,
        'mass number': 272,
        'mass': 272.14602 * u.u,
        'stable': False,
        'half-life': '200# ms',
    },

    'Ds-273': {
        'atomic number': 110,
        'mass number': 273,
        'mass': 273.14856 * u.u,
        'stable': False,
        'half-life': 0.00024 * u.s,
    },

    'Ds-274': {
        'atomic number': 110,
        'mass number': 274,
        'mass': 274.14941 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Ds-275': {
        'atomic number': 110,
        'mass number': 275,
        'mass': 275.15203 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Ds-276': {
        'atomic number': 110,
        'mass number': 276,
        'mass': 276.15303 * u.u,
        'stable': False,
        'half-life': '100# ms',
    },

    'Ds-277': {
        'atomic number': 110,
        'mass number': 277,
        'mass': 277.15591 * u.u,
        'stable': False,
        'half-life': 0.006 * u.s,
    },

    'Ds-278': {
        'atomic number': 110,
        'mass number': 278,
        'mass': 278.15704 * u.u,
        'stable': False,
        'half-life': '270# ms',
    },

    'Ds-279': {
        'atomic number': 110,
        'mass number': 279,
        'mass': 279.1601 * u.u,
        'stable': False,
        'half-life': 0.21 * u.s,
    },

    'Ds-280': {
        'atomic number': 110,
        'mass number': 280,
        'mass': 280.16131 * u.u,
        'stable': False,
        'half-life': 11.0 * u.s,
    },

    'Ds-281': {
        'atomic number': 110,
        'mass number': 281,
        'mass': 281.16451 * u.u,
        'stable': False,
        'half-life': 14.0 * u.s,
    },

    'Rg-272': {
        'atomic number': 111,
        'mass number': 272,
        'mass': 272.15327 * u.u,
        'stable': False,
        'half-life': 0.0045 * u.s,
    },

    'Rg-273': {
        'atomic number': 111,
        'mass number': 273,
        'mass': 273.15313 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'Rg-274': {
        'atomic number': 111,
        'mass number': 274,
        'mass': 274.15525 * u.u,
        'stable': False,
        'half-life': 0.029 * u.s,
    },

    'Rg-275': {
        'atomic number': 111,
        'mass number': 275,
        'mass': 275.15594 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Rg-276': {
        'atomic number': 111,
        'mass number': 276,
        'mass': 276.15833 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Rg-277': {
        'atomic number': 111,
        'mass number': 277,
        'mass': 277.15907 * u.u,
        'stable': False,
        'half-life': '10# ms',
    },

    'Rg-278': {
        'atomic number': 111,
        'mass number': 278,
        'mass': 278.16149 * u.u,
        'stable': False,
        'half-life': 0.008 * u.s,
    },

    'Rg-279': {
        'atomic number': 111,
        'mass number': 279,
        'mass': 279.16272 * u.u,
        'stable': False,
        'half-life': 0.18 * u.s,
    },

    'Rg-280': {
        'atomic number': 111,
        'mass number': 280,
        'mass': 280.16514 * u.u,
        'stable': False,
        'half-life': 4.3 * u.s,
    },

    'Rg-281': {
        'atomic number': 111,
        'mass number': 281,
        'mass': 281.16636 * u.u,
        'stable': False,
        'half-life': 24.0 * u.s,
    },

    'Rg-282': {
        'atomic number': 111,
        'mass number': 282,
        'mass': 282.16912 * u.u,
        'stable': False,
        'half-life': 96.0 * u.s,
    },

    'Rg-283': {
        'atomic number': 111,
        'mass number': 283,
        'mass': 283.17054 * u.u,
        'stable': False,
        'half-life': '30# s',
    },

    'Cn-276': {
        'atomic number': 112,
        'mass number': 276,
        'mass': 276.16141 * u.u,
        'stable': False,
        'half-life': '100# us',
    },

    'Cn-277': {
        'atomic number': 112,
        'mass number': 277,
        'mass': 277.16364 * u.u,
        'stable': False,
        'half-life': 0.00085 * u.s,
    },

    'Cn-278': {
        'atomic number': 112,
        'mass number': 278,
        'mass': 278.16416 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'Cn-279': {
        'atomic number': 112,
        'mass number': 279,
        'mass': 279.16654 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Cn-280': {
        'atomic number': 112,
        'mass number': 280,
        'mass': 280.16715 * u.u,
        'stable': False,
        'half-life': '5# ms',
    },

    'Cn-281': {
        'atomic number': 112,
        'mass number': 281,
        'mass': 281.16975 * u.u,
        'stable': False,
        'half-life': 0.18 * u.s,
    },

    'Cn-282': {
        'atomic number': 112,
        'mass number': 282,
        'mass': 282.1705 * u.u,
        'stable': False,
        'half-life': 0.0009 * u.s,
    },

    'Cn-283': {
        'atomic number': 112,
        'mass number': 283,
        'mass': 283.17327 * u.u,
        'stable': False,
        'half-life': 4.1 * u.s,
    },

    'Cn-284': {
        'atomic number': 112,
        'mass number': 284,
        'mass': 284.17416 * u.u,
        'stable': False,
        'half-life': 0.104 * u.s,
    },

    'Cn-285': {
        'atomic number': 112,
        'mass number': 285,
        'mass': 285.17712 * u.u,
        'stable': False,
        'half-life': 32.0 * u.s,
    },

    'Nh-278': {
        'atomic number': 113,
        'mass number': 278,
        'mass': 278.17058 * u.u,
        'stable': False,
        'half-life': 0.00034 * u.s,
    },

    'Nh-279': {
        'atomic number': 113,
        'mass number': 279,
        'mass': 279.17095 * u.u,
        'stable': False,
    },

    'Nh-280': {
        'atomic number': 113,
        'mass number': 280,
        'mass': 280.17293 * u.u,
        'stable': False,
    },

    'Nh-281': {
        'atomic number': 113,
        'mass number': 281,
        'mass': 281.17348 * u.u,
        'stable': False,
    },

    'Nh-282': {
        'atomic number': 113,
        'mass number': 282,
        'mass': 282.17567 * u.u,
        'stable': False,
        'half-life': 0.073 * u.s,
    },

    'Nh-283': {
        'atomic number': 113,
        'mass number': 283,
        'mass': 283.17657 * u.u,
        'stable': False,
        'half-life': 0.1 * u.s,
    },

    'Nh-284': {
        'atomic number': 113,
        'mass number': 284,
        'mass': 284.17873 * u.u,
        'stable': False,
        'half-life': 0.48 * u.s,
    },

    'Nh-285': {
        'atomic number': 113,
        'mass number': 285,
        'mass': 285.17973 * u.u,
        'stable': False,
        'half-life': 5.5 * u.s,
    },

    'Nh-286': {
        'atomic number': 113,
        'mass number': 286,
        'mass': 286.18221 * u.u,
        'stable': False,
        'half-life': 19.6 * u.s,
    },

    'Nh-287': {
        'atomic number': 113,
        'mass number': 287,
        'mass': 287.18339 * u.u,
        'stable': False,
        'half-life': 5.5 * u.s,
    },

    'Fl-285': {
        'atomic number': 114,
        'mass number': 285,
        'mass': 285.18364 * u.u,
        'stable': False,
        'half-life': 0.21 * u.s,
    },

    'Fl-286': {
        'atomic number': 114,
        'mass number': 286,
        'mass': 286.18423 * u.u,
        'stable': False,
        'half-life': 0.14 * u.s,
    },

    'Fl-287': {
        'atomic number': 114,
        'mass number': 287,
        'mass': 287.18678 * u.u,
        'stable': False,
        'half-life': 0.52 * u.s,
    },

    'Fl-288': {
        'atomic number': 114,
        'mass number': 288,
        'mass': 288.18757 * u.u,
        'stable': False,
        'half-life': 0.75 * u.s,
    },

    'Fl-289': {
        'atomic number': 114,
        'mass number': 289,
        'mass': 289.19042 * u.u,
        'stable': False,
        'half-life': 2.4 * u.s,
    },

    'Mc-287': {
        'atomic number': 115,
        'mass number': 287,
        'mass': 287.1907 * u.u,
        'stable': False,
        'half-life': 0.037 * u.s,
    },

    'Mc-288': {
        'atomic number': 115,
        'mass number': 288,
        'mass': 288.19274 * u.u,
        'stable': False,
        'half-life': 0.164 * u.s,
    },

    'Mc-289': {
        'atomic number': 115,
        'mass number': 289,
        'mass': 289.19363 * u.u,
        'stable': False,
        'half-life': 0.33 * u.s,
    },

    'Mc-290': {
        'atomic number': 115,
        'mass number': 290,
        'mass': 290.19598 * u.u,
        'stable': False,
        'half-life': 0.65 * u.s,
    },

    'Mc-291': {
        'atomic number': 115,
        'mass number': 291,
        'mass': 291.19707 * u.u,
        'stable': False,
    },

    'Lv-289': {
        'atomic number': 116,
        'mass number': 289,
        'mass': 289.19816 * u.u,
        'stable': False,
        'half-life': '2# ms',
    },

    'Lv-290': {
        'atomic number': 116,
        'mass number': 290,
        'mass': 290.19864 * u.u,
        'stable': False,
        'half-life': 0.008 * u.s,
    },

    'Lv-291': {
        'atomic number': 116,
        'mass number': 291,
        'mass': 291.20108 * u.u,
        'stable': False,
        'half-life': 0.028 * u.s,
    },

    'Lv-292': {
        'atomic number': 116,
        'mass number': 292,
        'mass': 292.20174 * u.u,
        'stable': False,
        'half-life': 0.024 * u.s,
    },

    'Lv-293': {
        'atomic number': 116,
        'mass number': 293,
        'mass': 293.20449 * u.u,
        'stable': False,
        'half-life': 0.08 * u.s,
    },

    'Ts-291': {
        'atomic number': 117,
        'mass number': 291,
        'mass': 291.20553 * u.u,
        'stable': False,
    },

    'Ts-292': {
        'atomic number': 117,
        'mass number': 292,
        'mass': 292.20746 * u.u,
        'stable': False,
    },

    'Ts-293': {
        'atomic number': 117,
        'mass number': 293,
        'mass': 293.20824 * u.u,
        'stable': False,
        'half-life': 0.022 * u.s,
    },

    'Ts-294': {
        'atomic number': 117,
        'mass number': 294,
        'mass': 294.21046 * u.u,
        'stable': False,
        'half-life': 0.051 * u.s,
    },

    'Og-293': {
        'atomic number': 118,
        'mass number': 293,
        'mass': 293.21356 * u.u,
        'stable': False,
    },

    'Og-294': {
        'atomic number': 118,
        'mass number': 294,
        'mass': 294.21392 * u.u,
        'stable': False,
        'half-life': 0.0007 * u.s,
    },

    'Og-295': {
        'atomic number': 118,
        'mass number': 295,
        'mass': 295.21624 * u.u,
        'stable': False,
        'half-life': 0.181 * u.s,
    },

}
