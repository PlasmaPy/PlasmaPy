"""
Dictionaries containing basic atomic data.

The periodic tabla data is from: http://periodic.lanl.gov/index.shtml
"""

import collections
import astropy.units as u

_PeriodicTable = collections.namedtuple(
    "periodic_table", ['group', 'category', 'block', 'period']
)

_Elements = {

    "H": {
        "atomic number": 1,
        "atomic mass": 1.008 * u.u,
        "element name": "hydrogen",
        "period": 1,
        "group": 1,
        "block": "s",
        "category": "nonmetal",
    },

    "He": {
        "atomic number": 2,
        "atomic mass": 4.002602 * u.u,
        "element name": "helium",
        "period": 1,
        "group": 18,
        "block": "s",
        "category": "noble gas",
    },

    "Li": {
        "atomic number": 3,
        "atomic mass": 6.94 * u.u,
        "element name": "lithium",
        "period": 2,
        "group": 1,
        "block": "s",
        "category": "alkali metal",
    },

    "Be": {
        "atomic number": 4,
        "atomic mass": 9.0121831 * u.u,
        "element name": "beryllium",
        "period": 2,
        "group": 2,
        "block": "s",
        "category": "alkaline earth metal",
    },

    "B": {
        "atomic number": 5,
        "atomic mass": 10.806 * u.u,
        "element name": "boron",
        "period": 2,
        "group": 13,
        "block": "p",
        "category": "metalloid",
    },

    "C": {
        "atomic number": 6,
        "atomic mass": 12.011 * u.u,
        "element name": "carbon",
        "period": 2,
        "group": 14,
        "block": "p",
        "category": "nonmetal",
    },

    "N": {
        "atomic number": 7,
        "atomic mass": 14.007 * u.u,
        "element name": "nitrogen",
        "period": 2,
        "group": 15,
        "block": "p",
        "category": "nonmetal",
    },

    "O": {
        "atomic number": 8,
        "atomic mass": 15.999 * u.u,
        "element name": "oxygen",
        "period": 2,
        "group": 16,
        "block": "p",
        "category": "nonmetal",
    },

    "F": {
        "atomic number": 9,
        "atomic mass": 18.998403163 * u.u,
        "element name": "fluorine",
        "period": 2,
        "group": 17,
        "block": "p",
        "category": "halogen",
    },

    "Ne": {
        "atomic number": 10,
        "atomic mass": 20.1797 * u.u,
        "element name": "neon",
        "period": 2,
        "group": 18,
        "block": "p",
        "category": "noble gas",
    },

    "Na": {
        "atomic number": 11,
        "atomic mass": 22.98976928 * u.u,
        "element name": "sodium",
        "period": 3,
        "group": 1,
        "block": "s",
        "category": "alkali metal",
    },

    "Mg": {
        "atomic number": 12,
        "atomic mass": 24.305 * u.u,
        "element name": "magnesium",
        "period": 3,
        "group": 2,
        "block": "s",
        "category": "alkaline earth metal",
    },

    "Al": {
        "atomic number": 13,
        "atomic mass": 26.9815385 * u.u,
        "element name": "aluminium",
        "period": 3,
        "group": 13,
        "block": "p",
        "category": "post-transition metal",
    },

    "Si": {
        "atomic number": 14,
        "atomic mass": 28.085 * u.u,
        "element name": "silicon",
        "period": 3,
        "group": 14,
        "block": "p",
        "category": "metalloid",
    },

    "P": {
        "atomic number": 15,
        "atomic mass": 30.973761998 * u.u,
        "element name": "phosphorus",
        "period": 3,
        "group": 15,
        "block": "p",
        "category": "nonmetal",
    },

    "S": {
        "atomic number": 16,
        "symbol": "S",
        "atomic mass": 32.06 * u.u,
        "element name": "sulfur",
        "period": 3,
        "group": 16,
        "block": "p",
        "category": "nonmetal",
    },

    "Cl": {
        "atomic number": 17,
        "atomic mass": 35.45 * u.u,
        "element name": "chlorine",
        "period": 3,
        "group": 17,
        "block": "p",
        "category": "halogen",
    },

    "Ar": {
        "atomic number": 18,
        "atomic mass": 39.948 * u.u,
        "element name": "argon",
        "period": 3,
        "group": 18,
        "block": "p",
        "category": "noble gas",
    },

    "K": {
        "atomic number": 19,
        "atomic mass": 39.0983 * u.u,
        "element name": "potassium",
        "period": 4,
        "group": 1,
        "block": "s",
        "category": "alkali metal",
    },

    "Ca": {
        "atomic number": 20,
        "atomic mass": 40.078 * u.u,
        "element name": "calcium",
        "period": 4,
        "group": 2,
        "block": "s",
        "category": "alkaline earth metal",
    },

    "Sc": {
        "atomic number": 21,
        "atomic mass": 44.955908 * u.u,
        "element name": "scandium",
        "period": 4,
        "group": 3,
        "block": "d",
        "category": "transition metal",
    },

    "Ti": {
        "atomic number": 22,
        "atomic mass": 47.867 * u.u,
        "element name": "titanium",
        "period": 4,
        "group": 4,
        "block": "d",
        "category": "transition metal",
    },

    "V": {
        "atomic number": 23,
        "atomic mass": 50.9415 * u.u,
        "element name": "vanadium",
        "period": 4,
        "group": 5,
        "block": "d",
        "category": "transition metal",
    },

    "Cr": {
        "atomic number": 24,
        "atomic mass": 51.9961 * u.u,
        "element name": "chromium",
        "period": 4,
        "group": 6,
        "block": "d",
        "category": "transition metal",
    },

    "Mn": {
        "atomic number": 25,
        "atomic mass": 54.938044 * u.u,
        "element name": "manganese",
        "period": 4,
        "group": 7,
        "block": "d",
        "category": "transition metal",
    },

    "Fe": {
        "atomic number": 26,
        "atomic mass": 55.845 * u.u,
        "element name": "iron",
        "period": 4,
        "group": 8,
        "block": "d",
        "category": "transition metal",
    },

    "Co": {
        "atomic number": 27,
        "atomic mass": 58.933 * u.u,
        "element name": "cobalt",
        "period": 4,
        "group": 9,
        "block": "d",
        "category": "transition metal",
    },

    "Ni": {
        "atomic number": 28,
        "atomic mass": 58.6934 * u.u,
        "element name": "nickel",
        "period": 4,
        "group": 10,
        "block": "d",
        "category": "transition metal",
    },

    "Cu": {
        "atomic number": 29,
        "atomic mass": 63.546 * u.u,
        "element name": "copper",
        "period": 4,
        "group": 11,
        "block": "d",
        "category": "transition metal",
    },

    "Zn": {
        "atomic number": 30,
        "atomic mass": 65.38 * u.u,
        "element name": "zinc",
        "period": 4,
        "group": 12,
        "block": "d",
        "category": "transition metal",
    },

    "Ga": {
        "atomic number": 31,
        "atomic mass": 69.723 * u.u,
        "element name": "gallium",
        "period": 4,
        "group": 13,
        "block": "p",
        "category": "post-transition metal",
    },

    "Ge": {
        "atomic number": 32,
        "atomic mass": 72.630 * u.u,
        "element name": "germanium",
        "period": 4,
        "group": 14,
        "block": "p",
        "category": "metalloid",
    },

    "As": {
        "atomic number": 33,
        "atomic mass": 74.921595 * u.u,
        "element name": "arsenic",
        "period": 4,
        "group": 15,
        "block": "p",
        "category": "metalloid",
    },

    "Se": {
        "atomic number": 34,
        "atomic mass": 78.971 * u.u,
        "element name": "selenium",
        "period": 4,
        "group": 16,
        "block": "p",
        "category": "nonmetal",
    },

    "Br": {
        "atomic number": 35,
        "atomic mass": 79.904 * u.u,
        "element name": "bromine",
        "period": 4,
        "group": 17,
        "block": "p",
        "category": "halogen",
    },

    "Kr": {
        "atomic number": 36,
        "atomic mass": 83.798 * u.u,
        "element name": "krypton",
        "period": 4,
        "group": 18,
        "block": "p",
        "category": "noble gas",
    },

    "Rb": {
        "atomic number": 37,
        "atomic mass": 85.4678 * u.u,
        "element name": "rubidium",
        "period": 5,
        "group": 1,
        "block": "s",
        "category": "alkali metal",
    },

    "Sr": {
        "atomic number": 38,
        "atomic mass": 87.62 * u.u,
        "element name": "strontium",
        "period": 5,
        "group": 2,
        "block": "s",
        "category": "alkaline earth metal",
    },

    "Y": {
        "atomic number": 39,
        "atomic mass": 88.90584 * u.u,
        "element name": "yttrium",
        "period": 5,
        "group": 3,
        "block": "d",
        "category": "transition metal",
    },

    "Zr": {
        "atomic number": 40,
        "atomic mass": 91.224 * u.u,
        "element name": "zirconium",
        "period": 5,
        "group": 4,
        "block": "d",
        "category": "transition metal",
    },

    "Nb": {
        "atomic number": 41,
        "atomic mass": 92.90637 * u.u,
        "element name": "niobium",
        "period": 5,
        "group": 5,
        "block": "d",
        "category": "transition metal",
    },

    "Mo": {
        "atomic number": 42,
        "atomic mass": 95.95 * u.u,
        "element name": "molybdenum",
        "period": 5,
        "group": 6,
        "block": "d",
        "category": "transition metal",
    },

    "Tc": {
        "atomic number": 43,
        "element name": "technetium",
        "period": 5,
        "group": 7,
        "block": "d",
        "category": "transition metal",
    },

    "Ru": {
        "atomic number": 44,
        "atomic mass": 101.07 * u.u,
        "element name": "ruthenium",
        "period": 5,
        "group": 8,
        "block": "d",
        "category": "transition metal",
    },

    "Rh": {
        "atomic number": 45,
        "atomic mass": 102.90550 * u.u,
        "element name": "rhodium",
        "period": 5,
        "group": 9,
        "block": "d",
        "category": "transition metal",
    },

    "Pd": {
        "atomic number": 46,
        "atomic mass": 106.42 * u.u,
        "element name": "palladium",
        "period": 5,
        "group": 10,
        "block": "d",
        "category": "transition metal",
    },

    "Ag": {
        "atomic number": 47,
        "atomic mass": 107.8682 * u.u,
        "element name": "silver",
        "period": 5,
        "group": 11,
        "block": "d",
        "category": "transition metal",
    },

    "Cd": {
        "atomic number": 48,
        "atomic mass": 112.414 * u.u,
        "element name": "cadmium",
        "period": 5,
        "group": 12,
        "block": "d",
        "category": "transition metal",
    },

    "In": {
        "atomic number": 49,
        "atomic mass": 114.818 * u.u,
        "element name": "indium",
        "period": 5,
        "group": 13,
        "block": "p",
        "category": "post-transition metal",
    },

    "Sn": {
        "atomic number": 50,
        "atomic mass": 118.710 * u.u,
        "element name": "tin",
        "period": 5,
        "group": 14,
        "block": "p",
        "category": "post-transition metal",
    },

    "Sb": {
        "atomic number": 51,
        "atomic mass": 121.760 * u.u,
        "element name": "antimony",
        "period": 5,
        "group": 15,
        "block": "p",
        "category": "metalloid",
    },

    "Te": {
        "atomic number": 52,
        "atomic mass": 127.60 * u.u,
        "element name": "tellurium",
        "period": 5,
        "group": 16,
        "block": "p",
        "category": "metalloid",
    },

    "I": {
        "atomic number": 53,
        "atomic mass": 126.90447 * u.u,
        "element name": "iodine",
        "period": 5,
        "group": 17,
        "block": "p",
        "category": "halogen",
    },

    "Xe": {
        "atomic number": 54,
        "atomic mass": 131.293 * u.u,
        "element name": "xenon",
        "period": 5,
        "group": 18,
        "block": "p",
        "category": "noble gas",
    },

    "Cs": {
        "atomic number": 55,
        "atomic mass": 132.90545196 * u.u,
        "element name": "caesium",
        "period": 6,
        "group": 1,
        "block": "s",
        "category": "alkali metal",
    },

    "Ba": {
        "atomic number": 56,
        "atomic mass": 137.327 * u.u,
        "element name": "barium",
        "period": 6,
        "group": 2,
        "block": "s",
        "category": "alkaline earth metal",
    },

    "La": {
        "atomic number": 57,
        "atomic mass": 138.90547 * u.u,
        "element name": "lanthanum",
        "period": 6,
        "group": 3,
        "block": "d",
        "category": "lanthanide",
    },

    "Ce": {
        "atomic number": 58,
        "atomic mass": 140.116 * u.u,
        "element name": "cerium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Pr": {
        "atomic number": 59,
        "atomic mass": 140.90766 * u.u,
        "element name": "praseodymium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Nd": {
        "atomic number": 60,
        "atomic mass": 144.242 * u.u,
        "element name": "neodymium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Pm": {
        "atomic number": 61,
        "element name": "promethium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Sm": {
        "atomic number": 62,
        "atomic mass": 150.36 * u.u,
        "element name": "samarium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Eu": {
        "atomic number": 63,
        "atomic mass": 151.964 * u.u,
        "element name": "europium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Gd": {
        "atomic number": 64,
        "atomic mass": 157.25 * u.u,
        "element name": "gadolinium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Tb": {
        "atomic number": 65,
        "atomic mass": 158.92535 * u.u,
        "element name": "terbium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Dy": {
        "atomic number": 66,
        "atomic mass": 162.500 * u.u,
        "element name": "dysprosium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Ho": {
        "atomic number": 67,
        "atomic mass": 164.93033 * u.u,
        "element name": "holmium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Er": {
        "atomic number": 68,
        "atomic mass": 167.259 * u.u,
        "element name": "erbium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Tm": {
        "atomic number": 69,
        "atomic mass": 168.93422 * u.u,
        "element name": "thulium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Yb": {
        "atomic number": 70,
        "atomic mass": 173.045 * u.u,
        "element name": "ytterbium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Lu": {
        "atomic number": 71,
        "atomic mass": 174.9668 * u.u,
        "element name": "lutetium",
        "period": 6,
        "group": 3,
        "block": "f",
        "category": "lanthanide",
    },

    "Hf": {
        "atomic number": 72,
        "atomic mass": 178.49 * u.u,
        "element name": "hafnium",
        "period": 6,
        "group": 4,
        "block": "d",
        "category": "transition metal",
    },

    "Ta": {
        "atomic number": 73,
        "atomic mass": 180.94788 * u.u,
        "element name": "tantalum",
        "period": 6,
        "group": 5,
        "block": "d",
        "category": "transition metal",
    },

    "W": {
        "atomic number": 74,
        "atomic mass": 183.84 * u.u,
        "element name": "tungsten",
        "period": 6,
        "group": 6,
        "block": "d",
        "category": "transition metal",
    },

    "Re": {
        "atomic number": 75,
        "atomic mass": 186.207 * u.u,
        "element name": "rhenium",
        "period": 6,
        "group": 7,
        "block": "d",
        "category": "transition metal",
    },

    "Os": {
        "atomic number": 76,
        "atomic mass": 190.23 * u.u,
        "element name": "osmium",
        "period": 6,
        "group": 8,
        "block": "d",
        "category": "transition metal",
    },

    "Ir": {
        "atomic number": 77,
        "atomic mass": 192.217 * u.u,
        "element name": "iridium",
        "period": 6,
        "group": 9,
        "block": "d",
        "category": "transition metal",
    },

    "Pt": {
        "atomic number": 78,
        "atomic mass": 195.084 * u.u,
        "element name": "platinum",
        "period": 6,
        "group": 10,
        "block": "d",
        "category": "transition metal",
    },

    "Au": {
        "atomic number": 79,
        "atomic mass": 196.966569 * u.u,
        "element name": "gold",
        "period": 6,
        "group": 11,
        "block": "d",
        "category": "transition metal",
    },

    "Hg": {
        "atomic number": 80,
        "atomic mass": 200.592 * u.u,
        "element name": "mercury",
        "period": 6,
        "group": 12,
        "block": "d",
        "category": "transition metal",
    },

    "Tl": {
        "atomic number": 81,
        "atomic mass": 204.38 * u.u,
        "element name": "thallium",
        "period": 6,
        "group": 13,
        "block": "p",
        "category": "post-transition metal",
    },

    "Pb": {
        "atomic number": 82,
        "atomic mass": 207.2 * u.u,
        "element name": "lead",
        "period": 6,
        "group": 14,
        "block": "p",
        "category": "post-transition metal",
    },

    "Bi": {
        "atomic number": 83,
        "atomic mass": 208.98040 * u.u,
        "element name": "bismuth",
        "period": 6,
        "group": 15,
        "block": "p",
        "category": "post-transition metal",
    },

    "Po": {
        "atomic number": 84,
        "element name": "polonium",
        "period": 6,
        "group": 16,
        "block": "p",
        "category": "metalloid",
    },

    "At": {
        "atomic number": 85,
        "element name": "astatine",
        "period": 6,
        "group": 17,
        "block": "p",
        "category": "halogen",
    },

    "Rn": {
        "atomic number": 86,
        "element name": "radon",
        "period": 6,
        "group": 18,
        "block": "p",
        "category": "noble gas",
    },

    "Fr": {
        "atomic number": 87,
        "element name": "francium",
        "period": 7,
        "group": 1,
        "block": "s",
        "category": "alkali metal",
    },

    "Ra": {
        "atomic number": 88,
        "element name": "radium",
        "period": 7,
        "group": 2,
        "block": "s",
        "category": "alkaline earth metal",
    },

    "Ac": {
        "atomic number": 89,
        "element name": "actinium",
        "period": 7,
        "group": 3,
        "block": "d",
        "category": "actinide",
    },

    "Th": {
        "atomic number": 90,
        "atomic mass": 232.0377 * u.u,
        "element name": "thorium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Pa": {
        "atomic number": 91,
        "atomic mass": 231.03588 * u.u,
        "element name": "protactinium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "U": {
        "atomic number": 92,
        "atomic mass": 238.02891 * u.u,
        "element name": "uranium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Np": {
        "atomic number": 93,
        "element name": "neptunium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Pu": {
        "atomic number": 94,
        "element name": "plutonium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Am": {
        "atomic number": 95,
        "element name": "americium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Cm": {
        "atomic number": 96,
        "element name": "curium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Bk": {
        "atomic number": 97,
        "element name": "berkelium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Cf": {
        "atomic number": 98,
        "element name": "californium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Es": {
        "atomic number": 99,
        "element name": "einsteinium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Fm": {
        "atomic number": 100,
        "element name": "fermium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Md": {
        "atomic number": 101,
        "element name": "mendelevium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "No": {
        "atomic number": 102,
        "element name": "nobelium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Lr": {
        "atomic number": 103,
        "element name": "lawrencium",
        "period": 7,
        "group": 3,
        "block": "f",
        "category": "actinide",
    },

    "Rf": {
        "atomic number": 104,
        "element name": "rutherfordium",
        "period": 7,
        "group": 4,
        "block": "d",
        "category": "transition metal",
    },

    "Db": {
        "atomic number": 105,
        "element name": "dubnium",
        "period": 7,
        "group": 5,
        "block": "d",
        "category": "transition metal",
    },

    "Sg": {
        "atomic number": 106,
        "element name": "seaborgium",
        "period": 7,
        "group": 6,
        "block": "d",
        "category": "transition metal",
    },

    "Bh": {
        "atomic number": 107,
        "element name": "bohrium",
        "period": 7,
        "group": 7,
        "block": "d",
        "category": "transition metal",
    },

    "Hs": {
        "atomic number": 108,
        "element name": "hassium",
        "period": 7,
        "group": 8,
        "block": "d",
        "category": "transition metal",
    },

    "Mt": {
        "atomic number": 109,
        "element name": "meitnerium",
        "period": 7,
        "group": 9,
        "block": "d",
        "category": "transition metal",
    },

    "Ds": {
        "atomic number": 110,
        "element name": "darmstadtium",
        "period": 7,
        "group": 10,
        "block": "d",
        "category": "transition metal",
    },

    "Rg": {
        "atomic number": 111,
        "element name": "roentgenium",
        "period": 7,
        "group": 11,
        "block": "d",
        "category": "transition metal",
    },

    "Cn": {
        "atomic number": 112,
        "element name": "copernicium",
        "period": 7,
        "group": 12,
        "block": "d",
        "category": "transition metal",
    },

    "Nh": {
        "atomic number": 113,
        "element name": "nihonium",
        "period": 7,
        "group": 13,
        "block": "p",
        "category": "post-transition metal",
    },

    "Fl": {
        "atomic number": 114,
        "element name": "flerovium",
        "period": 7,
        "group": 14,
        "block": "p",
        "category": "post-transition metal",
    },

    "Mc": {
        "atomic number": 115,
        "element name": "moscovium",
        "period": 7,
        "group": 15,
        "block": "p",
        "category": "post-transition metal",
    },

    "Lv": {
        "atomic number": 116,
        "element name": "livermorium",
        "period": 7,
        "group": 16,
        "block": "p",
        "category": "post-transition metal",
    },

    "Ts": {
        "atomic number": 117,
        "element name": "tennessine",
        "period": 7,
        "group": 17,
        "block": "p",
        "category": "halogen",
    },

    "Og": {
        "atomic number": 118,
        "element name": "oganesson",
        "period": 7,
        "group": 18,
        "block": "p",
        "category": "noble gas",
    },

}

_atomic_numbers_to_symbols = {
    elemdict['atomic number']: symb for (symb, elemdict) in _Elements.items()
}

_element_names_to_symbols = {
    elemdict['element name']: symb for (symb, elemdict) in _Elements.items()
}
