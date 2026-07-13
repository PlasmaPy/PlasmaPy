"""Bosch-Hale and Nevins-Swain fusion reaction parameters.

References
----------
H.-S. Bosch and G. M. Hale, "Improved formulas for fusion cross-sections
and thermal reactivities," Nuclear Fusion, vol. 32, no. 4, pp. 611-631, 1992.

W. M. Nevins and R. Swain, "The thermonuclear fusion rate coefficient for
p-11B reactions," Nuclear Fusion, vol. 40, no. 4, pp. 865-872, 2000.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class _ReactionParameters:
    """Parameters for a single fusion reaction channel."""

    name: str
    reactant1: str
    reactant2: str
    products: tuple[str, ...]
    BG: float
    m1c2_keV: float
    m2c2_keV: float
    A1: float
    A2: float
    A3: float
    A4: float
    A5: float
    B1: float
    B2: float
    B3: float
    B4: float
    C1: float
    C2: float
    C3: float
    C4: float
    C5: float
    C6: float
    C7: float
    min_E_keV: float
    max_E_keV: float
    min_T_keV: float
    max_T_keV: float
    Q_keV: float
    fraction_charged: float
    identical_particles: bool
    Z1: int
    Z2: int
    high_energy_A1: float = 0.0
    high_energy_A2: float = 0.0
    high_energy_A3: float = 0.0
    high_energy_A4: float = 0.0
    high_energy_A5: float = 0.0
    high_energy_B1: float = 0.0
    high_energy_B2: float = 0.0
    high_energy_B3: float = 0.0
    high_energy_B4: float = 0.0
    high_energy_cutoff_keV: float = 0.0
    high_energy_max_keV: float = 0.0


DT = _ReactionParameters(
    name="T(d,n)4He",
    reactant1="D",
    reactant2="T",
    products=("He-4", "n"),
    BG=34.3827,
    m1c2_keV=1875612.94257,
    m2c2_keV=2808921.13298,
    A1=6.927e4,
    A2=7.454e8,
    A3=2.050e6,
    A4=5.2002e4,
    A5=0.0,
    B1=6.38e1,
    B2=-9.95e-1,
    B3=6.981e-5,
    B4=1.728e-4,
    C1=1.17302e-9,
    C2=1.51361e-2,
    C3=7.51886e-2,
    C4=4.60643e-3,
    C5=1.35000e-2,
    C6=-1.06750e-4,
    C7=1.36600e-5,
    min_E_keV=0.5,
    max_E_keV=550.0,
    min_T_keV=0.2,
    max_T_keV=100.0,
    Q_keV=17.6e3,
    fraction_charged=0.2,
    identical_particles=False,
    Z1=1,
    Z2=1,
    high_energy_A1=-1.4714e6,
    high_energy_A2=0.0,
    high_energy_A3=0.0,
    high_energy_A4=0.0,
    high_energy_A5=0.0,
    high_energy_B1=-8.4127e-3,
    high_energy_B2=4.7983e-6,
    high_energy_B3=-1.0748e-9,
    high_energy_B4=8.5184e-14,
    high_energy_cutoff_keV=530.0,
    high_energy_max_keV=4700.0,
)

DHe3 = _ReactionParameters(
    name="3He(d,p)4He",
    reactant1="D",
    reactant2="He-3",
    products=("He-4", "p"),
    BG=68.7508,
    m1c2_keV=1875612.94257,
    m2c2_keV=2808391.60743,
    A1=5.7501e6,
    A2=2.5226e3,
    A3=4.5566e1,
    A4=0.0,
    A5=0.0,
    B1=-3.1995e-3,
    B2=-8.5530e-6,
    B3=5.9014e-8,
    B4=0.0,
    C1=5.51036e-10,
    C2=6.41918e-3,
    C3=-2.02896e-3,
    C4=-1.91080e-5,
    C5=1.35776e-4,
    C6=0.0,
    C7=0.0,
    min_E_keV=0.3,
    max_E_keV=900.0,
    min_T_keV=0.5,
    max_T_keV=190.0,
    Q_keV=18.3e3,
    fraction_charged=1.0,
    identical_particles=False,
    Z1=2,
    Z2=1,
    high_energy_A1=-8.3993e5,
    high_energy_A2=0.0,
    high_energy_A3=0.0,
    high_energy_A4=0.0,
    high_energy_A5=0.0,
    high_energy_B1=-2.6830e-3,
    high_energy_B2=1.1633e-6,
    high_energy_B3=-2.1332e-10,
    high_energy_B4=1.4250e-14,
    high_energy_cutoff_keV=900.0,
    high_energy_max_keV=4800.0,
)

DDp = _ReactionParameters(
    name="D(d,p)T",
    reactant1="D",
    reactant2="D",
    products=("T", "p"),
    BG=31.3970,
    m1c2_keV=1875612.94257,
    m2c2_keV=1875612.94257,
    A1=5.5576e4,
    A2=2.1054e2,
    A3=-3.2638e-2,
    A4=1.4987e-6,
    A5=1.8181e-10,
    B1=0.0,
    B2=0.0,
    B3=0.0,
    B4=0.0,
    C1=5.65718e-12,
    C2=3.41267e-3,
    C3=1.99167e-3,
    C4=0.0,
    C5=1.05060e-5,
    C6=0.0,
    C7=0.0,
    min_E_keV=0.5,
    max_E_keV=5000.0,
    min_T_keV=0.2,
    max_T_keV=100.0,
    Q_keV=4.03e3,
    fraction_charged=1.0,
    identical_particles=True,
    Z1=1,
    Z2=1,
)

DDn = _ReactionParameters(
    name="D(d,n)3He",
    reactant1="D",
    reactant2="D",
    products=("He-3", "n"),
    BG=31.3970,
    m1c2_keV=1875612.94257,
    m2c2_keV=1875612.94257,
    A1=5.3701e4,
    A2=3.3027e2,
    A3=-1.2706e-1,
    A4=2.9327e-5,
    A5=-2.5151e-9,
    B1=0.0,
    B2=0.0,
    B3=0.0,
    B4=0.0,
    C1=5.433603e-12,
    C2=5.85778e-3,
    C3=7.68222e-3,
    C4=0.0,
    C5=-2.96400e-6,
    C6=0.0,
    C7=0.0,
    min_E_keV=0.5,
    max_E_keV=4900.0,
    min_T_keV=0.2,
    max_T_keV=100.0,
    Q_keV=3.27e3,
    fraction_charged=0.25,
    identical_particles=True,
    Z1=1,
    Z2=1,
)

pB11 = _ReactionParameters(
    name="11B(p,4He)4He4He",
    reactant1="p",
    reactant2="B-11",
    products=("He-4", "He-4", "He-4"),
    BG=0.0,
    m1c2_keV=938272.08816,
    m2c2_keV=859526.0,
    A1=0.0,
    A2=0.0,
    A3=0.0,
    A4=0.0,
    A5=0.0,
    B1=0.0,
    B2=0.0,
    B3=0.0,
    B4=0.0,
    C1=4.4467e-14,
    C2=-5.9357e-2,
    C3=2.0165e-1,
    C4=1.0404e-3,
    C5=2.7621e-3,
    C6=-9.1653e-6,
    C7=9.8305e-7,
    min_E_keV=0.5,
    max_E_keV=3500.0,
    min_T_keV=0.1,
    max_T_keV=500.0,
    Q_keV=8.7e3,
    fraction_charged=1.0,
    identical_particles=False,
    Z1=5,
    Z2=1,
)

_REACTIONS_BY_NAME: dict[str, _ReactionParameters] = {
    "T(d,n)4He": DT,
    "D(d,p)T": DDp,
    "D(d,n)3He": DDn,
    "3He(d,p)4He": DHe3,
    "11B(p,4He)4He4He": pB11,
}

_REACTION_STRING_MAP: dict[str, _ReactionParameters] = {
    "D + T --> He-4 + n": DT,
    "D + He-3 --> He-4 + p": DHe3,
    "D + D --> T + p": DDp,
    "D + D --> He-3 + n": DDn,
    "p + B-11 --> He-4 + He-4 + He-4": pB11,
    "D + T --> alpha + n": DT,
    "D + He-3 --> alpha + p": DHe3,
    "p + B-11 --> 3 alpha": pB11,
}


def _lookup_reaction(reaction: str | tuple[str, str]):
    """Look up reaction parameters.

    Parameters
    ----------
    reaction : str or tuple of str
        Either a reaction string like ``"D + T --> He-4 + n"``,
        a short name like ``"T(d,n)4He"``,
        or a tuple of reactants like ``("D", "T")``.

    Returns
    -------
    _ReactionParameters
    """
    if isinstance(reaction, str):
        if reaction in _REACTION_STRING_MAP:
            return _REACTION_STRING_MAP[reaction]
        if reaction in _REACTIONS_BY_NAME:
            return _REACTIONS_BY_NAME[reaction]
        raise ValueError(
            f"Unknown reaction: '{reaction}'. Supported reactions are: "
            f"D + T --> He-4 + n, D + He-3 --> He-4 + p, "
            f"D + D --> T + p, D + D --> He-3 + n, p + B-11 --> He-4 + He-4 + He-4, "
            f"or short names: T(d,n)4He, D(d,p)T, D(d,n)3He, 3He(d,p)4He, "
            f"11B(p,4He)4He4He."
        )
    if isinstance(reaction, (tuple, list)) and len(reaction) == 2:
        r1, r2 = reaction
        s1, s2 = sorted((str(r1), str(r2)))
        if (s1, s2) == ("D", "D"):
            raise ValueError(
                "D + D has two reaction branches: "
                "'D(d,p)T' (neutron-producing) and "
                "'D(d,n)3He' (proton-producing). "
                "Please specify the full reaction string "
                "or use the short name."
            )
        if (s1, s2) == ("B-11", "p"):
            return pB11
        if (s1, s2) == ("D", "He-3"):
            return DHe3
        if (s1, s2) == ("D", "T"):
            return DT
        raise ValueError(
            f"Unknown reaction: {r1} + {r2}. "
            f"Supported: D + T, D + He-3, D + D, p + B-11."
        )
    raise TypeError(
        "Reaction must be a string or a tuple of two reactant strings."
    )
