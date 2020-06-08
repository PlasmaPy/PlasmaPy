"""
Classes, sets, and dictionaries to store data and taxonomy
information for special particles.
"""

from typing import Dict, Set

import numpy as np
from astropy import constants as const
from astropy import units as u

from plasmapy.particles.elements import _PeriodicTable


class _ParticleZooClass:
    """
    Create an object with taxonomy information for special particles.

    The `~plasmapy.particles.special_particles._ParticleZooClass._taxonomy_dict`
    attribute contains the name of each classification (e.g.,
    ``'lepton'``, ``'baryon'``, ``'matter'``, etc.) as the keys and a
    set of particle symbol strings of the particles belonging to that
    classification.

    The attributes of this class provide sets of strings representing
    particles in the corresponding category.

    Examples
    --------
    >>> ParticleZoo = _ParticleZooClass()
    >>> 'e-' in ParticleZoo.leptons
    True
    >>> 'nu_e' in ParticleZoo.antineutrinos
    False
    >>> 'mu+' in ParticleZoo.antiparticles
    True

    """

    def __init__(self):

        leptons = {"e-", "mu-", "tau-", "nu_e", "nu_mu", "nu_tau"}
        antileptons = {"e+", "mu+", "tau+", "anti_nu_e", "anti_nu_mu", "anti_nu_tau"}
        baryons = {"p+", "n"}
        antibaryons = {"p-", "antineutron"}

        particles = leptons | baryons
        antiparticles = antileptons | antibaryons
        fermions = leptons | antileptons | baryons | antibaryons

        bosons = set()

        neutrinos = {lepton for lepton in leptons if "nu" in lepton}

        antineutrinos = {antilepton for antilepton in antileptons if "nu" in antilepton}

        self._taxonomy_dict = {
            "lepton": leptons,
            "antilepton": antileptons,
            "baryon": baryons,
            "antibaryon": antibaryons,
            "fermion": fermions,
            "boson": bosons,
            "neutrino": neutrinos,
            "antineutrinos": antineutrinos,
            "matter": particles,
            "antimatter": antiparticles,
        }

    @property
    def leptons(self) -> Set[str]:
        """Return all leptons, excluding antileptons."""
        return self._taxonomy_dict["lepton"]

    @property
    def antileptons(self) -> Set[str]:
        """Return all antileptons."""
        return self._taxonomy_dict["antilepton"]

    @property
    def baryons(self) -> Set[str]:
        """Return all baryons, excluding antibaryons."""
        return self._taxonomy_dict["baryon"]

    @property
    def antibaryons(self) -> Set[str]:
        """Return all antibaryons."""
        return self._taxonomy_dict["antibaryon"]

    @property
    def fermions(self) -> Set[str]:
        """Return all fermions."""
        return self._taxonomy_dict["fermion"]

    @property
    def bosons(self) -> Set[str]:
        """Return all bosons."""
        return self._taxonomy_dict["boson"]

    @property
    def neutrinos(self) -> Set[str]:
        """Return all neutrinos."""
        return self._taxonomy_dict["neutrino"]

    @property
    def antineutrinos(self) -> Set[str]:
        """Return all antineutrinos."""
        return self._taxonomy_dict["antineutrinos"]

    @property
    def particles(self) -> Set[str]:
        """Return all particles, excluding antiparticles."""
        return self._taxonomy_dict["matter"]

    @property
    def antiparticles(self) -> Set[str]:
        """Return all antiparticles."""
        return self._taxonomy_dict["antimatter"]

    @property
    def everything(self) -> Set[str]:
        """Return all particles and antiparticles."""
        return self._taxonomy_dict["matter"] | self._taxonomy_dict["antimatter"]


ParticleZoo = _ParticleZooClass()


def _create_Particles_dict() -> Dict[str, dict]:
    """
    Create a dictionary of dictionaries that contains physical
    information for particles and antiparticles that are not elements or
    ions.

    The keys of the top-level dictionary are the standard particle
    symbols. The values of the top-level dictionary are dictionaries for
    each particle or antiparticle with strings such as ``'name'``,
    ``'mass'``, and ``'spin'`` as the keys and the corresponding atomic
    properties as symbols.
    """

    symbols_and_names = [
        ("e-", "electron"),
        ("e+", "positron"),
        ("mu-", "muon"),
        ("mu+", "antimuon"),
        ("tau-", "tau"),
        ("tau+", "antitau"),
        ("nu_e", "electron neutrino"),
        ("anti_nu_e", "electron antineutrino"),
        ("nu_mu", "muon neutrino"),
        ("anti_nu_mu", "muon antineutrino"),
        ("nu_tau", "tau neutrino"),
        ("anti_nu_tau", "tau antineutrino"),
        ("p+", "proton"),
        ("p-", "antiproton"),
        ("n", "neutron"),
        ("antineutron", "antineutron"),
    ]

    Particles = dict()

    for particle in ParticleZoo.everything:
        Particles[particle] = dict()

    for symbol, name in symbols_and_names:
        Particles[symbol]["name"] = name

    for fermion in ParticleZoo.fermions:
        Particles[fermion]["spin"] = 0.5

    for boson in ParticleZoo.bosons:  # coverage: ignore
        Particles[boson]["spin"] = 0

    for lepton in ParticleZoo.leptons:
        Particles[lepton]["class"] = "lepton"
        Particles[lepton]["lepton number"] = 1
        Particles[lepton]["baryon number"] = 0
        if lepton not in ParticleZoo.neutrinos:
            Particles[lepton]["integer charge"] = -1
        else:
            Particles[lepton]["integer charge"] = 0

    for antilepton in ParticleZoo.antileptons:
        Particles[antilepton]["class"] = "antilepton"
        Particles[antilepton]["lepton number"] = -1
        Particles[antilepton]["baryon number"] = 0
        if antilepton not in ParticleZoo.antineutrinos:
            Particles[antilepton]["integer charge"] = 1
        else:
            Particles[antilepton]["integer charge"] = 0

    for baryon in ParticleZoo.baryons:
        Particles[baryon]["class"] = "baryon"
        Particles[baryon]["lepton number"] = 0
        Particles[baryon]["baryon number"] = 1

    for antibaryon in ParticleZoo.antibaryons:
        Particles[antibaryon]["class"] = "antibaryon"
        Particles[antibaryon]["lepton number"] = 0
        Particles[antibaryon]["baryon number"] = -1

    for particle in ParticleZoo.leptons | ParticleZoo.antileptons:
        if "e" in particle:
            Particles[particle]["generation"] = 1
        elif "mu" in particle:
            Particles[particle]["generation"] = 2
        elif "tau" in particle:
            Particles[particle]["generation"] = 3

    for particle in ParticleZoo.leptons | ParticleZoo.antileptons:
        if "nu" not in particle:
            if "e" in particle:
                Particles[particle]["mass"] = const.m_e
            elif "mu" in particle:
                Particles[particle]["mass"] = 1.883_531_594e-28 * u.kg
                Particles[particle]["half-life"] = 2.1969811e-6 * u.s
            elif "tau" in particle:
                Particles[particle]["mass"] = 3.167_47e-27 * u.kg
                Particles[particle]["half-life"] = 2.906e-13 * u.s

    # Setting the neutrino mass to None reminds us that, while neutrinos
    # are not massless, we only have upper limits on what the neutrino
    # mass actually is.

    for particle in ParticleZoo.neutrinos | ParticleZoo.antineutrinos:
        Particles[particle]["mass"] = None

    special_attributes = {
        "p+": {
            "element": "H",
            "atomic number": 1,
            "mass number": 1,
            "element name": "hydrogen",
            "isotope": "H-1",
            "ion": "p+",
            "mass": const.m_p,
            "integer charge": 1,
            "periodic table": _PeriodicTable(
                group=1, period=1, block="s", category="nonmetal"
            ),
        },
        "p-": {"mass": const.m_p, "integer charge": -1},
        "n": {"mass": const.m_n, "half-life": 881.5 * u.s, "integer charge": 0},
        "antineutron": {
            "mass": const.m_n,
            "half-life": 881.5 * u.s,
            "integer charge": 0,
        },
    }

    for particle in special_attributes.keys():
        Particles[particle] = {**special_attributes[particle], **Particles[particle]}

    for particle in ParticleZoo.everything:
        if "half-life" not in Particles[particle].keys():
            Particles[particle]["half-life"] = np.inf * u.s

    for particle in ParticleZoo.particles:
        Particles[particle]["antimatter"] = False

    for antiparticle in ParticleZoo.antiparticles:
        Particles[antiparticle]["antimatter"] = True

    return Particles


_Particles = _create_Particles_dict()

_special_ion_masses = {
    "p+": const.m_p,
    "D 1+": 3.343583719e-27 * u.kg,
    "T 1+": 5.007356665e-27 * u.kg,
}

_antiparticles = {
    "p+": "p-",
    "n": "antineutron",
    "e-": "e+",
    "mu-": "mu+",
    "tau-": "tau+",
    "nu_e": "anti_nu_e",
    "nu_mu": "anti_nu_mu",
    "nu_tau": "anti_nu_tau",
    "p-": "p+",
    "antineutron": "n",
    "e+": "e-",
    "mu+": "mu-",
    "tau+": "tau-",
    "anti_nu_e": "nu_e",
    "anti_nu_mu": "nu_mu",
    "anti_nu_tau": "nu_tau",
}

if __name__ == "__main__":  # coverage: ignore
    from pprint import pprint

    print("Particles:")
    pprint(_Particles)
