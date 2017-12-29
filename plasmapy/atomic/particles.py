from astropy import units as u, constants as const
import numpy as np

Particles = {

    "e-": {
        "name": "electron",
        "mass": const.m_e,
        "charge": -1,
        "spin": 1/2,
        "half-life": np.inf * u.s,
        "antimatter": False,
        "generation": 1,
    },

    "e+": {
        "name": "positron",
        "mass": const.m_e,
        "charge": 1,
        "spin": 1/2,
        "half-life": np.inf * u.s,
        "antimatter": True,
        "generation": 1,
    },

    "mu-": {
        "name": "muon",
        "mass": 1.883_531_594e-28 * u.kg,
        "charge": -1,
        "spin": 1/2,
        "half-life": 2.1969811e-6 * u.s,
        "antimatter": False,
        "generation": 2,
    },

    "mu+": {
        "name": "antimuon",
        "mass": 1.883_531_594e-28 * u.kg,
        "charge": 1,
        "spin": 1/2,
        "half-life": 2.1969811e-6 * u.s,
        "antimatter": True,
        "generation": 2,
    },

    "tau-": {
        "name": "tau",
        "mass": 3.167_47e-27 * u.kg,
        "charge": -1,
        "spin": 1/2,
        "half-life": 2.906e-13 * u.s,
        "antimatter": False,
        "generation": 3,
     },

    "tau+": {
        "name": "antitau",
        "mass": 3.167_47e-27 * u.kg,
        "charge": 1,
        "spin": 1/2,
        "half-life": 2.906e-13 * u.s,
        "antimatter": True,
        "generation": 3,
    },

    "nu_e": {
        "name": "electron neutrino",
        "charge": 0,
        "spin": 1/2,
        "half-life": np.inf * u.s,
        "antimatter": False,
        "generation": 1,
    },

    "anti_nu_e": {
        "name": "electron antineutrino",
        "charge": 0,
        "spin": 1/2,
        "half-life": np.inf * u.s,
        "antimatter": True
        "generation": 1
    },

    "nu_mu": {
        "name": "muon neutrino",
        "spin": 1/2,
        "half-life": np.inf * u.s,
        "antimatter": False,
        "generation": 2,
    },

    "anti_nu_mu": {
        "name": "muon antineutrino",
        "spin": 1/2,
        "half-life": np.inf * u.s,
        "antimatter": True,
        "generation": 2,
    },

    "nu_tau": {
        "name": "tau neutrino",
        "spin": 1/2,
        "half-life": np.inf * u.s,
        "antimatter": False,
        "generation": 3,
    },

    "anti_nu_tau": {
        "name": "tau antineutrino",
        "spin": 1/2,
        "half-life": np.inf * u.s,
        "antimatter": True,
        "generation": 3,
    },

}
