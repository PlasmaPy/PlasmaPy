from astropy import units as u, constants as const
import numpy as np

"""Next we go through the process of creating the Particles dictionary."""


def _create_Particles_dict():

    leptons = ['e-', 'mu-', 'tau-', 'nu_e', 'nu_mu', 'nu_tau']
    antileptons = ['e+', 'mu+', 'tau+', 'anti_nu_e',
                   'anti_nu_mu', 'anti_nu_tau']

    baryons = ['p', 'n']
    antibaryons = ['p-', 'antineutron']

    everything = leptons + antileptons + baryons + antibaryons

    particles = leptons + antileptons + baryons + antibaryons
    antiparticles = antileptons + antibaryons

    fermions = leptons + antileptons + baryons + antibaryons
    bosons = []

    neutrinos = [lepton for lepton in leptons if 'nu' in lepton]
    antineutrinos = [antilepton for antilepton in antileptons
                     if 'nu' in antilepton]

    symbols_and_names = [
        ('e-', 'electron'),
        ('e+', 'positron'),
        ('mu-', 'muon'),
        ('mu+', 'antimuon'),
        ('tau-', 'tau'),
        ('tau+', 'antitau'),
        ('nu_e', 'electron neutrino'),
        ('anti_nu_e', 'electron antineutrino'),
        ('nu_mu', 'muon neutrino'),
        ('anti_nu_mu', 'muon antineutrino'),
        ('nu_tau', 'tau neutrino'),
        ('anti_nu_tau', 'tau antineutrino'),
        ('p', 'proton'),
        ('p-', 'antiproton'),
        ('n', 'neutron'),
        ('antineutron', 'antineutron'),
    ]

    Particles = dict()

    for thing in everything:
        Particles[thing] = dict()

    for symbol, name in symbols_and_names:
        Particles[symbol]['name'] = name

    for fermion in fermions:
        Particles[fermion]['spin'] = 0.5

    for boson in bosons:
        Particles[boson]['spin'] = 0

    for lepton in leptons:
        Particles[lepton]['class'] = 'lepton'
        Particles[lepton]['lepton number'] = 1
        Particles[lepton]['baryon number'] = 0

    for antilepton in antileptons:
        Particles[antilepton]['class'] = 'antilepton'
        Particles[antilepton]['lepton number'] = -1
        Particles[antilepton]['baryon number'] = 0

    for baryon in baryons:
        Particles[baryon]['class'] = 'baryon'
        Particles[baryon]['lepton number'] = 0
        Particles[baryon]['baryon number'] = 1

    for antibaryon in antibaryons:
        Particles[antibaryon]['class'] = 'antibaryon'
        Particles[antibaryon]['lepton number'] = 0
        Particles[antibaryon]['baryon number'] = -1

    for thing in leptons + antileptons:
        if 'e' in thing:
            Particles[thing]['generation'] = 1
        elif 'mu' in thing:
            Particles[thing]['generation'] = 2
        elif 'tau' in thing:
            Particles[thing]['generation'] = 3

    for thing in leptons + antileptons:
        if 'nu' not in thing:
            if 'e' in thing:
                Particles[thing]['mass'] = const.m_e
            elif 'mu' in thing:
                Particles[thing]['mass'] = 1.883_531_594e-28 * u.kg
                Particles[thing]['half-life'] = 2.1969811e-6 * u.s
            elif 'tau' in thing:
                Particles[thing]['mass'] = 3.167_47e-27 * u.kg
                Particles[thing]['half-life'] = 2.906e-13 * u.s

    for thing in ['p', 'p-']:
        Particles[thing]['mass'] = const.m_p

    for thing in ['n', 'antineutron']:
        Particles[thing]['mass'] = const.m_n
        Particles[thing]['half-life'] = 881.5 * u.s

    for thing in everything:
        if 'half-life' not in Particles[thing].keys():
            Particles[thing]['half-life'] = np.inf * u.s

    return Particles

_Particles = _create_Particles_dict()
