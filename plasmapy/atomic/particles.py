import typing
from astropy import units as u, constants as const
import numpy as np

_leptons = ['e-', 'mu-', 'tau-', 'nu_e', 'nu_mu', 'nu_tau']
_antileptons = ['e+', 'mu+', 'tau+', 'anti_nu_e',
                'anti_nu_mu', 'anti_nu_tau']

_baryons = ['p', 'n']
_antibaryons = ['p-', 'antineutron']

_everything = _leptons + _antileptons + _baryons + _antibaryons

_particles = _leptons + _baryons
_antiparticles = _antileptons + _antibaryons

_fermions = _leptons + _antileptons + _baryons + _antibaryons
_bosons = []

_neutrinos = [lepton for lepton in _leptons if 'nu' in lepton]
_antineutrinos = [antilepton for antilepton in _antileptons
                  if 'nu' in antilepton]

_special_particles = [
    'n', 'antineutron', 'p-',
    'e-', 'e+', 'nu_e', 'anti_nu_e',
    'mu-', 'mu+', 'nu_mu', 'anti_nu_mu',
    'tau-', 'tau+', 'nu_tau', 'anti_nu_tau']


def _create_Particles_dict() -> typing.Dict[str, dict]:
    """Create a dictionary of dictionaries that contains physical
    information for particles and antiparticles that are not
    elements or ions.

    The keys of the top-level dictionary are the standard
    particle symbols. The values of the top-level dictionary
    are dictionaries for each particle or antiparticle with
    strings such as 'name', 'mass', and 'spin' as the keys
    and the corresponding atomic properties as symbols."""

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

    for thing in _everything:
        Particles[thing] = dict()

    for symbol, name in symbols_and_names:
        Particles[symbol]['name'] = name

    for fermion in _fermions:
        Particles[fermion]['spin'] = 0.5

    for boson in _bosons:
        Particles[boson]['spin'] = 0

    for lepton in _leptons:
        Particles[lepton]['class'] = 'lepton'
        Particles[lepton]['lepton number'] = 1
        Particles[lepton]['baryon number'] = 0
        if lepton not in _neutrinos:
            Particles[lepton]['charge'] = -1
        else:
            Particles[lepton]['charge'] = 0

    for antilepton in _antileptons:
        Particles[antilepton]['class'] = 'antilepton'
        Particles[antilepton]['lepton number'] = -1
        Particles[antilepton]['baryon number'] = 0
        if antilepton not in _antineutrinos:
            Particles[antilepton]['charge'] = 1
        else:
            Particles[antilepton]['charge'] = 0

    for baryon in _baryons:
        Particles[baryon]['class'] = 'baryon'
        Particles[baryon]['lepton number'] = 0
        Particles[baryon]['baryon number'] = 1

    for antibaryon in _antibaryons:
        Particles[antibaryon]['class'] = 'antibaryon'
        Particles[antibaryon]['lepton number'] = 0
        Particles[antibaryon]['baryon number'] = -1

    for thing in _leptons + _antileptons:
        if 'e' in thing:
            Particles[thing]['generation'] = 1
        elif 'mu' in thing:
            Particles[thing]['generation'] = 2
        elif 'tau' in thing:
            Particles[thing]['generation'] = 3

    for thing in _leptons + _antileptons:
        if 'nu' not in thing:
            if 'e' in thing:
                Particles[thing]['mass'] = const.m_e
            elif 'mu' in thing:
                Particles[thing]['mass'] = 1.883_531_594e-28 * u.kg
                Particles[thing]['half-life'] = 2.1969811e-6 * u.s
            elif 'tau' in thing:
                Particles[thing]['mass'] = 3.167_47e-27 * u.kg
                Particles[thing]['half-life'] = 2.906e-13 * u.s

    for thing in _neutrinos + _antineutrinos:
        Particles[thing]['mass'] = None


    for thing in ['p', 'p-']:
        Particles[thing]['mass'] = const.m_p

    Particles['p']['charge'] = 1
    Particles['p-']['charge'] = -1

    for thing in ['n', 'antineutron']:
        Particles[thing]['mass'] = const.m_n
        Particles[thing]['half-life'] = 881.5 * u.s
        Particles[thing]['charge'] = 0

    for thing in _everything:
        if 'half-life' not in Particles[thing].keys():
            Particles[thing]['half-life'] = np.inf * u.s

    for particle in _particles:
        Particles[particle]['antimatter'] = False

    for antiparticle in _antiparticles:
        Particles[antiparticle]['antimatter'] = True

    for thing in _everything:
        try:
            Particles[thing]['half-life']
        except KeyError:
            Particles[thing]['half-life'] = np.inf * u.s

    return Particles


_Particles = _create_Particles_dict()


if __name__ == "__main__":  # coveralls: ignore
    from pprint import pprint
    print("Particles:")
    pprint(_Particles)
