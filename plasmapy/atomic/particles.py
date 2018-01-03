import typing
from astropy import units as u, constants as const
import numpy as np

# TODO: Extend docstrings in _create_Particles_dict and _create_aliases_dict
# TODO: Create an ion_symbol function
# TODO: Create a particle_symbol function
# TODO: Create a particle_mass function
# TODO: Create lepton_number, baryon_number, is_antimatter, is_lepton, etc.


def _create_Particles_dict() -> typing.Dict[str, dict]:
    """Create a dictionary of dictionaries that contains physical
    information for particles and antiparticles that are not
    elements or ions.

    The keys of the top-level dictionary are the standard
    particle symbols. The values of the top-level dictionary
    are dictionaries for each particle or antiparticle with
    strings such as 'name', 'mass', and 'spin' as the keys
    and the corresponding atomic properties as symbols."""

    leptons = ['e-', 'mu-', 'tau-', 'nu_e', 'nu_mu', 'nu_tau']
    antileptons = ['e+', 'mu+', 'tau+', 'anti_nu_e',
                   'anti_nu_mu', 'anti_nu_tau']

    baryons = ['p', 'n']
    antibaryons = ['p-', 'antineutron']

    everything = leptons + antileptons + baryons + antibaryons

    particles = leptons + baryons
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
        if lepton not in neutrinos:
            Particles[lepton]['charge'] = -1
        else:
            Particles[lepton]['charge'] = 0

    for antilepton in antileptons:
        Particles[antilepton]['class'] = 'antilepton'
        Particles[antilepton]['lepton number'] = -1
        Particles[antilepton]['baryon number'] = 0
        if antilepton not in antineutrinos:
            Particles[antilepton]['charge'] = 1
        else:
            Particles[antilepton]['charge'] = 0

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

    Particles['p']['charge'] = 1
    Particles['p-']['charge'] = -1

    for thing in ['n', 'antineutron']:
        Particles[thing]['mass'] = const.m_n
        Particles[thing]['half-life'] = 881.5 * u.s
        Particles[thing]['charge'] = 0

    for thing in everything:
        if 'half-life' not in Particles[thing].keys():
            Particles[thing]['half-life'] = np.inf * u.s

    for particle in particles:
        Particles[particle]['antimatter'] = False

    for antiparticle in antiparticles:
        Particles[antiparticle]['antimatter'] = True

    return Particles


def _create_alias_dicts(Particles: dict) -> (typing.Dict[str, str],
                                             typing.Dict[str, str]):
    """Create dictionaries for case sensitive aliases and case
    insensitive aliases of special particles and antiparticles.

    The keys of these dictionaries are the aliases, and the values
    are the corresponding standardized symbol for the particle or
    antiparticle."""

    case_sensitive_aliases = {}
    case_insensitive_aliases = {}

    for symbol in Particles.keys():
        name = Particles[symbol]['name']
        case_insensitive_aliases[name.lower()] = symbol

    case_sensitive_aliases_for_a_symbol = [
        (['beta-'], 'e-'),
        (['beta+'], 'e+'),
        (['p+'], 'p'),
        (['n-1'], 'n'),
        (['H-2'], 'D'),
        (['H-2+', 'H-2 1+', 'H-2 +1'], 'D 1+'),
        (['H-3+', 'H-3 1+', 'H-3 +1'], 'T 1+'),
    ]

    case_insensitive_aliases_for_a_symbol = [
        (['antielectron'], 'e+'),
        (['muon-'], 'mu-'),
        (['muon+'], 'mu+'),
        (['tau particle'], 'tau-'),
        (['protium'], 'H-1'),
        (['protium+', 'protium 1+', 'protium +1'], 'p'),
        (['deuterium', 'hydrogen-2'], 'D'),
        (['deuteron', 'deuterium+', 'deuterium 1+', 'deuterium +1', 'D+'],
         'D 1+'),
        (['tritium', 'hydrogen-3'], 'T'),
        (['triton', 'tritium+', 'tritium 1+', 'tritium +1'], 'T 1+'),
        (['alpha'], 'He-4 2+'),
    ]

    for aliases, symbol in case_sensitive_aliases_for_a_symbol:
        for alias in aliases:
            case_sensitive_aliases[alias] = symbol

    for aliases, symbol in case_insensitive_aliases_for_a_symbol:
        for alias in aliases:
            case_insensitive_aliases[alias.lower()] = symbol

    alias_keys = list(case_insensitive_aliases.keys())

    for alias in alias_keys:
        if 'anti' in alias and 'anti-' not in alias:
            symbol = case_insensitive_aliases[alias].lower()
            new_alias = alias.replace('anti', 'anti-')
            case_insensitive_aliases[new_alias] = symbol

    return case_sensitive_aliases, case_insensitive_aliases


_Particles = _create_Particles_dict()

_case_sensitive_aliases, _case_insensitive_aliases = \
    _create_alias_dicts(_Particles)


def _get_standard_symbol(alias: typing.Union[str, int]) -> str:
    """Returns the standard symbol for a particle or antiparticle
    when the argument is a valid alias.  If the argument is not a
    valid alias, then this function returns the original argument
    (which will usually be a string but may be an int representing
    atomic number)."""

    if not isinstance(alias, str):
        return alias

    if alias in _case_sensitive_aliases.keys():
        return _case_sensitive_aliases[alias]
    elif alias.lower() in _case_insensitive_aliases.keys():
        return _case_insensitive_aliases[alias.lower()]
    else:
        return alias


if __name__ == "__main__":  # coveralls: ignore
    from pprint import pprint
    print("Case insensitive aliases:")
    pprint(_case_insensitive_aliases)
    print(20*"=")
    print("Case sensitive aliases:")
    pprint(_case_sensitive_aliases)
    print(20*"=")
    print("Particles:")
    pprint(_Particles)
