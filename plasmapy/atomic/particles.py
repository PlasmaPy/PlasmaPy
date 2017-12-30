from astropy import units as u, constants as const
import numpy as np

# TODO: Add tests for _Particles and _Aliases
# TODO: Add type hints for a dictionary of dictionaries of whatever
# TODO: Extend docstrings in _create_Particles_dict and _create_aliases_dict
# TODO: Create function to take an alias and return the standardized symbol
# TODO: Create an ion_symbol function
# TODO: Create a particle_symbol function
# TODO: Create a particle_mass function
# TODO: Create lepton_number, baryon_number, is_antimatter, is_lepton, etc.


def _create_Particles_dict():
    """Create a dictionary that will contain physical information for
    all of the particles and antiparticles that we may care about."""

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


def _create_aliases_dict(Particles):
    """Create a dictionary to contain all of the horrible, horrible
     aliases used for different particles, antiparticles, isotopes,
     and ions."""

    aliases = {}

    for symbol in Particles.keys():
        aliases[symbol] = {'case sensitive': [], 'case insensitive': []}
        name = Particles[symbol]['name']
        aliases[symbol]['case insensitive'].append(name)
        if ' ' in name:
            aliases[symbol]['case insensitive'].append(name.replace(' ', '_'))

    aliases['e-']['case sensitive'].append('beta-')

    aliases['e+']['case sensitive'].append('beta+')
    aliases['e+']['case insensitive'].append('antielectron')

    aliases['mu-']['case insensitive'].append('muon-')
    aliases['mu+']['case insensitive'].append('muon+')

    aliases['p']['case sensitive'].append('p+')

    for symbol in ['D', 'D 1+', 'T', 'T 1+', 'He-4 2+']:
        aliases[symbol] = {'case sensitive': [], 'case insensitive': []}

    for symbol, mass_numb in [('D', 2), ('T', 3)]:
        aliases[symbol]['case sensitive'].append(f"H-{mass_numb}")
        aliases[symbol]['case insensitive'].append(f"hydrogen-{mass_numb}")

        ion = f"{symbol} 1+"
        aliases[ion]['case sensitive'].extend(
            [f"{symbol}+", f"{symbol} +1", f"H-{mass_numb}+",
             f"H-{mass_numb} 1+", f"H-{mass_numb} +1"])

        aliases[ion]['case insensitive'].extend(
            [f"hydrogen-{mass_numb}+", f"hydrogen-{mass_numb} 1+",
             f"hydrogen-{mass_numb} +1"])

    aliases['D']['case insensitive'].append('deuterium')
    aliases['T']['case insensitive'].append('tritium')

    aliases['D 1+']['case insensitive'].extend(
        ['deuteron', 'deuterium+', 'deuterium 1+', 'deuterium +1'])
    aliases['T 1+']['case insensitive'].extend(
        ['triton', 'tritium+', 'tritium 1+', 'tritium +1'])

    aliases['He-4 2+']['case insensitive'].extend(
        ['alpha', 'helium-4++', 'helium-4 2+', 'helium-4 +2'])
    aliases['He-4 2+']['case sensitive'].extend(
        ['He-4++', 'He-4 +2'])

    aliases['n']['case sensitive'].append('n-1')
    aliases['n']['case insensitive'].append('n0')

    for symbol in aliases.keys():
        for alias in aliases[symbol]['case insensitive']:
            # The second clause in the next line is to avoid an infinite loop
            if 'anti' in alias and 'anti-' not in alias:
                new_item = alias.replace('anti', 'anti-')
                aliases[symbol]['case insensitive'].append(new_item)

    return aliases


_Particles = _create_Particles_dict()
_aliases = _create_aliases_dict(_Particles)
