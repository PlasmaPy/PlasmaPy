r"""Functions that deal with string representations of atomic symbols
and numbers."""

import re
import warnings
from typing import (Optional, Any, Tuple)
from .parsing import _dealias_particle_aliases
from .particle_class import Particle
from .particle_input import particle_input
from ..utils import (AtomicWarning, InvalidParticleError)

# The @particle_input decorator takes the inputs for a function or
# method and passes through the corresponding instance of the Particle
# class via the argument annotated with Particle.  If the argument is
# named element, isotope, or ion; then this decorator will raise an
# InvalidElementError, InvalidIsotopeError, or InvalidIonError,
# respectively.  The Particle class constructor will raise an
# InvalidParticleError if the input does not correspond to a valid
# particle.


@particle_input
def atomic_symbol(element: Particle) -> str:
    r"""Returns the atomic symbol.

    Parameters
    ----------

    element: `str`, `int`, or `~plasmapy.atomic.Particle`
        A string representing an element, isotope, or ion; or an
        integer or string representing an atomic number.

    Returns
    -------

    symbol: `str`
        The atomic symbol of the element, isotope, or ion.

    Raises
    ------

    `~plasmapy.utils.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `TypeError`
        If the argument is not a string or integer.

    See also
    --------

    `~plasmapy.atomic.element_name` : returns the name of an element.

    `~plasmapy.atomic. isotope_symbol` : returns the isotope symbol
        instead of the atomic symbol.

    `~plasmapy.atomic.ion_symbol` : returns the ion symbol instead of
        the atomic symbol.

    `~plasmapy.atomic.particle_symbol` : returns the symbol of any valid
        particle.

    Notes
    -----

    This function returns the symbol of the element rather than the
    symbol of an isotope.  For example, `'deuterium'`, `'T'`, or
    `'hydrogen-2'` will yield `'H'`; `'alpha'` will yield `'He'`; and
    `'iron-56'` or `'Fe-56'` will yield `'Fe'`.

    This function is case insensitive when there is no ambiguity
    associated with case.  However, this function will return `'H'` for
    hydrogen for lower case `'p'` but capital `'P'` if the argument is
    `'P'` for phosphorus.  This function will return `'N'` for nitrogen
    if the argument is capital `'N'`, but will not accept lower case
    `'n'` for neutrons.

    Examples
    --------

    >>> atomic_symbol('helium')
    'He'
    >>> atomic_symbol(42)
    'Mo'
    >>> atomic_symbol('D')
    'H'
    >>> atomic_symbol('C-13')
    'C'
    >>> atomic_symbol('alpha')
    'He'
    >>> atomic_symbol('79')
    'Au'
    >>> atomic_symbol('N')  # Nitrogen
    'N'
    >>> atomic_symbol('P'), atomic_symbol('p')  # Phosphorus, proton
    ('P', 'H')

    """
    return element.element


@particle_input
def isotope_symbol(isotope: Particle, mass_numb: int = None) -> str:
    r"""Returns the symbol representing an isotope.

    Parameters
    ----------

    isotope: `str`, `int`, or `~plasmapy.atomic.Particle`
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    mass_numb: `int` or `str` (optional)
        An integer or string representing the mass number of the
        isotope.

    Returns
    -------

    symbol: `str`
        The isotopic symbol. The result will generally be returned as
        something like `'He-4'` or `'Au-197'`, but will return `'D'` for
        deuterium and `'T'` for tritium.

    Raises
    ------

    `~plasmapy.utils.InvalidIsotopeError`
        If the argument is a valid particle but not a valid isotope.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle
        or contradictory information is provided.

    `TypeError`
        If the argument is not a `str`, `int`, or `Particle`.

    Warns
    -----

    `~plasmapy.utils.AtomicWarning`
        If redundant isotope information is provided.

    See also
    --------

    `~plasmapy.atomic.atomic_symbol` : returns the atomic symbol instead
        of the isotope symbol

    `~plasmapy.atomic.ion_symbol` : returns the ion symbol instead of
        the isotope symbol.

    `~plasmapy.atomic.particle_symbol` : returns the symbol of any valid
        particle.


    Notes
    -----

    This function returns the symbol of the element rather than the
    symbol of an isotope.  For example, `'deuterium'`, `'T'`, or
    `'hydrogen-2'` will yield `'H'`; `'alpha'` will yield `'He'`; and
    `'iron-56'` or `'Fe-56'` will yield `'Fe'`.

    Examples
    --------
    >>> isotope_symbol('He', 4)
    'He-4'
    >>> isotope_symbol(79, 197)
    'Au-197'
    >>> isotope_symbol('hydrogen-2')
    'D'
    >>> isotope_symbol('carbon-13')
    'C-13'
    >>> isotope_symbol('alpha')
    'He-4'

    """
    return isotope.isotope


@particle_input
def ion_symbol(ion: Particle, mass_numb: int = None, Z: int = None) -> str:
    r"""Returns the symbol representing an ion.

    Parameters
    ----------

    ion: `int`, `str`, or `~plasmapy.atomic.Particle`
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    mass_numb: `int` or `str` (optional)
        An integer or string representing the mass number of the ion.

    Z: `int` or `str` (optional)
        An integer or string representing the integer charge of the ion.

    Returns
    -------

    symbol: `str`
        The ion symbol. The result will generally be returned as
        something like `'He-4 2+'`, `'D 1+'`, or `'p+'`.

    Raises
    ------

    `~plasmapy.utils.InvalidIonError`
        If the arguments correspond to a valid particle but not a valid
        ion.

    `~plasmapy.utils.InvalidParticleError`
        If arguments do not correspond to a valid particle or
        contradictory information is provided.

    `TypeError`
        If ion is not a string, integer, or Particle; or if either of
        `mass_numb` or `Z` is not an integer or a string representing an
        integer.

    Warns
    -----

    `~plasmapy.utils.AtomicWarning`
        If redundant mass number or charge information is provided.

    See also
    --------

    `~plasmapy.atomic.atomic_symbol` : returns the atomic symbol instead
        of the ion symbol

    `~plasmapy.atomic.isotope_symbol` : returns the isotope symbol
        instead of the ion symbol

    `~plasmapy.atomic.particle_symbol` : returns the symbol of any valid
        particle

    Notes
    -----

    This function returns the symbol of the element rather than the
    symbol of an isotope.  For example, `'deuterium'`, `'T'`, or
    `'hydrogen-2'` will yield `'H'`; `'alpha'` will yield `'He'`; and
    `'iron-56'` or `'Fe-56'` will yield `'Fe'`.

    Examples
    --------
    >>> ion_symbol('alpha')
    'He-4 2+'
    >>> ion_symbol(79, mass_numb=197, Z=12)
    'Au-197 12+'
    >>> ion_symbol('proton')
    'p+'
    >>> ion_symbol('D', Z=1)
    'D 1+'
    """
    return ion.ion


@particle_input
def particle_symbol(particle: Particle, mass_numb: int = None,
                    Z: int = None) -> str:
    r"""Returns the symbol of a particle.

    Parameters
    ----------

    particle: `int`, `str`, or `~plasmapy.atomic.Particle`
        A string representing a particle, element, isotope, or ion or an
        integer representing an atomic number

    mass_numb: `int` or `str` (optional)
        An integer or string representing the mass number of an isotope.

    Z: `int` or `str` (optional)
        An integer or string representing the integer charge of an ion.

    Returns
    -------

    symbol: `str`
        The particle symbol, containing charge and mass number
        information when available. The result will generally be
        returned as something like `'e-'`, `'Fe'`, `'He-4 2+'`,
        `'D'`, `'n'`, `'mu-'`, or `'p+'`.

    Raises
    ------

    `~plasmapy.utils.InvalidParticleError`
        If arguments do not correspond to a valid particle or
        contradictory information is provided.

    `TypeError`
        If ion is not a string, integer, or Particle; or if either of
        `mass_numb` or `Z` is not an integer or a string representing an
        integer.

    Warns
    -----

    `~plasmapy.utils.AtomicWarning`
        If redundant mass number or charge information is provided.

    See also
    --------

    `~plasmapy.atomic.atomic_symbol` : returns the atomic symbol instead

    `~plasmapy.atomic.isotope_symbol` : returns the isotope symbol
        instead

    `~plasmapy.atomic.ion_symbol` : returns the ion symbol instead

    Examples
    --------
    >>> particle_symbol('electron')
    'e-'
    >>> particle_symbol('proton')
    'p+'
    >>> particle_symbol('alpha')
    'He-4 2+'
    >>> particle_symbol('H-1', Z=-1)
    'H-1 1-'
    """
    return particle.particle


@particle_input
def element_name(element: Particle) -> str:
    r"""Returns the name of an element.

    Parameters
    ----------

    argument : `str`, `int`, or `~plasmapy.atomic.Particle`
        A string representing an element, isotope, or ion or an
        integer representing an atomic number

    Returns
    -------

    name : `str`
        The name of the element.

    Raises
    ------

    `~plasmapy.utils.InvalidElementError`
        If the argument is a valid particle but not a valid element.

    `~plasmapy.utils.InvalidParticleError`
        If the argument does not correspond to a valid particle.

    `TypeError`
        If the argument is not a string or integer.

    See also
    --------

    `~plasmapy.atomic.atomic_symbol` : returns the atomic symbol

    `~plasmapy.atomic.isotope_symbol` : returns the isotope symbol

    `~plasmapy.atomic.particle_symbol` : returns the symbol of any valid
        particle

    Examples
    --------

    >>> element_name("H")
    'hydrogen'
    >>> element_name("T")
    'hydrogen'
    >>> element_name("alpha")
    'helium'
    >>> element_name(42)
    'molybdenum'
    >>> element_name("C-12")
    'carbon'

    """
    return element.element_name
