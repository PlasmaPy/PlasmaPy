"""
This module contains functionality for calculating the timescales
for a range of configurations.
"""

__all__ = ["Hellinger", "Hellinger_2009", "Hellinger_2010", "Hellinger_2016"]

import astropy.units
import astropy.units as u
import numpy as np

from astropy.constants.si import eps0
from math import pi as pi
from mpmath import hyper2d
from scipy.stats import gausshyper as gh

from plasmapy.formulary.collisions import coulomb
from plasmapy.particles import Particle, ParticleList
from plasmapy.utils.decorators import validate_quantities


class validate:

    # Validate n_i argument
    def n_i(
        n_i: u.m**-3,
    ):
        if not isinstance(n_i, astropy.units.Quantity):
            raise TypeError(
                "Argument 'n_i' must be an astropy.units.Quantity, "
                f"instead got type of {type(n_i)}."
            )
        n_i = n_i.squeeze()
        if n_i.ndim != 0:
            raise ValueError(
                "Argument 'n_i' must be single value and not an array of"
                f" shape {n_i.shape}."
            )
        elif not isinstance(n_i.value, (int, float)):
            raise TypeError(
                "Argument 'n_i' must be an integer or float, received "
                f"{n_i} of type {type(n_i)}."
            )
        elif not n_i.value > 0:
            raise ValueError(
                f"Argument 'n_i' must be an positive argument, received "
                f"{n_i} of type {type(n_i)}."
            )

        return n_i

    # Validate ions argument
    def ions(
        ions: (Particle, Particle),
    ):
        if not isinstance(ions, (list, tuple, ParticleList)):
            ions = [ions]
        ions = ParticleList(ions)

        if not all(
            failed := [ion.is_ion and abs(ion.charge_number) > 0 for ion in ions]
        ):
            raise ValueError(
                "Particle(s) passed to 'ions' must be a charged"
                " ion. The following particle(s) is(are) not allowed "
                f"{[ion for ion, fail in zip(ions, failed) if not fail]}"
            )
        if len(ions) != 2:
            raise ValueError(
                f"Argument 'ions' can only take 2 inputs, received {ions}"
                f"with {len(ions)} inputs. Please try again."
            )
        return ions

    # Validate any speed argument
    def speeds(
        speeds: u.m / u.s,
    ):

        if len(speeds) != 2:
            raise ValueError(
                "Argument 'speeds' can only take 2 inputs, received "
                f"{speeds} with {speeds.ndim} inputs."
            )
        return speeds

    # Validate any temperature argument
    @validate_quantities(
        T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    )
    def temp(
        T: u.K,
    ):
        if not isinstance(T, astropy.units.Quantity):
            raise TypeError(
                "Argument 'T' must be an 'astropy.units.Quantity' with "
                f"a unit of K. Instead got type of {type(T)}."
            )
        elif T.shape != ():
            raise ValueError(
                "Argument 'T' must be single value and not an array of"
                f" shape {T.shape}."
            )
        elif not isinstance(T.value, (int, float)):
            raise TypeError(
                f"Argument 'T' must be an integer or float, received {T} "
                f"with type of {type(T)}."
            )
        elif not T.value > 0:
            raise ValueError(
                f"Argument 'T' must be a positive argument, received "
                f"{T} of type {type(T)}."
            )
        return T



def Hellinger(
    inputs,
    method,
):
    r"""


    Parameters
    ----------
    method: int or float
        An integer or float representing the year of the desired method.
        Current options are ``'2009'``, ``'2010'`` and ``'2016'``.
        Please see the note's section below for additional details.

    kwargs: dict
        Arguments required for the specified method, each method
        requires differing inputs. See notes section below for
        additinoal details.

    Returns
    -------


    Raises
    ------


    Notes
    -----

    +----------+------------------+
    |  Method  |     Function     |
    +----------+------------------+
    |   2009   | `Hellinger_2009` |
    |   2010   | `Hellinger_2010` |
    |   2016   | `Hellinger_2016` |
    +----------+------------------+

    species s on species t.

    Example
    -------

    """

    functions = {2009: Hellinger_2009, 2010: Hellinger_2010, 2016: Hellinger_2016}

    if not isinstance(method, (float, int)):
        raise TypeError(
            "Argument 'method' must be of type float or integer, "
            f"instead got type of {type(method)}."
        )
    elif method not in functions.keys():
        raise ValueError(
            f"Argument 'method' is not a valid entry, received {method} "
            f" and valid methods are {functions.keys()}. Please try again."
        )

    return functions[method](**inputs)


def Hellinger_2009(
    T: u.K,
    n_i: u.m**-3,
    ions: (Particle, Particle),
    par_speeds: u.m / u.s,
):
    r"""
    Compute the collisional timescale as presented by :cite:t:`hellinger:2009`.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        The scalar temperature magnitude in units convertible to K.

    n_i : `~astropy.units.Quantity`
        Ion number density in units convertible to m\ :sup:`-3`.  Must
        be single value and should be the ion of prime interest.

    ions :  a `list` of length 2 containing :term:`particle-like` objects
        A list of length 2 with an instance of the
        :term:`particle-like` object representing the ion species in
        each entry. (e.g., ``"p"`` for protons, ``"D+"`` for deuterium,
         `["p", ``D+``]).

    par_speeds : a `list` of length 2 containing `~astropy.units.Quantity` objects
        A list of length 2 with an `~astropy.units.Quantity` representing
        the PARALLEL velocity with units of  in each entry. (e.g [
        500 * u.m / u.s, 745 * u.m / u.s]).

    Returns
    -------
    :math:`\tau` : `~astropy.units.Quantity`
        The collisional timescale in units of seconds.

    Raises
    ______
    `TypeError`
        If applicable arguments are not instances of
        `~astropy.units.Quantity` or cannot be converted into one.

    `ValueError`
        Number of particles in ``ions`` is not 2 or the input values
        are not valid particles

    `ValueError`
        If ``n_i`` or ``T`` is negative or not a single value.

    `TypeError`
        If ``n_i`` or ``T`` is not of type integer or float.

    `ValueError`
        Number of parallel speeds in``par_speeds`` is not 2.

    `TypeError`
        If the parallel speeds in ``par_speeds`` is not of type
        integer or float

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the
        expected units.

    Notes
    -----
    Assuming a dominance of small angle deflections,

    we assume that all considered species have bi-Maxwellian velocity distribtion functions with a mean velocity parallel to the ambient magnetic field

    Example
    _______

    """

    # Validate arguments argument
    T = validate.temp(T)
    n_i = validate.n_i(n_i)
    ions = validate.ions(ions)
    par_speeds = validate.speeds(par_speeds)

    v_par = np.sqrt((par_speeds[0].value ** 2 + par_speeds[1].value ** 2) / 2)

    a = (ions[0].charge.value ** 2) * (ions[1].charge.value ** 2) * n_i.value

    b = (
        (12 * (pi**1.5))
        * (ions[0].mass.value * ions[1].mass.value)
        * (eps0**2)
        * (v_par**3)
    )

    c = coulomb.Coulomb_logarithm(T, n_i, ions)

    return ((a / b.value) * c) / u.s


def Hellinger_2010(
    T_par: u.K,
    T_perp: u.K,
    n_i: u.m**-3,
    ions: (Particle, Particle),
    par_speeds: u.m / u.s,
):
    r"""
    Compute the collisional timescale as presented by :cite:t:`hellinger:2010`.

    Parameters
    ----------
    T_par : `~astropy.units.Quantity`
        The parallel magnitude of the temperature in units convertible
         to K.

    T_perp : `~astropy.units.Quantity`
        The perpendicular magnitude of the temperature in units
        convertible to K.

    n_i : `~astropy.units.Quantity`
        Ion number density in units convertible to m\ :sup:`-3`.  Must
        be single value and should be the ion of prime interest.

    ions :  a `list` of length 2 containing :term:`particle-like` objects
        A list of length 2 with an instance of the
        :term:`particle-like` object representing the ion species in
        each entry. (e.g., ``"p"`` for protons, ``"D+"`` for deuterium,
         `["p", ``D+``]).

    par_speeds

    Returns
    -------

    Raises
    ______

    Notes
    _____
    The Hellinger_2009 timescale is defined by :cite:t:`hellinger:2009`,
     as denoted in equation 14. The equation assumes the dominance of
    two-particle small-angle collisions and that the distribution
    function is exactly bi-Maxwellian.

    .. math::
        \nu_{\alpha \beta} = \frac{q^{2}_{\alpha}q^{2}_{\beta}n_{\beta}}{12\pi^{3/2}\epsilon_{0}^{2}m_{\alpha}m_{\beta} v}

    index s and t denotes different species,
    The collisional variation in the distribution function of species s
    is given by a sum of terms giving the scattering on all species t in the form8

    :math:`q` is the charge of the respective species, :math:`n_{\gamma}` is the
    particle density of the target ion.




    """

    # Validate other arguments
    T_par = validate.temp(T_par)
    T_perp = validate.temp(T_perp)
    n_i = validate.n_i(n_i)
    ions = validate.ions(ions)
    par_speeds = validate.speeds(par_speeds)

    if T_par == 0:
        raise ValueError(
            "Argument 'T_par' must be a non zero value, received a "
            f"value of {T_par}. Please try again."
        )
    else:
        T = (2 * T_perp + T_par) / 3
        return (
            Hellinger_2009(T, n_i, ions, par_speeds)
            * 3
            / 5
            * gh(a=2, b=1.5, c=7 / 2, x=(1 - (T_perp / T_par)))
        )


def Hellinger_2016(
    T_par: (u.K, u.K),
    T_perp: (u.K, u.K),
    n_i: u.m**-3,
    ions: (Particle, Particle),
    par_speeds: (u.m / u.s, u.m / u.s),
    perp_speeds: (u.m / u.s, u.m / u.s),
):
    r"""
    Compute the collisional timescale as presented by :cite:t:`hellinger:2016`.

    Parameters
    ----------
    T_par

    T_perp

    n_i : `~astropy.units.Quantity`
        Ion number density in units convertible to m\ :sup:`-3`.  Must
        be single value and should be the ion of prime interest.

    ions :  a `list` of length 2 containing :term:`particle-like` objects
        A list of length 2 with an instance of the
        :term:`particle-like` object representing the ion species in
        each entry. (e.g., ``"p"`` for protons, ``"D+"`` for deuterium,
         `["p", ``D+``]).

    par_speeds

    perp_speeds

    Returns
    -------

    Raises
    ______

    Notes
    _____
    assume a homogeneous plasma consisting of species with bi-Maxwellian velocity distribution functions
    We estimated the importance of these parameters assuming isotropic populations


    """

    # Validate arguments
    T_par = validate.temp(T_par)
    T_perp = validate.temp(T_perp)
    n_i = validate.n_i(n_i)
    ions = validate.ions(ions)
    par_speeds = validate.speeds(par_speeds)
    perp_speeds = validate.speeds(perp_speeds)

    # Check for divide by zero error with t_par
    if T_par == 0:
        raise ValueError(
            "Argument 'T_par' must be a non zero value, received a "
            f"value of {T_par}. Please try again."
        )
    else:
        T = (2 * T_perp + T_par) / 3

        vstpar = np.sqrt((par_speeds[0] ** 2 + par_speeds[1] ** 2) / 2)

        Ast = (ions.mass[0] * T_perp[1] + ions.mass[1] * T_perp[0]) / (
            ions.mass[0] * T_par[1] + ions.mass[1] * T_par[0]
        )

        vs = (par_speeds[0] ** 2 + 2 * perp_speeds[0] ** 2) / 3
        vt = (par_speeds[1] ** 2 + 2 * perp_speeds[1] ** 2) / 3

        vst = vs - vt

        return Hellinger_2009(T, n_i, ions, par_speeds) * hyper2d(
            1, 1.5, 2.5, 1 - Ast, Ast * (vst**2 / 4 * vstpar**2)
        )
