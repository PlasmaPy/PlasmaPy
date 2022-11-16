"""
Module contains the functionality for computing different
timescales commonly associated with collisions.
"""


"""
This module contains functionality for calculating the timescales
for a range of configurations.
"""

__all__ = ["Hellinger", "USSR"]

import numpy as np
from math import pi as pi, gamma as gamma, factorial as fact
import astropy.units as u

from astropy.constants.si import eps0
from scipy.stats import gausshyper as gh

from plasmapy.formulary.collisions import coulomb
from plasmapy.particles import Particle, ParticleList
from plasmapy.utils.decorators import validate_quantities


def validate(
        n_i: u.m**-3,
        ions: (Particle, Particle),
        par_speeds: (u.m/u.s, u.m/u.s),
):
    # Validate ions argument
    if not isinstance(ions, (list, tuple, ParticleList)):
        ions = [ions]
    ions = ParticleList(ions)

    if not all(failed := [ion.is_ion and abs(ion.charge_number) > 0 for ion in ions]):
        raise ValueError(
            "Particle(s) passed to 'ions' must be a charged"
            " ion. The following particle(s) is(are) not allowed "
            f"{[ion for ion, fail in zip(ions, failed) if not fail]}"
        )

    # Validate ions dimension
    if len(ions) != 2:
        raise ValueError(
            f"Argument 'ions' can only take 2 inputs, received {ions}"
            f"with {len(ions)} inputs. Please try again."
        )

    # Validate n_i argument
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

    # Validate par_speeds argument
    if len(par_speeds) != 2:
        raise ValueError(
            "Argument 'par_speeds' can only take 2 inputs, received "
            f"{par_speeds} with {len(par_speeds)} inputs."
        )
    else:
        for j in range(2):
            if not isinstance(par_speeds[j].value, (float, int)):
                raise TypeError(
                    f"Argument {par_speeds[j].value} is of incorrect type, "
                    f"type int or float require and got {type(par_speeds[j].value)}"
                )


def validate_temp(
        T: u.K,
):
    # Validate temperature argument
    if T.shape != ():
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
        version,
):
    r"""
    Compute the collisional timescale as presented by :cite:t:`hellinger:2009`.

    Parameters
    ----------
    T : `~astropy.units.Quantity`
        The scalar temperature magnitude in units convertible to K.
        Applicable to version '2009'.

    T_par : `~astropy.units.Quantity`
        The parallel magnitude of the temperature in units convertible
         to K. Applicable to version '2010' and '2016'.

    T_perp : `~astropy.units.Quantity`
        The perpendicular magnitude of the temperature in units
        convertible to K. Applicable to version '2010' and '2016'.

    n_i : `~astropy.units.Quantity`
        Ion number density in units convertible to m\ :sup:`-3`.  Must
        be single value and should be the ion of prime interest.

    ions :  a `list` of length 2 containing :term:`particle-like` objects
        A list of length 2 with an instance of the
        :term:`particle-like` object representing the ion species in
        each entry. (e.g., ``"p"`` for protons, ``"D+"`` for deuterium,
         `["p", ``D+``]).

    par_speeds : a `list` of length 2 containing :term:`particle-like` objects
        A list of length 2 with an `~astropy.units.Quantity` representing
        the PARALLEL velocity with units of  in each entry. (e.g [
        500 * u.m / u.s, 745 * u.m / u.s]).

    Returns
    -------
    :math:`\tau` : `~astropy.units.Quantity`
        The collisional timescale in units of seconds.

    Raises
    ------
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

    species s on species t.

    Example
    -------
    >>> from astropy import units as u
    >>> from plasmapy.particles import Particle
    >>> from plasmapy.formulary.collisions.timescales import Hellinger
    >>> inputs = {
    ...     "T": 8.3e-9 * u.T,
    ...     "n_i": 4.0e5 * u.m**-3,
    ...     "ions": [Particle("H+"), Particle("He+")],
    ...     "par_speeds":  [500, 750] * u.m /u.s,
    ... }
    >>>
    <Quantity 1 / s>
    """

    valid_versions = [2009, 2010, 2016]
    valid_functions = [Hellinger_2009, Hellinger_2010, Hellinger_2016]

    if not isinstance(version, (float, int)):
        raise TypeError(
            "Argument 'version' must be of type float or integer, "
            f"instead got type of {type(version)}."
        )
    elif version not in valid_versions:
        raise ValueError(
            "Argument 'version' is not a valid entry, valid entries "
            f"are {valid_versions}. Please try again."
        )

    return valid_functions[valid_versions.index(version)](*inputs)



def Hellinger_2009(
        T: u.K,
        n_i: u.m**-3,
        ions: (Particle, Particle),
        par_speeds: (u.m/u.s, u.m/u.s),
):
    # Validate temperature aguments
    T = validate_temp(T)

    # Validate other arguments argument
    n_i, ions, par_speeds = validate(n_i, ions, par_speeds)

    v_par = np.sqrt((par_speeds[0].value ** 2 + par_speeds[1].value ** 2) / 2)

    a = ((ions[0].charge.value ** 2) * (ions[1].charge.value ** 2) * n_i.value)

    b = (12 * (pi ** 1.5)) * (ions[0].mass.value * ions[1].mass.value) * (eps0 ** 2) * (v_par ** 3)

    c = coulomb.Coulomb_logarithm(T, n_i, ions)

    return ((a / b.value) * c) / u.s


def Hellinger_2010(
        T_par: u.K,
        T_perp: u.K,
        n_i: u.m ** -3,
        ions: (Particle, Particle),
        par_speeds: (u.m / u.s, u.m / u.s)
):

    # Validate t_par and t_perp
    T_perp = validate_temp(T_perp)
    T_par = validate_temp(T_par)

    #Validate other arguments
    n_i, ions, par_speeds = validate(n_i, ions, par_speeds)

    if T_par == 0:
        raise ValueError(
            "Argument T_par must be a non zero value, please try again."
        )
    else:
        T = (2*T_perp + T_par)/3
        return Hellinger.CoulombCollisionsBiMaxwellian(T, n_i, ions, par_speeds) * 3 / 5 * gh(a=2, b=1.5, c=7 / 2, x=(1 - (T_perp / T_par)))


def Hellinger_2016(
    T_par: (u.K, u.K),
    T_perp: (u.K, u.K),
    n_i: u.m ** -3,
    ions: (Particle, Particle),
    par_speeds: (u.m / u.s, u.m / u.s)
):
    # Validate temperature arguments
    for arg in (T_par, T_perp):
        if arg.shape != 2:
            raise ValueError(
                f"Argument {arg} must be single value and not an array of"
                f" shape {arg.shape}."
            )
        for i in range(2):
            arg[i] = validate_temp(arg[i])

    # Validate other arguments
    n_i, ions, par_speeds = validate(n_i, ions, par_speeds)

    if T_par == 0:
        raise ValueError(
            "Argument T_par must be a non zero value, please try again."
        )
    else:
        T = (2 * T_perp + T_par) / 3
        return Hellinger.CoulombCollisionsBiMaxwellian(T, n_i, ions, par_speeds) * gh(a=2, b=1.5, c=7 / 2,
                                                                                              x=(1 - (T_perp / T_par)))


    return