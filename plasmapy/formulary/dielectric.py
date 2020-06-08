"""Functions to calculate plasma dielectric parameters"""
__all__ = [
    "cold_plasma_permittivity_SDP",
    "cold_plasma_permittivity_LRP",
    "permittivity_1D_Maxwellian",
]

from collections import namedtuple

import numpy as np
from astropy import units as u
from numpy import pi

from plasmapy.formulary import parameters
from plasmapy.formulary.dispersionfunction import plasma_dispersion_func_deriv
from plasmapy.utils.decorators import validate_quantities

r"""
Values should be returned as a `~astropy.units.Quantity` in SI units.
"""

StixTensorElements = namedtuple("StixTensorElements", ["sum", "difference", "plasma"])
RotatingTensorElements = namedtuple(
    "RotatingTensorElements", ["left", "right", "plasma"]
)


@validate_quantities(B={"can_be_negative": False}, omega={"can_be_negative": False})
def cold_plasma_permittivity_SDP(B: u.T, species, n, omega: u.rad / u.s):
    r"""
    Magnetized Cold Plasma Dielectric Permittivity Tensor Elements.

    Elements (S, D, P) are given in the "Stix" frame, ie. with B // z.

    The :math:`\exp(-i \omega t)` time-harmonic convention is assumed.

    Parameters
    ----------
    B : ~astropy.units.Quantity
        Magnetic field magnitude in units convertible to tesla.

    species : list of str
        List of the plasma particle species
        e.g.: ['e', 'D+'] or ['e', 'D+', 'He+'].

    n : list of ~astropy.units.Quantity
        `list` of species density in units convertible to per cubic meter
        The order of the species densities should follow species.

    omega : ~astropy.units.Quantity
        Electromagnetic wave frequency in rad/s.

    Returns
    -------
    sum : ~astropy.units.Quantity
        S ("Sum") dielectric tensor element.

    difference : ~astropy.units.Quantity
        D ("Difference") dielectric tensor element.

    plasma : ~astropy.units.Quantity
        P ("Plasma") dielectric tensor element.

    Notes
    -----
    The dielectric permittivity tensor is expressed in the Stix frame with
    the :math:`\exp(-i \omega t)` time-harmonic convention as
    :math:`\varepsilon = \varepsilon_0 A`, with :math:`A` being

    .. math::

        \varepsilon = \varepsilon_0 \left(\begin{matrix}  S & -i D & 0 \\
                              +i D & S & 0 \\
                              0 & 0 & P \end{matrix}\right)

    where:

    .. math::
        S = 1 - \sum_s \frac{\omega_{p,s}^2}{\omega^2 - \Omega_{c,s}^2}

        D = \sum_s \frac{\Omega_{c,s}}{\omega}
            \frac{\omega_{p,s}^2}{\omega^2 - \Omega_{c,s}^2}

        P = 1 - \sum_s \frac{\omega_{p,s}^2}{\omega^2}

    where :math:`\omega_{p,s}` is the plasma frequency and
    :math:`\Omega_{c,s}` is the signed version of the cyclotron frequency
    for the species :math:`s`.

    References
    ----------
    - T.H. Stix, Waves in Plasma, 1992.

    Examples
    --------
    >>> from astropy import units as u
    >>> from numpy import pi
    >>> B = 2*u.T
    >>> species = ['e', 'D+']
    >>> n = [1e18*u.m**-3, 1e18*u.m**-3]
    >>> omega = 3.7e9*(2*pi)*(u.rad/u.s)
    >>> permittivity = S, D, P = cold_plasma_permittivity_SDP(B, species, n, omega)
    >>> S
    <Quantity 1.02422...>
    >>> permittivity.sum   # namedtuple-style access
    <Quantity 1.02422...>
    >>> D
    <Quantity 0.39089...>
    >>> P
    <Quantity -4.8903...>
    """
    S, D, P = 1, 0, 1

    for s, n_s in zip(species, n):
        omega_c = parameters.gyrofrequency(B=B, particle=s, signed=True)
        omega_p = parameters.plasma_frequency(n=n_s, particle=s)

        S += -(omega_p ** 2) / (omega ** 2 - omega_c ** 2)
        D += omega_c / omega * omega_p ** 2 / (omega ** 2 - omega_c ** 2)
        P += -(omega_p ** 2) / omega ** 2
    return StixTensorElements(S, D, P)


@validate_quantities(B={"can_be_negative": False}, omega={"can_be_negative": False})
def cold_plasma_permittivity_LRP(B: u.T, species, n, omega: u.rad / u.s):
    r"""
    Magnetized Cold Plasma Dielectric Permittivity Tensor Elements.

    Elements (L, R, P) are given in the "rotating" basis, ie. in the basis
    :math:`(\mathbf{u}_{+}, \mathbf{u}_{-}, \mathbf{u}_z)`,
    where the tensor is diagonal and with B // z.

    The :math:`\exp(-i \omega t)` time-harmonic convention is assumed.

    Parameters
    ----------
    B : ~astropy.units.Quantity
        Magnetic field magnitude in units convertible to tesla.

    species : list of str
        The plasma particle species (e.g.: `['e', 'D+']` or
        `['e', 'D+', 'He+']`.

    n : list of ~astropy.units.Quantity
        `list` of species density in units convertible to per cubic meter.
        The order of the species densities should follow species.

    omega : ~astropy.units.Quantity
        Electromagnetic wave frequency in rad/s.

    Returns
    -------
    left : ~astropy.units.Quantity
        L ("Left") Left-handed circularly polarization tensor element.

    right : ~astropy.units.Quantity
        R ("Right") Right-handed circularly polarization tensor element.

    plasma : ~astropy.units.Quantity
        P ("Plasma") dielectric tensor element.

    Notes
    -----
    In the rotating frame defined by
    :math:`(\mathbf{u}_{+}, \mathbf{u}_{-}, \mathbf{u}_z)`
    with :math:`\mathbf{u}_{\pm}=(\mathbf{u}_x \pm \mathbf{u}_y)/\sqrt{2}`,
    the dielectric tensor takes a diagonal form with elements L, R, P with:

    .. math::
        L = 1 - \sum_s
                \frac{\omega_{p,s}^2}{\omega\left(\omega - \Omega_{c,s}\right)}

        R = 1 - \sum_s
                \frac{\omega_{p,s}^2}{\omega\left(\omega + \Omega_{c,s}\right)}

        P = 1 - \sum_s \frac{\omega_{p,s}^2}{\omega^2}

    where :math:`\omega_{p,s}` is the plasma frequency and
    :math:`\Omega_{c,s}` is the signed version of the cyclotron frequency
    for the species :math:`s`.

    References
    ----------
    - T.H. Stix, Waves in Plasma, 1992.

    Examples
    --------
    >>> from astropy import units as u
    >>> from numpy import pi
    >>> B = 2*u.T
    >>> species = ['e', 'D+']
    >>> n = [1e18*u.m**-3, 1e18*u.m**-3]
    >>> omega = 3.7e9*(2*pi)*(u.rad/u.s)
    >>> L, R, P = permittivity = cold_plasma_permittivity_LRP(B, species, n, omega)
    >>> L
    <Quantity 0.63333...>
    >>> permittivity.left    # namedtuple-style access
    <Quantity 0.63333...>
    >>> R
    <Quantity 1.41512...>
    >>> P
    <Quantity -4.8903...>
    """
    L, R, P = 1, 1, 1

    for s, n_s in zip(species, n):
        omega_c = parameters.gyrofrequency(B=B, particle=s, signed=True)
        omega_p = parameters.plasma_frequency(n=n_s, particle=s)

        L += -(omega_p ** 2) / (omega * (omega - omega_c))
        R += -(omega_p ** 2) / (omega * (omega + omega_c))
        P += -(omega_p ** 2) / omega ** 2
    return RotatingTensorElements(L, R, P)


@validate_quantities(
    kWave={"none_shall_pass": True}, validations_on_return={"can_be_complex": True}
)
def permittivity_1D_Maxwellian(
    omega: u.rad / u.s,
    kWave: u.rad / u.m,
    T: u.K,
    n: u.m ** -3,
    particle,
    z_mean: u.dimensionless_unscaled = None,
) -> u.dimensionless_unscaled:
    r"""
    The classical dielectric permittivity for a 1D Maxwellian plasma. This
    function can calculate both the ion and electron permittivities. No
    additional effects are considered (e.g. magnetic fields, relativistic
    effects, strongly coupled regime, etc.)

    Parameters
    ----------
    omega : ~astropy.units.Quantity
        The frequency in rad/s of the electromagnetic wave propagating
        through the plasma.

    kWave : ~astropy.units.Quantity
        The corresponding wavenumber, in rad/m, of the electromagnetic wave
        propagating through the plasma. This is often modulated by the
        dispersion of the plasma or by relativistic effects. See em_wave.py
        for ways to calculate this.

    T : ~astropy.units.Quantity
        The plasma temperature - this can be either the electron or the ion
        temperature, but should be consistent with density and particle.

    n : ~astropy.units.Quantity
        The plasma density - this can be either the electron or the ion
        density, but should be consistent with temperature and particle.

    particle : str
        The plasma particle species.

    z_mean : str
        The average ionization of the plasma. This is only required for
        calculating the ion permittivity.

    Returns
    -------
    chi : ~astropy.units.Quantity
        The ion or the electron dielectric permittivity of the plasma.
        This is a dimensionless quantity.

    Notes
    -----
    The dielectric permittivities for a Maxwellian plasma are described
    by the following equations [1]_

    .. math::
        \chi_e(k, \omega) = - \frac{\alpha_e^2}{2} Z'(x_e)

        \chi_i(k, \omega) = - \frac{\alpha_i^2}{2}\frac{Z}{} Z'(x_i)

        \alpha = \frac{\omega_p}{k v_{Th}}

        x = \frac{\omega}{k v_{Th}}

    :math:`chi_e` and :math:`chi_i` are the electron and ion permittivities
    respectively. :math:`Z'` is the derivative of the plasma dispersion
    function. :math:`\alpha` is the scattering parameter which delineates
    the difference between the collective and non-collective Thomson
    scattering regimes. :math:`x` is the dimensionless phase velocity
    of the EM wave propagating through the plasma.

    References
    ----------
    .. [1] J. Sheffield, D. Froula, S. H. Glenzer, and N. C. Luhmann Jr,
       Plasma scattering of electromagnetic radiation: theory and measurement
       techniques. Chapter 5 Pg 106 (Academic press, 2010).

    Example
    -------
    >>> from astropy import units as u
    >>> from numpy import pi
    >>> from astropy.constants import c
    >>> T = 30 * 11600 * u.K
    >>> n = 1e18 * u.cm**-3
    >>> particle = 'Ne'
    >>> z_mean = 8 * u.dimensionless_unscaled
    >>> vTh = parameters.thermal_speed(T, particle, method="most_probable")
    >>> omega = 5.635e14 * 2 * pi * u.rad / u.s
    >>> kWave = omega / vTh
    >>> permittivity_1D_Maxwellian(omega, kWave, T, n, particle, z_mean)
    <Quantity -6.72809...e-08+5.76037...e-07j>
    """
    # thermal velocity
    vTh = parameters.thermal_speed(T=T, particle=particle, method="most_probable")
    # plasma frequency
    wp = parameters.plasma_frequency(n=n, particle=particle, z_mean=z_mean)
    # scattering parameter alpha.
    # explicitly removing factor of sqrt(2) to be consistent with Froula
    alpha = np.sqrt(2) * (wp / (kWave * vTh)).to(u.dimensionless_unscaled)
    # The dimensionless phase velocity of the propagating EM wave.
    zeta = (omega / (kWave * vTh)).to(u.dimensionless_unscaled)
    chi = alpha ** 2 * (-1 / 2) * plasma_dispersion_func_deriv(zeta.value)
    return chi
