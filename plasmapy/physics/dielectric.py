"""Functions to calculate plasma dielectric parameters"""

from astropy import units as u
from plasmapy import utils
from plasmapy.physics import parameters
from plasmapy.constants import (pi, m_e, c, mu0, e, eps0)

r"""
Values should be returned as a `~astropy.units.Quantity` in SI units.
"""

__all__ = ['cold_plasma_permittivity_SDP',
           'cold_plasma_permittivity_LRP']

@utils.check_quantity({
    'B': {'units': u.T, 'can_be_negative': False},
    'omega': {'units': u.rad / u.s, 'can_be_negative': False},
})
def cold_plasma_permittivity_SDP(B, species, n, omega):
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
    S : ~astropy.units.Quantity
        S ("Sum") dielectric tensor element.

    D : ~astropy.units.Quantity
        D ("Difference") dielectric tensor element.

    P : ~astropy.units.Quantity
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
    >>> from plasmapy.constants import pi
    >>> B = 2*u.T
    >>> species = ['e', 'D+']
    >>> n = [1e18*u.m**-3, 1e18*u.m**-3]
    >>> omega = 3.7e9*(2*pi)*(u.rad/u.s)
    >>> S, D, P = cold_plasma_permittivity_SDP(B, species, n, omega)
    >>> S
    <Quantity 1.02422902>
    >>> D
    <Quantity 0.39089352>
    >>> P
    <Quantity -4.8903104>
    """
    S, D, P = 1, 0, 1

    for s, n_s in zip(species, n):
        omega_c = parameters.gyrofrequency(B=B, particle=s, signed=True)
        omega_p = parameters.plasma_frequency(n=n_s, particle=s)

        S += - omega_p ** 2 / (omega ** 2 - omega_c ** 2)
        D += omega_c / omega * omega_p ** 2 / (omega ** 2 - omega_c ** 2)
        P += - omega_p ** 2 / omega ** 2
    return S, D, P


def cold_plasma_permittivity_LRP(B, species, n, omega):
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
    L : ~astropy.units.Quantity
        L ("Left") Left-handed circularly polarization tensor element.

    R : ~astropy.units.Quantity
        R ("Right") Right-handed circularly polarization tensor element.

    P : ~astropy.units.Quantity
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
    >>> from plasmapy.constants import pi
    >>> B = 2*u.T
    >>> species = ['e', 'D+']
    >>> n = [1e18*u.m**-3, 1e18*u.m**-3]
    >>> omega = 3.7e9*(2*pi)*(u.rad/u.s)
    >>> L, R, P = cold_plasma_permittivity_LRP(B, species, n, omega)
    >>> L
    <Quantity 0.63333549>
    >>> R
    <Quantity 1.41512254>
    >>> P
    <Quantity -4.8903104>
    """
    L, R, P = 1, 1, 1

    for s, n_s in zip(species, n):
        omega_c = parameters.gyrofrequency(B=B, particle=s, signed=True)
        omega_p = parameters.plasma_frequency(n=n_s, particle=s)

        L += - omega_p ** 2 / (omega * (omega - omega_c))
        R += - omega_p ** 2 / (omega * (omega + omega_c))
        P += - omega_p ** 2 / omega ** 2
    return L, R, P
