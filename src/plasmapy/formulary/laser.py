"""
Functions for calculating quantities associated with laser pulses.

.. attention::

   |expect-api-changes|
"""

__all__ = [
    "em_angular_frequency",
    "electric_field_amplitude",
    "intensity",
    "normalized_vector_potential",
    "Gaussian_beam_waist_radius",
    "Gaussian_power",
    "Gaussian_Rayleigh_length",
    "Gaussian_spot_size_FWHM",
    "em_wavelength",
]
__aliases__ = [
    "a0_",
    "omega_",
    "E0_",
    "I_",
    "w0_",
]

import astropy.units as u
import numpy as np
from astropy.constants.si import c, e, eps0, m_e

from plasmapy.utils.decorators import validate_quantities

__all__ += __aliases__


@validate_quantities(
    intensity={"can_be_negative": False},
)
def electric_field_amplitude(
    intensity: u.Quantity[u.watt / u.m**2],
) -> u.Quantity[u.V / u.m]:
    r"""
    Calculate the electric field amplitude :math:`E_0` from the intensity :math:`I`
    of a laser.

    The electric field amplitude of an electromagnetic plane wave in vacuum
    is calculated using:

    .. math::
        E_0=\sqrt{\frac{2I}{c ε_0}},

    where :math:`c` is the speed of light and
    :math:`ε_0` is the permittivity of free space.

    **Aliases:** `E0_`

    Parameters
    ----------
    intensity : `~astropy.units.Quantity`
        Intensity of the laser pulse (convertible to W / m\ :sup:`2`).

    Returns
    -------
    E : `~astropy.units.Quantity`
        Maximum electric field amplitude for the intensity provided.

    Notes
    -----
    For details, see :cite:t:`ling:2016`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> electric_field_amplitude(1e-3 * u.watt / u.m**2)  # Electric Field Amplitude
    <Quantity 0.8680211 V / m>
    """

    E = np.sqrt((2 * intensity) / (c * eps0))
    return E.to(u.V / u.m)


E0_ = electric_field_amplitude
"""Alias to `~plasmapy.formulary.laser.electric_field_amplitude`."""


@validate_quantities(
    electric_field_amplitude={"can_be_negative": False},
)
def intensity(
    electric_field_amplitude: u.Quantity[u.V / u.m],
) -> u.Quantity[u.watt / u.m**2]:
    r"""
    Calculate the intensity :math:`I` of a laser from the
    electric field amplitude :math:`E_0`.

    The intensity of an electromagnetic plane wave in vacuum
    is calculated using:

    .. math::
        I=\frac{1}{2} c ε_0 E_0^2,

    where :math:`c` is the speed of light and
    :math:`ε_0` is the permittivity of free space.

    **Aliases:** `I_`

    Parameters
    ----------
    electric_field_amplitude: `~astropy.units.Quantity`
        Electric field amplitude of an electromagnetic plane wave
        (convertible to V / m).

    Returns
    -------
    Int : `~astropy.units.Quantity`
        Intensity for the electric field amplitude provided.

    Notes
    -----
    For details, see :cite:t:`ling:2016`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> intensity(0.8680211 * u.V / u.m)
    <Quantity 0.001 W / m2>
    """

    return (1 / 2) * c * eps0 * electric_field_amplitude**2


I_ = intensity
"""Alias to `~plasmapy.formulary.laser.intensity`."""


@validate_quantities(
    angular_frequency={"can_be_negative": False},
)
def em_wavelength(angular_frequency: u.Quantity[u.rad / u.s]) -> u.Quantity[u.m]:
    r"""
    Calculate the wavelength of a laser :math:`\lambda` given the
    the angular frequency :math:`\omega`.

    The wavelength of an electromagnetic wave :math:`\lambda`
    is calculated using:

    .. math::
        \lambda = \frac{2\pi c}{\omega},

    where :math:`\omega` is the angular frequency of the beam.

    Parameters
    ----------
    angular_frequency : `~astropy.units.Quantity`
        angular frequency of the laser beam (convertible to rad / s).

    Returns
    -------
    wavelength : `~astropy.units.Quantity`
        The wavelength of a laser with a given angular frequency.

    Notes
    -----
    For details, see :cite:t:`ling:2016`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> em_wavelength(2.354307546e15 * u.rad / u.s)
    <Quantity 8e-07 m>
    """
    return (2 * np.pi * c / angular_frequency).to(
        u.m, equivalencies=u.dimensionless_angles()
    )


@validate_quantities(
    wavelength={"can_be_negative": False},
)
def em_angular_frequency(wavelength: u.Quantity[u.m]) -> u.Quantity[u.rad / u.s]:
    r"""
    Calculate the angular frequency :math:`\omega` of a laser given the
    the wavelength of the beam :math:`\lambda`.

    The angular frequency of a wave :math:`\omega`
    can be calculated using:

    .. math::
        \omega = \frac{2\pi c}{\lambda},

    where :math:`\lambda` is the wavelength of the beam and :math:`c` is the speed of light.

    **Aliases:** `omega_`

    Parameters
    ----------
    wavelength : `~astropy.units.Quantity`
        Wavelength of the laser beam (convertible to m).

    Returns
    -------
    omega_0 : `~astropy.units.Quantity`
        The angular frequency of a laser with a given wavelength.

    Notes
    -----
    For details, see :cite:t:`ling:2016`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> em_angular_frequency(800 * u.nm)
    <Quantity 2.35456446e+15 rad / s>
    """
    return ((c / wavelength) * 2 * np.pi).to(
        u.rad / u.s, equivalencies=u.dimensionless_angles()
    )


omega_ = em_angular_frequency
"""Alias to `~plasmapy.formulary.laser.em_angular_frequency`."""


@validate_quantities(
    intensity={"can_be_negative": False},
    wavelength={"can_be_negative": False},
)
def normalized_vector_potential(
    intensity: u.Quantity[u.watt / u.m**2],
    wavelength: u.Quantity[u.m],
) -> float | np.ndarray:
    r"""
    Calculate the normalized vector potential :math:`a_0` from the intensity :math:`I`
    and the wavelength :math:`\lambda`.

    The normalized vector potential is also known as a dimensionless oscillation amplitude,
    quiver velocity, or normalized quiver momentum of an electron.

    The normalized vector potential of a laser
    is calculated using:

    .. math::
        a_0= \frac{e E_0}{m_e \omega c}=\frac{e \lambda}{m_e \pi} \sqrt{\frac{I} {2 \epsilon_0 c^5}},

    where :math:`e` is the fudamental charge,
    :math:`E_0` is the electric field amplitude,
    :math:`m_e` is the mass of an electron,
    :math:`\omega` is the angular frequency of the electromagnetic wave,
    :math:`c` is the speed of light,
    :math:`\lambda` is the wavelength,
    :math:`I` is the intensity of the elecromagnetic wave, and
    :math:`\epsilon_0` is the permitivity of free space.

    **Aliases:** `a0_`

    Parameters
    ----------
    intensity : `~astropy.units.Quantity`
        Intensity of the laser pulse (convertible to W / m\ :sup:`2`).

    wavelength : `~astropy.units.Quantity`
        Wavelength of the laser (convertible to m).

    Returns
    -------
    a_0 : float | numpy.ndarray
        The normalized vector potential of a plasma given the intensity and wavelength of
        a laser.

    Notes
    -----
    For details, see :cite:t:`gibbon:2016`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> normalized_vector_potential(1e18 * u.watt / u.cm**2, 1 * u.um)
    np.float64(0.8549297...)
    """

    a0 = (e * wavelength * np.sqrt(intensity / (2 * eps0 * c**5))) / (m_e * np.pi)
    return a0.to(u.dimensionless_unscaled).value  # type: ignore[no-any-return]


a0_ = normalized_vector_potential
"""Alias to `~plasmapy.formulary.laser.normalized_vector_potential`."""


@validate_quantities(
    intensity={"can_be_negative": False},
    beam_waist_radius={"can_be_negative": False},
)
def Gaussian_power(
    intensity: u.Quantity[u.watt / u.m**2],
    beam_waist_radius: u.Quantity[u.m],
) -> u.Quantity[u.Watt]:
    r"""
    Calculate the total power of a Gaussian beam :math:`P_0` from the intensity :math:`I`
    and the beam waist radius :math:`w_0`.

    The total power of a Gaussian beam
    is calculated using:

    .. math::
        P_0=\frac{1}{2}I_0 \pi w_0^2,

    where :math:`w_0` is the beam waist radius and
    :math:`I_0` is the intensity of the beam.

    Parameters
    ----------
    intensity : `~astropy.units.Quantity`
        Intensity of the laser pulse (convertible to W / m\ :sup:`2`).
    beam_waist_radius : `~astropy.units.Quantity`
        Beam waist of the laser pulse (convertible to m).

    Returns
    -------
    P_0 : `~astropy.units.Quantity`
        The total power of the Gaussian beam for the given intensity
        and spot size.

    Notes
    -----
    For details, see :wikipedia:`Gaussian beam`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> Gaussian_power(1e18 * u.watt / u.cm**2, 1 * u.um)  # Total beam power
    <Quantity 1.57079633e+10 W>
    """

    return (1 / 2) * intensity * np.pi * beam_waist_radius**2


@validate_quantities(
    spot_size_FWHM={"can_be_negative": False},
)
def Gaussian_beam_waist_radius(spot_size_FWHM: u.Quantity[u.m]) -> u.Quantity[u.m]:
    r"""
    Calculate the beam waist radius :math:`w_0` for the intensity profile
    of a Gaussian beam given the Full Width at Half Maximum spot size :math:`FWHM`.

    At focus, :math:`w_0` is the transverse distance from the center of the beam
    to where the intensity drops by a factor of :math:`1/e^2`.

    The beam waist radius of a Gaussian beam
    is calculated using:

    .. math::
        w_0=\frac{FWHM}{\sqrt{2 \ln {2}}},

    where :math:`FWHM` is the full width at half maximum spot size of the beam.

    **Aliases:** `w0_`

    Parameters
    ----------
    spot_size_FWHM : `~astropy.units.Quantity`
        Full Width at Half Maximum spot size of the Gaussian beam (convertible to m).

    Returns
    -------
    w_0 : `~astropy.units.Quantity`
        The beam waist radius of the Gaussian beam for the given FWHM spot size.

    See Also
    --------
    Gaussian_spot_size_FWHM

    Notes
    -----
    For details, see :wikipedia:`Gaussian beam`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> Gaussian_beam_waist_radius(8.242 * u.um)
    <Quantity 7e-6 m>
    """
    return spot_size_FWHM / np.sqrt(2 * np.log(2))


w0_ = Gaussian_beam_waist_radius
"""Alias to `~plasmapy.formulary.laser.Gaussian_beam_waist_radius`."""


@validate_quantities(
    beam_waist_radius={"can_be_negative": False},
)
def Gaussian_spot_size_FWHM(beam_waist_radius: u.Quantity[u.m]) -> u.Quantity[u.m]:
    r"""
    Calculate the Full Width at Half Maximum spot size :math:`FWHM` at focus given the
    beam waist radius of a Gaussian beam :math:`w_0`.

    The FWHM spot size of a Gaussian beam
    is calculated using:

    .. math::
        FWHM = w_0 \sqrt{2 \ln {2}},

    where :math:`w_0` is the beam waist radius of the beam.

    Parameters
    ----------
    beam_waist_radius : `~astropy.units.Quantity`
        beam waist radius of the Gaussian beam (convertible to nm).

    Returns
    -------
    FWHM : `~astropy.units.Quantity`
        The FWHM spot size of the Gaussian beam for the given beam waist.

    See Also
    --------
    Gaussian_beam_waist_radius

    Notes
    -----
    For details, see :wikipedia:`Gaussian beam`\ .

    Examples
    --------
    >>> import astropy.units as u
    >>> Gaussian_spot_size_FWHM(7 * u.um)
    <Quantity 8.242e-6 m>
    """
    return beam_waist_radius * np.sqrt(2 * np.log(2))


# use kwargs
@validate_quantities(
    wavelength={"can_be_negative": False},
    beam_waist_radius={"can_be_negative": False},
)
def Gaussian_Rayleigh_length(
    wavelength: u.Quantity[u.m],
    beam_waist_radius: u.Quantity[u.m],
) -> u.Quantity[u.m]:
    r"""
    Calculate the Rayleigh length :math:`z_R` from the beam waist radius :math:`w_0`
    and the wavelength :math:`\lambda`.

    The Rayleigh length of a Gaussian beam
    is calculated using:

    .. math::
        z_R=\frac{\pi w_0 ^2}{\lambda},

    where :math:`w_0` is the beam waist and
    :math:`\lambda` is the wavelength of the beam.

    Parameters
    ----------
    wavelength : `~astropy.units.Quantity`
        Wavelength of the laser pulse (convertible to m).
    beam_waist_radius : `~astropy.units.Quantity`
        Beam waist of the laser pulse (convertible to m).

    Returns
    -------
    z_R : `~astropy.units.Quantity`
        The Rayleigh length of the Gaussian beam for the given wavelength
        and beam waist.

    See Also
    --------
    Gaussian_beam_waist_radius
    Gaussian_spot_size_FWHM

    Notes
    -----
    For details, see :wikipedia:`Rayleigh length`.

    Examples
    --------
    >>> import astropy.units as u
    >>> Gaussian_Rayleigh_length(800 * u.nm, 1 * u.um)
    <Quantity 3.927e-06 m>
    """

    return (np.pi * beam_waist_radius**2) / wavelength
