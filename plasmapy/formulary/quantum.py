"""
Functions for quantum parameters, including electron degenerate
gases and warm dense matter.

"""
__all__ = [
    "chemical_potential",
    "deBroglie_wavelength",
    "Ef_",
    "Fermi_energy",
    "lambdaDB_",
    "lambdaDB_th_",
    "Thomas_Fermi_length",
    "thermal_deBroglie_wavelength",
    "Wigner_Seitz_radius",
]

import numpy as np
from astropy import units as u
from astropy.constants.si import c, e, eps0, h, hbar, k_B, m_e

from plasmapy import particles
from plasmapy.formulary import mathematics
from plasmapy.formulary.relativity import Lorentz_factor
from plasmapy.utils import RelativityError
from plasmapy.utils.decorators import validate_quantities


# TODO: Use @check_relativistic and @particle_input
@validate_quantities(
    V={"can_be_negative": True}, validations_on_return={"can_be_negative": False}
)
def deBroglie_wavelength(V: u.m / u.s, particle) -> u.m:
    r"""
    Calculates the de Broglie wavelength.

    **Aliases:** `lambdaDB_`

    Parameters
    ----------
    V : ~astropy.units.Quantity
        Particle velocity in units convertible to meters per second.

    particle : str or ~astropy.units.Quantity
        Representation of the particle species (e.g., `'e'`, `'p'`, `'D+'`,
        or `'He-4 1+'`, or the particle mass in units convertible to
        kilograms.

    Returns
    -------
    lambda_dB : ~astropy.units.Quantity
        The de Broglie wavelength in units of meters.

    Raises
    ------
    TypeError
        The velocity is not a `~astropy.units.Quantity` and cannot be
        converted into a ~astropy.units.Quantity.

    ~astropy.units.UnitConversionError
        If the velocity is not in appropriate units.

    ~plasmapy.utils.RelativityError
        If the magnitude of `V` is faster than the speed of light.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed

    Notes
    -----
    The de Broglie wavelength is given by

    .. math::

        \lambda_{dB} = \frac{h}{p} = \frac{h}{\gamma m V}

    where :math:`h` is the Planck constant, :math:`p` is the
    relativistic momentum of the particle, :math:`gamma` is the
    Lorentz factor, :math:`m` is the particle's mass, and :math:`V` is the
    particle's velocity.

    Examples
    --------
    >>> from astropy import units as u
    >>> velocity = 1.4e7 * u.m / u.s
    >>> deBroglie_wavelength(velocity, 'e')
    <Quantity 5.18997095e-11 m>
    >>> deBroglie_wavelength(V = 0 * u.m / u.s, particle = 'D+')
    <Quantity inf m>
    """

    V = np.abs(V)

    if np.any(V >= c):
        raise RelativityError(
            "Velocity input in deBroglie_wavelength cannot "
            "be greater than or equal to the speed of "
            "light."
        )

    if not isinstance(particle, u.Quantity):
        try:
            # TODO: Replace with more general routine!
            m = particles.particle_mass(particle)
        except Exception:
            raise ValueError("Unable to find particle mass.")
    else:
        try:
            m = particle.to(u.kg)
        except Exception:
            raise u.UnitConversionError(
                "The second argument for deBroglie"
                " wavelength must be either a "
                "representation of a particle or a"
                " Quantity with units of mass."
            )

    if V.size > 1:

        lambda_dBr = np.ones(V.shape) * np.inf * u.m
        indices = V.value != 0
        lambda_dBr[indices] = h / (m * V[indices] * Lorentz_factor(V[indices]))

    else:

        if V == 0 * u.m / u.s:
            lambda_dBr = np.inf * u.m
        else:
            lambda_dBr = h / (Lorentz_factor(V) * m * V)

    return lambda_dBr


lambdaDB_ = deBroglie_wavelength
""" Alias to :func:`deBroglie_wavelength`. """


@validate_quantities(
    T_e={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    validations_on_return={"can_be_negative": False},
)
def thermal_deBroglie_wavelength(T_e: u.K) -> u.m:
    r"""
    Calculate the thermal deBroglie wavelength for electrons.

    **Aliases:** `lambdaDB_th_`

    Parameters
    ----------
    T_e: ~astropy.units.Quantity
        Electron temperature.

    Returns
    -------
    lambda_dbTh: ~astropy.units.Quantity
        The thermal deBroglie wavelength for electrons in meters.

    Raises
    ------
    TypeError
        If argument is not a `~astropy.units.Quantity`.

    ~astropy.units.UnitConversionError
        If argument is in incorrect units.

    ValueError
        If argument contains invalid values.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.

    Notes
    -----
    The thermal deBroglie wavelength is approximately the average deBroglie
    wavelength for electrons in an ideal gas and is given by

    .. math::

       \lambda_{dbTh} = \frac{h}{\sqrt{2 \pi m_e k_B T_e}}

    Example
    -------
    >>> from astropy import units as u
    >>> thermal_deBroglie_wavelength(1 * u.eV)
    <Quantity 6.9193675e-10 m>
    """
    lambda_dbTh = h / np.sqrt(2 * np.pi * m_e * k_B * T_e)
    return lambda_dbTh


lambdaDB_th_ = thermal_deBroglie_wavelength
""" Alias to :func:`thermal_deBroglie_wavelength`. """


@validate_quantities(
    n_e={"can_be_negative": False}, validations_on_return={"can_be_negative": False}
)
def Fermi_energy(n_e: u.m ** -3) -> u.J:
    r"""
    Calculate the kinetic energy in a degenerate electron gas.

    **Aliases:** `Ef_`

    Parameters
    ----------
    n_e : ~astropy.units.Quantity
        Electron number density.

    Returns
    -------
    energy_F : ~astropy.units.Quantity
        The Fermi energy in Joules.

    Raises
    ------
    TypeError
        If argument is not a `~astropy.units.Quantity`.

    ~astropy.units.UnitConversionError
        If argument is in incorrect units.

    ValueError
        If argument contains invalid values.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.

    Notes
    -----
    The Fermi energy is the kinetic energy in a degenerate electron gas
    and is given by

    .. math::

       E_F = \frac{\pi^2 \hbar^2}{2 m_{e}}
       \left( \frac{3 n_{e}}{\pi} \right )^{2/3}

    This quantity is often used in place of thermal energy for analysis
    of cold, dense plasmas (e.g. warm dense matter, condensed matter).

    See also
    --------
    Thomas_Fermi_length

    Example
    -------
    >>> from astropy import units as u
    >>> Fermi_energy(1e23 * u.cm**-3)
    <Quantity 1.2586761e-18 J>
    """
    coeff = (np.pi * hbar) ** 2 / (2 * m_e)
    energy_F = coeff * (3 * n_e / np.pi) ** (2 / 3)
    return energy_F


Ef_ = Fermi_energy
""" Alias to :func:`Fermi_energy`. """


@validate_quantities(
    n_e={"can_be_negative": False}, validations_on_return={"can_be_negative": False}
)
def Thomas_Fermi_length(n_e: u.m ** -3) -> u.m:
    r"""
    Calculate the exponential scale length for charge screening
    for cold and dense plasmas.

    Parameters
    ----------
    n_e: ~astropy.units.Quantity
        Electron number density.

    Returns
    -------
    lambda_TF: ~astropy.units.Quantity
        The Thomas-Fermi screening length in meters.

    Raises
    ------
    TypeError
        If argument is not a `~astropy.units.Quantity`.

    ~astropy.units.UnitConversionError
        If argument is in incorrect units.

    ValueError
        If argument contains invalid values.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.

    Notes
    -----
    The Thomas-Fermi screening length is the exponential scale length for
    charge screening and is given by

    .. math::

       \lambda_TF = \sqrt{\frac{2 \epsilon_0 E_F}{3 n_e e^2}}

    for an electron degenerate gas.

    This quantity is often used in place of the Debye length for analysis
    of cold, dense plasmas (e.g. warm dense matter, condensed matter).

    The electrical potential will drop by a factor of 1/e every Thomas-Fermi
    screening length.

    Plasmas will generally be quasineutral on length scales significantly
    larger than the Thomas-Fermi screening length.

    See also
    --------
    Fermi_energy
    plasmapy.formulary.Debye_length

    Example
    -------
    >>> from astropy import units as u
    >>> Thomas_Fermi_length(1e23 * u.cm**-3)
    <Quantity 5.37991409e-11 m>

    """
    energy_F = Fermi_energy(n_e)
    lambda_TF = np.sqrt(2 * eps0 * energy_F / (3 * n_e * e ** 2))
    return lambda_TF


@validate_quantities(
    n={"can_be_negative": False}, validations_on_return={"can_be_negative": False}
)
def Wigner_Seitz_radius(n: u.m ** -3) -> u.m:
    r"""
    Calculate the Wigner-Seitz radius, which approximates the inter-
    particle spacing. It is the radius of a sphere whose volume is
    equal to the mean volume per atom in a solid. This parameter is
    often used to calculate the coupling parameter.
    When ion density is used, this is the ion sphere radius, i.e., the
    space occupied by a single ion with no other ions in that space. Higher
    density means less space for each ion, so the radius is smaller.

    Parameters
    ----------
    n: ~astropy.units.Quantity
        Particle number density.

    Returns
    -------
    radius: ~astropy.units.Quantity
        The Wigner-Seitz radius in meters.

    Raises
    ------
    TypeError
        If argument is not a ~astropy.units.Quantity.

    ~astropy.units.UnitConversionError
        If argument is in incorrect units.

    ValueError
        If argument contains invalid values.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.

    Notes
    -----
    The Wigner-Seitz radius approximates the interparticle spacing.
    It is the radius of a sphere whose volume is equal to the mean
    volume per atom in a solid:

    .. math::
        r = \left(\frac{3}{4 \pi n}\right)^{1/3}

    See also
    --------
    Fermi_energy

    Example
    -------
    >>> from astropy import units as u
    >>> Wigner_Seitz_radius(1e29 * u.m**-3)
    <Quantity 1.33650462e-10 m>

    """
    radius = (3 / (4 * np.pi * n)) ** (1 / 3)
    return radius


# TODO: remove NotImplementedError and 'doctest: +SKIP' when the following issues are addressed...
#       https://github.com/PlasmaPy/PlasmaPy/issues/726
#       https://github.com/astropy/astropy/issues/9721
@validate_quantities(
    n_e={"can_be_negative": False},
    T={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def chemical_potential(n_e: u.m ** -3, T: u.K) -> u.dimensionless_unscaled:
    r"""
    Calculate the ideal chemical potential.

    Parameters
    ----------
    n_e: ~astropy.units.Quantity
        Electron number density.

    T : ~astropy.units.Quantity
        The temperature.

    Returns
    -------
    beta_mu: ~astropy.units.Quantity
        The dimensionless ideal chemical potential. That is the ratio of
        the ideal chemical potential to the thermal energy.

    Raises
    ------
    TypeError
        If argument is not a `~astropy.units.Quantity`.

    ~astropy.units.UnitConversionError
        If argument is in incorrect units.

    ValueError
        If argument contains invalid values.

    Warns
    -----
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.

    Notes
    -----
    The ideal chemical potential is given by [1]_:

    .. math::
        \chi_a = I_{1/2}(\beta \mu_a^{ideal})

    where :math:`\chi` is the degeneracy parameter, :math:`I_{1/2}` is the
    Fermi integral with order 1/2, :math:`\beta` is the inverse thermal
    energy :math:`\beta = 1/(k_B T)`, and :math:`\mu_a^{ideal}`
    is the ideal chemical potential.

    The definition for the ideal chemical potential is implicit, so it must
    be obtained numerically by solving for the Fermi integral for values
    of chemical potential approaching the degeneracy parameter. Since values
    returned from the Fermi_integral are complex, a nonlinear
    Levenberg-Marquardt least squares method is used to iteratively approach
    a value of :math:`\mu` which minimizes
    :math:`I_{1/2}(\beta \mu_a^{ideal}) - \chi_a`

    This function returns :math:`\beta \mu^{ideal}` the dimensionless
    ideal chemical potential.

    Warning: at present this function is limited to relatively small
    arguments due to limitations in the `~mpmath` package's implementation
    of `~mpmath.polylog`, which PlasmaPy uses in calculating the Fermi
    integral.

    References
    ----------
    .. [1] Bonitz, Michael. Quantum kinetic theory. Stuttgart: Teubner, 1998.

    Example
    -------
    >>> from astropy import units as u
    >>> chemical_potential(n_e=1e21*u.cm**-3,T=11000*u.K)  # doctest: +SKIP
    <Quantity 2.00039985e-12>

    """

    raise NotImplementedError(
        "This function has been temporarily disabled due to a bug.\n"
        "Please refer to https://github.com/PlasmaPy/PlasmaPy/issues/726 \n"
        "and https://github.com/astropy/astropy/issues/9721 "
        "for progress in fixing it."
    )
    # deBroglie wavelength
    lambdaDB = thermal_deBroglie_wavelength(T)
    # degeneracy parameter
    degen = (n_e * lambdaDB ** 3).to(u.dimensionless_unscaled)

    def residual(params, data, eps_data):
        """Residual function for fitting parameters to Fermi_integral."""
        alpha = params["alpha"].value
        # note that alpha = mu / (k_B * T)
        model = mathematics.Fermi_integral(alpha, 0.5)
        complexResidue = (data - model) / eps_data
        return complexResidue.view(np.float)

    # setting parameters for fitting along with bounds
    alphaGuess = 1 * u.dimensionless_unscaled
    try:
        from lmfit import minimize, Parameters
    except (ImportError, ModuleNotFoundError) as e:
        from plasmapy.optional_deps import lmfit_import_error

        raise lmfit_import_error from e

    params = Parameters()
    params.add("alpha", value=alphaGuess, min=0.0)
    # calling minimize function from lmfit to fit by minimizing the residual
    data = np.array([degen])  # result of Fermi_integral - degen should be zero
    eps_data = np.array([1e-15])  # numerical error
    minFit = minimize(residual, params, args=(data, eps_data))
    beta_mu = minFit.params["alpha"].value * u.dimensionless_unscaled
    return beta_mu


# TODO: decorate with validate_quantities
# TODO: remove NotImplementedError and 'doctest: +SKIP' when the following issues are addressed...
#       https://github.com/PlasmaPy/PlasmaPy/issues/726
#       https://github.com/astropy/astropy/issues/9721
def _chemical_potential_interp(n_e, T):
    r"""
    Fitting formula for interpolating chemical potential between classical
    and quantum regimes.

    See [1]_, [2]_ for more information.

    Parameters
    ----------
    n_e: ~astropy.units.Quantity
        Electron number density.

    T : ~astropy.units.Quantity
        Temperature in units of temperature or energy.

    Returns
    -------
    beta_mu: ~astropy.units.Quantity
        The dimensionless chemical potential, which is a ratio of
        chemical potential energy to thermal kinetic energy.

    Raises
    ------
    TypeError
        If argument is not a ~astropy.units.Quantity.

    ~astropy.units.UnitConversionError
        If argument is in incorrect units.

    ValueError
        If argument contains invalid values.

    Warnings
    --------
    ~astropy.units.UnitsWarning
        If units are not provided, SI units are assumed.

    Notes
    -----
    The ideal chemical potential is given by [1]_:

    .. math::
        \frac{\mu}{k_B T_e} = - \frac{3}{2} \ln \Theta + \ln
        \frac{4}{3 \sqrt{\pi}} +
        \frac{A \Theta^{-b - 1} + B \Theta^{-(b + 1) / 2}}{1 + A \Theta^{-b}}

    where

    .. math::
        \Theta = \frac{k_B T_e}{E_F}

    is the degeneracy parameter, comparing the thermal energy to the Fermi
    energy, and the coefficients for the fitting formula
    are A=0.25945, B=0.0072, b=0.858.

    References
    ----------
    .. [1] Ichimaru, Statistical Plasma Physics Addison-Wesley,
       Reading, MA, 1991.

    .. [2] Gregori, G., et al. "Theoretical model of x-ray scattering as a
       dense matter probe." Physical Review E 67.2 (2003): 026412.

    Example
    -------
    >>> from astropy import units as u
    >>> _chemical_potential_interp(n_e=1e23*u.cm**-3, T=11000*u.K)  # doctest: +SKIP
    <Quantity 8.17649>

    """
    raise NotImplementedError(
        "This function has been temporarily disabled due to a bug.\n"
        "Please refer to https://github.com/PlasmaPy/PlasmaPy/issues/726 \n"
        "and https://github.com/astropy/astropy/issues/9721 "
        "for progress in fixing it."
    )
    A = 0.25945
    B = 0.072
    b = 0.858
    theta = k_B * T / Fermi_energy(n_e)
    term1 = -3 / 2 * np.log(theta)
    term2 = np.log(4 / (3 * np.sqrt(np.pi)))
    term3num = A * theta ** (-b - 1) + B * theta ** (-(b + 1) / 2)
    term3den = 1 + A * theta ** (-b)
    term3 = term3num / term3den
    beta_mu = term1 + term2 + term3
    return beta_mu.to(u.dimensionless_unscaled)
