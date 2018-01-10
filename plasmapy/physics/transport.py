"""Functions to calculate transport coefficients."""

from astropy import units
import numpy as np
import plasmapy.atomic as atomic
from plasmapy import utils
from plasmapy.utils.checks import check_quantity, _check_relativistic
from plasmapy.utils.exceptions import PhysicsError, PhysicsWarning
from plasmapy.constants import (m_p, m_e, c, mu0, k_B, e, eps0, pi, h, hbar)
from ..atomic import (ion_mass, charge_state)
from .parameters import (Debye_length, Hall_parameter,
                         collision_rate_electron_ion, collision_rate_ion_ion)
from inspect import stack
from copy import copy
import warnings


@utils.check_quantity({"T": {"units": units.K, "can_be_negative": False},
                       "n_e": {"units": units.m**-3}
                       })
def Coulomb_logarithm(T, n_e, particles, V=None):
    r"""
    Estimates the Coulomb logarithm.

    Parameters
    ----------

    T : Quantity
        Temperature in units of temperature or energy per particle,
        which is assumed to be equal for both the test particle and
        the target particle.

    n_e : Quantity
        The electron density in units convertible to per cubic meter.

    particles : tuple
        A tuple containing string representations of the test particle
        (listed first) and the target particle (listed second).

    V : Quantity, optional
        The relative velocity between particles.  If not provided,
        thermal velocity is assumed: :math:`\mu V^2 \sim 3 k_B T`
        where `mu` is the reduced mass.

    Returns
    -------
    lnLambda : float or numpy.ndarray
        An estimate of the Coulomb logarithm that is accurate to
        roughly its reciprocal.

    Raises
    ------
    ValueError
        If the mass or charge of either particle cannot be found, or
        any of the inputs contain incorrect values.

    UnitConversionError
        If the units on any of the inputs are incorrect.

    UserWarning
        If the inputted velocity is greater than 80% of the speed of
        light.

    TypeError
        If the n_e, T, or V are not Quantities.

    Notes
    -----
    The Coulomb logarithm is given by

    .. math::
        \ln{\Lambda} \equiv \ln\left( \frac{b_{max}}{b_{min}} \right)

    where :math:`b_{min}` and :math:`b_{max}` are the inner and outer
    impact parameters for Coulomb collisions [1]_.

    The outer impact parameter is given by the Debye length:
    :math:`b_{min} = \lambda_D` which is a function of electron
    temperature and electron density.  At distances greater than the
    Debye length, electric fields from other particles will be
    screened out due to electrons rearranging themselves.

    The choice of inner impact parameter is more controversial. There
    are two main possibilities.  The first possibility is that the
    inner impact parameter corresponds to a deflection angle of 90
    degrees.  The second possibility is that the inner impact
    parameter is a de Broglie wavelength, :math:`\lambda_B`
    corresponding to the reduced mass of the two particles and the
    relative velocity between collisions.  This function uses the
    standard practice of choosing the inner impact parameter to be the
    maximum of these two possibilities.  Some inconsistencies exist in
    the literature on how to define the inner impact parameter [2]_.

    Errors associated with the Coulomb logarithm are of order its
    inverse. If the Coulomb logarithm is of order unity, then the
    assumptions made in the standard analysis of Coulomb collisions
    are invalid.

    Examples
    --------
    >>> from astropy import units as u
    >>> n = 1e19*units.m**-3
    >>> T = 1e6*units.K
    >>> particles = ('e', 'p')
    >>> Coulomb_logarithm(T, n, particles)
    14.748259780491056
    >>> Coulomb_logarithm(T, n, particles, 1e6*u.m/u.s)
    11.363478214139432

    References
    ----------
    .. [1] Physics of Fully Ionized Gases, L. Spitzer (1962)

    .. [2] Comparison of Coulomb Collision Rates in the Plasma Physics
       and Magnetically Confined Fusion Literature, W. Fundamenski and
       O.E. Garcia, EFDA–JET–R(07)01
       (http://www.euro-fusionscipub.org/wp-content/uploads/2014/11/EFDR07001.pdf)

    """

    if not isinstance(particles, (list, tuple)) or len(particles) != 2:
        raise ValueError("The third input of Coulomb_logarithm must be a list "
                         "or tuple containing representations of two charged "
                         "particles.")

    masses = np.zeros(2) * units.kg
    charges = np.zeros(2) * units.C

    for particle, i in zip(particles, range(2)):

        try:
            masses[i] = atomic.ion_mass(particles[i])
        except Exception:
            raise ValueError("Unable to find mass of particle: " +
                             str(particles[i]) + " in Coulomb_logarithm.")

        try:
            charges[i] = np.abs(e * atomic.charge_state(particles[i]))
        except Exception:
            raise ValueError("Unable to find charge of particle: " +
                             str(particles[i]) + " in Coulomb_logarithm.")

    reduced_mass = masses[0] * masses[1] / (masses[0] + masses[1])

    # The outer impact parameter is the Debye length.  At distances
    # greater than the Debye length, the electrostatic potential of a
    # single particle is screened out by the electrostatic potentials
    # of other particles.  Past this distance, the electric fields of
    # individual particles do not affect each other much.  This
    # expression neglects screening by heavier ions.

    T = T.to(units.K, equivalencies=units.temperature_energy())

    b_max = Debye_length(T, n_e)

    # The choice of inner impact parameter is more controversial.
    # There are two broad possibilities, and the expressions in the
    # literature often differ by factors of order unity or by
    # interchanging the reduced mass with the test particle mass.

    # The relative velocity is a source of uncertainty.  It is
    # reasonable to make an assumption relating the thermal energy to
    # the kinetic energy: reduced_mass*velocity**2 is approximately
    # equal to 3*k_B*T.

    # If no relative velocity is inputted, then we make an assumption
    # that relates the thermal energy to the kinetic energy:
    # reduced_mass*velocity**2 is approximately equal to 3*k_B*T.

    if V is None:
        V = np.sqrt(3 * k_B * T / reduced_mass)
    utils._check_relativistic(V, 'Coulomb_logarithm', betafrac=0.8)

    # The first possibility is that the inner impact parameter
    # corresponds to a deflection of 90 degrees, which is valid when
    # classical effects dominate.

    b_perp = charges[0] * charges[1] / (4 * pi * eps0 * reduced_mass * V**2)

    # The second possibility is that the inner impact parameter is a
    # de Broglie wavelength.  There remains some ambiguity as to which
    # mass to choose to go into the de Broglie wavelength calculation.
    # Here we use the reduced mass, which will be of the same order as
    # mass of the smaller particle and thus the longer de Broglie
    # wavelength.

    b_deBroglie = hbar / (2 * reduced_mass * V)

    # Coulomb-style collisions will not happen for impact parameters
    # shorter than either of these two impact parameters, so we choose
    # the larger of these two possibilities.

    b_min = np.zeros_like(b_perp)

    for i in range(b_min.size):

        if b_perp.flat[i] > b_deBroglie.flat[i]:
            b_min.flat[i] = b_perp.flat[i]
        else:
            b_min.flat[i] = b_deBroglie.flat[i]

    # Now that we know how many approximations have to go into plasma
    # transport theory, we shall celebrate by returning the Coulomb
    # logarithm.

    ln_Lambda = np.log(b_max / b_min)
    ln_Lambda = ln_Lambda.to(units.dimensionless_unscaled).value

    return ln_Lambda


class classical_transport:
    r"""
    Braginskii-esque classical transport coefficients

    Notes
    -----
    Classical transport theory is derived by using kinetic theory to close the
    plasma two-fluid (electron and ion fluid) equations in the collisional
    limit. This function uses fitting functions from literature to calculate
    the transport coefficients, which are the resistivity, thermoelectric
    conductivity, thermal conductivity, and viscosity, which can be used
    to close the two-fluid equations.

    Note well the assumptions in the derivation of classical transport.
    Turbulent transport does not dominate, and the velocity distribution
    function is close to Maxwellian (no extremely strong gradients),
    which is equivalent to the following conditions:

    * collisional frequency >> gyrofrequency
    * collisional mean free path << gradient scale length along field
    * gyroradius << gradient scale length perpendicular to field

    When classical transport is not valid, e.g. due to the presence of strong
    gradients or turbulent transport, the transport is significantly increased
    by these other effects. Thus, classical transport often serves as a lower
    bound on the losses / transport encountered in a plasma.

    Parameters
    ----------
    T_e : Quantity
        Temperature in units of temperature or energy per particle

    n_e : Quantity
        The number density in units convertible to per cubic meter.

    T_i : Quantity
        Temperature in units of temperature or energy per particle

    n_i : Quantity
        The number density in units convertible to per cubic meter.

    ion_particle : string
        Representation of the ion species (e.g., 'p' for protons,
        'e' for electrons, 'D+' for deuterium, or 'He-4 +1' for singly
        ionized helium-4). If no charge state information is provided,
        then the particles are assumed to be singly charged.

    Z : integer or np.inf, optional
        The ion charge state. Overrides particle charge state if included.
        Different theories support different values of Z.
        For the original Braginskii model, Z can be any of [1,2,3,4,infinity].
        The Ji-Held model supports arbitrary Z. Average ionization state
        Zbar can be input using this input and model, but doing so may neglect
        effects caused by multiple ion populations.

    B : Quantity, optional
        The magnetic field strength in units convertible to Tesla. Defaults
        to zero.

    model: string
        Indication of whose formulation from literature to use. Allowed values
        are:
        'Braginskii',
        'Spitzer-Harm',
        'Epperlein-Haines' (not yet implemented),
        'Ji-Held'.
        See refs [1]_, [2]_, [3]_ and [4]_.

    field_orientation : string
        Either of 'parallel', 'par', 'perpendicular', 'perp', 'cross', or
        'all', indicating the cardinal orientation of the magnetic field with
        respect to the transport direction of interest. Note that 'perp' refers
        to transport perpendicular to the field direction in the direction of
        the temperature gradient, while 'cross' refers to the B X grad(T)
        direction. The option 'all' will return a numpy array of all three,
        np.array((par, perp, cross)).

    coulomb_log_ei: float or dimensionless Quantity, optional
        Force a particular value to be used for the Coulomb logarithm. If
        None, the PlasmaPy function Coulomb_Logarithm() will be used. Useful
        for comparing calculations.

    V_ei: Quantity, optional
       Supplied to coulomb_logarithm() function, not otherwise used.
       The relative velocity between particles.  If not provided,
       thermal velocity is assumed: :math:`\mu V^2 \sim 3 k_B T`
       where `mu` is the reduced mass.

    coulomb_log_ii: float or dimensionless Quantity, optional
        Force a particular value to be used for the Coulomb logarithm. If
        None, the PlasmaPy function Coulomb_Logarithm() will be used. Useful
        for comparing calculations.

    V_ii: Quantity, optional
       Supplied to coulomb_logarithm() function, not otherwise used.
       The relative velocity between particles.  If not provided,
       thermal velocity is assumed: :math:`\mu V^2 \sim 3 k_B T`
       where `mu` is the reduced mass.

    hall_e: float or dimensionless Quantity, optional
        Force a particular value to be used for the e Hall parameter. If
        None, the PlasmaPy function Hall_parameter() will be used. Useful
        for comparing calculations.

    hall_i: float or dimensionless Quantity, optional
        Force a particular value to be used for the ion Hall parameter. If
        None, the PlasmaPy function Hall_parameter() will be used. Useful
        for comparing calculations.

    mu: optional, float or dimensionless Quantity
        Ji-Held model only, may be used to include ion-electron effects
        on the ion transport coefficients. Defaults to m_e / m_i.
        Set to zero to disable these effects.

    theta: optional, float or dimensionless Quantity
        theta = T_e / T_i
        Ji-Held model only, may be used to include ion-electron effects
        on the ion transport coefficients. Defaults to T_e / T_i. Only
        has effect if mu is non-zero.


    Raises
    ------
    ValueError
        On incorrect or unknown values of arguments.
    plasmapy.utils.PhysicsError
        If input or calculated values for Coulomb logarithms are nonphysical.

    References
    ----------
    .. [1] S.I. Braginskii (1965)
    .. [2] Spitzer-Harm (1953)
    .. [3] Epperlein-Haines (1986)
    .. [4] Ji-Held (2013)
    """

    @utils.check_quantity({"T_e": {"units": units.K, "can_be_negative": False},
                           "n_e": {"units": units.m**-3},
                           "T_i": {"units": units.K, "can_be_negative": False},
                           "n_i": {"units": units.m**-3},
                           })
    def __init__(self, T_e, n_e, T_i, n_i, ion_particle, m_i=None, Z=None,
                 B=0.0 * units.T, model='Braginskii',
                 field_orientation='parallel',
                 coulomb_log_ei=None, V_ei=None,
                 coulomb_log_ii=None, V_ii=None,
                 hall_e=None, hall_i=None,
                 mu=None, theta=None):
        # check the model
        self.model = model.lower()  # string inputs should be case insensitive
        valid_models = ['braginskii',
                        'spitzer',
                        'spitzer-harm',
                        'ji-held',
                        ]
        is_valid_model = self.model in valid_models
        if not is_valid_model:
            raise ValueError(f"Unknown transport model '{self.model}'")

        # check the field orientation
        self.field_orientation = field_orientation.lower()
        valid_fields = ['parallel', 'par',
                        'perpendicular', 'perp',
                        'cross',
                        'all',
                        ]
        is_valid_field = self.field_orientation in valid_fields
        if not is_valid_field:
            raise ValueError(
                f"Unknown field orientation '{self.field_orientation}'")

        # density and temperature units have already been checked by decorator
        # so just convert
        self.T_e = T_e.to(units.K, equivalencies=units.temperature_energy())
        self.T_i = T_i.to(units.K, equivalencies=units.temperature_energy())
        self.n_e = n_e.to(units.m**-3)
        self.n_i = n_i.to(units.m**-3)

        # get ion mass and charge state
        if m_i is None:
            try:
                self.m_i = atomic.ion_mass(ion_particle)
            except Exception:
                raise ValueError("Unable to find mass of particle: " +
                                 str(ion_particle) + " in classical_transport")
        else:
            self.m_i = m_i.to(units.kg)
        if Z is None:
            try:
                self.Z = atomic.charge_state(ion_particle)
            except Exception:
                raise ValueError(
                    "Unable to find charge of particle: " +
                    str(ion_particle) +
                    " in classical_transport.")
        else:
            # red alert: the user has input a Z
            self.Z = Z
            # if it's not available for a particular model, they'll complain
            # later
            if Z < 0:
                raise ValueError("Z not allowed to be negative")

        # decide on the particle string for the electrons
        self.e_particle = 'e'
        self.ion_particle = ion_particle

        # save other arguments
        self.B = B
        self.V_ei = V_ei
        self.V_ii = V_ii

        # calculate Coulomb logs if not forced in input
        if coulomb_log_ei:
            self.coulomb_log_ei = coulomb_log_ei
        else:
            self.coulomb_log_ei = Coulomb_logarithm(
                T_e, n_e, [self.e_particle, ion_particle], V_ei)
        if self.coulomb_log_ei < 4:
            warnings.warn(
                f"coulomb_log_ei is {self.coulomb_log_ei}, you might "
                "have strong coupling effects",
                utils.PhysicsWarning)
        if self.coulomb_log_ei < 1:
            raise utils.PhysicsError(
                f"coulomb_log_ei is {self.coulomb_log_ei}, less than 1")

        if coulomb_log_ii is not None:
            self.coulomb_log_ii = coulomb_log_ii
        else:
            # TODO make comment below more clear?
            self.coulomb_log_ii = Coulomb_logarithm(T_i, n_e,  # not a typo, but worry
                                                    [ion_particle, ion_particle],
                                                    V_ii)

        if self.coulomb_log_ii < 4:
            warnings.warn(
                f"coulomb_log_ii is {self.coulomb_log_ii}, you might "
                "have strong coupling effects",
                utils.PhysicsWarning)
        if self.coulomb_log_ii < 1:
            raise utils.PhysicsError(
                f"coulomb_log_ii is {self.coulomb_log_ii}, less than 1")

        # calculate Hall parameters if not forced in input
        if hall_e is not None:
            self.hall_e = hall_e
        else:
            self.hall_e = Hall_parameter(
                n_e,
                T_e,
                B,
                self.e_particle,
                ion_particle,
                coulomb_log_ei,
                V_ei)
        if hall_i is not None:
            self.hall_i = hall_i
        else:
            self.hall_i = Hall_parameter(
                n_i, T_i, B, ion_particle, ion_particle, coulomb_log_ii, V_ii)
        # set up the ion non-dimensional coefficients for the Ji-Held model
        if mu is not None:
            self.mu = mu
        else:
            self.mu = m_e / self.m_i
        if theta is not None:
            self.theta = theta
        else:
            self.theta = self.T_e / self.T_i

    def resistivity(self):
        return resistivity(
            T_e=self.T_e,
            n_e=self.n_e,
            ion_particle=self.ion_particle,
            e_particle=self.e_particle,
            Z=self.Z,
            B=self.B,
            model=self.model,
            field_orientation=self.field_orientation,
            coulomb_log_ei=self.coulomb_log_ei,
            V_ei=self.V_ei,
            hall_e=self.hall_e,
        )

    def thermoelectric_conductivity(self):
        return thermoelectric_conductivity(
            T_e=self.T_e,
            n_e=self.n_e,
            e_particle=self.e_particle,
            Z=self.Z,
            B=self.B,
            model=self.model,
            field_orientation=self.field_orientation,
            coulomb_log_ei=self.coulomb_log_ei,
            V_ei=self.V_ei,
            hall_e=self.hall_e,
        )

    def ion_thermal_conductivity(self):
        return ion_thermal_conductivity(
            T_i=self.T_i,
            n_i=self.n_i,
            ion_particle=self.ion_particle,
            Z=self.Z,
            B=self.B,
            model=self.model,
            field_orientation=self.field_orientation,
            coulomb_log_ii=self.coulomb_log_ii,
            V_ii=self.V_ii,
            hall_i=self.hall_i,
            mu=self.mu,
            theta=self.theta,
        )

    def electron_thermal_conductivity(self):
        return electron_thermal_conductivity(
            T_e=self.T_e,
            n_e=self.n_e,
            ion_particle=self.ion_particle,
            e_particle=self.e_particle,
            Z=self.Z,
            B=self.B,
            model=self.model,
            field_orientation=self.field_orientation,
            coulomb_log_ei=self.coulomb_log_ei,
            V_ei=self.V_ei,
            hall_e=self.hall_e,
            mu=self.mu,
            theta=self.theta,
        )

    def ion_viscosity(self):
        return ion_viscosity(
            T_i=self.T_i,
            n_i=self.n_i,
            ion_particle=self.ion_particle,
            Z=self.Z,
            B=self.B,
            model=self.model,
            field_orientation=self.field_orientation,
            coulomb_log_ii=self.coulomb_log_ii,
            V_ii=self.V_ii,
            hall_i=self.hall_i,
            mu=self.mu,
            theta=self.theta,
        )

    def electron_viscosity(self):
        return electron_viscosity(
            T_e=self.T_e,
            n_e=self.n_e,
            ion_particle=self.ion_particle,
            e_particle=self.e_particle,
            Z=self.Z,
            B=self.B,
            model=self.model,
            field_orientation=self.field_orientation,
            coulomb_log_ei=self.coulomb_log_ei,
            V_ei=self.V_ei,
            hall_e=self.hall_e,
            mu=self.mu,
            theta=self.theta,
        )


def resistivity(T_e, n_e, ion_particle, e_particle, Z=None, B=0.0,
                model='Braginskii', field_orientation='parallel',
                coulomb_log_ei=None, V_ei=None, hall_e=None):
    r"""TODO"""
#    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())
#    n_e = n_e.to(units.m**-3)
#    model = model.lower()
#    field_orientation = field_orientation.lower()
#    e_particle = 'e'
#    if hall is None:
#        hall = Hall_parameter(n_e, T_e, B, e_particle, ion_particle,
#                              coulomb_log, V)
#    if Z is None:
#        Z = atomic.charge_state(ion_particle)
    alpha_hat = _nondim_resisitivity(hall_e, Z, e_particle, model,
                                     field_orientation)
    tau_e = 1 / collision_rate_electron_ion(T_e, n_e, ion_particle,
                                            coulomb_log_ei, V_ei)

#    e_cgs = e.value * c.value * 10
#    alpha = alpha_hat / (e_cgs**2 * n_e.value/1e6 * \
#        tau_e.value / (m_e.value*1000)) * 1/(4*np.pi*eps0.value)
#    return alpha * (units.ohm * units.m)
    alpha = alpha_hat / (n_e * e**2 * tau_e / m_e)
    return alpha.to(units.ohm * units.m)


def thermoelectric_conductivity(T_e, n_e, e_particle, Z=None,
                                B=0.0, model='Braginskii',
                                field_orientation='parallel',
                                coulomb_log_ei=None, V_ei=None, hall_e=None):
    r"""TODO"""
#    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())
#    n_e = n_e.to(units.m**-3)
#    model = model.lower()
#    field_orientation = field_orientation.lower()
#    e_particle = 'e'
#    if hall is None:
#        hall = Hall_parameter(n_e, T_e, B, e_particle, ion_particle,
#                              coulomb_log, V)
#    if Z is None:
#        Z = atomic.charge_state(ion_particle)
    beta_hat = _nondim_te_conductivity(hall_e, Z, e_particle, model,
                                       field_orientation)
    beta = beta_hat * units.s / units.s  # yay! already dimensionless
    return beta


def ion_thermal_conductivity(T_i, n_i, ion_particle, Z=None, B=0.0,
                             model='Braginskii', field_orientation='parallel',
                             coulomb_log_ii=None, V_ii=None, hall_i=None,
                             mu=None, theta=None):
    r"""Braginskii-esque two-fluid plasma thermal conductivities

    Notes
    -----
    See comments on the use of classical transport theory in the
    classical_transport function.

    """
#    T_i = T_i.to(units.K, equivalencies=units.temperature_energy())
#    n_i = n_i.to(units.m**-3)
#    model = model.lower()
#    field_orientation = field_orientation.lower()
#    if hall is None:
#        hall = Hall_parameter(n_i, T_i, B, ion_particle, ion_particle,
#                              coulomb_log, V)
#    if Z is None:
#        Z = atomic.charge_state(ion_particle)
    m_i = atomic.ion_mass(ion_particle)
    kappa_hat = _nondim_thermal_conductivity(hall_i, Z, ion_particle, model,
                                             field_orientation, mu, theta)
    tau_i = 1 / collision_rate_ion_ion(T_i, n_i, ion_particle,
                                       coulomb_log_ii, V_ii)
    kappa = kappa_hat * (n_i * k_B**2 * T_i * tau_i / m_i)
    return kappa.to(units.W / units.m / units.K)


def electron_thermal_conductivity(T_e, n_e, ion_particle, e_particle,
                                  Z=None, B=0.0, model='Braginskii',
                                  field_orientation='parallel',
                                  coulomb_log_ei=None, V_ei=None, hall_e=None,
                                  mu=None, theta=None):
    #    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())
    #    n_e = n_e.to(units.m**-3)
    #    model = model.lower()
    #    field_orientation = field_orientation.lower()
    #    e_particle = 'e'
    #    if hall is None:
    #        hall = Hall_parameter(n_e, T_e, B, e_particle,
    #                              ion_particle, coulomb_log, V)
    #    if Z is None:
    #        Z = atomic.charge_state(ion_particle)
    kappa_hat = _nondim_thermal_conductivity(hall_e, Z, e_particle, model,
                                             field_orientation, mu, theta)
    tau_e = 1 / collision_rate_electron_ion(T_e, n_e, ion_particle,
                                            coulomb_log_ei, V_ei)
    kappa = kappa_hat * (n_e * k_B**2 * T_e * tau_e / m_e)
    return kappa.to(units.W / units.m / units.K)


def ion_viscosity(T_i, n_i, ion_particle, Z=None, B=0.0, model='Braginskii',
                  field_orientation='parallel', coulomb_log_ii=None, V_ii=None,
                  hall_i=None, mu=None, theta=None):
    r"""TODO"""
#    T_i = T_i.to(units.K, equivalencies=units.temperature_energy())
#    n_i = n_i.to(units.m**-3)
#    model = model.lower()
#    field_orientation = field_orientation.lower()
#    if hall is None:
#        hall = Hall_parameter(n_i, T_i, B, ion_particle,
#                              ion_particle, coulomb_log, V)
#    if Z is None:
#        Z = atomic.charge_state(ion_particle)
    eta_hat = _nondim_viscosity(hall_i, Z, ion_particle, model,
                                field_orientation, mu, theta)
    tau_i = 1 / collision_rate_ion_ion(T_i, n_i, ion_particle,
                                       coulomb_log_ii, V_ii)
    if np.isclose(hall_i, 0, rtol=1e-8):
        eta1 = (eta_hat[0] * (n_i * k_B * T_i * tau_i),
                eta_hat[1] * (n_i * k_B * T_i * tau_i),
                eta_hat[2] * (n_i * k_B * T_i * tau_i),
                eta_hat[3] * (n_i * k_B * T_i * tau_i),
                eta_hat[4] * (n_i * k_B * T_i * tau_i))
    else:
        eta1 = (eta_hat[0] * (n_i * k_B * T_i * tau_i),
                eta_hat[1] * (n_i * k_B * T_i * tau_i) / hall_i**2,
                eta_hat[2] * (n_i * k_B * T_i * tau_i) / hall_i**2,
                eta_hat[3] * (n_i * k_B * T_i * tau_i) / hall_i,
                eta_hat[4] * (n_i * k_B * T_i * tau_i) / hall_i)
    if eta1[0].unit == eta1[2].unit and eta1[2].unit == eta1[4].unit:
        unit_val = eta1[0].unit
        eta = (np.array((eta1[0].value,
                         eta1[1].value,
                         eta1[2].value,
                         eta1[3].value,
                         eta1[4].value)) * unit_val).to(units.Pa * units.s)
    return eta


def electron_viscosity(T_e, n_e, ion_particle, e_particle, Z=None, B=0.0,
                       model='Braginskii', field_orientation='parallel',
                       coulomb_log_ei=None, V_ei=None, hall_e=None, mu=None,
                       theta=None):
    r"""TODO"""
#    T_e = T_e.to(units.K, equivalencies=units.temperature_energy())
#    n_e = n_e.to(units.m**-3)
#    model = model.lower()
#    field_orientation = field_orientation.lower()
#    e_particle = 'e'
#    if hall is None:
#        hall = Hall_parameter(n_e, T_e, B, e_particle,
#                              ion_particle, coulomb_log, V)
#    if Z is None:
#        Z = atomic.charge_state(ion_particle)
    eta_hat = _nondim_viscosity(hall_e, Z, e_particle, model,
                                field_orientation, mu, theta)
    tau_e = 1 / collision_rate_electron_ion(T_e, n_e, ion_particle,
                                            coulomb_log_ei, V_ei)
    if np.isclose(hall_e, 0, rtol=1e-8):
        eta1 = (eta_hat[0] * (n_e * k_B * T_e * tau_e),
                eta_hat[1] * (n_e * k_B * T_e * tau_e),
                eta_hat[2] * (n_e * k_B * T_e * tau_e),
                eta_hat[3] * (n_e * k_B * T_e * tau_e),
                eta_hat[4] * (n_e * k_B * T_e * tau_e))
    else:
        eta1 = (eta_hat[0] * (n_e * k_B * T_e * tau_e),
                eta_hat[1] * (n_e * k_B * T_e * tau_e) / hall_e**2,
                eta_hat[2] * (n_e * k_B * T_e * tau_e) / hall_e**2,
                eta_hat[3] * (n_e * k_B * T_e * tau_e) / hall_e,
                eta_hat[4] * (n_e * k_B * T_e * tau_e) / hall_e)
    if eta1[0].unit == eta1[2].unit and eta1[2].unit == eta1[4].unit:
        unit_val = eta1[0].unit
        eta = (np.array((eta1[0].value,
                         eta1[1].value,
                         eta1[2].value,
                         eta1[3].value,
                         eta1[4].value)) * unit_val).to(units.Pa * units.s)
    return eta


#
#
#
#
#
#
#
#                 no units allowed beyond this point!
#
#
#
#
#
#
#
#


def _nondim_thermal_conductivity(hall, Z, particle, model, field_orientation,
                                 mu=None, theta=None):
    """TODO"""
    from plasmapy.atomic.atomic import _is_electron
    if _is_electron(particle):
        if model == 'spitzer-harm' or model == 'spitzer':
            kappa_hat = _nondim_tc_e_spitzer(Z)
        elif model == 'braginskii':
            kappa_hat = _nondim_tc_e_braginskii(hall, Z, field_orientation)
        elif model == 'ji-held':
            kappa_hat = _nondim_tc_e_ji_held(hall, Z, field_orientation)
        else:
            raise ValueError(f"Unrecognized model '{model}' in "
                             "_nondim_thermal_conductivity")
    else:
        if model == 'braginskii':
            kappa_hat = _nondim_tc_i_braginskii(hall, field_orientation)
        elif model == 'ji-held':
            kappa_hat = _nondim_tc_i_ji_held(hall, Z, mu, theta,
                                             field_orientation)
        else:
            raise ValueError(f"Unrecognized model '{model}' in "
                             "_nondim_thermal_conductivity")
    return kappa_hat


def _nondim_viscosity(hall, Z, particle, model, field_orientation,
                      mu=None, theta=None):
    """TODO"""

    from plasmapy.atomic.atomic import _is_electron
    if _is_electron(particle):
        if model == 'braginskii':
            eta_hat = _nondim_visc_e_braginskii(hall, Z)
        elif model == 'ji-held':
            eta_hat = _nondim_visc_e_ji_held(hall, Z)
        else:
            raise ValueError(f"Unrecognized model '{model}' in "
                             "_nondim_viscosity")
    else:
        if model == 'braginskii':
            eta_hat = _nondim_visc_i_braginskii(hall)
        elif model == 'ji-held':
            eta_hat = _nondim_visc_i_ji_held(hall, Z, mu, theta)
        else:
            raise ValueError(f"Unrecognized model '{model}' in "
                             "_nondim_viscosity")
    return eta_hat


def _nondim_resisitivity(hall, Z, particle, model, field_orientation):
    """TODO"""

    if model == 'spitzer-harm' or model == 'spitzer':
        alpha_hat = _nondim_resist_spitzer(Z)
    elif model == 'braginskii':
        alpha_hat = _nondim_resist_braginskii(hall, Z, field_orientation)
    elif model == 'ji-held':
        alpha_hat = _nondim_resist_ji_held(hall, Z, field_orientation)
    else:
        raise ValueError(f"Unrecognized model '{model}' in "
                         "_nondim_resisitivity")
    return alpha_hat


def _nondim_te_conductivity(hall, Z, particle, model, field_orientation):
    """TODO"""

    if model == 'spitzer-harm' or model == 'spitzer':
        beta_hat = _nondim_tec_spitzer(Z)
    elif model == 'braginskii':
        beta_hat = _nondim_tec_braginskii(hall, Z, field_orientation)
    elif model == 'ji-held':
        beta_hat = _nondim_tec_ji_held(hall, Z, field_orientation)
    else:
        raise ValueError(f"Unrecognized model '{model}' in "
                         "_nondim_te_conductivity")
    return beta_hat


def check_Z(allowed_Z, Z):
    """TODO"""
    # first, determine if arbitrary Z values are allowed in the theory
    arbitrary_Z_allowed = False
    the_arbitrary_idx = np.nan
    for idx, allowed_Z_val in enumerate(allowed_Z):
        if allowed_Z_val == 'arbitrary':
            arbitrary_Z_allowed = True
            the_arbitrary_idx = idx
    # next, search the allowed_Z for a match to the current Z
    Z_idx = np.nan
    for idx, allowed_Z_val in enumerate(allowed_Z):
        if Z == allowed_Z_val:
            Z_idx = idx
    # at this point we have looped through allowed_Z and either found a match
    # or not. If we haven't found a match and arbitrary Z aren't allowed, break
    if np.isnan(Z_idx) and not arbitrary_Z_allowed:
        raise utils.PhysicsError(f"{Z} is not an allowed Z value")
    elif np.isnan(Z_idx):  # allowed arbitrary Z
        # return a Z_idx pointing to the 'arbitrary'
        Z_idx = the_arbitrary_idx
    else:  # allowed Z
        pass
    # we have got the Z_idx we want. return
    return Z_idx


def _get_spitzer_harm_coeffs(Z):
    allowed_Z = [1, 2, 4, 16, np.inf]
    Z_idx = check_Z(allowed_Z, Z)
    gamma_E = [0.5816, 0.6833, 0.7849, 0.9225, 1.0000]
    gamma_T = [0.2727, 0.4137, 0.5714, 0.8279, 1.0000]
    delta_E = [0.4652, 0.5787, 0.7043, 0.8870, 1.0000]
    delta_T = [0.2252, 0.3563, 0.5133, 0.7907, 1.0000]
    return (gamma_E[Z_idx], gamma_T[Z_idx], delta_E[Z_idx], delta_T[Z_idx])


def _nondim_tc_e_spitzer(Z):
    """TODO"""
    (gamma_E, gamma_T, delta_E, delta_T) = _get_spitzer_harm_coeffs(Z)
    kappa = (64 / np.pi) * delta_T * \
            (5 / 3 - (gamma_T * delta_E) / (delta_T * gamma_E))
    return kappa


def _nondim_resist_spitzer(Z):
    """TODO"""
    (gamma_E, gamma_T, delta_E, delta_T) = _get_spitzer_harm_coeffs(Z)
    alpha = (3 * np.pi / 32) * (1 / gamma_E)
#    alpha = 0.5064
    return alpha


def _nondim_tec_spitzer(Z):
    """TODO"""
    (gamma_E, gamma_T, delta_E, delta_T) = _get_spitzer_harm_coeffs(Z)
    beta = 5 / 2 * (8 / 5 * (delta_E / gamma_E) - 1)
#    beta = 0.703
    return beta


def _nondim_tc_e_braginskii(hall, Z, field_orientation):
    """TODO"""

    allowed_Z = [1, 2, 3, 4, np.inf]
    Z_idx = check_Z(allowed_Z, Z)

    delta_0 = [3.7703, 1.0465, 0.5814, 0.4106, 0.0961]
    delta_1 = [14.79, 10.80, 9.618, 9.055, 7.482]
    gamma_1_prime = [4.664, 3.957, 3.721, 3.604, 3.25]
    gamma_0_prime = [11.92, 5.118, 3.525, 2.841, 1.20]
    gamma_1_doubleprime = [2.500, 2.500, 2.500, 2.500, 2.500]
    gamma_0_doubleprime = [21.67, 15.37, 13.53, 12.65, 10.23]

    gamma_0 = gamma_0_prime[Z_idx] / delta_0[Z_idx]
    Delta = hall**4 + delta_1[Z_idx] * hall**2 + delta_0[Z_idx]
    kappa_par = gamma_0
    if field_orientation == 'parallel' or field_orientation == 'par':
        return kappa_par

    kappa_perp = (gamma_1_prime[Z_idx] * hall **
                  2 + gamma_0_prime[Z_idx]) / Delta
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return kappa_perp

    kappa_cross = (gamma_1_doubleprime[Z_idx] * hall**3 +
                   gamma_0_doubleprime[Z_idx] * hall) / Delta
    if field_orientation == 'cross':
        return kappa_cross

    if field_orientation == 'all':
        return np.array((kappa_par, kappa_perp, kappa_cross))


def _nondim_tc_i_braginskii(hall, field_orientation):
    """TODO"""

    kappa_par_coeff_0 = 3.906
    kappa_par = kappa_par_coeff_0
    if field_orientation == 'parallel' or field_orientation == 'par':
        return kappa_par

    kappa_perp_coeff_2 = 2.0
    kappa_perp_coeff_0 = 2.645
    delta_1 = 2.70
    delta_0 = 0.677
    Delta = hall**4 + delta_1 * hall**2 + delta_0
    kappa_perp = (kappa_perp_coeff_2 * hall**2 +
                  kappa_perp_coeff_0) / Delta
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return kappa_perp

    kappa_cross_coeff_3 = 2.5
    kappa_cross_coeff_1 = 4.65
    kappa_cross = (kappa_cross_coeff_3 * hall**3 +
                   kappa_cross_coeff_1 * hall) / Delta
    if field_orientation == 'cross':
        return kappa_cross

    if field_orientation == 'all':
        return np.array((kappa_par, kappa_perp, kappa_cross))


def _nondim_visc_e_braginskii(hall, Z):
    """TODO"""
    allowed_Z = [1]
    check_Z(allowed_Z, Z)
    eta_prime_0 = 0.733
    eta_doubleprime_2 = 2.05
    eta_doubleprime_0 = 8.50
    eta_tripleprime_2 = 1.0
    eta_tripleprime_0 = 7.91
    delta_1 = 13.8
    delta_0 = 11.6
    eta_0_e = eta_prime_0

    def f_eta_2(hall):
        Delta = hall**4 + delta_1 * hall**2 + delta_0
        return (eta_doubleprime_2 * hall**2 + eta_doubleprime_0) / Delta
    eta_2_e = f_eta_2(hall)
    eta_1_e = f_eta_2(2 * hall)

    def f_eta_4(hall):
        Delta = hall**4 + delta_1 * hall**2 + delta_0
        return (eta_tripleprime_2 * hall**3 +
                eta_tripleprime_0 * hall) / Delta
    eta_4_e = f_eta_4(hall)
    eta_3_e = f_eta_4(2 * hall)
    return np.array((eta_0_e, eta_1_e, eta_2_e, eta_3_e, eta_4_e))


def _nondim_visc_i_braginskii(hall):
    """TODO"""
    eta_prime_0 = 0.96
    eta_doubleprime_2 = 6 / 5
    eta_doubleprime_0 = 2.23
    eta_tripleprime_2 = 1.0
    eta_tripleprime_0 = 2.38
    delta_1 = 4.03
    delta_0 = 2.33
    eta_0_i = eta_prime_0

    def f_eta_2(hall):
        Delta = hall**4 + delta_1 * hall**2 + delta_0
        return (eta_doubleprime_2 * hall**2 + eta_doubleprime_0) / Delta
    eta_2_i = f_eta_2(hall)
    eta_1_i = f_eta_2(2 * hall)

    def f_eta_4(hall):
        Delta = hall**4 + delta_1 * hall**2 + delta_0
        return (eta_tripleprime_2 * hall**3 +
                eta_tripleprime_0 * hall) / Delta
    eta_4_i = f_eta_4(hall)
    eta_3_i = f_eta_4(2 * hall)
    return np.array((eta_0_i, eta_1_i, eta_2_i, eta_3_i, eta_4_i))


def _nondim_resist_braginskii(hall, Z, field_orientation):
    """TODO"""

    allowed_Z = [1, 2, 3, 4, np.inf]
    Z_idx = check_Z(allowed_Z, Z)

#    alpha_0 = 0.5129
    delta_0 = [3.7703, 1.0465, 0.5814, 0.4106, 0.0961]
    delta_1 = [14.79, 10.80, 9.618, 9.055, 7.482]
    alpha_1_prime = [6.416, 5.523, 5.226, 5.077, 4.63]
    alpha_0_prime = [1.837, 0.5956, 0.3515, 0.2566, 0.0678]
    alpha_1_doubleprime = [1.704, 1.704, 1.704, 1.704, 1.704]
    alpha_0_doubleprime = [0.7796, 0.3439, 0.2400, 0.1957, 0.0940]

    alpha_0 = 1 - alpha_0_prime[Z_idx] / delta_0[Z_idx]
    Delta = hall**4 + delta_1[Z_idx] * hall**2 + delta_0[Z_idx]
    alpha_par = alpha_0
    if field_orientation == 'parallel' or field_orientation == 'par':
        return alpha_par

    alpha_perp = (1 - (alpha_1_prime[Z_idx] * hall**2 +
                       alpha_0_prime[Z_idx]) / Delta)
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return alpha_perp

    alpha_cross = (alpha_1_doubleprime[Z_idx] * hall**3 +
                   alpha_0_doubleprime[Z_idx] * hall) / Delta
    if field_orientation == 'cross':
        return alpha_cross

    if field_orientation == 'all':
        return np.array((alpha_par, alpha_perp, alpha_cross))


def _nondim_tec_braginskii(hall, Z, field_orientation):
    """TODO"""

    allowed_Z = [1, 2, 3, 4, np.inf]
    Z_idx = check_Z(allowed_Z, Z)

    delta_0 = [3.7703, 1.0465, 0.5814, 0.4106, 0.0961]
    delta_1 = [14.79, 10.80, 9.618, 9.055, 7.482]
    beta_1_prime = [5.101, 4.450, 4.233, 4.124, 3.798]
    beta_0_prime = [2.681, 0.9473, 0.5905, 0.4478, 0.1461]
    beta_1_doubleprime = [1.5, 1.5, 1.5, 1.5, 1.5]
    beta_0_doubleprime = [3.053, 1.784, 1.442, 1.285, 0.877]

    Delta = hall**4 + delta_1[Z_idx] * hall**2 + delta_0[Z_idx]
    beta_0 = beta_0_prime[Z_idx] / delta_0[Z_idx]
#    beta_0 = 0.7110

    beta_par = beta_0
    if field_orientation == 'parallel' or field_orientation == 'par':
        return beta_par

    beta_perp = (beta_1_prime[Z_idx] * hall**2 +
                 beta_0_prime[Z_idx]) / Delta
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return beta_perp

    beta_cross = (beta_1_doubleprime[Z_idx] * hall**3 +
                  beta_0_doubleprime[Z_idx] * hall) / Delta
    if field_orientation == 'cross':
        return beta_cross

    if field_orientation == 'all':
        return np.array((beta_par, beta_perp, beta_cross))


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#               Abandon all hope, ye who enter here
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#


def _nondim_tc_e_ji_held(hall, Z, field_orientation):
    """TODO"""

    allowed_Z = [1, 2, 'arbitrary']
    Z_idx = check_Z(allowed_Z, Z)
    r = np.abs(Z * hall)

    def f_kappa_par_e(Z):
        numerator = 13.5 * Z**2 + 54.4 * Z + 25.2
        denominator = Z**3 + 8.35 * Z**2 + 15.2 * Z + 4.51
        return numerator / denominator

    def f_kappa_0(Z):
        numerator = 9.91 * Z**3 + 75.3 * Z**2 + 518 * Z + 333
        denominator = 1000
        return numerator / denominator

    def f_kappa_1(Z):
        numerator = 0.211 * Z**3 + 12.7 * Z**2 + 48.4 * Z + 6.45
        denominator = Z + 57.1
        return numerator / denominator

    def f_kappa_2(Z):
        numerator = 0.932 * Z**(7 / 3) + 0.135 * Z**2 + 12.3 * Z + 8.77
        denominator = Z + 4.84
        return numerator / denominator

    def f_kappa_3(Z):
        numerator = 0.246 * Z**3 + 2.65 * Z**2 - 92.8 * Z - 1.96
        denominator = Z**2 + 19.9 * Z + 35.3
        return numerator / denominator

    def f_kappa_4(Z):
        numerator = 2.76 * Z**(5 / 3) - 0.836 * Z**(2 / 3) - 0.0611
        denominator = Z - 0.214
        return numerator / denominator

    def f_k_0(Z):
        numerator = 0.0396 * Z**3 + 46.3 * Z + 176
        denominator = 1000
        return numerator / denominator

    def f_k_1(Z):
        numerator = 15.4 * Z**3 + 188 * Z**2 + 240 * Z + 35.3
        denominator = 1000 * Z + 397
        return numerator / denominator

    def f_k_2(Z):
        numerator = -0.159 * Z**2 - 12.5 * Z + 34.1
        denominator = Z**(2 / 3) + 0.741 * Z**(1 / 3) + 31.0
        return numerator / denominator

    def f_k_3(Z):
        numerator = 0.431 * Z**2 + 3.69 * Z + 0.0314
        denominator = Z + 3.62
        return numerator / denominator

    def f_k_4(Z):
        numerator = 0.0258 * Z**2 - 1.63 * Z + 0.711
        denominator = Z**(4 / 3) + 4.36 * Z**(2 / 3) + 2.75
        return numerator / denominator

    def f_k_5(Z):
        numerator = Z**3 + 11.9 * Z**2 + 28.8 * Z + 9.07
        denominator = 173 * Z + 133
        return numerator / denominator

    kappa_par_e = [3.204, 2.464, f_kappa_par_e(Z)]
    kappa_0 = [0.936, 1.749, f_kappa_0(Z)]
    kappa_1 = [1.166, 2.635, f_kappa_1(Z)]
    kappa_2 = [3.791, 5.644, f_kappa_2(Z)]
    kappa_3 = [-1.635, -2.212, f_kappa_3(Z)]
    kappa_4 = [2.370, 4.129, f_kappa_4(Z)]
    k_0 = [0.222, 0.269, f_k_0(Z)]
    k_1 = [0.343, 0.580, f_k_1(Z)]
    k_2 = [0.655, 0.252, f_k_2(Z)]
    k_3 = [0.899, 1.626, f_k_3(Z)]
    k_4 = [-0.110, -0.201, f_k_4(Z)]
    k_5 = [0.166, 0.255, f_k_5(Z)]

    kappa_par = kappa_par_e[Z_idx]
    if field_orientation == 'parallel' or field_orientation == 'par':
        return Z * kappa_par

    def f_kappa_perp(Z_idx):
        numerator = (13 / 4 * Z + np.sqrt(2)) * r + \
            kappa_0[Z_idx] * kappa_par_e[Z_idx]
        denominator = r**3 + \
            kappa_4[Z_idx] * r**(7 / 3) + \
            kappa_3[Z_idx] * r**(2) + \
            kappa_2[Z_idx] * r**(5 / 3) + \
            kappa_1[Z_idx] * r + \
            kappa_0[Z_idx]
        return numerator / denominator

    kappa_perp = f_kappa_perp(Z_idx)
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return Z * kappa_perp

    def f_kappa_cross(Z_idx):
        numerator = r * (5 / 2 * r + k_0[Z_idx] / k_5[Z_idx])
        denominator = r**3 + \
            k_4[Z_idx] * r**(7 / 3) + \
            k_3[Z_idx] * r**(2) + \
            k_2[Z_idx] * r**(5 / 3) + \
            k_1[Z_idx] * r + \
            k_0[Z_idx]
        return numerator / denominator

    kappa_cross = f_kappa_cross(Z_idx)
    if field_orientation == 'cross':
        return Z * kappa_cross

    if field_orientation == 'all':
        return np.array((Z * kappa_par, Z * kappa_perp, Z * kappa_cross))


def _nondim_resist_ji_held(hall, Z, field_orientation):
    """TODO"""

    allowed_Z = [1, 2, 'arbitrary']
    Z_idx = check_Z(allowed_Z, Z)
    r = np.abs(Z * hall)

    def f_alpha_par_e(Z):
        numerator = Z**(2 / 3)
        denominator = 1.46 * Z**(2 / 3) - 0.330 * Z**(1 / 3) + 0.888
        return 1 - numerator / denominator

    def f_alpha_0(Z):
        return 0.623 * Z**(5 / 3) - 2.61 * Z**(4 / 3) + 3.56 * Z + 0.557

    def f_alpha_1(Z):
        return 2.24 * Z**(2 / 3) - 1.11 * Z**(1 / 3) + 1.84

    def f_alpha_2(Z):
        return -0.0983 * Z**(1 / 3) + 0.0176

    def f_a_0(Z):
        return 0.0759 * Z**(8 / 3) + 0.897 * Z**2 + 2.06 * Z + 1.06

    def f_a_1(Z):
        return 2.18 * Z**(5 / 3) + 5.31 * Z + 3.73

    def f_a_2(Z):
        return 7.41 * Z + 1.11 * Z**(2 / 3) - 1.17

    def f_a_3(Z):
        return 3.89 * Z**(2 / 3) - 4.51 * Z**(1 / 3) + 6.76

    def f_a_4(Z):
        return 2.26 * Z**(1 / 3) + 0.281

    def f_a_5(Z):
        return 1.18 * Z**(5 / 3) - 1.03 * Z**(4 / 3) + 3.60 * Z + 1.32

    alpha_par_e = [0.504, 0.431, f_alpha_par_e(Z)]
    alpha_0 = [2.130, 3.078, f_alpha_0(Z)]
    alpha_1 = [2.970, 3.997, f_alpha_1(Z)]
    alpha_2 = [-0.081, -0.106, f_alpha_2(Z)]
    a_0 = [4.093, 9.250, f_a_0(Z)]
    a_1 = [11.22, 21.27, f_a_1(Z)]
    a_2 = [7.350, 15.41, f_a_2(Z)]
    a_3 = [6.140, 7.253, f_a_3(Z)]
    a_4 = [2.541, 3.128, f_a_4(Z)]
    a_5 = [5.070, 9.671, f_a_5(Z)]

    alpha_par = alpha_par_e[Z_idx]
    if field_orientation == 'parallel' or field_orientation == 'par':
        return alpha_par

    def f_alpha_perp(Z_idx):
        numerator = 1.46 * Z**(2 / 3) * r + \
            alpha_0[Z_idx] * (1 - alpha_par_e[Z_idx])
        denominator = r**(5 / 3) + \
            alpha_2[Z_idx] * r**(4 / 3) + \
            alpha_1[Z_idx] * r + \
            alpha_0[Z_idx]
        return 1 - numerator / denominator

    alpha_perp = f_alpha_perp(Z_idx)
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return alpha_perp

    def f_alpha_cross(Z_idx):
        numerator = Z**(2 / 3) * r * (2.53 * r + a_0[Z_idx] / a_5[Z_idx])
        denominator = r**(8 / 3) + \
            a_4[Z_idx] * r**(7 / 3) + \
            a_3[Z_idx] * r**(2) + \
            a_2[Z_idx] * r**(5 / 3) + \
            a_1[Z_idx] * r + \
            a_0[Z_idx]
        return numerator / denominator

    alpha_cross = f_alpha_cross(Z_idx)
    if field_orientation == 'cross':
        return alpha_cross

    if field_orientation == 'all':
        return np.array((alpha_par, alpha_perp, alpha_cross))


def _nondim_tec_ji_held(hall, Z, field_orientation):
    """TODO"""

    allowed_Z = [1, 2, 'arbitrary']
    Z_idx = check_Z(allowed_Z, Z)
    r = np.abs(Z * hall)

    def f_beta_par_e(Z):
        numerator = Z**(5 / 3)
        denominator = 0.693 * Z**(5 / 3) - 0.279 * Z**(4 / 3) + Z + 0.01
        return numerator / denominator

    def f_beta_0(Z):
        return 0.156 * Z**(8 / 3) + 0.994 * Z**2 + 3.21 * Z - 0.84

    def f_beta_1(Z):
        return 3.69 * Z**(5 / 3) + 3.77 * Z + 0.77

    def f_beta_2(Z):
        return 9.43 * Z + 4.22 * Z**(2 / 3) - 12.9 * Z**(1 / 3) + 4.56

    def f_beta_3(Z):
        return 2.70 * Z**(2 / 3) + 1.46 * Z**(1 / 3) - 0.17

    def f_beta_4(Z):
        return 2.58 * Z**(1 / 3) + 0.17

    def f_b_0(Z):
        numerator = 6.87 * Z**(3) + 78.2 * Z**2 + 623 * Z + 366
        denominator = 1000
        return numerator / denominator

    def f_b_1(Z):
        return 0.134 * Z**2 + 0.977 * Z + 0.17

    def f_b_2(Z):
        return 0.689 * Z**(4 / 3) - 0.377 * Z**(2 / 3) + \
            3.94 * Z**(1 / 3) + 0.644

    def f_b_3(Z):
        return -0.109 * Z + 1.33 * Z**(2 / 3) - 3.80 * Z**(1 / 3) + 0.289

    def f_b_4(Z):
        return 2.46 * Z**(2 / 3) + 0.522

    def f_b_5(Z):
        return 0.102 * Z**2 + 0.746 * Z + 0.072 * Z**(1 / 3) + 0.211

    beta_par_e = [0.702, 0.905, f_beta_par_e(Z)]
    beta_0 = [3.520, 10.55, f_beta_0(Z)]
    beta_1 = [8.230, 20.03, f_beta_1(Z)]
    beta_2 = [5.310, 13.87, f_beta_2(Z)]
    beta_3 = [3.990, 5.955, f_beta_3(Z)]
    beta_4 = [2.750, 3.421, f_beta_4(Z)]
    b_0 = [1.074, 1.980, f_b_0(Z)]
    b_1 = [1.281, 2.660, f_b_1(Z)]
    b_2 = [4.896, 6.746, f_b_2(Z)]
    b_3 = [-2.290, -2.605, f_b_3(Z)]
    b_4 = [2.982, 4.427, f_b_4(Z)]
    b_5 = [1.131, 2.202, f_b_5(Z)]

    beta_par = beta_par_e[Z_idx]
    if field_orientation == 'parallel' or field_orientation == 'par':
        return beta_par

    def f_beta_perp(Z_idx):
        numerator = 6.33 * Z**(5 / 3) * r + beta_0[Z_idx] * beta_par_e[Z_idx]
        denominator = r**(8 / 3) + \
            beta_4[Z_idx] * r**(7 / 3) + \
            beta_3[Z_idx] * r**(2) + \
            beta_2[Z_idx] * r**(5 / 3) + \
            beta_1[Z_idx] * r + \
            beta_0[Z_idx]
        return numerator / denominator

    beta_perp = f_beta_perp(Z_idx)
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return beta_perp

    def f_beta_cross(Z_idx):
        numerator = Z * r * (3 / 2 * r + b_0[Z_idx] / b_5[Z_idx])
        denominator = r**(3) + \
            b_4[Z_idx] * r**(7 / 3) + \
            b_3[Z_idx] * r**(2) + \
            b_2[Z_idx] * r**(5 / 3) + \
            b_1[Z_idx] * r + \
            b_0[Z_idx]
        return numerator / denominator

    beta_cross = f_beta_cross(Z_idx)
    if field_orientation == 'cross':
        return beta_cross

    if field_orientation == 'all':
        return np.array((beta_par, beta_perp, beta_cross))


def _nondim_visc_e_ji_held(hall, Z):
    """TODO"""

    allowed_Z = [1, 2, 'arbitrary']
    Z_idx = check_Z(allowed_Z, Z)
    r = np.abs(Z * hall)

    def f_eta_0_e(Z):
        return 1 / (0.55 * Z + 0.083 * Z**(1 / 3) + 0.732)

    def f_hprime_0(Z):
        return 0.0699 * Z**3 + 0.558 * Z**2 + 1.66 * Z + 1.06

    def f_hprime_1(Z):
        return 0.657 * Z**2 + 1.42 * Z + 0.416

    def f_hprime_2(Z):
        return -0.369 * Z**(4 / 3) + 0.379 * Z + 0.339 * Z**(1 / 3) + 2.17

    def f_hprime_3(Z):
        return 2.16 * Z - 0.657 * Z**(1 / 3) + 0.0347

    def f_hprime_4(Z):
        return -0.0703 * Z**(2 / 3) - 0.224 * Z**(1 / 3) + 0.333

    def f_h_0(Z):
        return 0.0473 * Z**3 + 0.323 * Z**2 + 0.951 * Z + 0.407

    def f_h_1(Z):
        return 0.171 * Z**2 + 0.523 * Z + 0.336

    def f_h_2(Z):
        return 0.362 * Z**(4 / 3) + 0.178 * Z + 1.06 * Z**(1 / 3) + 1.26

    def f_h_3(Z):
        return 0.599 * Z + 0.106 * Z**(2 / 3) - 0.444 * Z**(1 / 3) - 0.161

    def f_h_4(Z):
        return -0.16 * Z**(2 / 3) + 0.06 * Z**(1 / 3) + 0.232

    def f_h_5(Z):
        return 0.183 * Z**2 + 0.714 * Z + 0.0375 * Z**(1 / 3) + 0.47

    eta_0_e = [0.733, 0.516, f_eta_0_e(Z)]
    hprime_0 = [3.348, 7.171, f_hprime_0(Z)]
    hprime_1 = [2.493, 5.884, f_hprime_1(Z)]
    hprime_2 = [2.519, 2.425, f_hprime_2(Z)]
    hprime_3 = [1.538, 3.527, f_hprime_3(Z)]
    hprime_4 = [0.039, -0.061, f_hprime_4(Z)]
    h_0 = [1.728, 3.979, f_h_0(Z)]
    h_1 = [1.030, 2.066, f_h_1(Z)]
    h_2 = [2.860, 3.864, f_h_2(Z)]
    h_3 = [0.100, 0.646, f_h_3(Z)]
    h_4 = [0.132, 0.054, f_h_4(Z)]
    h_5 = [1.405, 2.677, f_h_5(Z)]

    eta_0 = eta_0_e[Z_idx]

    def f_eta_2(Z_idx, r):
        numerator = (6 / 5 * Z + 3 / 5 * np.sqrt(2)) * r + \
            hprime_0[Z_idx] * eta_0_e[Z_idx]
        denominator = r**(3) + \
            hprime_4[Z_idx] * r**(7 / 3) + \
            hprime_3[Z_idx] * r**(2) + \
            hprime_2[Z_idx] * r**(5 / 3) + \
            hprime_1[Z_idx] * r + \
            hprime_0[Z_idx]
        return numerator / denominator

    eta_2 = f_eta_2(Z_idx, r)

    eta_1 = f_eta_2(Z_idx, 2 * r)

    def f_eta_4(Z_idx, r):
        numerator = r * (r + h_0[Z_idx] / h_5[Z_idx])
        denominator = r**(3) + \
            h_4[Z_idx] * r**(7 / 3) + \
            h_3[Z_idx] * r**(2) + \
            h_2[Z_idx] * r**(5 / 3) + \
            h_1[Z_idx] * r + \
            h_0[Z_idx]
        return numerator / denominator

    eta_4 = f_eta_4(Z_idx, r)

    eta_3 = f_eta_4(Z_idx, 2 * r)

    return np.array((eta_0, eta_1, eta_2, eta_3, eta_4))


def _nondim_tc_i_ji_held(hall, Z, mu, theta, field_orientation):
    """TODO"""

#    mu = m_e / m_i
#    theta = T_e / T_i
    zeta = 1 / Z * np.sqrt(mu / theta)
    r = np.abs(hall / np.sqrt(2))

#    K = 2  # 2x2 moments, equivalent to original Braginskii
    K = 3  # 3x3 moments, more accurate

    if K == 3:
        Delta_par_i1 = 1 + 26.90 * zeta + 187.5 * zeta**2 + 346.9 * zeta**3
        kappa_par_i = (5.586 + 101.7 * zeta + 289.1 * zeta**2) / Delta_par_i1
    elif K == 2:
        Delta_par_i1 = 1 + 13.50 * zeta + 36.46 * zeta**2
        kappa_par_i = (5.524 + 30.38 * zeta) / Delta_par_i1
    if field_orientation == 'parallel' or field_orientation == 'par':
        return kappa_par_i / np.sqrt(2)

    if K == 3:
        Delta_perp_i1 = r**6 + \
            (3.635 + 29.15 * zeta + 83 * zeta**2) * r**4 + \
            (1.395 + 35.64 * zeta + 344.9 * zeta**2 +
             1345 * zeta**3 + 1891 * zeta**4) * r**2 + \
            0.09163 * Delta_par_i1**2
        kappa_perp_i = ((np.sqrt(2) + 15 / 2 * zeta) * r**4 +
                        (3.841 + 57.59 * zeta + 297.8 * zeta**2 +
                         555 * zeta**3) * r**2 +
                        0.09163 * kappa_par_i * (Delta_par_i1)**2
                        ) / Delta_perp_i1
    elif K == 2:
        Delta_perp_i1 = r**4 + \
            (1.352 + 12.49 * zeta + 34 * zeta**2) * r**2 + \
            0.1693 * Delta_par_i1**2
        kappa_perp_i = ((np.sqrt(2) + 15 / 2 * zeta) * r**2 +
                        0.1693 * kappa_par_i * (Delta_par_i1)**2
                        ) / Delta_perp_i1
    if field_orientation == 'perpendicular' or field_orientation == 'perp':
        return kappa_perp_i / np.sqrt(2)

    if K == 3:
        kappa_cross_i = r * (5 / 2 * r**4 +
                             (7.963 + 64.40 * zeta + 185 * zeta**2) * r**2 +
                             1.344 + 44.54 * zeta + 511.9 * zeta**2 +
                             2155 * zeta**3 + 3063 * zeta**4
                             ) / Delta_perp_i1
    elif K == 2:
        kappa_cross_i = r * (5 / 2 * r**2 +
                             2.323 + 22.73 * zeta + 62.5 * zeta**2
                             ) / Delta_perp_i1
    if field_orientation == 'cross':
        return kappa_cross_i / np.sqrt(2)

    if field_orientation == 'all':
        return np.array((kappa_par_i / np.sqrt(2), kappa_perp_i / np.sqrt(2),
                         kappa_cross_i / np.sqrt(2)))


def _nondim_visc_i_ji_held(hall, Z, mu, theta):

    zeta = 1 / Z * np.sqrt(mu / theta)
    r = np.abs(hall / np.sqrt(2))
    r13 = 2 * r

#    K = 2  # 2x2 moments, equivalent to original Braginskii
    K = 3  # 3x3 moments, more accurate

    if K == 3:
        Delta_par_i2 = 1 + 15.79 * zeta + 63.92 * zeta**2 + 71.69 * zeta**3
        Delta_perp_i2 = r**6 + \
            (4.391 + 26.69 * zeta + 56 * zeta**2) * r**4 + \
            (3.191 + 49.62 * zeta + 306.4 * zeta**2 +
             808.1 * zeta**3 + 784 * zeta**4) * r**2 + \
            0.4483 * Delta_par_i2**2
        eta_0_i = (1.365 + 16.75 * zeta + 35.84 * zeta**2) / Delta_par_i2
        eta_2_i = ((3 / 5 * np.sqrt(2) + 2 * zeta) * r**4 +
                   (2.680 + 25.98 * zeta + 90.71 * zeta**2 + 104 * zeta**3) * r**2 +
                   0.4483 * eta_0_i * (Delta_par_i2)**2
                   ) / Delta_perp_i2
        eta_1_i = ((3 / 5 * np.sqrt(2) + 2 * zeta) * r13**4 +
                   (2.680 + 25.98 * zeta + 90.71 * zeta**2 + 104 * zeta**3) * r13**2 +
                   0.4483 * eta_0_i * (Delta_par_i2)**2
                   ) / Delta_perp_i2
        eta_4_i = r * (r**4 +
                       (3.535 + 23.30 * zeta + 52 * zeta**2) * r**2 +
                       0.9538 + 21.81 * zeta + 174.2 * zeta**2 +
                       538.4 * zeta**3 + 576 * zeta**4
                       ) / Delta_perp_i2
        eta_3_i = r13 * (r13**4 +
                         (3.535 + 23.30 * zeta + 52 * zeta**2) * r13**2 +
                         0.9538 + 21.81 * zeta + 174.2 * zeta**2 +
                         538.4 * zeta**3 + 576 * zeta**4
                         ) / Delta_perp_i2
    elif K == 2:
        Delta_par_i2 = 1 + 7.164 * zeta + 10.49 * zeta**2
        Delta_perp_i2 = r**4 + \
            (2.023 + 11.68 * zeta + 20 * zeta**2) * r**2 + \
            0.5820 * Delta_par_i2**2
        eta_0_i = (1.357 + 5.243 * zeta) / Delta_par_i2
        eta_2_i = ((3 / 5 * np.sqrt(2) + 2 * zeta) * r**2 +
                   0.5820 * eta_0_i * (Delta_par_i2)**2
                   ) / Delta_perp_i2
        eta_1_i = ((3 / 5 * np.sqrt(2) + 2 * zeta) * r13**2 +
                   0.5820 * eta_0_i * (Delta_par_i2)**2
                   ) / Delta_perp_i2
        eta_4_i = r * (r**2 +
                       1.188 + 8.283 * zeta + 16 * zeta**2
                       ) / Delta_perp_i2
        eta_3_i = r13 * (r13**2 +
                         1.188 + 8.283 * zeta + 16 * zeta**2
                         ) / Delta_perp_i2
    return np.array((eta_0_i / np.sqrt(2), eta_1_i / np.sqrt(2),
                     eta_2_i / np.sqrt(2), eta_3_i / np.sqrt(2),
                     eta_4_i / np.sqrt(2)))

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#                      end of classical_transport
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
