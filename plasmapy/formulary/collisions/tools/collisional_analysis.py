"""
Module containing the Collisional Analysis formulation.
"""
__all__ = ["collisional_thermalization"]

import astropy.units as u
import numbers
import numpy as np

from plasmapy.particles import Particle, ParticleLike
from plasmapy.utils.decorators import validate_quantities


@validate_quantities(
    T_1={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    T_2={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def collisional_thermalization(
    r_0: u.au,
    r_n: u.au,
    n_1: u.cm**-3,
    n_2: u.cm**-3,
    v_1: u.m / u.s,
    T_1: u.K,
    T_2: u.K,
    ions: ParticleLike = [Particle("p+"), Particle("He-4++")],
    n_step: int = 1000,
    density_scale: float = -1.8,
    velocity_scale: float = -0.2,
    temperature_scale: float = -0.77,
):
    r"""
    Contains the functionality to calculate the collisional
    thermalization a plasma will undergo in transit, taken from
    :cite:t:`maruca:2013`. Coulomb collisions are soft, small angle
    deflections mediated by the electrostatic force and occur between
    constituent particles of a plasma :cite:t:`baumjohann:1997`.
    These collisions cause the plasma to thermalize over distance,
    i.e. the temperature of the plasma will approach thermal
    equilibrium. This function allows this thermalization to be
    modeled and predicts the temperature ratio, for different ion
    species within the plasma, at a different point in space.

    Parameters
    ----------
    r_0 : `~astropy.units.Quantity`
        Starting position of the plasma in units convertible
        to astronomical units.

    r_n : `~astropy.units.Quantity`
        Final position of the plasma in units convertible
        to astronomical units.

    n_1 : `~astropy.units.Quantity`
        The primary ion number density in units convertible
        to m\ :sup:`-3`.

    n_2 : `~astropy.units.Quantity`
        The secondary ion number density in units convertible
        to m\ :sup:`-3`.

    v_1 : `~astropy.units.Quantity`
        The primary ion speed in units convertible to ms\ :sup:`-1`.

    T_1 : `~astropy.units.Quantity`
        Temperature of the primary ion in units convertible to
        temperature K.

    T_2 : `~astropy.units.Quantity`
        Temperature of the secondary ion in units convertible to
        temperature K.

    ions : `Particle`
        Particle list containing two (2) particles, primary ion of
        interest is entered first, followed by the secondary ion.

    n_step : `int`
        The number of intervals used in solving a differential
        equation via the Euler method. Must be an int.

    density_scale : `float`
        The value used as the scaling parameter for the primary ion
        density, the default value is -1.8 and is taken from
        :cite:t:`hellinger:2011`.

    velocity_scale : `float`
        The value used as the scaling parameter for the primary ion
        velocity, the default value is -0.2 and is taken from
        :cite:t:`hellinger:2011`.

    temperature_scale : `float`
        The value used as the scaling parameter for the primary ion
        temperature, the default value is -0.77 and is taken from
        :cite:t:`hellinger:2011`.

    Returns
    -------
    theta : `float`
        The dimensionless ion temperature ratio prediction
        for the distance provided.

    Raises
    ------
    `TypeError`
        If applicable arguments are not instances of
        `~astropy.units.Quantity` or cannot be converted into one.

    ~astropy.units.UnitTypeError
        If applicable arguments do not have units convertible to the
        expected units.

    Notes
    -----
    The processes by which Coulomb collisions bring ion temperatures
    into local thermal equilibrium (LTE) has received considerable
    attention, :cite:t:`verscharen:2019`. The relative temperature
    between constituent plasma ion species is given as:

    .. math::

        \theta_{21} = \frac{T_{2}}{T_{1}} \, ,

    where :math:`T_{1}` and :math:`T_{2}` are the scalar temperatures
    for the primary ion of interest and the secondary ion respectively.
    The scalar temperature defined as:

    .. math::

        T_{\rm i} = \frac{2T_{{\rm i}, \perp} + T_{{\rm i}, \parallel}}{3} \, ,

    where :math:`T_{{\rm i}, \perp}` and :math:`T_{{\rm i}, \parallel}`
    are the temperature of the :math:`{\rm i}`-particles along the
    axes perpendicular and parallel to the ambient magnetic field.

    In order to determine how extensively an individual parcel of
    plasma has been processed by Coulomb collisions, :cite:t:`maruca:2013`
    introduced an approached called collisional analysis. Which seeks
    to quantify how collisions affect the plasma’s departures from
    LTE, the equation for collisional thermalization
    from :cite:t:`maruca:2013` is given below:

     .. math::

        \frac{d \theta_{21}}{dr} = A \left ( \frac{n_1}{v_1 T_1^{3/2}} \right ) \frac{\left( \mu_{1} \mu_{2}
        \right )^{1/2} Z_{1} Z_{2} \left( 1 - \theta_{21} \right ) \left(1 +
        \eta_{21}\theta_{21} \right )}{\left( \frac{\mu_{2}}{\mu_{1}} +
        \theta_{21} \right )^{3/2}} \lambda_{21}

    and

    .. math::

         \lambda_{21} = 9 + \ln \left| B \left ( \frac{T^{3}_{1}}{n_{1}} \right )^{1/2} \left( \frac{\theta_{21} +
         \frac{\mu_{2}}{\mu_1}}{Z_{1}Z_{2}(\mu_{1} + \mu_{2})} \right )
         \left( \frac{n_{2}Z_{2}^{2}}{n_{1}Z_{1}^{2}} + \theta_{21}
         \right)^{1/2}\right |

    With :math:`\eta = \frac{n_{2}}{n_{1}}`,
    :math:`\theta = \frac{T_{2}}{T_{1}}`, :math:`A = 2.60 \times 10^{7} \, {\rm cm}^{3} \, {\rm km} \, {\rm K}^{3/2} \, {\rm s}^{-1} \, {\rm au}^{-1}`,
    and :math:`B = 1 \, {\rm cm}^{-3/2}{\rm K}^{-3/2}`.

    The thermalization is from Coulomb collisions, which assumes
    “soft”, small-angle deflections mediated by the electrostatic
    force :cite:t:`baumjohann:1997`. It is assumed that there is no
    relative drift between the ion species and that it is a mixed ion
    collision, the Coulomb logarithm for a mixed ion collision is
    given by :cite:t:`nrlformulary:2019`.

    The density, velocity and temperature of the primary ion can be
    radially scaled, as seen below. The values for the scaling can be
    altered, though the default values are taken from
    :cite:t:`hellinger:2011`.

    .. math::

        n(r) \propto r^{-1.8}\ , \hspace{1cm} v_{r}(r) \propto r^{-0.2}\ ,
        \hspace{0.5cm} {\rm and} \hspace{0.5cm} T(r) \propto r^{-0.74}

    Application is primarily for the solar wind.

    Examples
    --------


    """

    # Validate number of ions
    if len(ions) != 2:
        raise ValueError(
            "Argument 'ions' can only take two input values. Instead "
            f"received {len(ions)} input values."
        )

    # Validate n_step argument
    if not isinstance(n_step, int):
        raise TypeError(
            "Argument 'n_step' is of incorrect type, type of "
            f"{type(n_step)} received instead. While 'n_step' must be "
            "of type int."
        )

    # Validate scaling arguments
    for arg in (density_scale, velocity_scale, temperature_scale):
        if not isinstance(arg, numbers.Real):
            raise TypeError(
                "Scaling argument is of incorrect type, type of "
                f"{type(arg)} received instead. Scaling argument "
                "should be of type float or int."
            )

    # Define differential equation function
    def df_eq(
        r_0,
        r_n,
        n_1,
        n_2,
        v_1,
        T_1,
        T_2,
        ions,
        n_step,
        density,
        velocity,
        temperature,
    ):
        # Initialize the alpha-proton charge and mass ratios.
        z_1 = ions[0].charge_number
        mu_1 = ions[0].mass_number

        z_2 = ions[1].charge_number
        mu_2 = ions[1].mass_number

        # Initialise.
        d_r = (r_0 - r_n) / (1.0 * n_step)
        d_r = d_r.value

        # Define constants
        A = 2.6 * 10**7 * (u.cm**3 * u.km * u.K**1.5) / (u.s * u.au)
        B = 1 / (u.cm * u.K) ** 1.5

        # Define Coulomb log for mixed ion collisions
        def lambda_ba(
            theta,
            T_1,
            n_1,
            n_2,
            z_1,
            z_2,
            mu_1,
            mu_2,
        ):
            a = np.sqrt(T_1**3 / n_1)
            b = (theta + mu_2 / mu_1) / (z_1 * z_2 * (mu_1 + mu_2))
            c = np.sqrt(n_2 * z_2**2 / n_1 * z_1**2 + theta)
            return 9 + np.log(B * a * b * c)

        # Loop.
        for i in range(n_step):
            r = r_n.value + ((i + 1) * d_r)

            eta = n_2 / n_1
            theta = T_2 / T_1

            l_ba = lambda_ba(theta, T_1, n_1, n_2, z_1, z_2, mu_1, mu_2)

            n_1 = n_1 * (r / r_n.value) ** density
            v_1 = v_1 * (r / r_n.value) ** velocity
            T_1 = T_1 * (r / r_n.value) ** temperature

            alpha = n_1 / (v_1 * (T_1**1.5))
            beta = (
                (np.sqrt(mu_1 * mu_2) * z_1 * z_2 * (1 - theta) * (1 + eta * theta))
                / (np.sqrt((mu_2 / mu_1) + theta) ** 3)
                * u.au
                / u.km
            )

            d_theta = d_r * u.m * alpha * l_ba * A * beta

            # print("dr, A, l_ba", d_r, A, l_ba, alpha, beta)

            print(theta, d_theta)
            theta = theta + d_theta

        return theta.value

    vars = [r_0, r_n, n_1, n_2, v_1, T_1, T_2]

    d_type = []
    for var in vars:
        if hasattr(var, "__len__"):
            d_type.append(True)
        else:
            d_type.append(False)

    var = all(i for i in d_type)

    if not var:
        return df_eq(
            r_0,
            r_n,
            n_1,
            n_2,
            v_1,
            T_1,
            T_2,
            ions,
            n_step,
            density_scale,
            velocity_scale,
            temperature_scale,
        )
    else:
        if all(len(vars[0]) == len(z) for z in vars[1:]):
            res = []
            for i in range(len(vars[0])):
                res.append(
                    df_eq(
                        r_0[i],
                        r_n[i],
                        n_1[i],
                        n_2[i],
                        v_1[i],
                        T_1[i],
                        T_2[i],
                        ions,
                        n_step,
                        density_scale,
                        velocity_scale,
                        temperature_scale,
                    )
                )
            return res

        else:
            raise ValueError(
                "Argument(s) are of unequal lengths, the following "
                "arguments should be of equal length: 'r_0', 'r_n', "
                "'n_1', 'n_2', 'v_1', 'T_1' and 'T_2'."
            )
