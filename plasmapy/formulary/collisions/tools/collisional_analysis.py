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
    Contains the functionality to calculate collisional
    thermalization a plasma will undergo in transit, taken from
    :cite:t:`maruca:2013`. Coulomb collisions, soft, small angle
    deflections mediated by the electrostatic force, occur between
    constituent particles of a plasma :cite:t:`baumjohann:1997`. This
    causes the plasma to thermalize over distance, i.e. the temperature
    of the plasma approaches thermal equilibrium. This function
    allows the thermalization of a plasma to be modeled and predicts
    the temperature ratio, for ions within the plasma, at a different
    point in space.

    Parameters
    ----------
    r_0 : `~astropy.units.Quantity`
        Starting position of the plasma in units convertible
        to meters or astronomical units.

    r_n : `~astropy.units.Quantity`
        Final position of the plasma in units of meters or
        astronomical units.

    n_1 : `~astropy.units.Quantity`
        The primary ion number density in units convertible to m\ :sup:`-3`.

    n_2 : `~astropy.units.Quantity`
        The secondary ion number density in units convertible to m\ :sup:`-3`.

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
        The value used for scaling parameter of the primary ion
        density, default value is -1.8 and is taken from
        :cite:t:`hellinger:2011`.

    velocity_scale : `float`
        The value used for scaling parameter of the primary ion
        velocity, default value is -0.2 and is taken from
        :cite:t:`hellinger:2011`.

    temperature_scale : `float`
        The value used for scaling parameter of the primary ion
        temperature, default value is -0.77 and is taken from
        :cite:t:`hellinger:2011`.

    Returns
    -------
    theta : `float`
        The dimensionless alpha-proton temperature ratio prediction
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
    Equation for collisional thermalization from :cite:t:`maruca:2013`,
    is given below:

     .. math::

        \frac{d \theta_{ba}}{dr} = k \frac{\left( \mu_{a} \mu_{b} Z_{a}
        Z_{b} \right )^{1/2} \left( 1 - \theta_{ba} \right ) \left(1 +
        \eta_{ba}\theta_{ba} \right )}{\left( \frac{\mu_{a}}{\mu_{b}} +
        \theta_{ba} \right )} \lambda_{ba}

    with

    .. math::

         \lambda_{ba} = 9 + \ln \left|  \left( \frac{\theta_{ba} +
         \frac{\mu_{b}}{\mu_a}}{Z_{a}Z_{b}(\mu_{a} + \mu_{b})} \right )
         \left( \frac{n_{b}Z_{b}^{2}}{n_{a}Z_{a}^{2}} + \theta_{ba}
         \right)^{1/2}\right |

    With :math:`\eta = \frac{n_{2}}{n_{1}}` and :math:`\theta = \frac{T_{2}}{T_{1}}`

    Applicable is primarily for the solar wind.


    It is assumed that there is no relative drift between the ion
    species and that it is a mixed ion collision. Thermalization is
    from Coulomb collisions, which assumes “soft”, small-angle
    deflections mediated by the electrostatic force.

    The density, velocity and temperature of the primary ion can be
    radially scaled, as seen below. The values for the scaling can be
    altered, though the default values are taken from
    :cite:t:`hellinger:2011`.

    .. math::

        n(r) \propto r^{-1.8}\ , \hspace{1cm} v_{r}(r) \propto r^{-0.2}\ , \hspace{0.5cm} {\rm and} \hspace{0.5cm} T(r) \propto r^{-0.74}


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

        # Define constant
        k = 0.174

        # Define Coulomb log for mixed ion collisions
        def lambda_ba(
            theta,
            n_1,
            n_2,
            z_1,
            z_2,
            mu_1,
            mu_2,
        ):
            return 9 + np.log(
                np.sqrt(n_2 * z_2**2 / theta * n_1 * z_1**2 + 1)
                * ((z_1 * z_2 * (mu_1 + mu_2)) / (theta + (mu_2 / mu_1)))
                * ((theta + mu_2 / mu_1) / (z_1 * z_2 * (mu_1 + mu_2)))
            )

        # Loop.
        for i in range(n_step):
            r = r_n.value + ((i + 1) * d_r)

            eta = n_2 / n_1
            theta = T_2 / T_1

            l_ba = lambda_ba(theta, n_1, n_2, z_1, z_2, mu_1, mu_2)

            n_1 = n_1 * (r / r_n) ** density
            v_1 = v_1 * (r / r_n) ** velocity
            T_1 = T_1 * (r / r_n) ** temperature

            d_theta = (
                d_r
                * l_ba
                * k
                * (np.sqrt(mu_1 * mu_2 * z_1 * z_2) * (1 - theta) * (1 + eta * theta))
                / (np.sqrt(mu_1 / mu_2 + theta) ** 3)
            )

            print(theta, d_theta)
            theta = theta + d_theta

        return theta

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
            L = 1
            for i in range(L):
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


import pandas as pd

obj = pd.read_pickle(
    r"/Users/elliotjohnson/GitHub/Collsionality/data/save/EA/solar_data.pkl"
)

key_list = list(obj.keys())
vars = {}

for key in key_list:
    vars[key] = list(obj[key].keys())

n_a = obj["proton"]["np1"] * u.cm**-3
n_b = obj["alpha"]["na"] * u.cm**-3
v_a = obj["proton"]["v_mag"] * u.m / u.s
T_a = obj["proton"]["Tperp1"] * u.K
T_b = obj["alpha"]["Trat"] * u.K
ions = [Particle("p+"), Particle("He-4++")]

l = len(obj["proton"]["np1"])
r_0 = [0.1] * l * u.au
r_n = [1.0] * l * u.au
n_step = 1

x = collisional_thermalization(r_0, r_n, n_a, n_b, v_a, T_a, T_b, ions, n_step)

import matplotlib.pyplot as plt

plt.hist(x, bins=50)
plt.show()
