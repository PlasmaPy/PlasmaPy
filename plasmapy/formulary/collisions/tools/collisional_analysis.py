"""
Module containing the Collisional Analysis formulation.
"""
__all__ = ["collisional_thermalization"]

import astropy.units as u
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
    desnity_scale: float = -1.8,
    velocity_scale: float = -0.2,
    temperature_scale: float = -0.77,
):
    r"""
    Coulomb collisions, soft, small angle deflections mediated by the
    electrostatic force, occur between constituent particles of a
    plasma. This can cause the plasma to thermalize over time, i.e.
    the temperature of the plasma approaches thermal equilibrium.
    This function allows the thermalization of the plasma to be
    modeled and predicts the temperature ratio, for ions within
    the plasma, at a different point in space.


    - particle thermalization
    - what is actual calculated
    - what each variable is and why its used, context



    Parameters
    ----------
    r_0 : `~astropy.units.Quantity`
        Starting position of the plasma in units of meters or
        astronomical units.

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

    !!!change to single vars
    scaling : `list`
        A list of values used for scaling parameters; density,
        velocity and temperature for the primary ion species.
        Defaults are taken from !!!Ref, and are -1.8, -0.2, -0.77
        respectively. The order for entry is density, velocity and
        temperature, DVT.


    Returns
    -------
    theta : `float`
        The dimensionless alpha-proton temperature ratio prediction
        for distance provided

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

    - how eta and theta are computed
    - applicable to what plasma, all mainly solar wind


    Big equation here

    assumptions
    - no relative drift
    - large angle deflection
    - mixed ion collisions assumption
    scalings,

    .. math::

        \frac{d \theta_{ba}}{dr} = k \frac{\left( \mu_{a} \mu_{b} Z_{a}
        Z_{b} \right )^{1/2} \left( 1 - \theta_{ba} \right ) \left(1 +
        \eta_{ba}\theta_{ba} \right )}{\left( \frac{\mu_{a}}{\mu_{b}} +
        \theta_{ba} \right )} \lambda_{ba}

    .. math::

         \lambda_{ba} = 9 + \ln \left|  \left( \frac{\theta_{ba} +
         \frac{\mu_{b}}{\mu_a}}{Z_{a}Z_{b}(\mu_{a} + \mu_{b})} \right )
         \left( \frac{n_{b}Z_{b}^{2}}{n_{a}Z_{a}^{2}} + \theta_{ba}
         \right)^{1/2}\right |

    Examples
    --------


    """

    # Validate n_step argument
    if not isinstance(n_step, int):
        raise TypeError(
            "Argument 'n_step' is of incorrect type, type of "
            f"{type(n_step)} received instead. While 'n_step' must be "
            "of type int."
        )

    # Validate scaling arguments
    for arg in (desnity_scale, velocity_scale, temperature_scale):
        if not isinstance(arg, (float, int)):
            raise TypeError(
                "Argument 'scaling' is of incorrect type, type of "
                f"{type(list)} received instead. While 'scaling' must be "
                "of type list. For additional information, please see "
                "the documentation."
            )


    def df_eq(
        r_0,
        r_n,
        n_a,
        n_b,
        v_a,
        T_a,
        T_b,
        ions,
        n_step,
        scaling_values = scaling,
    ):

        # Initialize the alpha-proton charge and mass ratios.
        z_a = ions[0].charge_number
        mu_a = ions[0].mass_number

        z_b = ions[1].charge_number
        mu_b = ions[1].mass_number


        # Initialise.
        d_r = (r_0 - r_n) / (1.0 * n_step)
        d_r = d_r.value


        # Define constant
        k = 0.174
        def lambda_ba(
                theta,
                n_a,
                n_b,
                z_a,
                z_b,
                mu_a,
                mu_b,

        ):

            return 9 #+ np.log(np.sqrt(n_b*z_b**2/theta*n_a*z_a**2 + 1)*((z_a*z_b*(mu_a+mu_b))/(theta + (mu_b/mu_a)))*((theta + mu_b/mu_a)/(z_a*z_b*(mu_a+mu_b))))


        # Loop.
        for i in range(n_step):
            r = r_n.value + ((i + 1) * d_r)

            eta = n_a/n_b
            theta = T_a / T_b

            l_ba = lambda_ba(theta, n_a, n_b, z_a, z_b, mu_a, mu_b)

            n_a = n_a * (r / r_n) ** scaling_values[0]
            v_a = v_a * (r / r_n) ** scaling_values[1]
            T_a = T_a * (r / r_n) ** scaling_values[2]

            d_theta = d_r*l_ba*k*(np.sqrt(mu_a*mu_b*z_a*z_b)*(1 - theta)*(1 + eta*theta))/(np.sqrt(mu_a / mu_b + theta) ** 3)

            print(theta, d_theta)
            theta = theta + d_theta

        return theta


    vars = [n_a, n_b, v_a, T_a, T_b]

    var = np.ndarray

    # isinstance(val, numbers.Integral)
    # numbers.Real
    # hasattr(“__len__”)


    if var == int:
        return df_eq(r_0, r_n, n_a, n_b, v_a, T_a, T_b, ions, n_step)
    elif var == np.ndarray:
        if all(len(vars[0]) == len(l) for l in vars[1:]):
            res = []
            L = 1
            for i in range(L):
                res.append(
                    df_eq(
                        r_0[i], r_n[i], n_a[i], n_b[i], v_a[i], T_a[i], T_b[i], ions, n_step
                    )
                )
            return res

        else:
            raise ValueError(
                "Argument(s) are of unequal lengths, the following "
                "arguments should be of equal length: 'n_a', 'n_b', "
                "'v_a', 'T_a', 'T_b'."
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
