"""
Module containing the Collisional Analysis formulation.
"""
__all__ = ["collisional_thermalization"]

import astropy.units as u
import numpy as np

from plasmapy.particles import Particle, ParticleLike
from plasmapy.utils.decorators import validate_quantities


@validate_quantities(
    T_a={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    T_b={"can_be_negative": False, "equivalencies": u.temperature_energy()},
)
def collisional_thermalization(
    r_0: u.au,
    r_n: u.au,
    n_a: u.cm**-3,
    n_b: u.cm**-3,
    v_a: u.m / u.s,
    T_a: u.K,
    T_b: u.K,
    ions: ParticleLike = [Particle("p+"), Particle("He-4++")],
    n_step: int = 1000,
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

    n_a : `~astropy.units.Quantity`
        The primary ion number density in units convertible to m\ :sup:`-3`.

    v_a : `~astropy.units.Quantity`
        The primary ion speed in units convertible to ms\ :sup:`-1`.

    T_a : `~astropy.units.Quantity`
        Temperature of the primary ion in units convertible to
        temperature K.

    T_b : `~astropy.units.Quantity`
        Temperature of the secondary ion in units convertible to
        temperature K.

    ions : `Particle`
        Temperature

    n_step : `int`
        Temperature


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
    scalings,


    Examples
    --------


    """

    # Validate n_step argument
    if not isinstance(n_step, int):
        raise TypeError(
            "Argument 'n_step' is of incorrect type, type of "
            f"{type(n_step)} received. While 'n_step' must be "
            "of type int."
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
        density_scale: float = -1.8,
        velocity_scale: float = -0.2,
        temperature_scale: float = -0.77,
    ):

        # Strip units
        r_0 = r_0.value
        r_n = r_n.value
        n_a = n_a.value
        n_b = n_b.value
        v_a = v_a.value
        T_a = T_a.value
        T_b = T_b.value

        # Initialize the alpha-proton charge and mass ratios.
        z_a = ions[0].charge_number
        mu_a = ions[0].mass_number

        z_b = ions[1].charge_number
        mu_b = ions[1].mass_number

        # Initialise.
        d_r = (r_0 - r_n) / (1.0 * n_step)

        # Loop.
        for i in range(n_step):
            r = r_n + ((i + 1) * d_r)

            theta = T_a / T_b

            n_a = n_a * (r / r_n) ** density_scale
            v_a = v_a * (r / r_n) ** velocity_scale
            T_a = T_a * (r / r_n) ** temperature_scale
            print(d_r)
            d_theta = (
                (
                    T_b * dT_dt(mu_a, mu_b, z_a, z_b, n_a, n_b, T_a, T_b)
                    - T_a * dT_dt(mu_b, mu_a, z_b, z_a, n_b, n_a, T_b, T_a)
                )
                / (v_a * T_b**2)
            ) * d_r
            print(theta, d_theta)
            theta = theta + d_theta

        return theta

    def c_log(mu_a, mu_b, Z_a, Z_b, n_a, n_b, T_a, T_b):
        alpha = Z_a * Z_b * (mu_a + mu_b)
        beta = mu_a * T_b + mu_b * T_a

        def v_therm(n, Z, T):
            return (n * Z**2) / T

        return 9 + np.log(
            (alpha / beta) * np.sqrt(v_therm(n_a, Z_a, T_a) + v_therm(n_b, Z_b, T_b))
        )

    def dT_dt(mu_a, mu_b, Z_a, Z_b, n_a, n_b, T_a, T_b):
        alpha = np.sqrt(mu_a * mu_b) * ((Z_a * Z_b) ** 2) * n_b
        beta = (mu_a * T_b + mu_b * T_a) ** (1.5)
        charlie = T_b - T_a
        delta = c_log(mu_a, mu_b, Z_a, Z_b, n_a, n_b, T_a, T_b)
        print(alpha, beta, charlie, delta)
        return (0.174) * (alpha / beta) * charlie * delta

    vars = [n_a, n_b, v_a, T_a, T_b]

    var = np.ndarray

    if var == int:
        return df_eq(r_0, r_n, n_a, n_b, v_a, T_a, T_b, ions, n_step)
    elif var == np.ndarray:
        if all(len(vars[0]) == len(l) for l in vars[1:]):
            res = []
            L = 10
            for i in range(L):
                res.append(
                    df_eq(
                        r_0, r_n, n_a[i], n_b[i], v_a[i], T_a[i], T_b[i], ions, n_step
                    )
                )
            return res

        else:
            raise ValueError(
                "Argument(s) are of unequal lengths, the following "
                "arguments should be of equal length: 'n_a', 'n_b', "
                "'v_a', 'T_a', 'T_b'."
            )

    return


import pandas as pd

obj = pd.read_pickle(
    r"/Users/elliotjohnson/GitHub/Collsionality/data/save/EA/solar_data.pkl"
)

key_list = list(obj.keys())
vars = {}

for key in key_list:
    vars[key] = list(obj[key].keys())


r_0 = 0.1 * u.au
r_n = 1.0 * u.au
n_a = obj["proton"]["np1"] * u.cm**-3
n_b = obj["alpha"]["na"] * u.cm**-3
v_a = obj["proton"]["v_mag"] * u.m / u.s
T_a = obj["proton"]["Tperp1"] * u.K
T_b = obj["alpha"]["Trat"] * u.K
ions = [Particle("p+"), Particle("He-4++")]
n_step: int = 10

x = collisional_thermalization(r_0, r_n, n_a, n_b, v_a, T_a, T_b, ions, n_step)


import matplotlib.pyplot as plt

plt.hist(x, bins="auto")
plt.show()
