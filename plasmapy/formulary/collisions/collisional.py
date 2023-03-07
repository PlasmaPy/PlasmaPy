"""
Module containing the Collisional Analysis formulation.
"""
__all__ = [
    "collisional_thermalization"
]

import astropy.units as u
import numpy as np

from plasmapy.particles import Particle, ParticleLike
from plasmapy.utils.decorators import validate_quantities

@validate_quantities(
    T_a={"can_be_negative": False, "equivalencies": u.temperature_energy()},
    T_b={"can_be_negative": False, "equivalencies": u.temperature_energy()}
)
def collisional_thermalization(
    r_0: u.m,
    r_n: u.m,
    n_a: u.cm**-3,
    n_b: u.cm**-3,
    v_a: u.m / u.s,
    T_a: u.K,
    T_b: u.K,
    ions: ParticleLike = [Particle("p+"), Particle("He-4++")],
    n_step: int = 1000
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



    def sub_function(
        density_scale: float = -1.8,
        velocity_scale: float = -0.2,
        temperature_scale: float = -0.77

    ):
        # Initialize the alpha-proton charge and mass ratios.
        z_a = ions[0].charge_number
        mu_a = ions[0].mass_number

        z_b = ions[1].charge_number
        mu_b = ions[1].mass_number

        # Initialise.
        d_r = (r_0 - r_n) / (1. * n_step)

        # Loop.
        for i in range(n_step):

            r = r_n + ((i + 1) * d_r)

            eta = n_a/n_b
            theta = T_a/T_b

            n_a = n_a * (r / r_n) ** density_scale
            v_a = v_a * (r / r_n) ** velocity_scale
            T_a = T_a * (r / r_n) ** temperature_scale


            d_theta = ((T_b*dT_dt(mu_a, mu_b, z_a, z_b, n_a, n_b, T_a, T_b) - T_a*dT_dt(mu_b, mu_a, z_b, z_a, n_b, n_a, T_b, T_a))/(v_a*T_b**2))*d_r

            theta = theta + d_theta
        return theta


    def c_log(mu_a, mu_b, Z_a, Z_b, n_a, n_b, T_a, T_b):
        alpha = Z_a*Z_b*(mu_a + mu_b)
        beta = mu_a*T_b + mu_b*T_a
        def v_therm(n, Z, T):
            return (n*Z**2)/T

        return 9 + np.log((alpha/beta)*np.sqrt(v_therm(n_a, Z_a, T_a) + v_therm(n_b, Z_b, T_b)))

    def dT_dt(mu_a, mu_b, Z_a, Z_b, n_a, n_b, T_a, T_b):
        alpha = (np.sqrt(mu_a*mu_b)*((Z_a*Z_b)**2)*n_b)
        beta = (mu_a*T_b + mu_b*T_a)**(1.5)
        charlie = (T_b - T_a)
        delta = c_log(mu_a, mu_b, Z_a, Z_b, n_a, n_b, T_a, T_b)

        return (0.174)*(alpha/beta)*charlie*delta

    return 


print("Hello")