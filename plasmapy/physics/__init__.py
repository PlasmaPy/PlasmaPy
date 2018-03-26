r"""
Values should be returned as an Astropy Quantity in SI units.

If a quantity has several names, then the function name should be the
one that provides the most physical insight into what the quantity
represents.  For example, 'gyrofrequency' indicates gyration, while
Larmor frequency indicates that this frequency is somehow related to a
human (or perhaps a cat?) named Larmor.  Similarly, using omega_ce as
a function name for this quantity will make the code less readable to
people who are unfamiliar with the notation or use a different symbol.

The docstrings for plasma parameter methods should describe the
physics associated with these quantities in ways that are
understandable to students who are taking their first course in plasma
physics while still being useful to experienced plasma physicists.

In many cases, units are enough to tell what field a quantity
represents.  The following line is an example.

>>> Alfven_speed(5*u.T, 8e-7*u.kg/u.m**3)

"""
from .parameters import (Alfven_speed,
                         ion_sound_speed,
                         thermal_speed,
                         kappa_thermal_speed,
                         gyrofrequency,
                         gyroradius,
                         plasma_frequency,
                         Debye_length,
                         Debye_number,
                         inertial_length,
                         magnetic_pressure,
                         magnetic_energy_density,
                         upper_hybrid_frequency,
                         lower_hybrid_frequency)

from .dielectric import (cold_plasma_permittivity_LRP,
                         cold_plasma_permittivity_SDP)

from .distribution import (Maxwellian_1D,
                           Maxwellian_velocity_3D,
                           Maxwellian_speed_1D,
                           Maxwellian_speed_3D,
                           kappa_velocity_3D,
                           kappa_velocity_1D)

from .quantum import (deBroglie_wavelength,
                      thermal_deBroglie_wavelength,
                      Fermi_energy,
                      Thomas_Fermi_length,
                      Wigner_Seitz_radius,
                      chemical_potential,
                      chemical_potential_interp)

from .relativity import Lorentz_factor

from .transport.collisions import (Coulomb_logarithm,
                                   b_perp,
                                   impact_parameter,
                                   collision_frequency,
                                   mean_free_path,
                                   mobility,
                                   Knudsen_number,
                                   coupling_parameter)
from plasmapy.physics.transport.braginskii import classical_transport
from . import transport
