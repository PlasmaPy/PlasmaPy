# 'physics' is a tentative name for this subpackage.  Another
# possibility is 'plasma'.  The organization is to be decided by v0.1.

from .parameters import Alfven_speed, ion_sound_speed, thermal_speed, kappa_thermal_speed, gyrofrequency, gyroradius, plasma_frequency, Debye_length, Debye_number, inertial_length, magnetic_pressure, magnetic_energy_density, upper_hybrid_frequency, lower_hybrid_frequency
from .quantum import deBroglie_wavelength, thermal_deBroglie_wavelength, Fermi_energy, Thomas_Fermi_length
from .relativity import Lorentz_factor
from .transport import Coulomb_logarithm, classical_transport
from .distribution import Maxwellian_1D, Maxwellian_velocity_3D, Maxwellian_speed_1D, Maxwellian_speed_3D,kappa_velocity_3D, kappa_velocity_1D
from .dielectric import cold_plasma_permittivity_LRP, cold_plasma_permittivity_SDP
