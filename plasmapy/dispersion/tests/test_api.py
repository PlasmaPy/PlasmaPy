import astropy.units as u
from astropy import constants
import pytest

@pytest.mark.xfail(strict=True)
def test_dispersion_api():
    """
    This is a broad overview of the initial ideas we had on the sample API for the dispersion solver
    as per https://github.com/PlasmaPy/plasmapy-meeting-notes/blob/master/2019-09-24-design-meeting.md

    It is marked as xfail - ideally a later PR will cause this 
    """

    p = Particle('p')

    # class Species(Particle):
    #     def __init__(self, particle_name, temperature, anisotropy, density, drift):
    #         ...
        
    from plasmapy.atomic import CustomParticle # probably the best place for it but negotiable
    custom_ion = CustomParticle(mass=25*electron_mass, charge=2 * constants.e.si)

    from plasmapy.dispersion import Species # place negotiable? better name?

    species0 = Species(custom_ion, 
        anisotropy = 1/3.1415,
        # or solve for anisotropy from T_perpendicular, T_parallel, but for now - it's fine to supply it
        density = 1e23 * u.m**-3,
        # in general, example: drift = [0, 0, 2] * (u.m/u.s), but for now don't worry about direction
        drift = 2 * (u.m/u.s),
    )

    # physics initialization 
    from plasmapy.dispersion import PhysicalInputs #, DispersionSolver
    dispersion_input = PhysicalInputs(species = [species0],
        magnetic_field = u.Quantity([0, 0, 5], u.T),
        # kwargs:
        wavenumber = np.linspace(0, 10, 1000) / u.m,
        # wavenumber = np.array([0.1, 10, 3,]) / u.m,
        angle = np.linspace(0, np.pi, 1000),
        # v_Alfven/c needs to be controlled - magnetic_field, density give us v
        # TODO what do we need for the relativistic dispersion?
        # TODO Ask Daniel Verscharen and other astrophysicists
        )    

    # class DispersionSolver:
    #     solvers = {'kinetic': KineticDispersionSolver,
    #                'two_fluid': tfst}
        
    #     def __init__(self, solver_type, **solver_kwargs):
    #         # doc issues but easy to extend
        
    from plasmapy.dispersion import KineticDispersionSolver
    kds = KineticDispersionSolver(dispersion_input,
                                  number_of_iterations = 10000,         
                                  #Threshold for the determinant of the dispersion tensor:
                                  # If det D <= det_D_threshold, the Newton iteration will 
                                  # be stopped
                                  det_D_threshold = 0.0001,  # TODO put reasonable float value here
                                  # Maximum of sum in Bessel functions (both regular 
                                  # and modified)
                                  # can be very low (e.g., 3) for quasi-parallel propagation
                                  nmax=1000,
                                  # If I_n is less than this value, higher n are neglected:
                                  Bessel_zero=1, # was 1.d-50, I'm unsure of what this syntax represents
                                  # Initial guess for the frequency (complex number):
                                  initial_guess=0.01+0.0001j,
                                  # Amplitude mode: If 1, then ampl = delta B_x / B0. 
                                  # If 2, then ampl = delta B_y / B0.
                                  # If 3, then ampl = delta B_z / B0
                                  # ampl_mode=3,
                                  # Amplitude for the calculation of polarization
                                  # properties:
                                  # ampl=1.d0
                                  amplBETTERNAME = (0, 0, 1)
                                  )
                     
    from plasmapy.dispersion import FluidDispersionSolver
    fds = FluidDispersionSolver(dispersion_input)
                     
    # def two_fluid_single_theta(beta=0.6, v_Alfven=1., me/mi=0.000545, theta=0., kmin=1e-2, kmax=10., npoints=200,wrt='n')
    # c/wpi=1, V_alfven=1, wci=1.
    # 
    # how to handle multifluid dispersion solvers?

