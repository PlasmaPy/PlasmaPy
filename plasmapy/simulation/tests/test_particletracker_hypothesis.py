try: 
    from hypothesis import given, settings
    from hypothesis.extra.numpy import arrays
except ImportError:
    import pytest
    pytest.skip("Optional hypothesis test")

from astropy import units as u
from plasmapy.classes.sources.analyticalplasma import AnalyticalPlasma
from plasmapy.simulation import ParticleTracker
from plasmapy.physics.parameters import gyrofrequency
import numpy as np

def magnetic_field(r):
    return u.Quantity([[4, 0, 0]]*len(r), u.T)

# precomputed for efficiency
E_unit = u.V / u.m
def electric_field(r):
    return u.Quantity([[0, 2, 0]]*len(r), E_unit)

N = 20
@given(arrays(dtype='float', shape=(N, 2)))
@settings(deadline=None)
def test_no_underflow(velocity):
    plasma = AnalyticalPlasma(magnetic_field, electric_field)
    freq = gyrofrequency(4 * u.T, 'p').to(u.Hz, equivalencies=u.dimensionless_angles())
    gyroperiod = (1/freq).to(u.s)
    steps_to_gyroperiod = 10
    timestep = gyroperiod / steps_to_gyroperiod
    number_steps = 5 * steps_to_gyroperiod * int(2 * np.pi)

    trajectory = ParticleTracker(plasma, 'p', N, 1, timestep/100, number_steps)
    trajectory._v[:, :2] = velocity
    trajectory.run()
