from plasmapy.classes.plasma_base import BasePlasma
import astropy.units as u

class Plasma3DSpatialOnly(BasePlasma):
    def electron_temperature(self, r):
        raise NotImplementedError
    def ion_temperature(self, r):
        raise NotImplementedError
    def electron_density(self, r):
        raise NotImplementedError
    def ion_density(self, r):
        raise NotImplementedError


class AnalyticTwoFluidPlasma(Plasma3DSpatialOnly):
    def __init__(self, electron_density,
                     electron_temperature,
                     ion_density,
                     ion_temperature,
                     electron_velocity = lambda x: u.Quantity([0, 0, 0], unit=u.m/u.s),
                     ion_velocity = lambda x: u.Quantity([0, 0, 0], unit=u.m/u.s),
                     magnetic_field = lambda x: u.Quantity([0, 0, 0], unit=u.T),
                     electric_field = lambda x: u.Quantity([0, 0, 0], unit=u.V/u.m)
                     ):
        self.electron_density = electron_density
        self.electron_temperature = electron_temperature
        self.ion_density = ion_density
        self.ion_temperature = ion_temperature
        self.electron_velocity = electron_velocity
        self.ion_velocity = ion_velocity
        self.electric_field = electric_field
        self.magnetic_field = magnetic_field

