from plasmapy.classes.plasma_base import BasePlasma

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
                     magnetic_field,
                     ):
        self.electron_density = electron_density
        self.electron_temperature = electron_temperature
        self.ion_density = ion_density
        self.ion_temperature = ion_temperature
        self.magnetic_field = magnetic_field

