"""Classes that describes normalizations for different systems of equations."""

import abc
import astropy.units as u

class AbstractNormalizations(abc.ABC):

    @abc.abstractmethod
    def magnetic_field(self) -> u.T:
        raise NotImplementedError

    @abc.abstractmethod
    def length(self) -> u.m:
        raise NotImplementedError

    @abc.abstractmethod
    def time(self) -> u.s:
        raise NotImplementedError

    @abc.abstractmethod
    def velocity(self) -> u.m / u.s:
        raise NotImplementedError

    @abc.abstractmethod
    def current_density(self) -> u.A / u.m ** 2:
        raise NotImplementedError

    @abc.abstractmethod
    def pressure(self) -> u.Pa:
        raise NotImplementedError

    @abc.abstractmethod
    def electric_field(self) -> u.V / u.m:
        raise NotImplementedError
