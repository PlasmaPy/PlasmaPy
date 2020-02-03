from abc import ABC, abstractmethod
from typing import NoReturn, Union
from numbers import Real, Integral
from astropy import units as u


class AbstractSimulation(ABC):

    @abstractmethod
    def check_inputs(self) -> NoReturn:
        """
        Check the inputted simulation parameters.  Raise an exception or
        issue a warning if any potential or actual problems are found.
        """

        pass

    @abstractmethod
    def summarize(self) -> str:
        """
        Return a string containing a summary of the simulation parameters
        and status.
        """

        pass

    @abstractmethod
    def initialize(self) -> NoReturn:
        """Prepare the simulation to be run."""

        pass

    @abstractmethod
    def simulate(self) -> NoReturn:
        """Perform the actual simulation."""

        pass

    @abstractmethod
    def finalize(self) -> NoReturn:
        """Perform the steps to close the simulation and """

        pass


class AbstractTimeDependentSimulation(AbstractSimulation):

    @abstractmethod
    @property
    def start_time(self) -> Union[Real, u.Quantity]:
        """The start time of the simulation."""

        pass

    @abstractmethod
    @property
    def end_time(self) -> Union[Real, u.Quantity]:
        """The end time of the simulation."""

        pass


    @abstractmethod
    @property
    def time_step(self) -> Union[Real, u.Quantity]:
        """Return the time step of the simulation."""

        pass

    @abstractmethod
    def time_advance(self) -> NoReturn:
        """Perform one step of the time advance."""

        pass

    @abstractmethod
    @property
    def step_number(self) -> Integral:
        """The current step number of the simulation."""

        pass


