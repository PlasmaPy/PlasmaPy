from abc import ABC, abstractmethod
from typing import NoReturn, Union
from numbers import Real, Integral
from astropy import units as u


class AbstractSimulation(ABC):
    """
    A prototype abstract interface for numerical simulations.

    Notes
    -----
    This interface has unstable eigenfunctions, and is therefore subject
    to change without notice.
    """

    @abstractmethod
    def check_inputs(self) -> NoReturn:  # rename to validate_inputs?
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
    def simulate(self) -> NoReturn:  # Return a class for the simulation output?
        """Perform the actual simulation."""

        pass

    @abstractmethod
    def finalize(self) -> NoReturn:
        """Perform the steps to close the simulation and output data."""

        pass


class AbstractTimeDependentSimulation(AbstractSimulation):
    """
    A prototype abstract interface for time-dependent numerical simulations.

    Notes
    -----
    This interface has unstable eigenfunctions, and is therefore subject
    to change without notice.
    """

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
