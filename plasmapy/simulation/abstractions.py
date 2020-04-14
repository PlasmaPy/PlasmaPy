"""Abstract classes for numerical simulations."""

__all__ = ["AbstractSimulation", "AbstractTimeDependentSimulation"]

from abc import ABC, abstractmethod
from typing import NoReturn


class AbstractSimulation(ABC):
    """
    A prototype abstract interface for numerical simulations.

    Notes
    -----
    This interface is incomplete and unstable, and is thus subject to
    change at any time.
    """

    @abstractmethod
    def summarize(self):
        """
        Print out a summary of the simulation parameters and status.
        """

        pass

    @abstractmethod
    def initialize(self) -> NoReturn:
        """Prepare the simulation to be run."""

        pass

    @abstractmethod
    def simulate(self):
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
    This interface is incomplete and unstable, and is thus subject to
    change at any time.
    """

    pass
