import numpy as np
from .particle_class import Particle
from .particle_input import particle_input
from ..utils import AtomicError

class IonizationState:
    """
    Describe the ionization state of a single element.
    """
    @particle_input(require='element', exclude=['isotope', 'ion'])
    def __init__(self,
                 element: Particle,
                 distribution=None,
                 T_e=None,  # not yet implemented!
                 ):
        """
        Initialize the ionization state distribution.

        Notes
        -----
        The method to calculate collisional ionization equilibrium has
        not yet been implemented and will raise a `NotImplementedError`.

        """
        self._element = element

        if distribution is None:
            raise NotImplementedError(
                "IonizationState instances are not yet able to "
                "calculate collisional ionization equilibrium.")

        if isinstance(distribution, (tuple, list)):
            distribution = np.array(distribution)
        elif not isinstance(distribution, np.ndarray):
            raise TypeError

        if len(distribution) != element.atomic_number + 1:
            raise AtomicError("Incorrect number of ionization states.")

        ionsum = np.sum(distribution)
        if not np.isclose(ionsum, 1.0):
            raise ValueError(
                f"The sum of the ionization states of {self.element} "
                f"is {ionsum}, which does not approximately equal one.")

        self._ionization_state = distribution

    def __getitem__(self, value):
        if isinstance(value, int) and 0 <= value <= self.atomic_number:
            return self.ionization_state[value]
        elif isinstance(value, Particle):
            if value.element == self.element:
                return self.ionization_state[value.integer_charge]

        # TODO: Add slicing

    def __str__(self):
        return f"IonizationStates of {self.element}"

    def __eq__(self, other):
        ...



    @property
    def atomic_number(self):
        return self._element.atomic_number

    @property
    def nstates(self):
        return self.atomic_number + 1

    @property
    def element(self):
        """
        Return the `~plasmapy.atomic.Particle` instance corresponding to
        the element.
        """
        return self._element.element

    @property
    def ionization_state(self):
        return self._ionization_state

#    @ionization_state.getter
#    def ionization_state(self):
#        return self.