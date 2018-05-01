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
                 T_e=None,
                 atol=1e-16,
                 ):
        """
        Initialize the ionization state distribution.

        Parameters
        ----------
        element: str, int, or ~plasmapy.atomic.Particle
            A `str` or `~plasmapy.atomic.Particle` instance representing
            an element, or an `int` representing the atomic number of an
            element.

        distribution: ~numpy.ndarray, list, or tuple
            The ionization fractions of an element.  This parameter must
            have a length equal to one more than the atomic number, and
            must sum to one within an absolute tolerance of `atol`.

        T_e: ~astropy.units.Quantity, optional
            The electron temperature.

        atol: float, optional
            The absolute tolerance of how much the sum of the ionization
            states may differ from one.

        Notes
        -----
        The method to calculate collisional ionization equilibrium has
        not yet been implemented and will raise a `NotImplementedError`.

        """
        self._element = element

        if distribution is None or T_e is not None:
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
        if not np.isclose(ionsum, 1.0, atol=atol):
            raise ValueError(
                f"The sum of the ionization states of {self.element} "
                f"is {ionsum}, which does not approximately equal one.")

        self._ionization_state = distribution

    def __getitem__(self, value):
        """Return the ionization fraction(s)."""
        if isinstance(value, slice):
            return self._ionization_state[value.start:value.stop:value.step]
        elif isinstance(value, int) and 0 <= value <= self.atomic_number:
            return self.ionization_state[value]
        elif isinstance(value, Particle):
            if value.element == self.element:
                return self.ionization_state[value.integer_charge]


    def __str__(self):
        return f"IonizationStates of {self.element}"

    def __repr__(self):
        return f"IonizationStates of {self.element}"

    def __eq__(self, other):
        raise NotImplementedError

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

    def set_collisional_equilibrium(self, T_e=None):
        raise NotImplementedError
