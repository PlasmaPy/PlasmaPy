import pytest
import numpy as np
from ..ionization_states import IonizationState
from ...utils import AtomicError, RunTestError, InvalidIsotopeError
from ...atomic import (
    atomic_number,
    atomic_symbol,
    isotope_symbol,
    Particle,
)

test_cases = {
    'Li': {
        'particle': 'Li',
        'ionic_fractions': np.array([0.4, 0.3, 0.2, 0.1]),
        'tol': 1e-15,
    },

    'Li ground state': {
        'particle': 'Li',
        'ionic_fractions': np.array([1, 0, 0, 0], dtype=np.int64),
        'tol': 1e-15,
    },

    'H': {
        'particle': 'H',
        'ionic_fractions': [0.6, 0.4],
        'tol': 1e-8,
    },

    'H acceptable error': {
        'particle': 'H',
        'ionic_fractions': [0.6, 0.400_000_001],
        'tol': 1e-8,
    },

}


class Test_IonizationState:
    """Test instances of IonizationState."""

    @classmethod
    def setup_class(cls):
        "Set up the class and test instantantiation."
        cls.instances = {}
        for test_case in test_cases.keys():
            try:
                cls.instances[test_case] = IonizationState(**test_cases[test_case])
            except Exception as exc:
                raise RunTestError(
                    f"Unable to create IonizationState instance for "
                    f"test case {test_case}.")

    def test_equality1(self):
        assert self.instances['Li'] == self.instances['Li'], \
            "Identical IonizationState instances are not equal."

    def test_equality2(self):
        assert self.instances['H'] == self.instances['H acceptable error'], \
            ("Two IonizationState instances that are approximately the "
             "same to within the tolerance are not testing as equal.")

    def test_inequality(self):
        assert self.instances['Li ground state'] != self.instances['Li'], \
            "Different IonizationState instances are equal."

    def test_equality_error(self):
        """
        Test that comparisons of IonizationState instances for
        different elements
        """
        with pytest.raises(AtomicError):
            self.instances['Li'] == self.instances['H']

    @pytest.mark.parametrize('test_name', test_cases.keys())
    def test_iteration(self, test_name: str):
        """Test that IonizationState instances iterate impeccably."""
        try:
            states = [state for state in self.instances[test_name]]
        except Exception:
            raise AtomicError(
                f"Unable to perform iteration for {test_name}.")

        try:
            integer_charges = [state.integer_charge for state in states]
            ionic_fractions = np.array([state.ionic_fraction for state in states])
            ionic_symbols = [state.ionic_symbol for state in states]
        except Exception:
            raise AtomicError("An attribute may be misnamed or missing.")

        try:
            base_symbol = isotope_symbol(ionic_symbols[0])
        except InvalidIsotopeError:
            base_symbol = atomic_symbol(ionic_symbols[0])
        finally:
            atomic_numb = atomic_number(ionic_symbols[1])

        errors = []

        expected_charges = np.arange(atomic_numb + 1)
        if not np.all(integer_charges == expected_charges):
            errors.append(
                f"The resulting integer charges are {integer_charges}, "
                f"which are not equal to the expected integer charges, "
                f"which are {expected_charges}.")

        expected_fracs = test_cases[test_name]['ionic_fractions']
        if not np.allclose(ionic_fractions, expected_fracs):
            errors.append(
                f"The resulting ionic fractions are {ionic_fractions}, "
                f"which are not equal to the expected ionic fractions "
                f"of {expected_fracs}.")

        expected_particles = [Particle(base_symbol, Z=charge) for charge in integer_charges]
        expected_symbols = [particle.ionic_symbol for particle in expected_particles]
        if not ionic_symbols == expected_symbols:
            errors.append(
                f"The resulting ionic symbols are {ionic_symbols}, "
                f"which are not equal to the expected ionic symbols of "
                f"{expected_symbols}.")

        if errors:
            errors.insert(0, (
                f"The test of IonizationState named '{test_name}' has "
                f"resulted in the following errors when attempting to "
                f"iterate."))
            errmsg = " ".join(errors)
            raise AtomicError(errmsg)


